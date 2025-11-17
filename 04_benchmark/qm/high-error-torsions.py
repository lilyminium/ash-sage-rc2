from collections import defaultdict
import typing
import pathlib

import click
import tqdm

import numpy as np
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.compute as pc
import pyarrow.parquet as pq

from openff.toolkit import Molecule, ForceField, unit

def compute_angle_difference(mm_value, qm_value):
    difference = mm_value - qm_value
    if difference > 180:
        difference -= 360
    if difference < -180:
        difference += 360
    return difference


@click.command()
@click.option(
    "--input-directory",
    "-i",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    default="topology-values",
    help="Directory containing topology values to label.",
)
@click.option(
    "--output-directory",
    "-o",
    type=click.Path(exists=False, file_okay=False, dir_okay=True),
    default="torsion-error-with-parameters/v3-k100",
    help="Directory to write output files to.",
)
@click.option(
    "--input-labels-directory",
    "-l",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    default="topology-labels/fb-fit-v3-single-mean-k100_unconstrained/",
    help="Directory containing topology labels.",
)
@click.option(
    "--ff",
    "--forcefield",
    "forcefields",
    multiple=True,
    type=str,
    help="Force field(s) to label. If not provided, all force fields in the input directory will be labeled.",
)
def main(
    input_directory: str = "topology-values",
    output_directory: str = "topology-labeled-by-parameter",
    input_labels_directory: str = "topology-labels/fb-fit-v3-single-mean-k100_unconstrained/",
    forcefields: list | None = None,
):
    dataset = ds.dataset(input_directory).filter(
        pc.field("topology") == "ProperTorsions"
    )
    labels = ds.dataset(input_labels_directory).filter(
        pc.field("topology") == "ProperTorsions"
    )

    # label torsions by central bond
    ALL_LABELS = defaultdict(lambda: defaultdict(list))
    for row in labels.to_table().to_pylist():
        atom_indices = tuple(row["atom_indices"])
        sorted_central_bond = tuple(sorted(atom_indices[1:3]))
        ALL_LABELS[row["mapped_smiles"]][sorted_central_bond].append(row["parameter_id"])

    # filter for FFs
    if not forcefields:
        forcefields = dataset.to_table(
            columns=["method"]
        ).to_pydict()["method"]
        forcefields = sorted(set(forcefields) - {"qm"})

    qm_pylist = dataset.filter(
        pc.field("method") == "qm"
    ).to_table().to_pylist()
    qm_pylist_by_smiles = defaultdict(lambda: defaultdict(dict))
    for row in qm_pylist:
        sorted_central_bond = tuple(sorted(row["atom_indices"][1:3]))
        qm_pylist_by_smiles[row["qcarchive_id"]][sorted_central_bond][tuple(row["atom_indices"])] = row["value"]

    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    # label molecules
    for ff in forcefields:

        subset = dataset.filter(
            pc.field("method") == ff
        )
        entries = []
        ffdf = subset.to_table().to_pandas()
        central_bonds = []
        all_atom_indices = []
        for _, row in ffdf.iterrows():
            atom_indices = tuple(row["atom_indices"])
            all_atom_indices.append(atom_indices)
            central_bond = tuple(sorted(atom_indices[1:3]))
            central_bonds.append(central_bond)
        ffdf["central_bond"] = central_bonds
        ffdf["all_atom_indices"] = all_atom_indices

        for (central_bond, qcarchive_id), cbdf in tqdm.tqdm(
            ffdf.groupby(["central_bond", "qcarchive_id"]),
            desc=f"Computing high-error torsions for {ff}",
        ):
            cbdf = cbdf.sort_values("all_atom_indices")
            mapped_smi = cbdf.iloc[0]["mapped_smiles"]
            parameter_labels_for_molecule = ALL_LABELS[mapped_smi][central_bond]
            assert len(parameter_labels_for_molecule) > 0, f"No parameter labels found for {mapped_smi} central bond {central_bond}: {ALL_LABELS[mapped_smi].keys()}"
            qm_values_for_molecule = qm_pylist_by_smiles[qcarchive_id][central_bond]
            mm_values = cbdf["value"].values
            qm_values = np.zeros_like(mm_values)
            for i, (_, row) in enumerate(cbdf.iterrows()):
                atom_indices = tuple(row["atom_indices"])
                qm_value = qm_values_for_molecule[atom_indices]
                qm_values[i] = qm_value
            differences = np.array([
                compute_angle_difference(mm, qm)
                for mm, qm in zip(mm_values, qm_values)
            ])
            mean_signed_error = differences.mean()
            mean_unsigned_error = np.abs(differences).mean()
            rmse = np.sqrt((differences ** 2).mean())
            entry = {
                "qcarchive_id": qcarchive_id,
                "mapped_smiles": mapped_smi,
                "method": ff,
                "topology": "ProperTorsions",
                "atom_indices": list(central_bond),
                "parameter_ids": parameter_labels_for_molecule,
                "mean_signed_error": mean_signed_error,
                "mean_unsigned_error": mean_unsigned_error,
                "rmse": rmse,
            }
            entries.append(entry)
        
        table = pa.Table.from_pylist(entries)
        filename = output_directory / f"{ff}.parquet"
        pq.write_table(table, filename)
        print(f"Wrote {len(entries)} entries to {filename}")

if __name__ == "__main__":
    main()
