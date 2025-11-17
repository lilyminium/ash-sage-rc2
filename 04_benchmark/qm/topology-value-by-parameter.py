from collections import defaultdict
import typing
import pathlib

import click
import tqdm

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
    default="topology-labeled-by-parameter",
    help="Directory to write output files to.",
)
@click.option(
    "--topology-group",
    "-t",
    type=click.Choice(["Bonds", "Angles", "ProperTorsions", "ImproperTorsions"], case_sensitive=True),
    default="Bonds",
    help="The topology group to label.",
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
    topology_group: typing.Literal["Bonds", "Angles", "ProperTorsions", "ImproperTorsions"] = "Bonds",
    input_labels_directory: str = "topology-labels/fb-fit-v3-single-mean-k100_unconstrained/",
    forcefields: list | None = None,
):
    dataset = ds.dataset(input_directory).filter(
        pc.field("topology") == topology_group
    )
    labels = ds.dataset(input_labels_directory).filter(
        pc.field("topology") == topology_group
    )
    ALL_LABELS = defaultdict(dict)
    for row in labels.to_table().to_pylist():
        ALL_LABELS[row["mapped_smiles"]][tuple(row["atom_indices"])] = row["parameter_id"]

    # filter for FFs
    if not forcefields:
        forcefields = dataset.to_table(
            columns=["method"]
        ).to_pydict()["method"]
        forcefields = sorted(set(forcefields) - {"qm"})

    qm_pylist = dataset.filter(
        pc.field("method") == "qm"
    ).to_table().to_pylist()
    qm_pylist_by_smiles = defaultdict(dict)
    for row in qm_pylist:
        qm_pylist_by_smiles[row["qcarchive_id"]][tuple(row["atom_indices"])] = row["value"]

    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    # label molecules
    for ff in forcefields:
        ffdir = output_directory / ff
        ffdir.mkdir(parents=True, exist_ok=True)

        subset = dataset.filter(
            pc.field("method") == ff
        )
        entries = []
        for row in tqdm.tqdm(subset.to_table().to_pylist(), desc=f"Labeling {topology_group} for {ff}"):
            parameter_labels_for_molecule = ALL_LABELS[row["mapped_smiles"]]
            qm_values_for_molecule = qm_pylist_by_smiles[row["qcarchive_id"]]
            value = row["value"]
            parameter_id = parameter_labels_for_molecule[tuple(row["atom_indices"])]
            qm_value = qm_values_for_molecule[tuple(row["atom_indices"])]
            entry = dict(row)
            if topology_group == "Bonds":
                difference = value - qm_value
            else:
                difference = compute_angle_difference(value, qm_value)
            entry.update(
                {
                    "parameter_id": parameter_id,
                    "reference_value": qm_value,
                    "difference": difference,
                    "abs_difference": abs(difference),
                }
            )
            entries.append(entry)
        
        table = pa.Table.from_pylist(entries)
        filename = ffdir / f"{topology_group}.parquet"
        pq.write_table(table, filename)
        print(f"Wrote {len(entries)} entries to {filename}")

if __name__ == "__main__":
    main()
