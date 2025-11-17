import pathlib
import typing

import click
import tqdm
from loguru import logger

import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.compute as pc
import pyarrow.parquet as pq

from openff.toolkit import Molecule, ForceField, unit

UNITS = {
    "Bonds": unit.angstrom,
    "Angles": unit.degree,
    "ProperTorsions": unit.degree,
}


@click.command()
@click.option(
    "--input-directory",
    "-i",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    default="topology-values",
    help="Directory containing topology values to label.",
)
@click.option(
    "--forcefield",
    "-f",
    type=str,
    default="../forcefields/fb-fit-v3-single-mean-k100_unconstrained.offxml",
    help="Force field to use for labeling.",
)
@click.option(
    "--output-directory",
    "-o",
    type=click.Path(exists=False, file_okay=False, dir_okay=True),
    default="topology-labels",
    help="Directory to write output files to.",
)
def main(
    input_directory: str = "topology-values",
    forcefield: str = "../forcefields/fb-fit-v3-single-mean-k100_unconstrained.offxml",
    output_directory: str = "topology-labels",
):
    dataset = ds.dataset(input_directory).filter(
        pc.field("method") == "qm"
    )

    # load smiles
    mapped_smiles: list[str] = sorted(
        set(
            dataset.to_table(
                columns=["mapped_smiles"]
            ).to_pydict()["mapped_smiles"]
        ),
        key=len,
        reverse=True
    )

    # load force field
    ff = ForceField(forcefield)


    ff_name = pathlib.Path(forcefield).stem
    output_directory = pathlib.Path(output_directory) / ff_name
    output_directory.mkdir(parents=True, exist_ok=True)

    # label molecules
    for i, mapped_smi in enumerate(tqdm.tqdm(mapped_smiles)):
        if mapped_smi == "[H:12][C:1]1=[C:3]([S:10][C:4](=[C:2]1[H:13])[S:11](=[O:8])(=[O:9])[C:7]([H:21])([C:5]([H:15])([H:16])[H:17])[C:6]([H:18])([H:19])[H:20])[H:14]":
            print(i, mapped_smi)
        all_rows = []
        mol = Molecule.from_mapped_smiles(
            mapped_smi,
            allow_undefined_stereo=True
        )
        labels = ff.label_molecules(
            mol.to_topology(),
        )[0]

        for parameter_type in ["Bonds", "Angles", "ProperTorsions"]:
            param_labels = labels[parameter_type]
            for atom_indices, parameter in param_labels.items():
                entry = {
                    "mapped_smiles": mapped_smi,
                    "method": ff_name,
                    "topology": parameter_type,
                    "atom_indices": list(atom_indices),
                    "parameter_id": parameter.id,
                }
                all_rows.append(entry)

        # parameter_type = "ImproperTorsions"
        # param_labels = labels[parameter_type]
        # for atom_indices, parameter in param_labels.items():
        #     central = atom_indices[1]
        #     non_central_indices = [
        #         atom_indices[0], atom_indices[2], atom_indices[3]
        #     ]
        #     for i, j, k in [
        #         (0, 1, 2),
        #         (1, 2, 0),
        #         (2, 0, 1),
        #     ]:
        #         entry = {
        #             "mapped_smiles": mapped_smi,
        #             "method": ff_name,
        #             "topology": parameter_type,
        #             "atom_indices": [non_central_indices[i], central, non_central_indices[j], non_central_indices[k]],
        #             "parameter_id": parameter.id,
        #         }
        #         all_rows.append(entry)
                

        table = pa.Table.from_pylist(all_rows)
        output_path = output_directory / f"{i:06d}.parquet"
        pq.write_table(table, output_path)
        # logger.info(f"Wrote {len(all_rows)} labels to {output_path}")

if __name__ == "__main__":
    main()
