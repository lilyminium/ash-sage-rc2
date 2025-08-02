"""
This script runs through a set of tables containing optimization data
and labels each of them with the appropriate force field parameters.
It saves the labeled data in a specified output directory.

This step is needed for the data selection process.

Note: only bonds, angles, and proper and improper torsions are labeled.
"""

import logging
import sys

import collections
import pathlib
import tqdm
import click

import pyarrow as pa
import numpy as np
import pyarrow.compute as pc
import pyarrow.parquet as pq
import pyarrow.dataset as ds

from openff.toolkit import Molecule, ForceField


logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    stream=sys.stdout
)

def get_unique_smiles(input_file: str, column: str = "smiles") -> set[str]:
    # Load the table
    table = pq.read_table(input_file)
    logger.info(f"Loaded {table.num_rows} rows from {input_file}")

    # Get the unique smiles
    unique_smiles = set(pc.unique(table.column(column)))
    logger.info(f"Loaded {len(unique_smiles)} unique {column}")
    return unique_smiles


def label_table_with_forcefield(
    forcefield_file: str,
    input_file: str,
    output_directory: str,
    forcefield_name: str = None,
):
    forcefield = ForceField(forcefield_file)
    unique_smiles = get_unique_smiles(input_file, "cmiles")

    file_number = 0

    # Read existing data
    output_path = pathlib.Path(output_directory)
    if not output_path.exists():
        output_path.mkdir(parents=True, exist_ok=True)
    else:
        existing_dataset = ds.dataset(output_path)
        if existing_dataset.count_rows():
            subset = existing_dataset.filter(
                pc.field("forcefield") == forcefield_name
            )
            existing_smiles = set(
                subset.to_table(columns=["cmiles"]).to_pydict()["cmiles"]
            )
            logger.info(f"Loaded {len(existing_smiles)} existing cmiles")
            unique_smiles = unique_smiles - existing_smiles
            logger.info(f"New cmiles: {len(unique_smiles)}")
            file_number = len(existing_dataset.files)
    
    entries = []

    for smiles in tqdm.tqdm(unique_smiles):
        try:
            mol = Molecule.from_mapped_smiles(
                str(smiles),
                allow_undefined_stereo=True,
            )
        except:
            continue
        labels = forcefield.label_molecules(mol.to_topology())[0]

        for parameter_type in ["Bonds", "Angles", "ProperTorsions", "ImproperTorsions"]:
            parameter_indices = collections.defaultdict(list)
            for indices, parameter in labels[parameter_type].items():
                parameter_indices[parameter.id].append(list(indices))
            for k, v in parameter_indices.items():
                ix = np.array(v).flatten().tolist()
                entry = {
                    "cmiles": smiles,
                    "forcefield": forcefield_name,
                    "parameter_type": parameter_type,
                    "parameter_id": k,
                    "parameter_indices": ix,
                    "parameter_indices_str": "-".join([str(x) for x in ix]),
                }
                entries.append(entry)
    
    # Create a new table and save
    if entries:
        new_table = pa.Table.from_pylist(entries)
        new_filename = output_path / f"batch-{file_number:04d}.parquet"
        assert not new_filename.exists(), f"File {new_filename} already exists"
        pq.write_table(new_table, new_filename)
        logger.info(f"Wrote {len(entries)} rows to {new_filename}")


@click.command()
@click.option(
    "--input-directory",
    "-i",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    default="data/tables/optimization",
    help="Directory containing the input tables."
)
@click.option(
    "--output-directory",
    "-o",
    type=click.Path(exists=False, file_okay=False, dir_okay=True),
    default="parameters/valence",
    help="Directory to save the labeled tables."
)
@click.option(
    "--forcefield-file",
    "-ff",
    type=click.Path(exists=True, dir_okay=False),
    default="openff_unconstrained-2.2.1.offxml",
    help="Path to the forcefield file."
)
@click.option(
    "--forcefield-name",
    "-ffn",
    type=str,
    default=None,
    help="Name of the forcefield."
)
def main(
    input_directory: str = "data/tables/optimization",
    output_directory: str = "parameters/valence",
    forcefield_file: str = "openff_unconstrained-2.2.1.offxml",
    forcefield_name: str = None,
):
    table_files = sorted(pathlib.Path(input_directory).glob("*.parquet"))
    logger.info(f"Found {len(table_files)} input files")

    if forcefield_name is None:
        forcefield_name = pathlib.Path(forcefield_file).stem

    output_directory = pathlib.Path(output_directory) / forcefield_name
    output_directory.mkdir(parents=True, exist_ok=True)

    logger.info(f"Using forcefield {forcefield_file}")
    logger.info(f"Output directory: {output_directory}")

    for table_file in tqdm.tqdm(table_files, desc="Processing files"):
        label_table_with_forcefield(
            forcefield_file,
            table_file,
            output_directory,
            forcefield_name=forcefield_name,
        )


if __name__ == "__main__":
    main()
