"""
This script runs through a set of tables containing torsion drive data
and labels each of them with the appropriate force field parameters.
It saves the labeled data in a specified output directory.

This step is needed for the data selection process.
"""

import logging
import sys

import click
import pathlib
import tqdm

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

def label_torsion_table_with_forcefield(
    forcefield_file: str,
    input_file: str,
    output_directory: str,
    forcefield_name: str = None,
):
    forcefield = ForceField(forcefield_file)

    table = pq.read_table(input_file)
    logger.info(f"Loaded {table.num_rows} rows from {input_file}")

    df = table.to_pandas()
    unique_torsion_rows = []

    # we need to split up potential 2D torsions
    for _, row in df.iterrows():
        cmiles = row["cmiles"]
        torsiondrive_id = row["id"]
        torsions = list(row["dihedral"])
        while torsions:
            unique_torsion_rows.append({
                "cmiles": cmiles,
                "id": torsiondrive_id,
                "dihedral": torsions[:4]
            })
            torsions = torsions[4:]

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
            
            unique_torsion_rows = [
                row
                for row in unique_torsion_rows
                if row["cmiles"] not in existing_smiles
            ]
            logger.info(f"New rows: {len(unique_torsion_rows)}")
            file_number = len(existing_dataset.files)
    
    entries = []

    for row in tqdm.tqdm(unique_torsion_rows, desc="Processing rows"):
        try:
            mol = Molecule.from_mapped_smiles(
                str(row["cmiles"]),
                allow_undefined_stereo=True,
            )
        except:
            continue
        labels = forcefield.label_molecules(mol.to_topology())[0]
        torsion_labels = labels["ProperTorsions"]
        dih = tuple(row["dihedral"])
        if dih not in torsion_labels:
            dih = tuple(reversed(dih))
        
        try:
            torsion_id = torsion_labels[dih].id
        except KeyError:
            # logger.warning(f"Could not find torsion label for {dih}")
            continue
        # get all torsions running through central bond
        sorted_central = sorted(dih[1:3])
        all_parameter_ids = []
        for ix, parameter in torsion_labels.items():
            if sorted(ix[1:3]) == sorted_central:
                all_parameter_ids.append(parameter.id)
        all_parameter_ids = sorted(set(all_parameter_ids))

        entry = {
            "id": row["id"],
            "cmiles": row["cmiles"],
            "forcefield": forcefield_name,
            "parameter_id": torsion_id,
            "dihedral": list(dih),
            "all_parameter_ids": all_parameter_ids,
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
    default="data/tables/torsiondrive",
    help="Directory containing the input tables."
)
@click.option(
    "--output-directory",
    "-o",
    type=click.Path(exists=False, file_okay=False, dir_okay=True),
    default="parameters/torsions",
    help="Directory to save the labeled parameters."
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
    input_directory: str = "data/tables/torsiondrive",
    output_directory: str = "parameters/torsions",
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
        label_torsion_table_with_forcefield(
            forcefield_file,
            table_file,
            output_directory,
            forcefield_name=forcefield_name,
        )

if __name__ == "__main__":
    main()
