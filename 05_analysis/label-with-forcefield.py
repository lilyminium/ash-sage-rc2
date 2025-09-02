import pathlib
import sys

import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq

import tqdm
import click
from rdkit import Chem
from loguru import logger
from openff.toolkit import Molecule, ForceField
from yammbs.checkmol import analyze_functional_groups

logger.remove()
logger.add(sys.stdout)

def sanitize_smiles(smiles: str) -> str:
    """Sanitize SMILES"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Chem.MolToSmiles(mol)
    return ""

@click.command()
@click.option(
    "--input",
    "-i",
    "input_path",
    type=click.Path(exists=True, dir_okay=True, file_okay=True),
    help=(
        "Input path. This should reference either "
        "a directory for a PyArrow dataset "
        "or a file for a Pandas dataframe."
    )
)
@click.option(
    "--output",
    "-o",
    "output_path",
    type=click.Path(dir_okay=False, file_okay=True),
    help=(
        "Output path. This should reference a file for a Parquet dataset."
    )
)
@click.option(
    "--smiles-col",
    "-s",
    "smiles_columns",
    multiple=True,
    type=str,
    help=(
        "SMILES columns to label. This should reference one or more columns "
        "in the input dataset that contain SMILES strings."
    )
)
@click.option(
    "--forcefield",
    "-ff",
    "forcefield",
    type=str,
    help=(
        "Path to the forcefield file to use for labeling."
    )
)
@click.option(
    "--forcefield-name",
    "-ffn",
    "forcefield_name",
    type=str,
    help=(
        "Name of the forcefield to use for labeling."
    )
)
def main(
    input_path: str,
    output_path: str,
    forcefield: str,
    forcefield_name: str,
    smiles_columns: list[str] = [],
):
    forcefield = ForceField(forcefield)

    input_path = pathlib.Path(input_path)
    
    # read in pyarrow dataset
    if input_path.is_dir():
        input_dataset = ds.dataset(input_path)
        input_df = input_dataset.to_table(
            columns=smiles_columns
        ).to_pandas()
    # else read in pandas csv
    elif input_path.is_file():
        input_df = pd.read_csv(input_path)

    # get smiles to label
    all_smiles = set()
    for col in smiles_columns:
        all_smiles.update(set(input_df[col].dropna().unique()))

    logger.info(f"Found {len(all_smiles)} unique smiles to label")

    all_rows = []
    for smiles in tqdm.tqdm(sorted(all_smiles), desc="Labeling SMILES"):
        smiles = sanitize_smiles(smiles)
        mol = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
        labels = forcefield.label_molecules(mol.to_topology())[0]

        PARAMETER_TYPES = ["Bonds", "Angles", "ProperTorsions", "ImproperTorsions", "vdW"]
        for parameter_type in PARAMETER_TYPES:
            parameter_values = labels[parameter_type]
            parameter_ids = sorted(set([p.id for p in parameter_values.values()]))
            for parameter_id in parameter_ids:
                row = {
                    "smiles": smiles,
                    "forcefield": forcefield_name,
                    "parameter_type": parameter_type,
                    "parameter_id": parameter_id,
                }
                all_rows.append(row)

    table = pa.Table.from_pylist(all_rows)

    output_path = pathlib.Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    pq.write_table(table, output_path)
    logger.info(f"Wrote {len(all_rows)} rows to {output_path}")


if __name__ == "__main__":
    main()
