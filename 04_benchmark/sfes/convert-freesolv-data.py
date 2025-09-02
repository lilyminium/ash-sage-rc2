"""
Convert FreeSolv data to a reference dataframe.
"""

import pathlib
import sys

import pandas as pd
import click
import tqdm
from rdkit import Chem
from openff.units import unit

from loguru import logger

logger.remove()
logger.add(sys.stdout)

def sanitize_smiles(smi: str) -> str:
    """Sanitize SMILES"""
    return Chem.MolToSmiles(Chem.MolFromSmiles(smi))


@click.command(help=__doc__)
@click.option(
    "--input-database-file",
    "-i",
    type=str,
    default="input-data/freesolv-database.txt",
    help="Path to the input database file"
)
@click.option(
    "--output-database-file",
    "-o",
    type=str,
    default="results/freesolv/experiment.csv",
    help="Path to the output database file"
)
def main(
    input_database_file: str = "input-data/freesolv-database.txt",
    output_database_file: str = "results/freesolv/experiment.csv"
):
    df = pd.read_csv(
        input_database_file, sep="; ",
        skiprows=2,
        engine="python",
    )
    logger.info(f"Read {len(df)} entries from {input_database_file}")
    rows = []
    counter = 1
    for _, row in tqdm.tqdm(df.iterrows(), desc="Converting rows"):
        try:
            solute = sanitize_smiles(row["SMILES"])
            solvent = sanitize_smiles("O")
        except KeyError:
            continue
    
        experimental_dg = row["experimental value (kcal/mol)"]

        uncertainty = row["experimental uncertainty (kcal/mol)"]

        row = {
            "Id": row["# compound id (and file prefix)"],
            "Temperature (K)": 298.15,
            "Pressure (kPa)": 101.325,
            "Solute": solute,
            "Solvent": solvent,
            "Value (kcal / mol)": experimental_dg,
            "Uncertainty (kcal / mol)": uncertainty,
            "Method": "Experimental",
            "Dataset": "FreeSolv"
        }
        rows.append(row)
        counter += 1

    final_df =  pd.DataFrame(rows)

    output_database_file = pathlib.Path(output_database_file)
    output_database_file.parent.mkdir(exist_ok=True, parents=True)
    final_df.to_csv(output_database_file, index=False)
    logger.info(f"Saved {len(final_df)} entries to {output_database_file}")


if __name__ == "__main__":
    main()
