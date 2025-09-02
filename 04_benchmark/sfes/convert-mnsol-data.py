"""
Convert MNSol data to a reference dataframe.
"""

import json
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
    default="input-data/MNSol_alldata.txt",
    help="Path to the input database file"
)
@click.option(
    "--name-to-smiles-file",
    "-n",
    type=str,
    default="input-data/mnsol-name-to-smiles.json",
    help="Path to the name to SMILES mapping file"
)
@click.option(
    "--output-database-file",
    "-o",
    type=str,
    default="results/mnsol/experiment.csv",
    help="Path to the output database file"
)
def main(
    input_database_file: str = "input-data/MNSol_alldata.txt",
    name_to_smiles_file: str = "input-data/mnsol-name-to-smiles.json",
    output_database_file: str = "results/mnsol/experiment.csv"
):

    with open(name_to_smiles_file, "r") as f:
        name_to_smiles = json.load(f)

    df = pd.read_csv(input_database_file, sep="\t")
    logger.info(f"Read {len(df)} entries from {input_database_file}")
    rows = []
    counter = 1
    for _, row in tqdm.tqdm(df.iterrows(), desc="Converting rows"):
        if row.type != "abs":
            continue
        try:
            solute = sanitize_smiles(name_to_smiles[row.SoluteName])
            solvent = sanitize_smiles(name_to_smiles[row.Solvent])
        except KeyError:
            continue
    
        experimental_dg = row.DeltaGsolv
        uncertainty = 0.3
        row = {
            "Id": f"mnsol-{counter:04d}",
            "Temperature (K)": 298,
            "Pressure (kPa)": 100,
            "Solute": solute,
            "Solvent": solvent,
            "Value (kcal / mol)": experimental_dg,
            "Uncertainty (kcal / mol)": uncertainty,
            "Method": "Experimental",
            "Dataset": "MNSol",
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
