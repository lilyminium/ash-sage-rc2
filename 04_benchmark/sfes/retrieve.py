import pathlib
import sys
import os
import warnings
from loguru import logger

from alchemiscale import AlchemiscaleClient, Scope

import click
import tqdm
import pandas as pd
from rdkit import Chem

from openff.units import unit
from openfe import Transformation

logger.remove()
logger.add(sys.stdout)

def sanitize_smiles(smi: str) -> str:
    """Sanitize SMILES"""
    return Chem.MolToSmiles(Chem.MolFromSmiles(smi))


def get_alchemiscale_client() -> AlchemiscaleClient:
    """Get client from environmental variables"""
    user_id = os.environ["ALCHEMISCALE_ID"]
    user_key = os.environ["ALCHEMISCALE_KEY"]
    asc = AlchemiscaleClient("https://api.alchemiscale.org", user_id, user_key)
    return asc


def get_results_with_transformations(
    asc: AlchemiscaleClient,
    scope_key: str,
) -> dict[Transformation, dict]:
    """
    Retrieve results with transformations from Alchemiscale client

    Parameters
    ----------
    asc: AlchemiscaleClient
        Alchemiscale client
    scope_key: str
        Scope key
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        results = asc.get_network_results(scope_key)
        results_with_transformation = {asc.get_transformation(k): v for k, v in results.items()}
    return results_with_transformation


@click.command()
@click.option(
    "--forcefield-name",
    "-ff",
    type=str,
    help="Force field name to use"
)
@click.option(
    "--scope-key-file",
    "-i",
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Path to the input scope key file"
)
@click.option(
    "--dataset-name",
    "-d",
    type=str,
    help="Name of the dataset"
)
@click.option(
    "--output-file",
    "-o",
    type=click.Path(exists=False, dir_okay=False, writable=True),
    help="Path to the output CSV file"
)
def main(
    forcefield_name: str,
    scope_key_file: str,
    dataset_name: str,
    output_file: str,
):
    with open(scope_key_file, "r") as f:
        scope_key = f.read().strip()

    logger.info(f"Loaded scope key from {scope_key_file}: {scope_key}")

    asc = get_alchemiscale_client()
    results_by_transformation = get_results_with_transformations(asc, scope_key)
    rows = []

    for transformation, result in tqdm.tqdm(results_by_transformation.items()):
        # check everything is as it should be
        vacuum_ffs = set(
            transformation.protocol.settings.vacuum_forcefield_settings.forcefields
        )
        solvent_ffs = set(
            transformation.protocol.settings.solvent_forcefield_settings.forcefields
        )
        assert vacuum_ffs.issubset(solvent_ffs)
        assert len(solvent_ffs) <= 2

        # save data in same format
        solvent = sanitize_smiles(transformation.stateA.components["solvent"].smiles)
        solute = sanitize_smiles(transformation.stateA.components["solute"].smiles)
        value = result.get_estimate().m_as(unit.kilocalories_per_mole)
        error = result.get_uncertainty().m_as(unit.kilocalories_per_mole)

        row = {
            "Id": transformation.name,
            "Temperature (K)": 298.15,
            "Pressure (kPa)": 101.325,
            "Solute": solute,
            "Solvent": solvent,
            "Value (kcal / mol)": value,
            "Uncertainty (kcal / mol)": error,
            "Method": forcefield_name,
            "Dataset": dataset_name
        }
        rows.append(row)

    df = pd.DataFrame(rows)

    output_file = pathlib.Path(output_file)
    output_file.parent.mkdir(exist_ok=True, parents=True)
    df.to_csv(output_file, index=False)
    logger.info(f"Saved {len(df)} entries to {output_file}")


if __name__ == "__main__":
    main()
