"""
Download ThermoML dataset initially.

This saves the file "thermoml.csv".
This file contains the initial ThermoML data downloaded from the NIST website.
It filters out any properties with multi-component SMILES or malformed data.
It also hashes the physical properties for unique identification.
It is used as the starting point for further processing and filtering in subsequent scripts.
"""

import pathlib
import logging
import time

import click
import tqdm
import pandas as pd

from openff.evaluator.datasets.datasets import PhysicalPropertyDataSet
from openff.evaluator.plugins import register_default_plugins
from openff.evaluator.datasets.curation.components import thermoml

register_default_plugins()

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


def hash_physical_property(physprop) -> str:
    state = physprop.thermodynamic_state
    source = physprop.source
    attrs = (
        str(physprop.substance),
        physprop.phase.value,
        (
            str(state.temperature),
            str(state.pressure),
        ),
        (
            str(physprop.value),
            str(physprop.uncertainty),
        ),
        (
            source.doi,
            source.reference,
        )
    )
    return str(hash(attrs))

@click.command
@click.option(
    "--output-directory",
    "-o",
    default="input",
    help="The directory to save the output CSV files",
)
@click.option(
    "--n-processes",
    "-np",
    default=1,
    help="The number of processes to use for loading the data",
)
def main(
    output_directory: str = "output",
    n_processes: int = 1,
):
    now = time.time()
    logger.info(f"Starting at {time.ctime(now)}")

    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(exist_ok=True, parents=True)

    df = thermoml.ImportThermoMLData.apply(
        pd.DataFrame(),
        thermoml.ImportThermoMLDataSchema(
            cache_file_name="input/initial-thermoml.csv"
        ),
        n_processes
    )

    logger.info(f"Downloaded {len(df)} initial data points")

    # filter for mis-formed/multi-component SMILES
    filtered_properties = []

    smiles_columns = [x for x in df.columns if x.startswith("Component ")]
    for _, row in tqdm.tqdm(df.iterrows(), desc="Filtering", total=len(df)):
        smiles = [row[col] for col in smiles_columns if pd.notna(row[col])]
        # filter for multi-component smiles -- not caught in first step
        if any("." in smi for smi in smiles):
            continue

        filtered_properties.append(row)

    df = pd.DataFrame(filtered_properties)

    logger.info(f"Starting to convert to Evaluator dataset at {time.ctime(time.time())}")
    ds = PhysicalPropertyDataSet.from_pandas(df)
    logger.info(f"Finished converting to Evaluator dataset at {time.ctime(time.time())}")
    for physprop in tqdm.tqdm(ds.properties, desc="Hashing"):
        physprop.id = hash_physical_property(physprop)

    df = ds.to_pandas()
    
    output_path = output_directory / "thermoml.csv"
    df.to_csv(output_path)
    logger.info(
        f"Loaded {len(df)} data points "
        f"and saved to {output_path}"
    )

    logger.info(f"Finished at {time.ctime(time.time())}")
    logger.info(f"Elapsed time: {time.time() - now} seconds")


if __name__ == "__main__":
    main()
