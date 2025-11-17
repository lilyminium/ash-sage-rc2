"""
Gather separately-calculated benchmarks into a single CSV with reference values.
\b
This script writes four files:
    - `all-benchmarks.csv`
    - `summary-benchmarks.csv`
    - `summary-benchmarks-Density.csv`
    - `summary-benchmarks-EnthalpyOfMixing.csv`

\b
The all-benchmarks.csv file contains the columns:
    - id (str): the physical property ID
    - index (int): the index of the property in the dataset
    - type (str): the type of the physical property
    - substance (str): the substance the property was measured for
    - forcefield (str): the name of the force field used
    - replicate (int): the replicate number
    - value (float): the calculated value of the property in default units (g/mL for densities, kJ/mol for enthalpies of mixing)
    - uncertainty (float): the uncertainty of the calculated value in default units (g/mL for densities, kJ/mol for enthalpies of mixing)
    - reference_value (float): the reference value of the property in default units (g/mL for densities, kJ/mol for enthalpies of mixing)
    - reference_uncertainty (float): the uncertainty of the reference value in default units (g/mL for densities, kJ/mol for enthalpies of mixing)

\b
The summary-benchmarks.csv file contains the columns:
    - id (str): the physical property ID
    - type (str): the type of the physical property
    - substance (str): the substance the property was measured for
    - forcefield (str): the name of the force field used
    - mean (float): the mean value of the property in default units (g/mL for densities, kJ/mol for enthalpies of mixing)
    - uncertainty (float): the standard deviation of the value in default units (g/mL for densities, kJ/mol for enthalpies of mixing)
    - reference_value (float): the reference value of the property in default units (g/mL for densities, kJ/mol for enthalpies of mixing)
    - reference_uncertainty (float): the uncertainty of the reference value in default units (g/mL for densities, kJ/mol for enthalpies of mixing)
    - n_replicates (int): the number of replicates used to calculate the mean and uncertainty
"""

import pathlib
import sys

import click

import pandas as pd
from loguru import logger

logger.remove()
logger.add(sys.stdout)


@click.command()
@click.option(
    "--input-directory",
    "-i",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    default="training",
    help="Directory containing the input CSV files to combine.",
)
@click.option(
    "--output-directory",
    "-o",
    type=click.Path(file_okay=False, dir_okay=True),
    default="output",
    help="Directory to write the output CSV file to.",
)
def main(
    input_directory: str = "results/training",
    output_directory: str = "output/training",
):
    input_directory = pathlib.Path(input_directory)
    csv_files = sorted(input_directory.glob("all/*.csv"))

    dfs = []
    for csv_file in csv_files:
        df = pd.read_csv(csv_file)
        dfs.append(df)

    df = pd.concat(dfs, ignore_index=True)
    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    csv_file = output_directory / "all-benchmarks.csv"
    df.to_csv(csv_file)
    logger.info(f"Wrote {len(df)} entries to {csv_file}")

    csv_files = sorted(input_directory.glob("summary/*.csv"))
    dfs = []
    for csv_file in csv_files:
        df = pd.read_csv(csv_file)
        dfs.append(df)
    mean_values = pd.concat(dfs, ignore_index=True)
    output_file = output_directory / "summary-benchmarks.csv"
    mean_values.to_csv(output_file)
    logger.info(f"Wrote {len(mean_values)} summary entries to {output_file}")

    # break up into density and enthalpies of mixing
    for type_, subdf in mean_values.groupby("type"):
        output_file = output_directory / f"summary-benchmarks-{type_}.csv"
        subdf.to_csv(output_file)
        logger.info(f"Wrote {len(subdf)} summary entries to {output_file}")

if __name__ == "__main__":
    main()
