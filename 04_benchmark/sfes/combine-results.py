"""
Combine all CSVs of SFEs in an input directory into one DataFrame with reference values.

This script reads in all CSVs of an input directory. Reference values are labelled with the `reference_name`.
Combined data is saved to the specified output file.
"""

import pathlib
import sys

import click

import pandas as pd
from loguru import logger

logger.remove()
logger.add(sys.stdout)


@click.command(help=__doc__)
@click.option(
    "--input-directory",
    "-i",
    type=click.Path(exists=True, dir_okay=True, readable=True),
    help="Path to the input directory"
)
@click.option(
    "--output-file",
    "-o",
    type=click.Path(exists=False, dir_okay=False, writable=True),
    help="Path to the output CSV file"
)
@click.option(
    "--reference-name",
    "-r",
    type=str,
    default="Experimental",
    help="Name of the reference method"
)
def main(
    input_directory: str,
    output_file: str,
    reference_name: str = "Experimental",
):
    input_directory = pathlib.Path(input_directory)
    csv_files = sorted(input_directory.glob("*.csv"))
    logger.info(f"Found {len(csv_files)} CSV files in {input_directory}")

    dfs = []
    for csv_file in csv_files:
        df_ = pd.read_csv(csv_file)
        dfs.append(df_)

    combined_df = pd.concat(dfs, ignore_index=True)

    # get experimental
    reference = combined_df[combined_df["Method"] == reference_name]
    logger.info(f"Found {len(reference)} rows for reference method '{reference_name}'")

    results_df = combined_df[combined_df["Method"] != reference_name]
    # get reference results for each row
    rows = []
    for _, row in results_df.iterrows():
        solute = row["Solute"]
        solvent = row["Solvent"]
        ref_row = reference[(reference["Solute"] == solute) & (reference["Solvent"] == solvent)]
        assert len(ref_row) == 1, f"Expected one reference row for {solute} in {solvent}, found {len(ref_row)}"
        ref_value = ref_row["Value (kcal / mol)"].values[0]
        ref_error = ref_row["Uncertainty (kcal / mol)"].values[0]
        row = dict(row)
        row["Reference Value (kcal / mol)"] = ref_value
        row["Reference Uncertainty (kcal / mol)"] = ref_error
        row["Id"] = ref_row["Id"].values[0]
        rows.append(row)


    final_df = pd.DataFrame(rows)
    unique_methods = sorted(final_df["Method"].unique())
    n_unique_methods = len(unique_methods)
    logger.info(f"Combined results for {n_unique_methods} unique methods: {', '.join(unique_methods)}")

    output_file = pathlib.Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    final_df.to_csv(output_file, index=False)
    logger.info(f"Combined results (n={len(final_df)}) saved to {output_file}")


if __name__ == "__main__":
    main()
