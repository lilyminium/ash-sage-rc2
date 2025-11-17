"""
This script reads in a statistics dataframe and compares a reference
and target force field. It ranks groups by size of meaningful difference
for each statistic.

This writes out a file of significant changes, with columns:
    - force field column (usually 'FF') (str): force field
    - group (str): group name
    - mle_diff (float): difference in maximum likelihood estimate (target - reference)
    - Change (str): 'Improved' or 'Worsened' depending on whether the target force field
      significantly improved or worsened over the reference force field.
"""

from collections import defaultdict
import pathlib
import sys
import typing

import click


from loguru import logger

import numpy as np
import pandas as pd
import pyarrow.compute as pc
import pyarrow.dataset as ds

logger.remove()
logger.add(sys.stdout)


@click.command(help=__doc__)
@click.option(
    "--input",
    "-i",
    "input_files",
    multiple=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help=(
        "Input files. These should be CSV files containing statistics data."
    )
)
@click.option(
    "--output",
    "-o",
    "output_file",
    type=click.Path(exists=False, dir_okay=False, file_okay=True, writable=True),
    default="output/mm-qm-differences-labelled.csv",
    help="Output CSV file.",
    show_default=True,
)
@click.option(
    "--reference-force-field",
    "-ref",
    type=str,
    default="Sage 2.2.1",
    help="Reference force field to compare against.",
)
@click.option(
    "--target-force-field",
    "-tgt",
    type=str,
    default="v1-k100",
    help="Target force field to compare.",
)
@click.option(
    "--force-field-col",
    "-cff",
    type=str,
    default="FF",
    help="Column name for force field in the output CSV file.",
)
@click.option(
    "--statistic",
    "-stat",
    type=click.Choice(["RMSE", "MUE", "MSE"], case_sensitive=False),
    default="RMSE",
    help="Statistic to compare.",
)
def main(
    input_files: str,
    output_file: str,
    reference_force_field: str = "Sage 2.2.1",
    target_force_field: str = "v1-k100",
    force_field_col: str = "FF",
    statistic: typing.Literal["RMSE", "MUE", "MSE"] = "RMSE",
):
    dfs = []
    for input_file in input_files:
        df = pd.read_csv(input_file)
        dfs.append(df)
    df = pd.concat(dfs, ignore_index=True)
    logger.info(f"Loaded to {len(df)} rows")

    df = df[df[force_field_col].isin([reference_force_field, target_force_field])]
    logger.info(f"Filtered to {len(df)} rows with given force fields")

    df = df[df.stat == statistic]
    logger.info(f"Filtered to {len(df)} rows with statistic {statistic}")

    logger.info(f"Comparing {target_force_field} to {reference_force_field}")

    reference_df = df[df[force_field_col] == reference_force_field]
    target_df = df[df[force_field_col] == target_force_field]
    assert len(reference_df) == len(target_df)

    # join on group
    merged_df = pd.merge(
        reference_df,
        target_df,
        on="group",
        suffixes=("_ref", "_tgt"),
    )
    logger.info(f"Merged to {len(merged_df)} rows")
    assert len(merged_df) == len(reference_df)

    merged_df["mle_diff"] = merged_df["mle_tgt"] - merged_df["mle_ref"]
    merged_df = merged_df.sort_values(by="mle_diff", ascending=True)

    # compute cases where target is better than reference (assuming lower is better)
    sig_improved = merged_df[merged_df["high_tgt"] < merged_df["low_ref"]]
    sig_improved["Change"] = "Improved"
    sig_worsened = merged_df[merged_df["low_tgt"] > merged_df["high_ref"]]
    sig_worsened["Change"] = "Worsened"
    changes = pd.concat([sig_improved, sig_worsened], ignore_index=True).sort_values(
        by="mle_diff"
    )

    logger.info(
        f"Found {len(changes)} significant changes ({len(sig_improved)} improved, {len(sig_worsened)} worsened)"
    )
    logger.info(f"Improved: {', '.join(sig_improved['group'].tolist())}")
    logger.info(f"Worsened: {', '.join(sig_worsened['group'].tolist())}")

    output_file = pathlib.Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    changes.to_csv(output_file, index=False)
    logger.info(f"Wrote changes to {output_file}")

if __name__ == "__main__":
    main()
