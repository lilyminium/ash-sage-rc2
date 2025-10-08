"""
Plot a distribution of changes in target property.
"""

import click
import pathlib
import sys

import tqdm
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from loguru import logger

logger.remove()
logger.add(sys.stdout)


@click.command(help=__doc__)
@click.option(
    "--input",
    "-i",
    "input_files",
    multiple=True,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    required=True,
    help="Input CSV files to process.",
)
@click.option(
    "--output",
    "-o",
    "output_file",
    default="output/change-in-target.csv",
    type=click.Path(dir_okay=False, file_okay=True, writable=True),
    help="Path to save the output CSV file with computed changes.",
)
@click.option(
    "--image-file",
    "-im",
    default="images/change-in-target.png",
    type=click.Path(dir_okay=False, file_okay=True, writable=True),
    help="Path to save the output image file (e.g., PNG).",
)
@click.option(
    "--target-forcefield",
    "-tgt",
    type=str,
    default="v3-k100",
    help="Name of the target force field to compare against the reference.",
)
@click.option(
    "--reference-forcefield",
    "-ref",
    type=str,
    default="Sage 2.2.1",
    help="Name of the reference force field.",
)
@click.option(
    "--forcefield-col",
    "-cff",
    type=str,
    default="FF",
    help="Name of the column containing force field names.",
)
@click.option(
    "--value-col",
    "-cv",
    type=str,
    default="value",
    help="Name of the column containing property values.",
)
@click.option(
    "--id-col",
    "-id",
    type=str,
    default="id",
    help="Name of the column containing unique identifiers for each data point.",
)
@click.option(
    "--smiles-col",
    "-cs",
    "smiles_cols",
    type=str,
    multiple=True,
    help="Columns containing SMILES strings for searching.",
)
@click.option(
    "--property-name",
    "-pn",
    type=str,
    default="Density",
    help="Name of the property to analyze.",
)
@click.option(
    "--property-units",
    "-pu",
    type=str,
    default="g/mL",
    help="Units of the property to analyze.",
)
def main(
    input_files: list[str],
    output_file: str,
    image_file: str,
    target_forcefield: str = "v3-k100",
    reference_forcefield: str = "Sage 2.2.1",
    forcefield_col: str = "FF",
    value_col: str = "value",
    id_col: str = "id",
    smiles_cols: list[str] = [],
    property_name: str = "Density",
    property_units: str = "g/mL",
):
    """
    Plot a distribution of changes in target property.

    Parameters
    ----------
    input_files : list[str]
        List of input CSV files containing the data.
    output_file : str
        Path to save the output CSV file with computed changes.
    image_file : str
        Path to save the output image file (e.g., PNG).
    target_forcefield : str, optional
        Name of the target force field to compare against the reference, by default "v3-k100".
    reference_forcefield : str, optional
        Name of the reference force field, by default "Sage 2.2.1".
    forcefield_col : str, optional
        Name of the column containing force field names, by default "FF".
    value_col : str, optional
        Name of the column containing property values, by default "value".
    id_col : str, optional
        Name of the column containing unique identifiers for each data point, by default "id".
    smiles_cols : list[str], optional
        List of columns containing SMILES strings for searching, by default [].
    """
    # Load and concatenate all input files
    dataframes = []
    for input_file in input_files:
        df = pd.read_csv(input_file)
        logger.info(f"Loaded {len(df)} rows from {input_file}")
        dataframes.append(df)
    data = pd.concat(dataframes, ignore_index=True)
    logger.info(f"Concatenated dataframe has {len(data)} rows")

    # Filter data for target and reference force fields
    target_data = data[data[forcefield_col] == target_forcefield]
    reference_data = data[data[forcefield_col] == reference_forcefield]

    logger.info(
        f"Target force field '{target_forcefield}' has {len(target_data)} rows"
    )
    logger.info(
        f"Reference force field '{reference_forcefield}' has {len(reference_data)} rows"
    )

    # pivot
    df = df[df[forcefield_col].isin([target_forcefield, reference_forcefield])]

    smiles_cols = list(smiles_cols)
    wide = df.pivot(
        index=[id_col] + smiles_cols,
        columns=forcefield_col,
        values=value_col,
    ).reset_index()
    # drop rows with missing data
    wide = wide.dropna(subset=[target_forcefield, reference_forcefield])
    logger.info(f"After pivoting, {len(wide)} rows have both target and reference data")

    # compute difference
    difference_col = f"{target_forcefield} - {reference_forcefield}"
    wide[difference_col] = wide[target_forcefield] - wide[reference_forcefield]


    # add column looking for water ("O") in any smiles column
    if smiles_cols:
        wide["water"] = wide[smiles_cols].apply(
            lambda row: any("O" in str(smiles) for smiles in row), axis=1
        )
        n_with_water = wide["water"].sum()
        logger.info(f"Found {n_with_water} rows with water in any of {smiles_cols}")

    # save data
    output_file = pathlib.Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    wide.to_csv(output_file, index=False)
    logger.info(f"Wrote data with {len(wide)} rows to {output_file}")

    # plot
    ax = sns.histplot(data=wide, x=difference_col, hue="water" if smiles_cols else None, kde=True)
    ax.axvline(0, ls="--", color="k")

    # label median difference
    median_diff = np.median(wide[difference_col])
    ax.axvline(median_diff, ls=":", color="r")
    # set top left
    ax.text(
        0.05,
        0.95,
        f"Median change: {median_diff:.3f} {property_units}",
        transform=ax.transAxes,
        verticalalignment="top",
        color="r",
    )

    ax.set_title(difference_col)
    ax.set_xlabel(f"Change in {property_name} ({property_units})")
    plt.tight_layout()

    image_file = pathlib.Path(image_file)
    image_file.parent.mkdir(parents=True, exist_ok=True)

    plt.savefig(image_file, dpi=300)
    logger.info(f"Saved image to {image_file}")


if __name__ == "__main__":
    main()
