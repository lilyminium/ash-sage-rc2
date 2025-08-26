"""
This script compares the performance of the first (0th) and last (15th) iterations of a force field
on training targets that are physical properties.
Note, the 15th iteration does not necessarily correspond to the final iteration of the force field.

Output: images/iter_00-15.png
"""
import sys

import click
import pandas as pd
import numpy as np
import seaborn as sns
from loguru import logger
from matplotlib import pyplot as plt

from utils import get_limits

logger.remove()
logger.add(sys.stdout)

@click.command()
@click.option(
    "--input-file",
    "-i",
    default="output/training-per-iteration.csv",
    type=click.Path(exists=True, dir_okay=False),
    help="Path to the input CSV file containing training data per iteration."
)
@click.option(
    "--output-file",
    "-o",
    default="images/iter_00-15.png",
    type=click.Path(exists=False, dir_okay=False),
    help="Path to save the output image comparing iterations."
)
def main(
    input_file: str = "output/training-per-iteration.csv",
    output_file: str = "images/iter_00-15.png"
):
    df = pd.read_csv(input_file, index_col=0)

    ff_cols = [col for col in df.columns if "iter" in col]
    df["Id"] = df.index
    id_cols = ["Reference", "Id", "Property type"]
    other_cols = [x for x in df.columns if x not in id_cols + ff_cols]
    melted = df.melt(
        id_vars=id_cols + other_cols,
        value_vars=ff_cols,
        value_name="Value",
        var_name="Force field",
    )

    row = "Property type"
    col = "Force field"
    y = "Value"

    g = sns.FacetGrid(
        data=melted[melted["Force field"].isin(["iter_0000", "iter_0015"])],
        col=col,
        row=row,
        sharex=False,
        sharey=False,
        height=3,
        aspect=1,
        margin_titles=True
    )
    g.map(sns.scatterplot, "Reference", y, s=3)
    g.set_titles(col_template="{col_name}", row_template="{row_name}")

    for (row_name, col_name), ax in g.axes_dict.items():
        subdf = melted[
            (melted[row] == row_name)
            & (melted[col] == col_name)
        ]
        obs = subdf[y].values
        ref = subdf["Reference"].values
        rmse = ((obs - ref) ** 2).mean() ** 0.5
        
        limits = get_limits(subdf, y)
        ax.set_xlim(limits)
        ax.set_ylim(limits)
        ax.text(0.1, 0.9, f"RMSE={rmse:.2e}", transform=ax.transAxes)
        ax.plot(limits, limits, ls="--", color="gray", lw=1)

        if row_name == "Density":
            ax.set_ylabel("Computed [g/mL]")
            ax.set_xlabel("Reference [g/mL]")
        else:
            ax.set_ylabel("Computed [kJ/mol]")
            ax.set_xlabel("Reference [kJ/mol]")

    plt.tight_layout()

    g.savefig(output_file, dpi=300)
    logger.info(f"Saved figure to {output_file}")



if __name__ == "__main__":
    main()
