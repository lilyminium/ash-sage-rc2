"""
Plot aggregated statistics with error bars.

This returns a FacetGrid plot with two panels:
* Error (RMSE, MUE, MSE)
* Correlation (rho, R2)
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

PANEL_MAP = {
    "RMSE": "Error",
    "MUE": "Error",
    "MSE": "Error",
    "rho": "Correlation",
    "R2": "Correlation",
}


def draw_errorbars(
    data,
    x: str = "xpos",
    y: str = "mle",
    xerr_low: str | None = None,
    xerr_high: str | None = None,
    yerr_low: str | None = "low",
    yerr_high: str | None = "high",
    fmt: str = "o",
    capsize: int = 3,
    **kwargs,
):
    """Draw error bars for each row in the dataframe"""
    ax = plt.gca()

    xerr, yerr = None, None
    if xerr_low is not None or xerr_high is not None:
        if xerr_low is None or xerr_high is None:
            raise ValueError("Both xerr_low and xerr_high must be provided")

    if yerr_low is not None or yerr_high is not None:
        if yerr_low is None or yerr_high is None:
            raise ValueError("Both yerr_low and yerr_high must be provided")

    for _, row in data.iterrows():
        if xerr_low is not None:
            xerr = np.array(
                [[max(row[x] - row[xerr_low], 0)], [max(row[xerr_high] - row[x], 0)]]
            )
        if yerr_low is not None:
            yerr = np.array(
                [[max(row[y] - row[yerr_low], 0)], [max(row[yerr_high] - row[y], 0)]]
            )
        ax.errorbar(
            row[x], row[y], xerr=xerr, yerr=yerr, fmt=fmt, capsize=capsize, **kwargs
        )


def plot_facetgrid(
    df: pd.DataFrame,
    force_field_order: list[str],
    force_field_col: str = "FF",
    dodge_spacing: float = 0.1,
    height: float = 4,
    aspect: float = 1,
    unit_str: str = "[kJ/mol]",
    title: str | None = None,
):
    df = df.copy()

    # Assign panels and preserve order
    df["Panel"] = df["stat"].map(PANEL_MAP)
    STAT_ORDER = list(PANEL_MAP)
    df["stat"] = pd.Categorical(df["stat"], categories=STAT_ORDER, ordered=True)

    # Assign numeric x positions for dodge
    x_lookup = {cat: i for i, cat in enumerate(STAT_ORDER)}
    df["xbase"] = df["stat"].map(x_lookup).astype(float)

    # Compute dodge offsets by FF
    offsets = {
        ff: (i - (len(force_field_order) - 1) / 2) * dodge_spacing
        for i, ff in enumerate(force_field_order)
    }
    df["xpos"] = df["xbase"] + df[force_field_col].map(offsets).astype(float)

    # plot facetgrid
    g = sns.FacetGrid(
        df,
        col="Panel",
        col_order=["Error", "Correlation"],
        hue=force_field_col,
        height=height,
        aspect=aspect,
        sharey=False,
        sharex=False,
        legend_out=True,
        hue_order=force_field_order,
    )
    g.map_dataframe(draw_errorbars)

    for panel, ax in g.axes_dict.items():
        ax.set_xticks([x_lookup[s] for s in STAT_ORDER if PANEL_MAP[s] == panel])
        ax.set_xticklabels([s for s in STAT_ORDER if PANEL_MAP[s] == panel])
        ax.set_xlabel("")
        if panel == "Error":
            ax.set_ylabel(f"{panel} {unit_str}")
            ax.axhline(0, ls="--", color="gray", lw=1)
        else:
            ax.set_ylabel(panel)

    g.set_titles("")
    if title is not None:
        g.figure.suptitle(title)
    plt.tight_layout()

    g.add_legend(title=force_field_col)
    return g


@click.command(help=__doc__)
@click.option(
    "--input-file",
    "-i",
    "input_files",
    multiple=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    help="Path to the input/s CSV file/s.",
    required=False,
)
@click.option(
    "--input-file-directory",
    "-id",
    type=click.Path(exists=True, dir_okay=True, file_okay=False),
    default=None,
    required=False,
    help="Path to the directory containing input CSV files.",
)
@click.option(
    "--image-directory",
    "-im",
    type=click.Path(exists=False, dir_okay=True, file_okay=False),
    default="images/aggregated-statistics",
    help="Path to the output directory.",
)
@click.option(
    "--force-field-order",
    "-ff",
    multiple=True,
    type=str,
    required=True,
    help="Order of force fields to include in the plot.",
)
@click.option(
    "--force-field-col",
    "-cff",
    type=str,
    default="FF",
    help="Column name for force field in the input CSV file.",
)
@click.option(
    "--dodge-spacing",
    "-s",
    type=float,
    default=0.1,
    help="Spacing between dodged points.",
)
@click.option(
    "--height",
    "-h",
    type=float,
    default=4,
    help="Height of each facet.",
)
@click.option(
    "--aspect",
    "-a",
    type=float,
    default=1,
    help="Aspect ratio of each facet.",
)
@click.option(
    "--unit-str",
    "-u",
    type=str,
    default="[kJ/mol]",
    help="Unit string to display on the y-axis.",
)
def main(
    input_files: list[str],
    image_directory: str,
    force_field_order: list[str],
    input_file_directory: str | None = None,
    force_field_col: str = "FF",
    dodge_spacing: float = 0.1,
    height: float = 4,
    aspect: float = 1,
    unit_str: str = "[kJ/mol]",
):
    dfs = []
    logger.info(f"{len(input_files)} input files specified")
    if input_file_directory is not None:
        input_file_directory = pathlib.Path(input_file_directory)
        additional_files = list(input_file_directory.glob("*.csv"))
        logger.info(f"Adding {len(additional_files)} files from {input_file_directory}")
        input_files = list(input_files) + [str(f) for f in additional_files]


    for input_file in input_files:
        df = pd.read_csv(input_file)
        logger.info(f"Read {len(df)} rows from {input_file}")
        dfs.append(df)
    df = pd.concat(dfs, ignore_index=True)

    df = df[df[force_field_col].isin(force_field_order)]
    logger.info(f"Filtered to {len(df)} rows with given force fields")

    image_directory = pathlib.Path(image_directory)
    image_directory.mkdir(parents=True, exist_ok=True)

    # plot for each group
    for group, subdf in tqdm.tqdm(df.groupby("group")):
        n_samples = subdf["n"].values[0]
        title = f"{group} (n={n_samples})"
        g = plot_facetgrid(
            subdf,
            force_field_order=force_field_order,
            force_field_col=force_field_col,
            dodge_spacing=dodge_spacing,
            height=height,
            aspect=aspect,
            unit_str=unit_str,
            title=title,
        )
        output_file = image_directory / f"stats-{group}.png"
        plt.savefig(output_file, dpi=300)
        logger.info(f"Saved plot to {output_file}")
        plt.close()


if __name__ == "__main__":
    main()
