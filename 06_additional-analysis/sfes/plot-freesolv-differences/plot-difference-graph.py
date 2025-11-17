"""
This script generates difference graphs comparing two methods from a CSV file.
It groups the data by a specified column (e.g., solvent) and creates separate plots for each group.
The plots are saved in a specified output directory.
\b
The CSV file is expected to have the following columns:
    - Id: Unique identifier for each measurement
    - Method: Name of the method used
    - Value (kcal / mol): Measured value by the method
    - Reference Value (kcal / mol): Reference value for comparison
    - Solvent: Name of the solvent used
    - Solute: Name of the solute used
"""

import pathlib
import sys

import click
import tqdm

from loguru import logger
import pandas as pd
from matplotlib import pyplot as plt


logger.remove()
logger.add(sys.stdout)


def plot(
    target_df: pd.DataFrame,
    reference_col: str,
    target_col: str,
    plot_by: str,
    other_col: str,
    width: float = 8,
    markersize: int = 50,
    xlim: float = 5,
):
    """Plot a difference graph comparing two methods."""
    df_plot = target_df.sort_values([plot_by, other_col], kind="stable")

    # Build y-positions with gaps between solvent segments
    step = 0.5  # vertical step between items
    ypos = {}
    y = 0.0

    for plot_mol, g in df_plot.groupby(plot_by, sort=False):
        for other_mol in g[other_col]:
            ypos[(plot_mol, other_mol)] = y
            y += step

    # Create the plot
    n_items = len(df_plot)
    height = max(2, 0.1 * n_items + 0.3 * len(ypos))  # heuristic sizing
    fig, ax = plt.subplots(figsize=(width, height))

    for _, row in df_plot.iterrows():
        key = (row[plot_by], row[other_col])
        y_i = ypos[key]

        # endpoints: circle at new method, diamond at reference
        ax.scatter(
            [row[target_col]],
            [y_i],
            s=markersize,
            marker="o",
            color="orange",
            zorder=10,
        )
        ax.scatter(
            [row[reference_col]],
            [y_i],
            s=markersize,
            marker="D",
            color="#015480",
            zorder=10,
        )

        # line
        ax.plot(
            [row[target_col], row[reference_col]],
            [y_i, y_i],
            linewidth=1.5,
            color="#3E424A",
            zorder=1,
        )

    # Y ticks: show y-axis names in plotted order
    yticks = []
    yticklabels = []
    for _, row in df_plot.iterrows():
        key = (row[plot_by], row[other_col])
        yticks.append(ypos[key])
        yticklabels.append(row[other_col])

    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels)
    ax.set_ylim(min(yticks) - 1, max(yticks) + 1)
    # hard-code this
    ax.set_xlim(-xlim, xlim)

    # Axis labels & title
    ax.set_xlabel("Computed - Reference SFE (kcal / mol)")
    ax.set_ylabel(other_col)

    # Add a faint vertical reference line at zero (optional but handy)
    ax.axvline(0, linewidth=1, linestyle="--", alpha=1, color="gray", zorder=0)
    for i in range(1, 5):
        ax.axvline(i, linewidth=1, linestyle="--", alpha=0.3, color="gray", zorder=0)
        ax.axvline(-i, linewidth=1, linestyle="--", alpha=0.3, color="gray", zorder=0)

    # Legend explaining markers
    legend_handles = [
        plt.Line2D(
            [0], [0], marker="o", linestyle="None", color="orange", label=target_col
        ),
        plt.Line2D(
            [0], [0], marker="D", linestyle="None", color="#015480", label=reference_col
        ),
    ]
    ax.legend(handles=legend_handles, title="Endpoints", frameon=False, loc="best")
    return ax


@click.command(help=__doc__)
@click.option(
    "--input-file",
    "-i",
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Input CSV file with MNsol results.",
    default="../../../04_benchmark/sfes/output/mnsol.csv",
)
@click.option(
    "--output-directory",
    "-o",
    type=click.Path(file_okay=False, dir_okay=True, writable=True),
    help="Output directory to save plots.",
    default="images/mnsol-by-solvent",
)
@click.option(
    "--reference-method",
    "-r",
    type=str,
    help="Reference method name (as in CSV).",
    default="Sage 2.2.1 + ELF10",
)
@click.option(
    "--target-method",
    "-t",
    type=str,
    help="Target method name (as in CSV).",
    default="v3-k100",
)
@click.option(
    "--plot-by",
    "-cp",
    type=click.Choice(["Solvent", "Solute"]),
    help="Column name to group molecules by (e.g., 'Solvent').",
    default="Solvent",
)
@click.option(
    "--other-col",
    "-co",
    type=click.Choice(["Solvent", "Solute"]),
    help="Column name for the other variable plotted on y-axis (e.g., 'Solute').",
    default="Solute",
)
@click.option(
    "--width",
    "-w",
    type=float,
    help="Width of the plot in inches.",
    default=6,
)
@click.option(
    "--markersize",
    "-m",
    type=int,
    help="Size of the markers in the plot.",
    default=50,
)
@click.option(
    "--xlim",
    "-x",
    type=float,
    help="X-axis limit (symmetric around zero).",
    default=5,
)
def main(
    input_file: str = "../../../04_benchmark/sfes/output/mnsol.csv",
    output_directory: str = "images/mnsol-by-solvent",
    reference_method: str = "Sage 2.2.1 + ELF10",
    target_method: str = "v3-k100",
    plot_by: str = "Solvent",
    other_col: str = "Solute",
    width: float = 6,
    markersize: int = 50,
    xlim: float = 5,
):
    df = pd.read_csv(input_file)
    df["Difference (kcal / mol)"] = (
        df["Value (kcal / mol)"] - df["Reference Value (kcal / mol)"]
    )

    rows_reference = {
        row["Id"]: row for _, row in df[df.Method == reference_method].iterrows()
    }

    # Filter to target method and add reference differences
    target_df = pd.DataFrame(df[df.Method == target_method])
    difference_values = []
    for _, row in target_df.iterrows():
        reference_row = rows_reference[row["Id"]]
        difference_values.append(reference_row["Difference (kcal / mol)"])

    target_df[f"{target_method}"] = target_df["Difference (kcal / mol)"]
    target_df[f"{reference_method}"] = difference_values

    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    logger.info(f"Saving plots to {output_directory}")
    i = 0
    for molecule, subdf in tqdm.tqdm(target_df.groupby(plot_by)):
        logger.info(f"Plotting {plot_by}: {molecule} with {len(subdf)} items")
        ax = plot(
            subdf,
            reference_col=reference_method,
            target_col=target_method,
            plot_by=plot_by,
            other_col=other_col,
            width=width,
            markersize=markersize,
            xlim=xlim,
        )
        ax.set_title(f"{plot_by}: {molecule}")
        plt.tight_layout()
        # can't just use molecule smiles in filename, case insensitive
        output_file = output_directory / f"{plot_by}-{molecule}-{i}.png"
        plt.savefig(output_file, dpi=300)
        logger.info(f"Saved plot to {output_file}")
        plt.close()
        i += 1


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
