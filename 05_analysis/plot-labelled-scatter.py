"""
Plot labelled scatter plots and calculate statistics.
\b
This script saves a figure for each group plotted,
with force fields separated per-column on a FacetGrid.
It also saves a summary CSV file with the computed statistics.
The CSV contains the columns:
    - force field column (usually 'FF') (str): force field
    - stat (str): statistic in question (RMSE, MUE, rho, or R2)
    - n (int): the number of observations
    - mle (float): the MLE value
    - mean (float): the mean value
    - low (float): the lower bound of the confidence interval
    - high (float): the upper bound of the confidence interval
    - group (str): the functional group
"""

import click
import pathlib
import sys

import tqdm
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from loguru import logger
from cinnabar.stats import bootstrap_statistic
from plot_helpers import get_limits

logger.remove()
logger.add(sys.stdout)


def plot_errorbar(
    data,
    x_col: str = None,
    y_col: str = None,
    xerr_col: str = None,
    yerr_col: str = None,
    color: str = None,
    markersize: int = 1,
    linewidth: float = 1,
):
    """Plot points with error bars"""
    ax = plt.gca()

    x = data[x_col]
    y = data[y_col]
    xerr = data[xerr_col] if xerr_col else None
    yerr = data[yerr_col] if yerr_col else None

    ax.errorbar(
        x,
        y,
        xerr=xerr,
        yerr=yerr,
        linestyle="",
        color=color,
        marker=".",
        markersize=markersize,
        lw=linewidth,
    )


def compute_stats(
    df: pd.DataFrame,
    forcefield_col: str = "FF",
    x_col: str = None,
    y_col: str = None,
    xerr_col: str = None,
    yerr_col: str = None,
) -> pd.DataFrame:
    stat_rows = []

    for ff, subdf in df.groupby(by=forcefield_col):
        y = subdf[y_col].values
        yerr = subdf[yerr_col].values
        x = subdf[x_col].values
        xerr = subdf[xerr_col].values

        for i, stat in enumerate(["RMSE", "MUE", "R2", "rho"], 1):
            s = bootstrap_statistic(
                x,
                y,
                xerr,
                yerr,
                statistic=stat,
            )
            stat_row = {forcefield_col: ff, "stat": stat, "n": len(subdf)}
            stat_row.update(s)
            stat_rows.append(stat_row)
    return pd.DataFrame(stat_rows)


def scatterplot_group_with_stats(
    df: pd.DataFrame,
    group: str = "All",
    forcefield_col: str = "FF",
    forcefield_order: list[str] = None,
    height: float = 4,
    aspect: float = 0.9,
    x_col: str = None,
    y_col: str = None,
    xerr_col: str = None,
    yerr_col: str = None,
    plot_hue: bool = False,
    stat_fontsize: float = 11,
    unit_str: str = "",
    stat_df: pd.DataFrame | None = None,
) -> tuple[sns.FacetGrid, pd.DataFrame]:
    """
    Create a Seaborn FacetGrid with scatter plots and stats

    Parameters
    ----------
    df : pd.DataFrame
        The input data frame containing the data to plot.
    group : str
        The name of the group to plot.
    forcefield_col : str
        The name of the column containing the force field information.
    forcefield_order : list[str]
        The order of the force fields to plot.
    height : float
        The height of the plots.
    aspect : float
        The aspect ratio of the plots.
    x_col : str
        The name of the column containing the x-axis data.
    y_col : str
        The name of the column containing the y-axis data.
    xerr_col : str
        The name of the column containing the x-axis error data.
    yerr_col : str
        The name of the column containing the y-axis error data.
    plot_hue : bool
        Whether to plot the hue based on the group.
    stat_fontsize : float
        The font size for the statistics text.
    unit_str : str
        The unit string to display on the axes.
    stat_df : pd.DataFrame | None
        The data frame containing the statistics to display.
        If this is pre-computed, this is used for statistics.
        If this is not provided, this is re-computed.

    Returns
    -------
    tuple[sns.FacetGrid, pd.DataFrame]
        The created FacetGrid and the statistics data frame.
    """
    extra_kwargs = {}
    subdf = df[df[group]]
    if forcefield_order:
        extra_kwargs["col_order"] = forcefield_order
    if plot_hue:
        extra_kwargs["hue"] = group
    else:
        df = subdf

    g = sns.FacetGrid(
        data=df,
        col=forcefield_col,
        **extra_kwargs,
        sharex=False,
        sharey=False,
        margin_titles=True,
        height=height,
        aspect=aspect,
    )
    g.map_dataframe(
        plot_errorbar, x_col=x_col, y_col=y_col, xerr_col=xerr_col, yerr_col=yerr_col
    )
    if plot_hue:
        g.add_legend()

    if stat_df is None:
        stat_df = compute_stats(
            subdf, x_col=x_col, y_col=y_col, xerr_col=xerr_col, yerr_col=yerr_col
        )
        stat_df["group"] = group

    limits = get_limits(df, x_col, y_col)

    for col_name, ax in g.axes_dict.items():
        ax.set_xlim(limits)
        ax.set_ylim(limits)
        ax.plot(limits, limits, ls="--", color="gray", lw=1)

        rows = stat_df[stat_df[forcefield_col] == col_name]
        assert len(rows) == 4, f"Expected 4 stats for {col_name}, got {len(rows)}"
        stat_str = f"{col_name}\n"
        for _, row in rows.iterrows():
            stat_str += (
                f"{row['stat']}={row['mle']:.2e} "
                f"[95%: {row['low']:.2e}, {row['high']:.2e}] \n"
            )
        ax.set_title(stat_str, fontsize=stat_fontsize)
        ax.set_ylabel(f"Computed {unit_str}")
        ax.set_xlabel(f"Reference {unit_str}")

    n_samples = stat_df["n"].values[0]
    g.figure.suptitle(f"{group} (n={n_samples})")

    return g, stat_df


@click.command(help=__doc__)
@click.option(
    "--input",
    "-i",
    "input_file",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    help="Input CSV file containing the data to plot.",
)
@click.option(
    "--output",
    "-o",
    "output_file",
    type=click.Path(exists=False, file_okay=True, dir_okay=False),
    help="Output CSV file for saving computed statistics.",
)
@click.option(
    "--images",
    "-im",
    "image_directory",
    type=click.Path(exists=False, file_okay=False, dir_okay=True),
    help="Output directory for saving plots.",
)
@click.option("--x-col", "-cx", "x_col", type=str, help="Column name for x-axis.")
@click.option("--y-col", "-cy", "y_col", type=str, help="Column name for y-axis.")
@click.option(
    "--xerr-col",
    "-cxe",
    "xerr_col",
    type=str,
    required=False,
    default=None,
    help="Column name for x-error.",
)
@click.option(
    "--yerr-col",
    "-cye",
    "yerr_col",
    type=str,
    required=False,
    default=None,
    help="Column name for y-error.",
)
@click.option(
    "--forcefield-col",
    "-cff",
    "forcefield_col",
    type=str,
    default="FF",
    help="Column name for forcefield.",
)
@click.option(
    "--forcefield",
    "-ff",
    "forcefield_order",
    multiple=True,
    type=str,
    help="Order of forcefields to plot.",
)
@click.option(
    "--height", "-h", "height", type=float, default=4, help="Height of the plots."
)
@click.option(
    "--aspect",
    "-a",
    "aspect",
    type=float,
    default=0.9,
    help="Aspect ratio of the plots.",
)
@click.option(
    "--unit-str",
    "-u",
    "unit_str",
    type=str,
    default="",
    help="Unit string for the axes.",
)
@click.option(
    "--plot-all-groups",
    "plot_all_groups",
    is_flag=True,
    default=False,
    help="Plot all groups.",
)
@click.option(
    "--plot-groups",
    "-g",
    "plot_groups",
    multiple=True,
    type=str,
    help="Specific groups to plot.",
)
@click.option(
    "--plot-hue",
    "plot_hue",
    is_flag=True,
    default=False,
    help=(
        "Use hue for plotting. "
        "If True, the points will be colored by the hue variable."
        "If False, only the plots in the subset will be plotted."
    ),
)
@click.option(
    "--n-min-entries",
    "-n",
    "n_min_entries",
    type=int,
    default=10,
    help="Minimum number of entries required to plot a group.",
)
def main(
    input_file: str,
    output_file: str,
    image_directory: str,
    x_col: str = None,
    y_col: str = None,
    xerr_col: str = None,
    yerr_col: str = None,
    forcefield_col: str = "FF",
    forcefield_order: list[str] = None,
    height: float = 4,
    aspect: float = 0.9,
    plot_all_groups: bool = False,
    plot_groups: list[str] = [],
    plot_hue: bool = False,
    stat_fontsize: float = 11,
    unit_str: str = "",
    n_min_entries: int = 10,
):
    df = pd.read_csv(input_file)
    df["All"] = True

    # get groups to plot
    groups_to_plot = list(plot_groups)
    if plot_all_groups:
        # use 'FF' or forcefield_col as start point
        index = list(df.columns).index(forcefield_col)
        groups_to_plot = df.columns[index + 1 :].tolist()

    # move 'All' to the start
    groups_to_plot.remove("All")
    groups_to_plot.insert(0, "All")

    image_directory = pathlib.Path(image_directory)
    image_directory.mkdir(parents=True, exist_ok=True)

    output_file = pathlib.Path(output_file)
    if output_file.exists():
        stat_df = pd.read_csv(output_file)
    else:
        stat_df = None

    stat_dfs = []
    first_ff = list(df[forcefield_col].unique())[0]
    ffdf = df[df[forcefield_col] == first_ff]
    for group in tqdm.tqdm(groups_to_plot, desc="Plotting groups"):
        # skip <n_min_entries entries
        n_entries = len(ffdf[ffdf[group]])
        if n_entries < n_min_entries:
            logger.info(f"Skipping group {group} with {n_entries} entries (< {n_min_entries})")
            continue

        group_stat_df = None
        if stat_df is not None:
            group_stat_df = stat_df[stat_df["group"] == group]

        g, stats_df = scatterplot_group_with_stats(
            df,
            group=group,
            forcefield_col=forcefield_col,
            forcefield_order=forcefield_order,
            height=height,
            aspect=aspect,
            x_col=x_col,
            y_col=y_col,
            xerr_col=xerr_col,
            yerr_col=yerr_col,
            plot_hue=plot_hue,
            stat_fontsize=stat_fontsize,
            unit_str=unit_str,
            stat_df=group_stat_df,
        )
        plt.tight_layout()

        subfile = f"scatter-{group}.png"
        if plot_hue:
            subfile = f"scatter-{group}-hue.png"
        imgfile = image_directory / subfile
        g.savefig(imgfile, dpi=300)
        logger.info(f"Saved scatter plot for group {group} to {imgfile}")
        plt.close()

        stat_dfs.append(stats_df)

    # Concatenate all stats DataFrames
    output_file = pathlib.Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    if not (output_file.is_file() and output_file.exists()):
        all_stats_df = pd.concat(stat_dfs, ignore_index=True)
        all_stats_df.to_csv(output_file, index=False)
        logger.info(f"Saved combined stats to {output_file}")


if __name__ == "__main__":
    main()
