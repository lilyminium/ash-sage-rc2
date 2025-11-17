import pathlib
import matplotlib as mpl
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

OPENFF_BLUE = "#015480"
OPENFF_LIGHT_BLUE = "#2F9ED2"
OPENFF_ORANGE = "#F08521"
OPENFF_RED = "#F03A21"
OPENFF_GRAY = "#3E424A"

COLORS = {
    "blue": OPENFF_BLUE,
    "cyan": OPENFF_LIGHT_BLUE,
    "orange": OPENFF_ORANGE,
    "red": OPENFF_RED,
    "gray": OPENFF_GRAY
}

sns.set_context("paper")
sns.set_style("ticks")
mpl.rcParams['font.sans-serif'] = ["muli"]

import pyarrow.dataset as ds
import pyarrow.compute as pc
import seaborn as sns
from matplotlib import pyplot as plt

PANEL_MAP = {
    "RMSE": "Error",
    "MUE": "Error",
    "MSE": "Error",
    "rho": "Correlation",
    "R2": "Correlation",
}

FORCEFIELD_ORDER = ["Sage 2.2.1 + ELF10", "Sage 2.2.1 + AshGC", "Sage 2.3.0"]


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


def read_into_df(
    directory: pathlib.Path
):
    dfs = []
    FFS = {
        "v3-k100": "Sage 2.3.0"
    }
    directory = pathlib.Path(directory)
    for fn in [
        # directory / "Sage-2.1.0.csv",
        directory / "Sage-2.2.1-+-ELF10.csv",
        directory / "Sage-2.2.1-+-AshGC.csv",
        directory / "v3-k100.csv"
    ]:
        df = pd.read_csv(fn)
        df["FF"] = [FFS.get(x, x) for x in df.FF.values]
        dfs.append(df)

    df = pd.concat(dfs).reset_index()
    return df

def print_latex_stats(df):
    for _, row in df.iterrows():
        print(f"{row['FF']} | {row['stat']}={row['mle']:.3f}_{{{row['low']:.3f}}}^{{{row['high']:.3f}}}")
