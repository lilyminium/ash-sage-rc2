import pathlib
import matplotlib as mpl
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

sns.set_context("paper")
sns.set_style("ticks")
mpl.rcParams['font.sans-serif'] = ["muli"]

FORCEFIELD_ORDER = ["Sage 2.2.1", "Sage 2.2.1 + AshGC", "Sage 2.3.0"]


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


def main():
    densities = pd.read_csv("processed_density_validation_stats.csv")
    densities["Dataset"] = "Density$_{}$" + f"(n={densities.n.values[0]})"
    dhmixes = pd.read_csv("processed_dhmix_validation_stats.csv")
    dhmixes["Dataset"] = "∆H$_{mix}$ " + f"(n={dhmixes.n.values[0]})"
    

    df = pd.concat([densities, dhmixes])
    
    STAT_ORDER = ["RMSE", "MSE"]
    df["stat"] = pd.Categorical(df["stat"], categories=STAT_ORDER, ordered=True)

    # Assign numeric x positions for dodge
    x_lookup = {cat: i for i, cat in enumerate(STAT_ORDER)}
    df["xbase"] = df["stat"].map(x_lookup).astype(float)
    dodge_spacing = 0.2

    # Compute dodge offsets by FF
    offsets = {
        ff: (i - (len(FORCEFIELD_ORDER) - 1) / 2) * dodge_spacing
        for i, ff in enumerate(FORCEFIELD_ORDER)
    }
    df["xpos"] = df["xbase"] + df["FF"].map(offsets).astype(float)

    g = sns.FacetGrid(
        data=df,
        col="Dataset",
        hue="FF",
        hue_order=FORCEFIELD_ORDER,
        sharex=False,
        sharey=False,
        legend_out=True,
        height=2.,
        aspect=0.8
    )
    g.map_dataframe(draw_errorbars)

    for dataset, ax in g.axes_dict.items():
        ax.set_xticks([x_lookup[s] for s in STAT_ORDER])
        ax.set_xticklabels([s for s in STAT_ORDER])
        ax.set_xlabel("")
        ax.axhline(0, ls="--", color="gray")
        # ax.set_ylim(-0.6, 1.7)
        ax.set_xlim(-0.5, 1.5)
        ax.set_ylabel("")

    ax1, ax2 = g.axes.flatten()
    ax1.set_ylabel("Error [g/mL]")
    ax2.set_ylabel("Error [kJ/mol]")

    plt.tight_layout()
    g.add_legend()
    g.set_titles(col_template="{col_name}")
    plt.savefig("images/combined-phys-props.png", dpi=300)
    print("Saved figure to images/combined-phys-props.png")

if __name__ == "__main__":
    main()