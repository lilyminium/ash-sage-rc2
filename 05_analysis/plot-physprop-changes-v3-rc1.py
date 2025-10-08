import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

PROPS = {
    "density": {"label": "Density", "units": "g/mL"},
    "dhmix": {"label": "∆Hmix", "units": "kJ/mol"},
}

PALETTE = {
    "Sage 2.3.0rc1": "tab:green",
    "v3-k100": "tab:purple",
}


def main():
    for prop, info in PROPS.items():
        v3 = pd.read_csv(f"output/phys-prop/{prop}-all/v3-k100-vs-Sage 2.2.1.csv")
        rc1 = pd.read_csv(f"output/phys-prop/{prop}-all/Sage 2.3.0rc1-vs-Sage 2.2.1.csv")
        rc1["Difference from 2.2.1"] = rc1["Sage 2.3.0rc1 - Sage 2.2.1"]
        rc1["FF"] = "Sage 2.3.0rc1"
        v3["Difference from 2.2.1"] = v3["v3-k100 - Sage 2.2.1"]
        v3["FF"] = "v3-k100"
        df = pd.concat([rc1, v3])

        ax = sns.histplot(
            data=df,
            x="Difference from 2.2.1",
            hue="FF",
            palette=PALETTE,
            kde=True,
            alpha=0.5,
            lw=1
        )
        ax.set_title(f"{info['label']} (n={len(df['id'].unique())})")
        ax.axvline(0, ls="--", color="k")

        # add medians
        for ff in df["FF"].unique():
            median_diff = df[df["FF"] == ff]["Difference from 2.2.1"].median()
            color = PALETTE.get(ff, "r")
            ax.axvline(median_diff, ls=":", color=color)
            # set top left
            ax.text(
                0.05,
                0.95 - 0.05 * list(df["FF"].unique()).index(ff),
                f"{ff} median: {median_diff:.3f} {info['units']}",
                transform=ax.transAxes,
                verticalalignment="top",
                color=color,
            )

        ax.set_xlabel(f"{info['label']} change from 2.2.1 ({info['units']})")
        plt.tight_layout()
        plt.savefig(f"images/compare-refit/{prop}-changes/change-in-{prop}-v3-rc1.png", dpi=300)
        plt.close()
        print(f"Saved images/compare-refit/{prop}-changes/change-in-{prop}-v3-rc1.png")


if __name__ == "__main__":
    main()
