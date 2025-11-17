import pandas as pd
import numpy as np
import pyarrow.compute as pc
import pyarrow.dataset as ds

import seaborn as sns
import tqdm
from matplotlib import pyplot as plt

import argparse
from collections import namedtuple
import pathlib
import matplotlib as mpl
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import yaml

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

def main():
    df = pd.read_csv("processed_bulk_qm_properties.csv")
    df_ = pd.read_csv("processed_improper_torsions.csv")
    df = pd.concat([df, df_], axis=0, ignore_index=True)

    ff_order = ["Sage 2.1.0", "Sage 2.2.1", "Sage 2.3.0",]

    g = sns.catplot(
        kind="box",
        data=df,
        col="Property",
        col_order=[
            "rmsd", "tfd", "|ddE|",
            # "Bonds ICRMSD", "Angles ICRMSD", "ProperTorsions ICRMSD",
            # "ImproperTorsions ICRMSD",
        ],
        col_wrap=3,
        hue="FF",
        y="FF",
        x="Value",
        order=ff_order,
        hue_order=ff_order,
        sharex=False,
        height=1.5,
        aspect=1.5,
        fliersize=2,
        linewidth=0.5,
        whis=(5, 95),
        flierprops={"markeredgewidth": 0.5},
    )
    axes = g.axes.flatten()
    axes[0].set_xlabel("RMSD [$\AA$]")
    axes[1].set_xlabel("TFD")
    axes[2].set_xlabel("|ddE| [kcal/mol]")
    g.set_titles("")

    plt.savefig("images/qm-box.png", dpi=300)
    plt.close()

    g = sns.catplot(
        kind="box",
        data=df,
        col="Property",
        col_order=[
            # "rmsd", "tfd", "|ddE|",
            "Bonds ICRMSD", "Angles ICRMSD", "ProperTorsions ICRMSD",
            "ImproperTorsions ICRMSD",
        ],
        col_wrap=2,
        hue="FF",
        y="FF",
        x="Value",
        order=ff_order,
        hue_order=ff_order,
        sharex=False,
        height=1.5,
        aspect=1.5,
        fliersize=2,
        linewidth=0.5,
        whis=(5, 95),
        flierprops={"markeredgewidth": 0.5},
    )
    axes = g.axes.flatten()
    axes[0].set_xlabel("Bond ICRMSD [$\AA$]")
    axes[1].set_xlabel("Angle ICRMSD [$\degree$]")
    axes[2].set_xlabel("Proper torsion ICRMSD [$\degree$]")
    axes[3].set_xlabel("Proper torsion ICRMSD [$\degree$]")
    g.set_titles("")

    plt.savefig("images/qm-box-icrmsd.png", dpi=300)
    plt.close()
    # print("Saved figure to images/qm-box.png")

if __name__ == "__main__":
    main()
