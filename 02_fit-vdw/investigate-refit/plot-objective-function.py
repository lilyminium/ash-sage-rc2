import pathlib
import re
import sys

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import click
from loguru import logger

logger.remove()
logger.add(sys.stdout)

NUM  = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
PATTERN = re.compile(
    rf"(?ms)"
    rf"\s*#\|\s*Iteration\s+(?P<iteration>\d+):\s*Evaluating\s+objective\s+function\s*\|#.*?"
    rf"#\|\s*Objective\s+Function\s+Breakdown\s*\|#.*?"
    rf"phys-prop\s+{NUM}\s+{NUM}\s+(?P<phys_prop>{NUM})[^\n]*\n"
    rf"Regularization\s+{NUM}\s+{NUM}\s+(?P<regularization>{NUM})[^\n]*\n"
    rf"Total[^\n]*?(?P<total>{NUM})[^\n]*"
)

PALETTE = {
    "phys_prop": "tab:blue",
    "regularization": "tab:orange",
    "total": "tab:green",
}

def main(
    refit_directory: str = "../refit",
):
    fb_files = sorted(
        pathlib.Path(refit_directory).glob("force_balance*.log")
    )
    logger.info(f"Found {len(fb_files)} force balance log files")

    all_rows = []
    for file in fb_files:
        contents = file.read_text()
        contents = re.sub(r"\x1b\[[0-9;]*m", "", contents)
        rows = [m.groupdict() for m in PATTERN.finditer(contents)]
        for row in rows:
            row["iteration"] = int(row["iteration"])
            row["phys_prop"] = float(row["phys_prop"])
            row["regularization"] = float(row["regularization"])
            row["total"] = float(row["total"])
            all_rows.append(row)

    df = pd.DataFrame(all_rows).sort_values("iteration")

    melted = df.melt(id_vars="iteration", value_name="Contribution", var_name="Target")
    ax = sns.lineplot(
        data=melted,
        x="iteration",
        y="Contribution",
        hue="Target",
        palette=PALETTE,
    )
    min_iter = df[df.total == df.total.min()].iteration.values[0]
    min_physprop_iter = df[df.phys_prop == df.phys_prop.min()].iteration.values[0]
    ax.set_xlabel("Iteration")
    ax.axvline(min_iter, ls="--", color=PALETTE["total"])
    ax.axvline(min_physprop_iter, ls="--", color=PALETTE["phys_prop"])
    plt.savefig("images/objective-function.png", dpi=300)
    logger.info("Saved objective function plot to images/objective-function.png")
    plt.close()

if __name__ == "__main__":
    main()
