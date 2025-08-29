"""
This script plots the percentage changes in epsilon and rmin_half vdW parameters between
two force field files.

It saves an output CSV for each parameterized parameter (marked by an `_parameterize` attribute).
The CSV is saved in the output directory under the name `parameter-changes.csv`.
It contains the columns `smirks`, `SMIRKS`, `initial_epsilon`, `final_epsilon`, `initial_rmin_half`, `final_rmin_half`, `% epsilon`, and `% rmin_half`.
The difference between the `smirks` and `SMIRKS` columns is that the SMARTS specifying
electronegative atoms (`'#7,#8,#9,#16,#17,#35'`) is replaced with the short-hand `ENA`
for better visual rendering.

Images are saved in the specified image directory under the names
`percent-changes-epsilon.png` and `percent-changes-rmin_half.png`.
Several attributes of the images (e.g. xlim) are hardcoded in the script
and may need to be customized in future applications.
"""

import pathlib
import sys

from loguru import logger
from openff.toolkit import Molecule, ForceField, unit
import pandas as pd
import matplotlib.pyplot as plt

import click

logger.remove()
logger.add(sys.stdout)

def plot_percent(df, col="% epsilon", xlabel="$\epsilon$, ", lim=(-20, 20)):
    fig, ax = plt.subplots(figsize=(7, 5))
    bars = ax.barh(
        df.SMIRKS.values,
        df[col].values,
    )

    # Add a vertical line at x=0 for the origin
    ax.axvline(0, color='black', linewidth=1, ls="--")

    # Add value labels
    for bar in bars:
        width = bar.get_width()
        x = width + (1 if width > 0 else -1)
        y = bar.get_y() + bar.get_height() / 2
        ha ='left' if width > 0 else 'right'
        
        if abs(width) > lim[1]:
            inc = (0.1 * (lim[1] - lim[0]))
            x = lim[1] - inc
            if width < 0:
                x = lim[0] + inc
            ax.text(x, y, "//  ", va="center", ha=ha, fontsize=18)
            y -= bar.get_height()
        ax.text(x, y, f'{width:.1f} %', va='center', ha=ha)

    # Formatting
    ax.set_xlim(*lim)
    ax.set_xlabel('% change ' + xlabel)
    plt.tight_layout()
    return ax


@click.command(help=__doc__)
@click.option(
    "--input-forcefield",
    "-if",
    default="openff-2.2.1.offxml",
    help="The input force field file to compare",
)
@click.option(
    "--output-forcefield",
    "-of",
    default="forcefield/force-field.offxml",
    help="The output force field file to compare",
)
@click.option(
    "--output-directory",
    "-o",
    default="output",
    help="The directory to save the output CSV file",
)
@click.option(
    "--image-directory",
    "-im",
    default="images",
    help="Directory to save the high-resolution images of parameter changes",
)
def main(
    input_forcefield: str = "openff-2.2.1.offxml",
    output_forcefield: str = "../refit/result/optimize/force-field.offxml",
    output_directory: str = "output",
    image_directory: str = "images",
):
    old = ForceField(input_forcefield)
    old_vdw_handler = old.get_parameter_handler("vdW")

    new = ForceField(output_forcefield, allow_cosmetic_attributes=True)
    new_vdw_handler = new.get_parameter_handler("vdW")

    data = []
    for new_param in new_vdw_handler.parameters:
        if hasattr(new_param, "_parameterize"):
            old_param = old_vdw_handler[new_param.smirks]
            row = {
                "smirks": new_param.smirks,
                "initial_epsilon": old_param.epsilon.m_as(unit.kilocalories_per_mole),
                "final_epsilon": new_param.epsilon.m_as(unit.kilocalories_per_mole),
                "initial_rmin_half": old_param.rmin_half.m_as(unit.angstrom),
                "final_rmin_half": new_param.rmin_half.m_as(unit.angstrom),
            }
            data.append(row)

    df = pd.DataFrame(data)

    df["% epsilon"] = 100 * (df.final_epsilon - df.initial_epsilon) / df.initial_epsilon
    df["% rmin_half"] = 100 * (df.final_rmin_half - df.initial_rmin_half) / df.initial_rmin_half
    df["SMIRKS"] = [
        x.replace('#7,#8,#9,#16,#17,#35', 'ENA')
        for x in df.smirks.values
    ]

    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    output_file = output_directory / "parameter-changes.csv"
    df.to_csv(output_file, index=False)
    logger.info(f"Saved parameter changes to {output_file}")

    # Plotting
    image_directory = pathlib.Path(image_directory)
    image_directory.mkdir(parents=True, exist_ok=True)

    ax = plot_percent(df, col="% epsilon", xlabel="$\epsilon$", lim=(-70, 70))
    imgfile = image_directory / "percent-changes-epsilon.png"
    plt.savefig(imgfile, dpi=300)
    logger.info(f"Saved epsilon percent changes plot to {imgfile}")
    plt.close()

    ax = plot_percent(df, col="% rmin_half", xlabel="$r_{min}/2$", lim=(-10, 10))
    imgfile = image_directory / "percent-changes-rmin-half.png"
    plt.savefig(imgfile, dpi=300)
    logger.info(f"Saved rmin_half percent changes plot to {imgfile}")
    plt.close()
    

if __name__ == "__main__":
    main()
