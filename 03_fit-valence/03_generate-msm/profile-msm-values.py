"""
Plot MSM distributions of force field parameters.

This script reads in an MSM data directory and a force field JSON file
and generates distribution plots for each parameter type in the force field.
Plots are saved to an output directory organized by parameter type, by parameter ID.
"""

import click
import json
import sys
import pathlib

import tqdm
from loguru import logger

import seaborn as sns
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

from openff.units import unit

logger.remove()
logger.add(sys.stdout)


@click.command(help=__doc__)
@click.option(
    "--input-file",
    "-i",
    type=str,
    default="msm-ff/ff-v4.json",
    help="Input force field JSON file.",
)
@click.option(
    "--output-directory",
    "-o",
    type=str,
    default="images/distributions/ff-v4",
    help="Output directory to write plots.",
)
def main(
    input_file: str = "msm-ff/ff-v0.json",
    output_directory: str = "images/distributions/ff-v0"
):
    with open(input_file, "r") as f:
        data = json.load(f)

    kj_per_mol_per_nm2 = unit.kilojoule_per_mole / (unit.nanometer ** 2)
    kcal_per_mol_per_a2 = unit.kilocalorie_per_mole / (unit.angstrom ** 2)
    kj_per_mol_per_rad2 = unit.kilojoule_per_mole / (unit.radian ** 2)
    kcal_per_mol_per_rad2 = unit.kilocalorie_per_mole / (unit.radian ** 2)

    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    
    for parameter_type, subdict in data.items():
        output_subdirectory = output_directory / parameter_type.lower()
        output_subdirectory.mkdir(parents=True, exist_ok=True)

        if parameter_type == "Bonds":
            xlabel = "Equilibrium length (Å)"
            ylabel = "Force constant (kcal/(mol·Å²))"
            xfactor = (1 * unit.nanometer).to(unit.angstrom).m
            yfactor = (1 * kj_per_mol_per_nm2).to(kcal_per_mol_per_a2).m

        elif parameter_type == "Angles":
            xlabel = "Equilibrium angle (degrees)"
            ylabel = "Force constant (kcal/(mol·rad²))"
            xfactor = (1 * unit.radian).to(unit.degree).m
            yfactor = (1 * kj_per_mol_per_rad2).to(kcal_per_mol_per_rad2).m

        else:
            raise ValueError(f"Unknown parameter type {parameter_type}")
        
        for parameter_id, parameter_data in tqdm.tqdm(subdict.items()):
            smirks = parameter_data["smirks"]
            eq_values = np.array(parameter_data["eq"])
            k_values = np.array(parameter_data["k"])

            df = pd.DataFrame({
                xlabel: eq_values * xfactor,
                ylabel: k_values * yfactor,
            })
            ax = sns.kdeplot(
                data=df,
                x=xlabel,
                y=ylabel,
                fill=True,
            )
            ax.set_title(f"{parameter_id} (n={len(df)})\n{smirks}")
            plt.tight_layout()
            filename = output_subdirectory / f"{parameter_id}.png"
            plt.savefig(filename, dpi=300)
            plt.close()
            logger.info(f"Wrote {filename}")


if __name__ == "__main__":
    main()
