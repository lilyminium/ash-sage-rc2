"""
Plot parameter values over fitting iterations.

This takes as input the refitting directory and assumes it contains
an `optimize.tmp` directory with the necessary files.
It also assumes that the reference force field is in `forcefield/force-field.offxml`.
"""
from collections import defaultdict
import pathlib
import sys

import click
import tqdm

import pandas as pd

from loguru import logger

from openff.toolkit import Molecule, ForceField
import seaborn as sns
import matplotlib.pyplot as plt


logger.remove()
logger.add(sys.stdout)

@click.command(help=__doc__)
@click.option(
    "--input-path",
    "-i",
    default="../refit",
    help="The directory containing the results files from the refit",
)
@click.option(
    "--output-file",
    "-o",
    default="images/parameters-over-time.png",
    help="The output image file to save the parameters over time",
)
def main(
    input_path: str = "../refit",
    output_file: str = "images/parameters-over-time.png"
):
    input_path = pathlib.Path(input_path)

    ff_directory = input_path / "optimize.tmp"
    ff_files = sorted(ff_directory.glob("phys-prop/iter*/force-field.offxml"))

    ref_ff_path = input_path / "forcefield" / "force-field.offxml"
    reference = ForceField(str(ref_ff_path.resolve()), allow_cosmetic_attributes=True)
    ref_vdw_handler = reference.get_parameter_handler("vdW")

    # get rmins and epsilons over iterations
    rmin_data = {}
    epsilon_data = {}
    epsilon_constrained_data = defaultdict(list)
    for prop in ref_vdw_handler.parameters:
        if not hasattr(prop, "_parameterize"):
            continue
        rmin_data[prop.id] = [prop.rmin_half.m]
        epsilon_data[prop.id] = [prop.epsilon.m]

    iter_cols = []
    for ff_file in tqdm.tqdm(ff_files):
        iter_cols.append(ff_file.parent.name)
        ff = ForceField(ff_file, allow_cosmetic_attributes=True)
        vdw_handler = ff.get_parameter_handler("vdW")
        for param in vdw_handler.parameters:
            if not hasattr(param, "_parameterize"):
                continue
            rmin_data[param.id].append(param.rmin_half.m)
            epsilon_data[param.id].append(param.epsilon.m)
            epsilon_constrained_data[param.id].append(param._constrained_epsilon.m)

    for v in epsilon_constrained_data.values():
        v.insert(0, v[0])

    rmin_df = pd.DataFrame.from_dict(
        rmin_data,
        orient="index",
        columns=["Reference"] + iter_cols
    ).reset_index(
        names=["Parameter"]
    ).melt(
        id_vars=["Parameter", "Reference"],
        value_vars=iter_cols,
        var_name="Iteration",
        value_name="Value",
    )
    rmin_df["Iteration"] = [int(x.split("_")[-1]) for x in rmin_df.Iteration.values]
    rmin_df["Attribute"] = "rmin_half"

    epsilon_df = pd.DataFrame.from_dict(
        epsilon_data,
        orient="index",
        columns=["Reference"] + iter_cols
    ).reset_index(
        names=["Parameter"]
    ).melt(
        id_vars=["Parameter", "Reference"],
        value_vars=iter_cols,
        var_name="Iteration",
        value_name="Value",
    )
    epsilon_df["Iteration"] = [int(x.split("_")[-1]) for x in epsilon_df.Iteration.values]
    epsilon_df["Attribute"] = "epsilon"

    constrained_epsilon_df = pd.DataFrame.from_dict(
        epsilon_constrained_data,
        orient="index",
        columns=["Reference"] + iter_cols
    ).reset_index(
        names=["Parameter"]
    ).melt(
        id_vars=["Parameter", "Reference"],
        value_vars=iter_cols,
        var_name="Iteration",
        value_name="Value",
    )
    constrained_epsilon_df["Iteration"] = [int(x.split("_")[-1]) for x in constrained_epsilon_df.Iteration.values]
    constrained_epsilon_df["Attribute"] = "constrained_epsilon"

    both_df = pd.concat([rmin_df, epsilon_df, constrained_epsilon_df], ignore_index=True)

    # map SMIRKS
    id_to_smirks = {
        p.id: p.smirks
        for p in ref_vdw_handler.parameters
        if hasattr(p, "_parameterize")
    }
    both_df["smirks"] = [id_to_smirks[x] for x in both_df.Parameter.values]
    both_df["SMIRKS"] = [
        x.replace('#7,#8,#9,#16,#17,#35', 'ENA')
        for x in both_df.smirks.values
    ]
    both_df.to_csv("output/parameters-over-iterations.csv", index=False)

    g = sns.FacetGrid(
        data=both_df,
        row="Parameter",
        col="Attribute",
        aspect=2.5,
        height=1.2,
        sharey=False,
        margin_titles=True
    )
    g.map(sns.lineplot, "Iteration", "Value")
    g.set_titles(col_template="{col_name}", row_template="{row_name}")
    plt.tight_layout()
    g.savefig(output_file, dpi=300)
    logger.info(f"Saved figure to {output_file}")


if __name__ == "__main__":
    main()
