import pathlib

import click
import tqdm

import pandas as pd

from openff.evaluator.client.client import RequestResult
from openff.evaluator.datasets.datasets import PhysicalPropertyDataSet
from openff.toolkit import Molecule, ForceField
import seaborn as sns
import matplotlib.pyplot as plt


@click.command()
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

    both_df = pd.concat([rmin_df, epsilon_df], ignore_index=True)

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


if __name__ == "__main__":
    main()
