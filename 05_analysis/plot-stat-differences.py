import pathlib
import sys

import numpy as np
import pandas as pd
from loguru import logger

import click
import seaborn as sns
from matplotlib import pyplot as plt

logger.remove()
logger.add(sys.stdout)

def plot_stats(df, output_directory, super_group: str):
    height = len(df.group.unique()) // 4
    for stat, subdf in df.groupby("stat"):
        fig, ax = plt.subplots(figsize=(5, height))

        # Encode group as numeric positions for y-axis
        groups = subdf["group"].unique()
        y_positions = {g: i for i, g in enumerate(groups)}
        for ff, subsubdf in subdf.groupby("FF"):
            # Plot each FF separately with distinct colours
            y = [y_positions[g] for g in subsubdf["group"]]
            xerr = [subsubdf["mle"] - subsubdf["low"], subsubdf["high"] - subsubdf["mle"]]
            ax.errorbar(
                subsubdf["mle"], y, xerr=xerr,
                fmt="o", capsize=4, label=ff
            )
        ax.set_yticks(np.arange(len(groups)))
        ax.set_yticklabels(groups)
        ax.set_xlabel(stat)
        ax.legend(title="FF")
        # plt.title(stat)
        plt.tight_layout()

        filename = output_directory / f"{super_group}_{stat}.png"
        plt.savefig(filename, dpi=300)
        logger.info(f"Saved plot to {filename}")
        plt.close()


@click.command()
@click.option(
    "--input",
    "-i",
    "input_file",
    type=click.Path(exists=True, readable=True, file_okay=True, dir_okay=False),
    default="output/phys-prop/density/stats.csv",
    help="Path to the input CSV file"
)
@click.option(
    "--output",
    "-o",
    "output_directory",
    type=click.Path(exists=False, writable=True, dir_okay=True),
    default="images/compare-smee/densities",
    help="Path to the output directory for images"
)
@click.option(
    "--min-n-samples",
    "-n",
    "min_n_samples",
    type=int,
    default=10,
    help="Minimum number of samples required"
)
@click.option(
    "--force-field-col",
    "-cff",
    "force_field_col",
    type=str,
    default="FF",
    help="Column name for force field"
)
@click.option(
    "--reference-force-field",
    "-r",
    "reference_force_field",
    type=str,
    default="Sage 2.2.1",
    help="Reference force field"
)
@click.option(
    "--target-force-field",
    "-t",
    "target_force_field",
    type=str,
    default="smee-v2",
    help="Target force field"
)
def main(
    output_directory: str = "images/compare-smee/densities",
    input_file: str = "output/phys-prop/density/stats.csv",
    min_n_samples: int = 10,
    force_field_col: str = "FF",
    reference_force_field: str = "Sage 2.2.1",
    target_force_field: str = "smee-v2"
):
    df = pd.read_csv(input_file)
    logger.info(f"Read {len(df)} entries from {input_file}")

    df = df[df["n"] >= min_n_samples]
    logger.info(f"Filtered to {len(df)} entries with n >= {min_n_samples}")


    reference_df = df[df[force_field_col] == reference_force_field]
    target_df = df[df[force_field_col] == target_force_field]

    logger.info(f"Reference group size: {len(reference_df)}")
    logger.info(f"Target group size: {len(target_df)}")

    # merge reference/target df on "group", "stat"
    merged_df = pd.merge(
        reference_df, target_df,
        on=["group", "stat"],
        suffixes=("_ref", "_target")
    )

    logger.info(f"Merged DataFrame size: {len(merged_df)}")

    merged_df["mle_diff"] = merged_df["mle_target"] - merged_df["mle_ref"]
    group_sorting = []
    for group in merged_df.sort_values("mle_diff")["group"].values:
        if group not in group_sorting:
            group_sorting.append(group)


    df = df[
        df[force_field_col].isin([reference_force_field, target_force_field])
    ]
    # sort groups according to group_sorting
    df = df.sort_values("group", key=lambda x: pd.Categorical(x, categories=group_sorting, ordered=True))

    # break up into checkmol/bonds/angles/propers/impropers
    bond_groups = [x for x in df.group.unique() if x.startswith("b")]
    angle_groups = [x for x in df.group.unique() if x.startswith("a")]
    improper_groups = [x for x in df.group.unique() if x.startswith("i")]
    proper_groups = [x for x in df.group.unique() if x.startswith("t")]
    vdw_groups = [x for x in df.group.unique() if x.startswith("n")]
    checkmol_groups = [x for x in df.group.unique() if x[0].lower() != x[0]]

    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)

    checkmol_df = df[df["group"].isin(checkmol_groups)]
    plot_stats(checkmol_df, output_directory, "checkmol")

    bond_df = df[df["group"].isin(bond_groups)]
    plot_stats(bond_df, output_directory, "bonds")

    angle_df = df[df["group"].isin(angle_groups)]
    plot_stats(angle_df, output_directory, "angles")

    # improper_df = df[df["group"].isin(improper_groups)]
    # plot_stats(improper_df, output_directory, "impropers")

    proper_df = df[df["group"].isin(proper_groups)]
    plot_stats(proper_df, output_directory, "propers")

    vdw_df = df[df["group"].isin(vdw_groups)]
    plot_stats(vdw_df, output_directory, "vdw")

if __name__ == "__main__":
    main()
