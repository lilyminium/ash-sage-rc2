import pathlib
import sys
import typing

import click

from loguru import logger
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

logger.remove()
logger.add(sys.stdout)


@click.command()
@click.option(
    "--input-file",
    "-i",
    type=click.Path(exists=True, dir_okay=False, readable=True),
    required=True,
    help="Path to the input CSV file with ratio representation"
)
@click.option(
    "--output-file",
    "-o",
    type=click.Path(dir_okay=False, writable=True),
    required=True,
    help="Path to the output PNG file"
)
@click.option(
    "--description",
    "-d",
    type=str,
    default="Densities > 0.05 g/mL",
    help="Description of error set for the plot label"
)
@click.option(
    "--group-type",
    "-g",
    type=click.Choice(["Checkmol", "Bonds", "Angles", "ProperTorsions", "ImproperTorsions", "vdW"]),
    default="Checkmol",
    help="Type of group to plot"
)
@click.option(
    "--force-field-order",
    "-ff",
    multiple=True,
    type=str,
)
@click.option(
    "--force-field-col",
    "-cff",
    type=str,
    default="FF",
)
@click.option(
    "--threshold",
    "-t",
    type=float,
    default=2,
    help="Threshold for ratio representation"
)
def main(
    input_file: str,
    output_file: str,
    description: str,
    group_type: typing.Literal["Checkmol", "Bonds", "Angles", "ProperTorsions", "ImproperTorsions", "vdW"],
    force_field_order: list[str],
    force_field_col: str = "FF",
    threshold: float = 2,
):
    df = pd.read_csv(input_file).sort_values("ratio_representation", ascending=False)
    logger.info(f"Read {len(df)} entries from {input_file}")

    # filter group type
    df = df[df["Type"] == group_type]
    logger.info(f"Filtered to {len(df)} entries of type {group_type}")

    # filter just to FFs
    df = df[df[force_field_col].isin(force_field_order)]
    logger.info(f"Filtered to {len(df)} entries for force fields {force_field_order}")

    over_threshold = df[df.ratio_representation > threshold]
    logger.info(f"Found {len(df)} entries with ratio_representation > {threshold}")

    # get all force fields for all groups
    matching_groups = over_threshold.Group.unique()
    logger.info(f"Found {len(matching_groups)} unique groups with ratio_representation > {threshold}")

    # filter to include groups for all FFs even if under threshold
    df = df[df.Group.isin(matching_groups)]
    logger.info(f"Filtered to {len(df)} entries in {len(matching_groups)} groups")

    # barplot
    ax = sns.barplot(
        data=df,
        y="Group",
        x="ratio_representation",
        hue=force_field_col,
        hue_order=force_field_order,
        legend=True,
    )
    ax.set_xlabel(f"Ratio of fraction of mols with {description}\n/ fraction of mols in entire set for each group")
    ax.set_ylabel("Group")
    ax.set_title(group_type)
    ax.axvline(1, color="gray", ls="--", lw=1)
    plt.tight_layout()

    output_file = pathlib.Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    plt.savefig(output_file, dpi=300)
    logger.info(f"Saved plot to {output_file}")


if __name__ == "__main__":
    main()
