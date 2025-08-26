"""
Plot ddE data from QM benchmark, using all-to-all RMSD data.

\b
This script saves four plots:
- dde.png: A plot of ddE values for each force field, limited to 12 kcal/mol.
- dde-close.png: A close-up of the ddE plot, limited to 5 kcal/mol.
- dde-step.png: A step plot of ddE values without a legend.
- dde-step-legend.png: A step plot of ddE values with a legend.
"""

import functools
import pathlib
import logging
import click


import pandas as pd
import pyarrow.compute as pc
import pyarrow.dataset as ds


from openff.toolkit import Molecule


import seaborn as sns
from matplotlib import pyplot as plt

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)


@functools.cache
def smiles_to_inchi(smiles: str) -> str:
    mol = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
    return mol.to_inchi(fixed_hydrogens=True)


def plot_ecdf(
    data: pd.DataFrame,
    palette: dict,
    xmax=1,
    xlabel: str = r"|ddE| ($\mathrm{kcal \, mol^{-1}}$)",
    imgfile: pathlib.Path = None,
):
    """
    Plot an ECDF of ddE values.
    This is a side-effecting function that saves the plot to a file.

    Parameters
    ----------
    data : pd.DataFrame
        The data to plot. Must contain columns "|ddE|", "Force field", and "method".
    palette : dict
        A dictionary mapping force field names to colors.
    xmax : float, optional
        The maximum x value to plot, by default 1.
    xlabel : str, optional
        The label for the x axis, by default r"|ddE| ($\mathrm{kcal \, mol^{-1}}$)".
    imgfile : pathlib.Path, optional
        The file to save the plot to.
    """

    fig, ax = plt.subplots(figsize=(8, 6))
    rmsd_ax = sns.ecdfplot(
        ax=ax,
        data=data,
        x="|ddE|",
        hue="Force field",
        palette=palette,
        stat="count",
    )
    rmsd_ax.set_xlim(0, xmax)
    rmsd_ax.set_xlabel(xlabel)
    plt.tight_layout()
    plt.savefig(imgfile, dpi=300)
    logger.info(f"Saved RMSD plot to {imgfile}")
    plt.close()


@click.command(help=__doc__)
@click.option(
    "--ff-name-and-stem",
    "-ff",
    "ff_name_and_stems",
    type=(str, str),
    multiple=True,
    help=(
        "Force fields to include and their stem names. "
        "The first argument should be the name of the force field in the plot, "
        "the second argument should be the stem name of the force field (e.g., 'openff_unconstrained-2.2.1'). "
        "This option can be specified multiple times to include multiple force fields."
    )
)
@click.option(
    "--input-directory",
    "-i",
    "input_directory",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    default="rmsd",
    help="Directory to read RMSD/TFD data from.",
)
@click.option(
    "--output-directory",
    "-o",
    "output_directory",
    type=click.Path(exists=False, file_okay=False, dir_okay=True),
    default="images",
    help="Directory to write output files to.",
)
def main(
    ff_name_and_stems: list[tuple[str, str]],
    input_directory: str = "ddEs",
    output_directory: str = "images",
):

    # load input ddE data
    dataset = ds.dataset(input_directory)
    unique_ffs = [
        ff_name_and_stem[0] for ff_name_and_stem in ff_name_and_stems
    ]
    dataset = dataset.filter(
        pc.field("Force field").isin(unique_ffs)
    )
    df = dataset.to_table().to_pandas()
    df["|ddE|"] = df["ddE"].abs()
    logger.info(f"Loaded {len(df)} ddE records from {input_directory}")

    # plot info about loaded data
    n_ddes = df.groupby("Force field").size().to_dict()
    count_str = ", ".join(
        [f"{ff}: {count}" for ff, count in n_ddes.items()]
    )
    logger.info(f"Counts of ddE entries per force field -- {count_str}")

    n_conformers = df.groupby("Force field")["n_conformers"].sum().to_dict()
    count_str = ", ".join(
        [f"{ff}: {count}" for ff, count in n_conformers.items()]
    )
    logger.info(f"Counts of conformers per force field -- {count_str}")


    # count loaded FFs
    n_unique_ffs = len(unique_ffs)
    logger.info(f"Found {n_unique_ffs} unique force fields in input {input_directory}")

    palette = sns.color_palette(n_colors=len(unique_ffs))
    PALETTE = {
        unique_ffs[i]: palette[i] for i in range(len(unique_ffs))
    }

    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)

    # plot ddE, close-up
    for xmax, filename in [(12, "dde.png"), (5, "dde-close.png")]:
        plot_ecdf(
            data=df,
            palette=PALETTE,
            xmax=xmax,
            xlabel=r"ddE ($\mathrm{kcal \, mol^{-1}}$)",
            imgfile=output_directory / filename,
        )

    # plot step plot
    for legend in [True, False]:
        fig, ax = plt.subplots(figsize=(6, 4))
        sns.histplot(
            ax=ax,
            data=df,
            x="ddE",
            hue="Force field",
            palette=PALETTE,
            binrange=(-10.25, 10.25),
            binwidth=0.5,
            element="step",
            fill=False,
            legend=legend,
        )
        ax.set_xlabel("ddE (kcal/mol)")
        plt.tight_layout()
        imgfile = output_directory / f"dde-step{'-legend' if legend else ''}.png"
        plt.savefig(imgfile, dpi=300)
        logger.info(f"Saved ddE step plot to {imgfile}")
        plt.close()



if __name__ == "__main__":
    main()
