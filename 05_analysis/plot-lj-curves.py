"""
Generate and plot Lennard-Jones curves for comparison between two force fields.
Only vdW parameters that differ between the two force fields are plotted.
"""

import click
import sys
import pathlib

import tqdm
import numpy as np

from loguru import logger
from matplotlib import pyplot as plt
from openff.toolkit import ForceField, unit

logger.remove()
logger.add(sys.stdout)

OPENFF_BLUE = "#015480"
OPENFF_ORANGE = "#F08521"

def calculate_lj(r: float | np.ndarray, sigma: float, epsilon: float) -> float | np.ndarray:
    """Calculate the Lennard-Jones potential at a given distance r."""
    r6 = (sigma / r) ** 6
    r12 = r6 ** 2
    return 4 * epsilon * (r12 - r6)


@click.command()
@click.option(
    "--reference-force-field",
    "-rf",
    default="openff-2.2.1.offxml",
    help="The reference force field file.",
)
@click.option(
    "--reference-name",
    "-rn",
    default="Sage 2.2.1",
    help="The name to use for the reference force field in plots.",
)
@click.option(
    "--target-force-field",
    "-tf",
    default="../04_benchmark/forcefields/fb-fit-v3-single-mean-k100.offxml",
    help="The target force field file to compare against the reference.",
)
@click.option(
    "--target-name",
    "-tn",
    default="v3-k100",
    help="The name to use for the target force field in plots.",
)
@click.option(
    "--image-directory",
    "-id",
    default="images/lj-curves/v3-vs-221",
    type=click.Path(exists=False, dir_okay=True, file_okay=False),
    help="Path to the directory to save the output images.",
)
def main(
    reference_force_field: str = "openff-2.2.1.offxml",
    reference_name: str = "Sage 2.2.1",
    target_force_field: str = "../04_benchmark/forcefields/fb-fit-v3-single-mean-k100.offxml",
    target_name: str = "v3-k100",
    image_directory: str = "images/lj-curves/v3-vs-221",
):
    image_directory = pathlib.Path(image_directory)
    image_directory.mkdir(parents=True, exist_ok=True)

    reference_ff = ForceField(reference_force_field)
    reference_vdw_handler = reference_ff.get_parameter_handler("vdW")
    target_ff = ForceField(target_force_field)
    target_vdw_handler = target_ff.get_parameter_handler("vdW")

    for ref_param in tqdm.tqdm(reference_vdw_handler.parameters):
        tgt_param = target_vdw_handler[ref_param.smirks]

        if (
            (ref_param.epsilon == tgt_param.epsilon)
            and (ref_param.rmin_half == tgt_param.rmin_half)
        ):
            logger.info(f"No change for {ref_param.smirks}")
            continue
        
        ref_epsilon = ref_param.epsilon.m_as(unit.kilocalories_per_mole)
        ref_sigma = ref_param.sigma.m_as(unit.angstrom)
        ref_rmin_half = ref_param.rmin_half.m_as(unit.angstrom)
        tgt_epsilon = tgt_param.epsilon.m_as(unit.kilocalories_per_mole)
        tgt_sigma = tgt_param.sigma.m_as(unit.angstrom)
        tgt_rmin_half = tgt_param.rmin_half.m_as(unit.angstrom)

        xmin = ref_sigma - 0.1
        xs = np.linspace(xmin, xmin + 2, 100)

        logger.info(f"Plotting {ref_param.smirks}: {ref_epsilon} -> {tgt_epsilon}, {ref_sigma} -> {tgt_sigma}")

        ref_ys = calculate_lj(xs, ref_sigma, ref_epsilon)
        tgt_ys = calculate_lj(xs, tgt_sigma, tgt_epsilon)

        fig, ax = plt.subplots(figsize=(5, 3))
        # plot curves
        ax.plot(xs, ref_ys, label=reference_name, color=OPENFF_BLUE)
        ax.plot(xs, tgt_ys, label=target_name, color=OPENFF_ORANGE)
        ax.axhline(0, ls="--", color="gray", lw=1)
        # plot rmin_half lines
        ylim = ax.get_ylim()
        ax.vlines(
            ref_rmin_half * 2,
            min(ref_ys), 0,
            ls="--",
            colors=OPENFF_BLUE,
        )
        ax.vlines(
            tgt_rmin_half * 2,
            min(tgt_ys), 0,
            ls="--",
            colors=OPENFF_ORANGE,
        )
        ax.set_xlabel("Distance (Å)")
        ax.set_ylabel("Lennard-Jones Potential (kcal/mol)")
        ax.set_title(f"{tgt_param.id}: {tgt_param.smirks.replace('#7,#8,#9,#16,#17,#35', 'ENA')}")
        ax.legend()
        plt.tight_layout()
        image_file = image_directory / f"{tgt_param.id}.png"
        plt.savefig(image_file, dpi=300)
        plt.close()
        logger.info(f"Saved image to {image_file}")


if __name__ == "__main__":
    main()
