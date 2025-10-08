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


def main(
    image_directory = "images/lj-curves-combined",
):
    image_directory = pathlib.Path(image_directory)
    image_directory.mkdir(parents=True, exist_ok=True)

    reference_ff = ForceField("openff-2.2.1.offxml")
    reference_vdw_handler = reference_ff.get_parameter_handler("vdW")
    rc1 = ForceField("../04_benchmark/forcefields/openff-2.3.0rc1.offxml")
    rc1_vdw_handler = rc1.get_parameter_handler("vdW")
    v3 = ForceField("../04_benchmark/forcefields/fb-fit-v3-single-mean-k100.offxml")
    v3_vdw_handler = v3.get_parameter_handler("vdW")

    for ref_param in tqdm.tqdm(reference_vdw_handler.parameters):
        v3_param = v3_vdw_handler[ref_param.smirks]
        rc1_param = rc1_vdw_handler[ref_param.smirks]

        if (
            (ref_param.epsilon == v3_param.epsilon)
            and (ref_param.rmin_half == v3_param.rmin_half)
        ):
            logger.info(f"No change for {ref_param.smirks}")
            continue
        
        ref_epsilon = ref_param.epsilon.m_as(unit.kilocalories_per_mole)
        ref_sigma = ref_param.sigma.m_as(unit.angstrom)
        ref_rmin_half = ref_param.rmin_half.m_as(unit.angstrom)
        v3_epsilon = v3_param.epsilon.m_as(unit.kilocalories_per_mole)
        v3_sigma = v3_param.sigma.m_as(unit.angstrom)
        v3_rmin_half = v3_param.rmin_half.m_as(unit.angstrom)
        rc1_epsilon = rc1_param.epsilon.m_as(unit.kilocalories_per_mole)
        rc1_sigma = rc1_param.sigma.m_as(unit.angstrom)
        rc1_rmin_half = rc1_param.rmin_half.m_as(unit.angstrom)

        xmin = ref_sigma - 0.1
        xs = np.linspace(xmin, xmin + 2, 100)

        logger.info(f"Plotting {ref_param.smirks}: {ref_epsilon} -> {v3_epsilon}, {ref_sigma} -> {v3_sigma}")

        ref_ys = calculate_lj(xs, ref_sigma, ref_epsilon)
        v3_ys = calculate_lj(xs, v3_sigma, v3_epsilon)
        rc1_ys = calculate_lj(xs, rc1_sigma, rc1_epsilon)

        fig, ax = plt.subplots(figsize=(5, 3))
        # plot curves
        ax.plot(xs, ref_ys, label="Sage 2.2.1", color="tab:blue")
        ax.plot(xs, rc1_ys, label="Sage 2.3.0-rc1", color="tab:green")
        ax.plot(xs, v3_ys, label="v3-k100", color="tab:orange")
        ax.axhline(0, ls="--", color="gray", lw=1)
        # plot rmin_half lines
        ylim = ax.get_ylim()
        ax.vlines(
            ref_rmin_half * 2,
            min(ref_ys), 0,
            ls="--",
            colors="tab:blue",
        )
        ax.vlines(
            rc1_rmin_half * 2,
            min(rc1_ys), 0,
            ls="--",
            colors="tab:green",
        )
        ax.vlines(
            v3_rmin_half * 2,
            min(v3_ys), 0,
            ls="--",
            colors="tab:orange",
        )
        ax.set_xlabel("Distance (Å)")
        ax.set_ylabel("Lennard-Jones Potential (kcal/mol)")
        ax.set_title(f"{ref_param.id}: {ref_param.smirks.replace('#7,#8,#9,#16,#17,#35', 'ENA')}")
        ax.legend()
        plt.tight_layout()
        image_file = image_directory / f"{ref_param.id}.png"
        plt.savefig(image_file, dpi=300)
        plt.close()
        logger.info(f"Saved image to {image_file}")


if __name__ == "__main__":
    main()
