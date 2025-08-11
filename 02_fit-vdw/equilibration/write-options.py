"""
This script generates a `options.json` file with equilibration options for the OpenFF Evaluator.

Most options are hardcoded, but can be modified as needed.
Previous ash-sage-rc1 script only allowed up to 12 ns equilibration,
but this version allows up to 200 ns for some slowly-equilibrating systems.
"""

import click

from openff.evaluator.properties import Density, EnthalpyOfMixing
from openff.evaluator.client import RequestOptions

from openff.evaluator.client import RequestOptions
from openff.evaluator.layers.equilibration import EquilibrationProperty
from openff.evaluator.utils.observables import ObservableType


@click.command()
@click.option(
    "--n-molecules",
    "-n",
    type=int,
    default=1000,
    help=(
        "Number of molecules per box."
    )
)
def main(
    n_molecules: int = 1000, 
):

    potential_energy = EquilibrationProperty()
    potential_energy.relative_tolerance = 0.05
    potential_energy.observable_type = ObservableType.PotentialEnergy
    potential_energy.n_uncorrelated_samples = 300

    density = EquilibrationProperty()
    density.relative_tolerance = 0.05
    density.observable_type = ObservableType.Density
    density.n_uncorrelated_samples = 300

    options = RequestOptions()
    options.calculation_layers = ["EquilibrationLayer"]
    density_schema = Density.default_equilibration_schema(
        n_molecules=n_molecules,
        error_tolerances=[potential_energy, density],
        # every iteration is 200 ps
        max_iterations=1000, # go up to 200 ns
        error_on_failure=False,
    )

    dhmix_schema = EnthalpyOfMixing.default_equilibration_schema(
        n_molecules=n_molecules,
        error_tolerances=[potential_energy, density],
        max_iterations=1000,
        error_on_failure=False,
    )

    # note: output frequency is every 10 ps.
    options.add_schema("EquilibrationLayer", "Density", density_schema)
    options.add_schema("EquilibrationLayer", "EnthalpyOfMixing", dhmix_schema)
    options.json("options.json", format=True) 

if __name__ == "__main__":
    main()

