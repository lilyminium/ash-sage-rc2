"""
This script generates a new force field for a valence fit by:
* removing the Constraints section from the OpenFF force field
* splitting torsion multiplicities into separate parameters
* modifying some torsions suggested by Cresset, and adding a new biaryl torsion
* discarding any cosmetic attributes from training vdW

No other changes are made.
"""

import click
from helper_functions import (
    remove_constraints,
    split_torsion_multiplicities,
    modify_cresset_torsions,
    print_number_parameters
)


@click.command()
@click.option(
    "--output",
    "-o",
    "output_path",
    default="output/initial-force-field-v1.offxml",
    type=click.Path(exists=False, dir_okay=False, file_okay=True),
    required=True,
    help="The path to the output force field file.",
)
@click.option(
    "--input",
    "-i",
    "input_path",
    type=str,
    default="openff_unconstrained-2.2.1.offxml",
    help="The path of the force field to download.",
    show_default=True,
)
def generate_force_field(
    output_path: str,
    input_path: str = "openff_unconstrained-2.2.1.offxml",
):
    from openff.toolkit import ForceField

    force_field = ForceField(input_path, allow_cosmetic_attributes=True)

    remove_constraints(force_field)

    torsion_handler = force_field.get_parameter_handler("ProperTorsions")

    # torsion multiplicity changes -- Brent
    split_torsion_multiplicities(torsion_handler)

    modify_cresset_torsions(torsion_handler)

    print_number_parameters(force_field)

    # Write out file
    force_field.to_file(output_path, discard_cosmetic_attributes=True)
    print(f"Force field written to {output_path}")


if __name__ == "__main__":
    generate_force_field()
