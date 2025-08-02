import click

def remove_constraints(force_field):
    constraint_handler = force_field.get_parameter_handler("Constraints")
    parameters = constraint_handler.get_parameter({"id": "c1"})
    if parameters:
        constraint_handler._parameters.remove(parameters[0])

@click.command()
@click.option(
    "--output",
    "-o",
    "output_path",
    default="output/initial-force-field-v0.offxml",
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

    # Write out file
    force_field.to_file(output_path, discard_cosmetic_attributes=True)
    print(f"Force field written to {output_path}")


if __name__ == "__main__":
    generate_force_field()
