import click
import re
import pathlib

from openff.interchange.drivers.all import get_openmm_energies
from openff.toolkit import ForceField, Molecule, Topology
from openff.toolkit.typing.engines.smirnoff import LibraryChargeHandler, vdWHandler, ConstraintHandler
from openff.units import Quantity, unit


@click.command()
@click.option(
    "--input",
    "-i",
    "input_forcefield",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    required=True,
    help="The path to the input force field file.",
)
@click.option(
    "--output",
    "-o",
    "output_forcefield",
    type=click.Path(exists=False, dir_okay=False, file_okay=True),
    required=True,
    help="The path to the output force field file.",
)
def main(
    input_forcefield: str = "../04_fit-forcefield/nor4/fb-fit/result/optimize/force-field.offxml",
    output_forcefield: str = "openff-2.2.1-rc1.offxml",
):

    # Delete cosmetic attributes
    force_field = ForceField(
        input_forcefield,
        allow_cosmetic_attributes=True
    )
    # remove all cosmetic attributes
    for _, handler in force_field._parameter_handlers.items():
        for parameter in handler.parameters:
            try:
                parameter.delete_cosmetic_attribute('parameterize')
            except AttributeError:
                pass


    # Add the new parameters, from:
    # Tang, K.T., Toennies, J.P. New combining rules for well parameters and shapes
    # of the van der Waals potential of mixed rare gas systems. Z Phys D - Atoms,
    # Molecules and Clusters 1, 91–101 (1986). https://doi.org/10.1007/BF01384663
    # epsilon is converted from Hartrees using 0.894e-3 Ha * 627.5094740631 kcal/mol / Ha = 0.560993470 kcal/mol
    # --> 0.561 kcal/mol preserving the initial number of significant figures
    xe_vdw_param = vdWHandler.vdWType(
        smirks="[#54:1]",
        epsilon=0.561 * unit.kilocalorie / unit.mole,
        sigma=4.363 * unit.angstrom,
        id="n36",
    )

    xe_charge_param = LibraryChargeHandler.LibraryChargeType(
        smirks="[#54:1]", charge1=0.0 * unit.elementary_charge, id="Xe"
    )

    try:
        force_field.get_parameter_handler("vdW").add_parameter(parameter=xe_vdw_param)
    except:
        pass
    try:
        force_field.get_parameter_handler("LibraryCharges").add_parameter(parameter=xe_charge_param)
    except:
        pass
    
    # Save unconstrained version
    output_forcefield = pathlib.Path(output_forcefield)
    stem = output_forcefield.stem
    unconstrained_path = output_forcefield.parent / f"{stem}_unconstrained.offxml"
    force_field.to_file(unconstrained_path)
    print(f"Saved unconstrained force field to {unconstrained_path}")

    # Add constraints
    constraint_handler = force_field.get_parameter_handler('Constraints')
    constraint_param = constraint_handler.ConstraintType(smirks="[#1:1]-[*:2]", id="c1")
    constraint_handler.add_parameter(parameter=constraint_param,before=0)
    force_field.to_file(output_forcefield)
    print(f"Saved constrained force field to {output_forcefield}")

    # Make sure new FF loads with most recent toolkit version
    test_uc = ForceField(unconstrained_path)
    test_c = ForceField(output_forcefield)

    # make sure that Xe is treated correctly
    topology = Topology.from_molecules(
        [
            Molecule.from_smiles("[Xe]"),
            Molecule.from_smiles("[Xe]"),
        ]
    )

    x = 5
    # Place them a few A apart
    topology.set_positions(
        Quantity(
            [
                [0.0, 0.0, 0.0],
                [
                    x,
                    0.0,
                    0.0,
                ],
            ],
            unit.angstrom,
        )
    )

    interchange = test_c.create_interchange(topology)

    get_openmm_energies(interchange)['Nonbonded'].magnitude


if __name__ == "__main__":
    main()
