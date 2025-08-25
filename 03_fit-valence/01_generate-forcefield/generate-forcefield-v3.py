"""
This script generates a new force field for a valence fit by:
* removing the Constraints section from the OpenFF force field
* splitting torsion multiplicities into separate parameters
* modifying some torsions suggested by Cresset, and adding a new biaryl torsion
* discarding any cosmetic attributes from training vdW
* fixing some problematic bonds and angles
"""

import copy
import click
from helper_functions import (
    remove_constraints,
    split_torsion_multiplicities,
    modify_cresset_torsions,
    print_number_parameters
)

def modify_r3_angles(angle_handler):
    # move r3 angles to the end
    previous = angle_handler.parameters[-1]
    r3_angles = [
        angle_handler.get_parameter({"id": param})[0]
        for param in ["a4", "a5", "a6"]
    ]
    for param in r3_angles:
        angle_handler.parameters.remove(param)
        angle_handler.add_parameter(parameter=param, after=previous.smirks)
        previous = param


def split_msm_angles(angle_handler):
    # add a non-C r4 internal angle
    a7 = angle_handler.get_parameter({"id": "a7"})[0]
    a7a = copy.deepcopy(a7)
    a7a.id = "a7a"
    a7a.smirks = "[r4:1]1-;@[r4:2]-;@[r4:3]~[*]~1"
    angle_handler.add_parameter(parameter=a7a)

    # r4 internal angle applies to external as well
    assert a7.smirks == "[#6r4:1]-;@[#6r4:2]-;@[#6r4:3]"
    a7.smirks = "[#6r4:1]1-;@[#6r4:2]-;@[#6r4:3]~[*]~1"
    # move to end
    angle_handler._parameters.remove(a7)
    angle_handler.add_parameter(
        parameter=a7,
        after=angle_handler.parameters[-1].smirks
    )

    # split a28 into H-O-X and X-O-X
    a28 = angle_handler.get_parameter({"id": "a28"})[0]
    assert a28.smirks == "[*:1]-[#8:2]-[*:3]"
    a28a = copy.deepcopy(a28)
    a28a.id = "a28a"
    a28a.smirks = "[#1:1]-[#8:2]-[*:3]"
    angle_handler.add_parameter(parameter=a28a, after=a28.smirks)

    # add O-S-O angle
    a32 = angle_handler.get_parameter({"id": "a32"})[0]
    a42 = copy.deepcopy(a32)
    a42.id = "a42"
    a42.smirks = "[#8X1:1]~[#16X4:2](~[#8X1,#7X2])~[#8X1:3]"
    angle_handler.add_parameter(parameter=a42, after=a32.smirks)

    # add bond to r3, r5 angles
    a3 = angle_handler.get_parameter({"id": "a3"})[0]
    assert a3.smirks == "[*;r3:1]1~;@[*;r3:2]~;@[*;r3:3]1"
    a3.smirks = "[*;r3:1]1~;@[*;r3:2]~;@[*;r3:3]~1"
    a41 = angle_handler.get_parameter({"id": "a41"})[0]
    assert a41.smirks == "[*;r5:1]1@[*;r5:2]@[*;r5:3]@[*;r5]@[*;r5]1"
    a41.smirks = "[*;r5:1]1@[*;r5:2]@[*;r5:3]@[*;r5]@[*;r5]~1"
    a41a = angle_handler.get_parameter({"id": "a41a"})[0]
    assert a41a.smirks == "[*;r5:1]1@[#16;r5:2]@[*;r5:3]@[*;r5]@[*;r5]1"
    a41a.smirks = "[*;r5:1]1@[#16;r5:2]@[*;r5:3]@[*;r5]@[*;r5]~1"


def split_msm_bonds(bond_handler):
    # split b8 into separate parameter for nitro groups
    b8 = bond_handler.get_parameter({"id": "b8"})[0]
    assert b8.smirks == "[#6X3:1]-[#7X3:2]"
    b8a = copy.deepcopy(b8)
    b8a.id = "b8a"
    b8a.smirks = "[#6X3:1]-[#7X3:2](~[#8X1])~[#8X1]"
    bond_handler.add_parameter(parameter=b8a, after=b8.smirks)


def fix_problematic_bonds(bond_handler):
    b57 = bond_handler.get_parameter({"id": "b57"})[0]
    assert b57.smirks == "[#16X4,#16X3:1]~[#7:2]"
    b57b = copy.deepcopy(b57)
    b57b.id = "b57b"
    b57b.smirks = "[#16X4,#16X3:1]~[#7+1:2]"
    bond_handler.add_parameter(parameter=b57b, after=b57.smirks)

    # split N-P-N bond
    b62 = bond_handler.get_parameter({"id": "b62"})[0]
    assert b62.smirks == "[#15:1]-[#7:2]"
    b62a = copy.deepcopy(b62)
    b62a.id = "b62a"
    b62a.smirks = "[#7,#8]-[#15X4:1]-[#7:2]"
    bond_handler.add_parameter(parameter=b62a, after=b62.smirks)


def fix_problematic_angles(angle_handler):
    # add carboxylate angle specifically
    a15 = angle_handler.get_parameter({"id": "a15"})[0]
    assert a15.smirks == "[#8X1:1]~[#6X3:2]~[#8:3]"
    a15a = copy.deepcopy(a15)
    a15a.id = "a15a"
    a15a.smirks = "[#8X1:1]~[#6X3:2]~[#8X1:3]"
    angle_handler.add_parameter(parameter=a15a, after=a15.smirks)

    # split out N-
    a18 = angle_handler.get_parameter({"id": "a18"})[0]
    a18b = copy.deepcopy(a18)
    a18b.id = "a18b"
    a18b.smirks = "[*:1]~[#7X2-1:2]~[*:3]"
    angle_handler.add_parameter(parameter=a18b, after=a18.smirks)

    # fix oddly inclusive ring exclusion
    a18a = angle_handler.get_parameter({"id": "a18a"})[0]
    assert a18a.smirks == "[*:1]@-[r!r6;#7X4,#7X3,#7X2-1:2]@-[*:3]"
    a18a.smirks = "[*:1]@-[r5;#7X4,#7X3,#7X2-1:2]@-[*:3]"

    # move it down near the end
    angle_handler._parameters.remove(a18a)
    a41a = angle_handler.get_parameter({"id": "a41a"})[0]
    angle_handler.add_parameter(
        parameter=a18a,
        after=a41a.smirks
    )

    # move a13 near end too
    a13 = angle_handler.get_parameter({"id": "a13"})[0]
    assert a13.smirks == "[*;r6:1]~;@[*;r5:2]~;@[*;r5;x2:3]"
    angle_handler._parameters.remove(a13)
    angle_handler.add_parameter(
        parameter=a13,
        after=a18a.smirks
    )

    # add new angle for 6X3 in r3
    a6 = angle_handler.get_parameter({"id": "a6"})[0]
    a4 = angle_handler.get_parameter({"id": "a4"})[0]
    a4a = copy.deepcopy(a4)
    assert a4a.smirks == "[*;r3:1]~;@[*;r3:2]~;!@[*:3]"
    a4a.smirks = "[*;r3:1]~;@[#6X3;r3:2]~;!@[*:3]"
    a4a.id = "a4a"
    # add last
    angle_handler.add_parameter(parameter=a4a, after=a6.smirks)

    # add P-C-P angle specifically
    a1 = angle_handler.get_parameter({"id": "a1"})[0]
    assert a1.smirks == "[*:1]~[#6X4:2]-[*:3]"
    a1a = copy.deepcopy(a1)
    a1a.id = "a1a"
    a1a.smirks = "[#15:1]~[#6X4:2]-[#15:3]"
    angle_handler.add_parameter(parameter=a1a, after=a1.smirks)

    # N=S~O differentiation
    a24 = angle_handler.get_parameter({"id": "a24"})[0]
    a24a = copy.deepcopy(a24)
    a24a.id = "a24a"
    a24a.smirks = "[#1:1]-[#7X2:2]~[#16X4:3]"
    angle_handler.add_parameter(parameter=a24a, after=a24.smirks)

    # split out amide
    a11 = angle_handler.get_parameter({"id": "a11"})[0]
    a11a = copy.deepcopy(a11)
    a11a.id = "a11a"
    a11a.smirks = "[#1:1]-[#6X3:2](=[#8])-[#7:3]"
    angle_handler.add_parameter(parameter=a11a, after=a11.smirks)

    # split a10
    a10 = angle_handler.get_parameter({"id": "a10"})[0]
    assert a10.smirks == "[*:1]~[#6X3:2]~[*:3]"
    a10a = copy.deepcopy(a10)
    a10a.id = "a10a"
    a10a.smirks = "[*:1]~[#6X3:2](=[#8])-[#8X2:3]"
    angle_handler.add_parameter(parameter=a10a, after=a10.smirks)

    # split a33 to be safe -- maybe should just modify smirks
    a33 = angle_handler.get_parameter({"id": "a33"})[0]
    assert a33.smirks == "[*:1]~[#16X3$(*~[#8X1,#7X2]):2]~[*:3]"
    a33a = copy.deepcopy(a33)
    a33a.id = "a33a"
    a33a.smirks = "[*:1]~[#16X3:2](~[#8X1,#7X2])~[*:3]"
    angle_handler.add_parameter(parameter=a33a, after=a33.smirks)

def move_r4_angles_to_end(angle_handler):
    final_order = [
        # "a13", "a13a", "a14", #r5
        "a7", "a8", "a9", #r4
    ]
    previous = angle_handler.get_parameter({"id": "a41a"})[0]
    angle_parameters = [
        angle_handler.get_parameter({"id": param})[0]
        for param in final_order
    ]
    for param in angle_parameters:
        angle_handler.parameters.remove(param)
        angle_handler.add_parameter(parameter=param, after=previous.smirks)
        previous = param




@click.command()
@click.option(
    "--output",
    "-o",
    "output_path",
    default="output/initial-force-field-v3.offxml",
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

    # angles
    angle_handler = force_field.get_parameter_handler("Angles")
    modify_r3_angles(angle_handler)

    split_msm_angles(angle_handler)

    bond_handler = force_field.get_parameter_handler("Bonds")
    split_msm_bonds(bond_handler)

    fix_problematic_angles(angle_handler)
    move_r4_angles_to_end(angle_handler)
    fix_problematic_bonds(bond_handler)

    print_number_parameters(force_field)


    # Write out file
    force_field.to_file(output_path, discard_cosmetic_attributes=True)
    print(f"Force field written to {output_path}")


if __name__ == "__main__":
    generate_force_field()
