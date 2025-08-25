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

    # split S-6X3-N
    b52 = bond_handler.get_parameter({"id": "b52"})[0]
    assert b52.smirks == "[#16X2,#16X1-1,#16X3+1:1]-[#6X3:2]"
    b52a = copy.deepcopy(b52)
    b52a.id = "b52a"
    b52a.smirks = "[#16X2,#16X1-1,#16X3+1:1]-[#6X3:2](=[*])-[#7,#8]"
    bond_handler.add_parameter(parameter=b52a, after=b52.smirks)

    # === new ===

    # split b10
    b10 = bond_handler.get_parameter({"id": "b10"})[0]
    assert b10.smirks == "[#6X3:1](=[#8X1+0])-[#7X3:2]"
    b10a = copy.deepcopy(b10)
    b10a.id = "b10a"
    b10a.smirks = "[#6X3:1](=[#8X1+0])-[#7X3R:2]"
    bond_handler.add_parameter(parameter=b10a, after=b10.smirks)

    # split b57 for sulfonamides
    b57c = copy.deepcopy(b57b)
    b57c.id = "b57c"
    b57c.smirks = "[#16X4:1](~[#8X1])(~[#8X1])-[#7:2]"
    bond_handler.add_parameter(parameter=b57c, after=b57b.smirks)

    # split b4
    b4 = bond_handler.get_parameter({"id": "b4"})[0]
    assert b4.smirks == "[#6X3:1]-[#6X3:2]"
    b4a = copy.deepcopy(b4)
    b4a.id = "b4a"
    b4a.smirks = "[#6X3R:1]-[#6X3R:2]"
    bond_handler.add_parameter(parameter=b4a, after=b4.smirks)

    # split b56
    b56 = bond_handler.get_parameter({"id": "b56"})[0]
    assert b56.smirks == "[#16X4,#16X3!+1:1]-[#6:2]"
    b56a = copy.deepcopy(b56)
    b56a.id = "b56a"
    b56a.smirks = "[#8X1]~[#16X4:1](~[#8X1])-[#6:2]"
    bond_handler.add_parameter(parameter=b56a, after=b56.smirks)

    # split b35
    b35 = bond_handler.get_parameter({"id": "b35"})[0]
    assert b35.smirks == "[#7X3:1]-[#7X2:2]"
    b35a = copy.deepcopy(b35)
    b35a.id = "b35a"
    b35a.smirks = "[#7X3R:1]@[#7X2R:2]"
    bond_handler.add_parameter(parameter=b35a, after=b35.smirks)

    # split b42 for nitro
    b42 = bond_handler.get_parameter({"id": "b42"})[0]
    assert b42.smirks == "[#7:1]~[#8X1:2]"
    b42a = copy.deepcopy(b42)
    b42a.id = "b42a"
    b42a.smirks = "[#7+1:1]~[#8X1:2](~[#8X1])"
    bond_handler.add_parameter(parameter=b42a, after=b42.smirks)

    # split b58
    b58 = bond_handler.get_parameter({"id": "b58"})[0]
    assert b58.smirks == "[#16X4,#16X3:1]-[#8X2:2]"
    b58a = copy.deepcopy(b58)
    b58a.id = "b58a"
    b58a.smirks = "[#8]~[#16X4:1](~[#8])(~[#8])-[#8X2:2]"
    bond_handler.add_parameter(parameter=b58a, after=b58.smirks)

    # split b53
    b53 = bond_handler.get_parameter({"id": "b53"})[0]
    assert b53.smirks == "[#16X2,#16X1-1:1]-[#7:2]"
    b53a = copy.deepcopy(b53)
    b53a.id = "b53a"
    b53a.smirks = "[#16X2R:1]@[#7R:2]"
    bond_handler.add_parameter(parameter=b53a, after=b53.smirks)



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

    # add new angle for N in two rings
    # unintuitively, this is in 3 rings because of the super-ring
    a18c = copy.deepcopy(a18a)
    a18c.id = "a18c"
    a18c.smirks = "[r5:1]@-[r5x3;#7X4,#7X3,#7X2-1:2]@-[r:3]"
    angle_handler.add_parameter(parameter=a18c, after=a18a.smirks)


    # add angle for atoms in 5-membered ring with external ring bond
    a18d = copy.deepcopy(a18a)
    a18d.id = "a18d"
    a18d.smirks = "[*:1]@~[r5;#7X3,#6X3:2]-;!@[*:3]"
    angle_handler.add_parameter(parameter=a18d, after=a18b.smirks)

    # add new angle for 6X3 in r3
    a6 = angle_handler.get_parameter({"id": "a6"})[0]
    a4 = angle_handler.get_parameter({"id": "a4"})[0]
    a4a = copy.deepcopy(a4)
    assert a4a.smirks == "[*;r3:1]~;@[*;r3:2]~;!@[*:3]"
    a4a.smirks = "[*;r3:1]~;@[#6X3;r3:2]~;!@[*:3]"
    a4a.id = "a4a"
    # add last
    angle_handler.add_parameter(parameter=a4a, after=a6.smirks)

    # add new N-P angle specifically
    a40 = angle_handler.get_parameter({"id": "a40"})[0]
    assert a40.smirks == "[*:1]~[#15:2]~[*:3]"
    a40a = copy.deepcopy(a40)
    a40a.id = "a40a"
    a40a.smirks = "[#7:1]-[#15X4:2]~[*:3]"
    angle_handler.add_parameter(parameter=a40a, after=a40.smirks)

    # add new P-N angle specifically
    a21 = angle_handler.get_parameter({"id": "a21"})[0]
    assert a21.smirks == "[#1:1]-[#7X3$(*~[#6X3,#6X2,#7X2+0]):2]-[*:3]"
    a21a = copy.deepcopy(a21)
    a21a.id = "a21a"
    a21a.smirks = "[#1:1]-[#7X3$(*~[#6X3,#6X2,#7X2+0]):2]-[#15:3]"
    angle_handler.add_parameter(parameter=a21a, after=a21.smirks)

    # add P-C-P angle specifically
    a1 = angle_handler.get_parameter({"id": "a1"})[0]
    assert a1.smirks == "[*:1]~[#6X4:2]-[*:3]"
    a1a = copy.deepcopy(a1)
    a1a.id = "a1a"
    a1a.smirks = "[#15:1]~[#6X4:2]-[#15:3]"
    angle_handler.add_parameter(parameter=a1a, after=a1.smirks)

    # sulfonamide angles
    a19 = angle_handler.get_parameter({"id": "a19"})[0]
    a19a = copy.deepcopy(a19)
    a19a.id = "a19a"
    a19a.smirks = "[#1:1]-[#7X3:2]-[#16X4:3](~[#8X1])~[#8X1]"
    angle_handler.add_parameter(parameter=a19a, after=a19.smirks)

    a31 = angle_handler.get_parameter({"id": "a31"})[0]
    assert a31.smirks == "[*:1]~[#16X4:2]~[*:3]"
    a31a = copy.deepcopy(a31)
    a31a.id = "a31a"
    a31a.smirks = "[#7X2:1]~[#16X4:2](~[#8X1])(~[#8X1])~[*:3]"
    angle_handler.add_parameter(parameter=a31a, after=a31.smirks)

    # differentiate angle in and not in ring
    a22 = angle_handler.get_parameter({"id": "a22"})[0]
    a22a = copy.deepcopy(a22)
    a22a.id = "a22a"
    a22a.smirks = "[*:1]~[#7X2+0;!R:2]~[*:3]"
    angle_handler.add_parameter(parameter=a22a, after=a22.smirks)

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

    # expand a23
    a23 = angle_handler.get_parameter({"id": "a23"})[0]
    assert a23.smirks == "[*:1]~[#7X2+0:2]~[#6X2:3](~[#16X1])"
    a23.smirks = "[*:1]~[#7X2+0:2]~[#6X2:3](~[#16X1,#8X1])"

    # split a39
    a39 = angle_handler.get_parameter({"id": "a39"})[0]
    assert a39.smirks == "[#6X3:1]-[#16X2:2]-[#1:3]"
    a39a = copy.deepcopy(a39)
    a39a.id = "a39a"
    a39a.smirks = "[#8X1,#16X1]=[#6X3:1]-[#16X2:2]-[#1:3]"
    angle_handler.add_parameter(parameter=a39a, after=a39.smirks)

    # split a26
    a26 = angle_handler.get_parameter({"id": "a26"})[0]
    assert a26.smirks == "[#8X1:1]~[#7X3:2]~[#8X1:3]"
    a26a = copy.deepcopy(a26)
    a26a.id = "a26a"
    a26a.smirks = "[#8X1:1]~[#7X3:2](~[#8,#7])~[#8X1:3]"
    angle_handler.add_parameter(parameter=a26a, after=a26.smirks)

    # === new ===

    # split a10
    a10 = angle_handler.get_parameter({"id": "a10"})[0]
    assert a10.smirks == "[*:1]~[#6X3:2]~[*:3]"
    a10a = copy.deepcopy(a10)
    a10a.id = "a10a"
    a10a.smirks = "[*:1]~[#6X3:2](=[#8])-[#8X2:3]"
    angle_handler.add_parameter(parameter=a10a, after=a10.smirks)

    a10b = copy.deepcopy(a10)
    a10b.id = "a10b"
    a10b.smirks = "[#7:1]-[#6X3:2]=[#8:3]"
    angle_handler.add_parameter(parameter=a10b, after=a10a.smirks)

    # split a21
    a21 = angle_handler.get_parameter({"id": "a21"})[0]
    assert a21.smirks == "[#1:1]-[#7X3$(*~[#6X3,#6X2,#7X2+0]):2]-[*:3]"
    a21a = copy.deepcopy(a21)
    a21a.id = "a21a"
    a21a.smirks = "[#1:1]-[#7X3$(*~[#6X3](~[#7])~[#7]):2]~[#6:3]"
    angle_handler.add_parameter(parameter=a21a, after=a21.smirks)

    a21b = copy.deepcopy(a21)
    a21b.id = "a21b"
    a21b.smirks = "[#1:1]-[#7X3$(*-[#6X3]):2]-[#1:3]"
    angle_handler.add_parameter(parameter=a21b, after=a21a.smirks)

    # split a1 again
    a1 = angle_handler.get_parameter({"id": "a1"})[0]
    assert a1.smirks == "[*:1]~[#6X4:2]-[*:3]"
    a1b = copy.deepcopy(a1)
    a1b.id = "a1b"
    a1b.smirks = "[R:1]@-[#6X4R:2]!@;-[*:3]"
    angle_handler.add_parameter(parameter=a1b, after=a1.smirks)

    # add new angle for 8X1~S-NH2
    a42 = angle_handler.get_parameter({"id": "a42"})[0]
    a43 = copy.deepcopy(a42)
    a43.id = "a43"
    a43.smirks = "[#7X3,#8X2:1]-[#16X4:2]~[#8X1:3](~[#8X1])~[#8]"
    angle_handler.add_parameter(parameter=a43, after=a42.smirks)

    # split a18 again... too much?
    a18a = angle_handler.get_parameter({"id": "a18a"})[0]
    a18e = copy.deepcopy(a18a)
    a18e.id = "a18e"
    a18e.smirks = "[r5;#7:1]@-[r5;#7X3:2]@-[r5;#7:3]"
    angle_handler.add_parameter(parameter=a18e, after=a18a.smirks)

    # split a33 to be safe -- maybe should just modify smirks
    a33 = angle_handler.get_parameter({"id": "a33"})[0]
    assert a33.smirks == "[*:1]~[#16X3$(*~[#8X1,#7X2]):2]~[*:3]"
    a33a = copy.deepcopy(a33)
    a33a.id = "a33a"
    a33a.smirks = "[*:1]~[#16X3:2](~[#8X1,#7X2])~[*:3]"
    angle_handler.add_parameter(parameter=a33a, after=a33.smirks)

    # split a15 again
    a15a = angle_handler.get_parameter({"id": "a15a"})[0]
    assert a15a.smirks == "[#8X1:1]~[#6X3:2]~[#8X1:3]"
    a15b = copy.deepcopy(a15a)
    a15b.id = "a15b"
    a15b.smirks = "[#8X1:1]~[#6X3:2](-[a])~[#8X1:3]"
    angle_handler.add_parameter(parameter=a15b, after=a15a.smirks)

    # split a40 again
    a40a = angle_handler.get_parameter({"id": "a40a"})[0]
    assert a40a.smirks == "[#7:1]-[#15X4:2]~[*:3]"
    a40b = copy.deepcopy(a40a)
    a40b.id = "a40b"
    a40b.smirks = "[#8:1]~[#15X4:2]~[#8,#16:3]"
    angle_handler.add_parameter(parameter=a40b, after=a40a.smirks)

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
    default="output/initial-force-field-v2.offxml",
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
