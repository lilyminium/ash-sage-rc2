import copy
import click

def remove_constraints(force_field):
    constraint_handler = force_field.get_parameter_handler("Constraints")
    parameters = constraint_handler.get_parameter({"id": "c1"})
    if parameters:
        constraint_handler._parameters.remove(parameters[0])


def split_torsion_multiplicities(torsion_handler):

    # t58
    t58 = torsion_handler.get_parameter({"id": "t58"})[0]
    assert t58.smirks == "[#6X4:1]-[#6X4:2]-[#7X4,#7X3:3]-[#6X4:4]"
    t58.smirks = "[#6X4:1]-[#6X4:2]-[#7X4:3]-[#6X4:4]"

    t58a = copy.deepcopy(t58)
    t58a.id = "t58a"
    t58a.smirks = "[#6X4:1]-[#6X4:2]-[#7X3:3]-[#6X4:4]"
    torsion_handler.add_parameter(parameter=t58a, after=t58.smirks)

    # t59
    t59 = torsion_handler.get_parameter({"id": "t59"})[0]
    assert t59.smirks == "[#1:1]-[#7X4,#7X3:2]-[#6X4;r3:3]-[*:4]"
    t59.smirks = "[#1:1]-[#7X4:2]-[#6X4;r3:3]-[*:4]"

    t59a = copy.deepcopy(t59)
    t59a.id = "t59a"
    t59a.smirks = "[#1:1]-[#7X3:2]-[#6X4;r3:3]-[*:4]"
    torsion_handler.add_parameter(parameter=t59a, after=t59.smirks)

    # t60
    t60 = torsion_handler.get_parameter({"id": "t60"})[0]
    assert t60.smirks == "[#1:1]-[#7X4,#7X3:2]-[#6X4;r3:3]-[#6X4;r3:4]"
    t60.smirks = "[#1:1]-[#7X4:2]-[#6X4;r3:3]-[#6X4;r3:4]"

    t60a = copy.deepcopy(t60)
    t60a.id = "t60a"
    t60a.smirks = "[#1:1]-[#7X3:2]-[#6X4;r3:3]-[#6X4;r3:4]"
    torsion_handler.add_parameter(parameter=t60a, after=t60.smirks)

    # t61
    t61 = torsion_handler.get_parameter({"id": "t61"})[0]
    assert t61.smirks == "[!#1:1]-[#7X4,#7X3:2]-[#6X4;r3:3]-[*:4]"
    t61.smirks = "[!#1:1]-[#7X4:2]-[#6X4;r3:3]-[*:4]"

    t61a = copy.deepcopy(t61)
    t61a.id = "t61a"
    t61a.smirks = "[!#1:1]-[#7X3:2]-[#6X4;r3:3]-[*:4]"
    torsion_handler.add_parameter(parameter=t61a, after=t61.smirks)

    # t62
    t62 = torsion_handler.get_parameter({"id": "t62"})[0]
    assert t62.smirks == "[!#1:1]-[#7X4,#7X3:2]-[#6X4;r3:3]-[#6X4;r3:4]"
    t62.smirks = "[!#1:1]-[#7X4:2]-[#6X4;r3:3]-[#6X4;r3:4]"

    t62a = copy.deepcopy(t62)
    t62a.id = "t62a"
    t62a.smirks = "[!#1:1]-[#7X3:2]-[#6X4;r3:3]-[#6X4;r3:4]"
    torsion_handler.add_parameter(parameter=t62a, after=t62.smirks)

    # t71
    t71 = torsion_handler.get_parameter({"id": "t71"})[0]
    assert t71.smirks == "[#6X3:1]=[#7X2,#7X3+1:2]-[#6X4:3]-[#1:4]"
    t71.smirks = "[#6X3:1]=[#7X2:2]-[#6X4:3]-[#1:4]"

    t71a = copy.deepcopy(t71)
    t71a.id = "t71a"
    t71a.smirks = "[#6X3:1]=[#7X3+1:2]-[#6X4:3]-[#1:4]"
    torsion_handler.add_parameter(parameter=t71a, after=t71.smirks)


    # t72
    t72 = torsion_handler.get_parameter({"id": "t72"})[0]
    assert t72.smirks == "[#6X3:1]=[#7X2,#7X3+1:2]-[#6X4:3]-[#6X3,#6X4:4]"
    t72.smirks = "[#6X3:1]=[#7X2:2]-[#6X4:3]-[#6X3,#6X4:4]"

    t72a = copy.deepcopy(t72)
    t72a.id = "t72a"
    t72a.smirks = "[#6X3:1]=[#7X3+1:2]-[#6X4:3]-[#6X3,#6X4:4]"
    torsion_handler.add_parameter(parameter=t72a, after=t72.smirks)

    # t73 split into t73 and t73a
    t73 = torsion_handler.get_parameter({"id": "t73"})[0]
    assert t73.smirks == "[*:1]~[#7X3,#7X2-1:2]-[#6X3:3]~[*:4]"
    t73.smirks = "[*:1]~[#7X3:2]-[#6X3:3]~[*:4]"

    t73a = copy.deepcopy(t73)
    t73a.id = "t73a"
    t73a.smirks = "[*:1]~[#7X2-1:2]-[#6X3:3]~[*:4]"
    torsion_handler.add_parameter(parameter=t73a, after=t73.smirks)

    # t74 split into t74 and t74a
    t74 = torsion_handler.get_parameter({"id": "t74"})[0]
    assert t74.smirks == "[*:1]~[#7X3,#7X2-1:2]-!@[#6X3:3]~[*:4]"
    t74.smirks = "[*:1]~[#7X3:2]-!@[#6X3:3]~[*:4]"

    t74a = copy.deepcopy(t74)
    t74a.id = "t74a"
    t74a.smirks = "[*:1]~[#7X2-1:2]-!@[#6X3:3]~[*:4]"
    torsion_handler.add_parameter(parameter=t74a, after=t74.smirks)

    # t82 split into t82 and t82b (t82a already exists)
    t82 = torsion_handler.get_parameter({"id": "t82"})[0]
    assert t82.smirks == "[*:1]=[#7X2,#7X3+1:2]-[#6X3:3]-[*:4]"
    t82.smirks = "[*:1]=[#7X2:2]-[#6X3:3]-[*:4]"

    t82b = copy.deepcopy(t82)
    t82b.id = "t82b"
    t82b.smirks = "[*:1]=[#7X3+1:2]-[#6X3:3]-[*:4]"
    torsion_handler.add_parameter(parameter=t82b, after=t82.smirks)
    
    # t83 split into t83 and t83ax
    t83 = torsion_handler.get_parameter({"id": "t83"})[0]
    assert t83.smirks == "[*:1]=[#7X2,#7X3+1:2]-[#6X3:3]=,:[*:4]"
    t83.smirks = "[*:1]=[#7X2:2]-[#6X3:3]=,:[*:4]"

    t83b = copy.deepcopy(t83)
    t83b.id = "t83b"
    t83b.smirks = "[*:1]=[#7X3+1:2]-[#6X3:3]=,:[*:4]"
    torsion_handler.add_parameter(parameter=t83b, after=t83.smirks)

    # -- new addition from Brent's
    # t84 split into t84 and t84a
    t84 = torsion_handler.get_parameter({"id": "t84"})[0]
    assert t84.smirks == "[*:1]~[#7X2,#7X3$(*~[#8X1]):2]:[#6X3:3]~[*:4]"
    t84.smirks = "[*:1]~[#7X2:2]:[#6X3:3]~[*:4]"

    t84b = copy.deepcopy(t84)
    t84b.id = "t84b"
    t84b.smirks = "[*:1]~[#7X3$(*~[#8X1]):2]:[#6X3:3]~[*:4]"
    torsion_handler.add_parameter(parameter=t84b, after=t84.smirks)

    # t115
    t115 = torsion_handler.get_parameter({"id": "t115"})[0]
    assert t115.smirks == "[*:1]-[#16X2,#16X3+1:2]-[#6:3]~[*:4]"
    t115.smirks = "[*:1]-[#16X2:2]-[#6:3]~[*:4]"

    t115a = copy.deepcopy(t115)
    t115a.id = "t115a"
    t115a.smirks = "[*:1]-[#16X3+1:2]-[#6:3]~[*:4]"
    torsion_handler.add_parameter(parameter=t115a, after=t115.smirks)

    # t116
    t116 = torsion_handler.get_parameter({"id": "t116"})[0]
    assert t116.smirks == "[*:1]-[#16X2,#16X3+1:2]-[#6:3]-[#1:4]"
    t116.smirks = "[*:1]-[#16X2:2]-[#6:3]-[#1:4]"

    t116a = copy.deepcopy(t116)
    t116a.id = "t116a"
    t116a.smirks = "[*:1]-[#16X3+1:2]-[#6:3]-[#1:4]"
    torsion_handler.add_parameter(parameter=t116a, after=t116.smirks)

    # t117
    t117 = torsion_handler.get_parameter({"id": "t117"})[0]
    assert t117.smirks == "[#6X3:1]-@[#16X2,#16X1-1,#16X3+1:2]-@[#6X3,#7X2;r5:3]=@[#6,#7;r5:4]"
    t117.smirks = "[#6X3:1]-@[#16X3+1:2]-@[#7X2;r5:3]=@[#6,#7;r5:4]"

    t117a = copy.deepcopy(t117)
    t117a.id = "t117a"
    t117a.smirks = "[#6X3:1]-@[#16X2:2]-@[#6X3;r5:3]=@[#6,#7;r5:4]"
    torsion_handler.add_parameter(parameter=t117a, after=t117.smirks)

    t117b = copy.deepcopy(t117)
    t117b.id = "t117b"
    t117b.smirks = "[#6X3:1]-@[#16X1-1:2]-@[#6X3;r5:3]=@[#6,#7;r5:4]"

    torsion_handler.add_parameter(parameter=t117b, after=t117a.smirks)

    t117c = copy.deepcopy(t117)
    t117c.id = "t117c"
    t117c.smirks = "[#6X3:1]-@[#16X3+1:2]-@[#6X3;r5:3]=@[#6,#7;r5:4]"
    torsion_handler.add_parameter(parameter=t117c, after=t117b.smirks)

    t117d = copy.deepcopy(t117)
    t117d.id = "t117d"
    t117d.smirks = "[#6X3:1]-@[#16X2:2]-@[#7X2;r5:3]=@[#6,#7;r5:4]"
    torsion_handler.add_parameter(parameter=t117d, after=t117c.smirks)

    t117e = copy.deepcopy(t117)
    t117e.id = "t117e"
    t117e.smirks = "[#6X3:1]-@[#16X1-1:2]-@[#7X2;r5:3]=@[#6,#7;r5:4]"
    torsion_handler.add_parameter(parameter=t117e, after=t117d.smirks)

    # t118
    t118 = torsion_handler.get_parameter({"id": "t118"})[0]
    assert t118.smirks == "[*:1]~[#16X4,#16X3!+1:2]-[#6X4:3]-[*:4]"
    t118.smirks = "[*:1]~[#16X4:2]-[#6X4:3]-[*:4]"

    t118a = copy.deepcopy(t118)
    t118a.id = "t118a"
    t118a.smirks = "[*:1]~[#16X3!+1:2]-[#6X4:3]-[*:4]"
    torsion_handler.add_parameter(parameter=t118a, after=t118.smirks)

    # t119
    t119 = torsion_handler.get_parameter({"id": "t119"})[0]
    assert t119.smirks == "[#6X4:1]-[#16X4,#16X3+0:2]-[#6X4:3]-[#1:4]"
    t119.smirks = "[#6X4:1]-[#16X4:2]-[#6X4:3]-[#1:4]"

    t119a = copy.deepcopy(t119)
    t119a.id = "t119a"
    t119a.smirks = "[#6X4:1]-[#16X3+0:2]-[#6X4:3]-[#1:4]"
    torsion_handler.add_parameter(parameter=t119a, after=t119.smirks)

    # t120
    t120 = torsion_handler.get_parameter({"id": "t120"})[0]
    assert t120.smirks == "[#6X4:1]-[#16X4,#16X3+0:2]-[#6X4:3]~[#6X4:4]"
    t120.smirks = "[#6X4:1]-[#16X4:2]-[#6X4:3]~[#6X4:4]"

    t120a = copy.deepcopy(t120)
    t120a.id = "t120a"
    t120a.smirks = "[#6X4:1]-[#16X3+0:2]-[#6X4:3]~[#6X4:4]"
    torsion_handler.add_parameter(parameter=t120a, after=t120.smirks)

    # t121
    t121 = torsion_handler.get_parameter({"id": "t121"})[0]
    assert t121.smirks == "[*:1]~[#16X4,#16X3+0:2]-[#6X3:3]~[*:4]"
    t121.smirks = "[*:1]~[#16X4:2]-[#6X3:3]~[*:4]"

    t121a = copy.deepcopy(t121)
    t121a.id = "t121a"
    t121a.smirks = "[*:1]~[#16X3+0:2]-[#6X3:3]~[*:4]"
    torsion_handler.add_parameter(parameter=t121a, after=t121.smirks)

    # t122
    t122 = torsion_handler.get_parameter({"id": "t122"})[0]
    assert t122.smirks == "[#6:1]-[#16X4,#16X3+0:2]-[#6X3:3]~[*:4]"
    t122.smirks = "[#6:1]-[#16X4:2]-[#6X3:3]~[*:4]"

    t122a = copy.deepcopy(t122)
    t122a.id = "t122a"
    t122a.smirks = "[#6:1]-[#16X3+0:2]-[#6X3:3]~[*:4]"
    torsion_handler.add_parameter(parameter=t122a, after=t122.smirks)

    # t130
    t130 = torsion_handler.get_parameter({"id": "t130"})[0]
    assert t130.smirks == "[*:1]~[#7X4,#7X3:2]-[#7X4,#7X3:3]~[*:4]"
    t130.smirks = "[*:1]~[#7X4:2]-[#7X4:3]~[*:4]"

    t130a = copy.deepcopy(t130)
    t130a.id = "t130a"
    t130a.smirks = "[*:1]~[#7X3:2]-[#7X4:3]~[*:4]"
    torsion_handler.add_parameter(parameter=t130a, after=t130.smirks)

    t130b = copy.deepcopy(t130)
    t130b.id = "t130b"
    t130b.smirks = "[*:1]~[#7X3:2]-[#7X3:3]~[*:4]"
    torsion_handler.add_parameter(parameter=t130b, after=t130a.smirks)

    t130c = copy.deepcopy(t130)
    t130c.id = "t130c"
    t130c.smirks = "[*:1]~[#7X4:2]-[#7X3:3]~[*:4]"
    torsion_handler.add_parameter(parameter=t130c, after=t130b.smirks)

    # t131
    t131 = torsion_handler.get_parameter({"id": "t131"})[0]
    assert t131.smirks == "[#1:1]-[#7X4,#7X3:2]-[#7X4,#7X3:3]-[#1:4]"
    t131.smirks = "[#1:1]-[#7X4:2]-[#7X4:3]-[#1:4]"

    t131a = copy.deepcopy(t131)
    t131a.id = "t131a"
    t131a.smirks = "[#1:1]-[#7X3:2]-[#7X4:3]-[#1:4]"
    torsion_handler.add_parameter(parameter=t131a, after=t131.smirks)

    t131b = copy.deepcopy(t131)
    t131b.id = "t131b"
    t131b.smirks = "[#1:1]-[#7X3:2]-[#7X3:3]-[#1:4]"
    torsion_handler.add_parameter(parameter=t131b, after=t131a.smirks)

    t131c = copy.deepcopy(t131)
    t131c.id = "t131c"
    t131c.smirks = "[#1:1]-[#7X4:2]-[#7X3:3]-[#1:4]"
    torsion_handler.add_parameter(parameter=t131c, after=t131b.smirks)

    # t132
    t132 = torsion_handler.get_parameter({"id": "t132"})[0]
    assert t132.smirks == "[#6X4:1]-[#7X4,#7X3:2]-[#7X4,#7X3:3]-[#1:4]"
    t132.smirks = "[#6X4:1]-[#7X4:2]-[#7X4:3]-[#1:4]"

    t132a = copy.deepcopy(t132)
    t132a.id = "t132a"
    t132a.smirks = "[#6X4:1]-[#7X3:2]-[#7X4:3]-[#1:4]"
    torsion_handler.add_parameter(parameter=t132a, after=t132.smirks)

    t132b = copy.deepcopy(t132)
    t132b.id = "t132b"
    t132b.smirks = "[#6X4:1]-[#7X3:2]-[#7X3:3]-[#1:4]"
    torsion_handler.add_parameter(parameter=t132b, after=t132a.smirks)

    t132c = copy.deepcopy(t132)
    t132c.id = "t132c"
    t132c.smirks = "[#6X4:1]-[#7X4:2]-[#7X3:3]-[#1:4]"
    torsion_handler.add_parameter(parameter=t132c, after=t132b.smirks)

    # t133
    t133 = torsion_handler.get_parameter({"id": "t133"})[0]
    assert t133.smirks == "[#6X4:1]-[#7X4,#7X3:2]-[#7X4,#7X3:3]-[#6X4:4]"
    t133.smirks = "[#6X4:1]-[#7X4:2]-[#7X4:3]-[#6X4:4]"

    t133a = copy.deepcopy(t133)
    t133a.id = "t133a"
    t133a.smirks = "[#6X4:1]-[#7X3:2]-[#7X4:3]-[#6X4:4]"
    torsion_handler.add_parameter(parameter=t133a, after=t133.smirks)

    t133b = copy.deepcopy(t133)
    t133b.id = "t133b"
    t133b.smirks = "[#6X4:1]-[#7X3:2]-[#7X3:3]-[#6X4:4]"
    torsion_handler.add_parameter(parameter=t133b, after=t133a.smirks)

    t133c = copy.deepcopy(t133)
    t133c.id = "t133c"
    t133c.smirks = "[#6X4:1]-[#7X4:2]-[#7X3:3]-[#6X4:4]"
    torsion_handler.add_parameter(parameter=t133c, after=t133b.smirks)

    # t134
    t134 = torsion_handler.get_parameter({"id": "t134"})[0]
    assert t134.smirks == "[*:1]-[#7X4,#7X3:2]-[#7X3$(*~[#6X3,#6X2]):3]~[*:4]"
    t134.smirks = "[*:1]-[#7X4:2]-[#7X3$(*~[#6X3,#6X2]):3]~[*:4]"

    t134a = copy.deepcopy(t134)
    t134a.id = "t134a"
    t134a.smirks = "[*:1]-[#7X3:2]-[#7X3$(*~[#6X3,#6X2]):3]~[*:4]"
    torsion_handler.add_parameter(parameter=t134a, after=t134.smirks)

    # 142
    t142 = torsion_handler.get_parameter({"id": "t142"})[0]
    assert t142.smirks == "[*:1]-[#16X2,#16X3+1:2]-[!#6:3]~[*:4]"
    t142.smirks = "[*:1]-[#16X2:2]-[!#6;X2:3]~[*:4]"

    t142a = copy.deepcopy(t142)
    t142a.id = "t142a"
    t142a.smirks = "[*:1]-[#16X3+1:2]-[!#6;X2:3]~[*:4]"
    torsion_handler.add_parameter(parameter=t142a, after=t142.smirks)

    t142b = copy.deepcopy(t142)
    t142b.id = "t142b"
    t142b.smirks = "[*:1]-[#16X2:2]-[!#6;X3:3]~[*:4]"
    torsion_handler.add_parameter(parameter=t142b, after=t142a.smirks)

    t142c = copy.deepcopy(t142)
    t142c.id = "t142c"
    t142c.smirks = "[*:1]-[#16X3+1:2]-[!#6;X3:3]~[*:4]"
    torsion_handler.add_parameter(parameter=t142c, after=t142b.smirks)

    t142d = copy.deepcopy(t142)
    t142d.id = "t142d"
    t142d.smirks = "[*:1]-[#16X2:2]-[!#6;X4:3]~[*:4]"
    torsion_handler.add_parameter(parameter=t142d, after=t142c.smirks)

    t142e = copy.deepcopy(t142)
    t142e.id = "t142e"
    t142e.smirks = "[*:1]-[#16X3+1:2]-[!#6;X4:3]~[*:4]"
    torsion_handler.add_parameter(parameter=t142e, after=t142d.smirks)

    # t143
    t143 = torsion_handler.get_parameter({"id": "t143"})[0]
    assert t143.smirks == "[*:1]~[#16X4,#16X3+0:2]-[#7:3]~[*:4]"
    t143.smirks = "[*:1]~[#16X4:2]-[#7X2:3]~[*:4]"

    t143a = copy.deepcopy(t143)
    t143a.id = "t143a"
    t143a.smirks = "[*:1]~[#16X3+0:2]-[#7X2:3]~[*:4]"
    torsion_handler.add_parameter(parameter=t143a, after=t143.smirks)

    t143b = copy.deepcopy(t143)
    t143b.id = "t143b"
    t143b.smirks = "[*:1]~[#16X4:2]-[#7X3:3]~[*:4]"
    torsion_handler.add_parameter(parameter=t143b, after=t143a.smirks)

    t143c = copy.deepcopy(t143)
    t143c.id = "t143c"
    t143c.smirks = "[*:1]~[#16X3+0:2]-[#7X3:3]~[*:4]"
    torsion_handler.add_parameter(parameter=t143c, after=t143b.smirks)

    t143d = copy.deepcopy(t143)
    t143d.id = "t143d"
    t143d.smirks = "[*:1]~[#16X4:2]-[#7X4:3]~[*:4]"
    torsion_handler.add_parameter(parameter=t143d, after=t143c.smirks)

    t143e = copy.deepcopy(t143)
    t143e.id = "t143e"
    t143e.smirks = "[*:1]~[#16X3+0:2]-[#7X4:3]~[*:4]"
    torsion_handler.add_parameter(parameter=t143e, after=t143d.smirks)

    # t144
    t144 = torsion_handler.get_parameter({"id": "t144"})[0]
    assert t144.smirks == "[#6X4:1]-[#16X4,#16X3+0:2]-[#7X4,#7X3:3]-[#1:4]"
    t144.smirks = "[#6X4:1]-[#16X4:2]-[#7X4:3]-[#1:4]"

    t144a = copy.deepcopy(t144)
    t144a.id = "t144a"
    t144a.smirks = "[#6X4:1]-[#16X3+0:2]-[#7X4:3]-[#1:4]"
    torsion_handler.add_parameter(parameter=t144a, after=t144.smirks)

    t144b = copy.deepcopy(t144)
    t144b.id = "t144b"
    t144b.smirks = "[#6X4:1]-[#16X4:2]-[#7X3:3]-[#1:4]"
    torsion_handler.add_parameter(parameter=t144b, after=t144a.smirks)

    t144c = copy.deepcopy(t144)
    t144c.id = "t144c"
    t144c.smirks = "[#6X4:1]-[#16X3+0:2]-[#7X3:3]-[#1:4]"
    torsion_handler.add_parameter(parameter=t144c, after=t144b.smirks)

    # t145
    t145 = torsion_handler.get_parameter({"id": "t145"})[0]
    assert t145.smirks == "[#6X3:1]-[#16X4,#16X3+0:2]-[#7X4,#7X3:3]-[#1:4]"
    t145.smirks = "[#6X3:1]-[#16X4:2]-[#7X4:3]-[#1:4]"

    t145a = copy.deepcopy(t145)
    t145a.id = "t145a"
    t145a.smirks = "[#6X3:1]-[#16X3+0:2]-[#7X4:3]-[#1:4]"
    torsion_handler.add_parameter(parameter=t145a, after=t145.smirks)

    t145b = copy.deepcopy(t145)
    t145b.id = "t145b"
    t145b.smirks = "[#6X3:1]-[#16X4:2]-[#7X3:3]-[#1:4]"
    torsion_handler.add_parameter(parameter=t145b, after=t145a.smirks)

    t145c = copy.deepcopy(t145)
    t145c.id = "t145c"
    t145c.smirks = "[#6X3:1]-[#16X3+0:2]-[#7X3:3]-[#1:4]"
    torsion_handler.add_parameter(parameter=t145c, after=t145b.smirks)

    # t146
    t146 = torsion_handler.get_parameter({"id": "t146"})[0]
    assert t146.smirks == "[#6X4:1]-[#16X4,#16X3+0:2]-[#7X4,#7X3:3]-[#6X4:4]"
    t146.smirks = "[#6X4:1]-[#16X4:2]-[#7X4:3]-[#6X4:4]"

    t146a = copy.deepcopy(t146)
    t146a.id = "t146a"
    t146a.smirks = "[#6X4:1]-[#16X3+0:2]-[#7X4:3]-[#6X4:4]"
    torsion_handler.add_parameter(parameter=t146a, after=t146.smirks)

    t146b = copy.deepcopy(t146)
    t146b.id = "t146b"
    t146b.smirks = "[#6X4:1]-[#16X4:2]-[#7X3:3]-[#6X4:4]"
    torsion_handler.add_parameter(parameter=t146b, after=t146a.smirks)

    t146c = copy.deepcopy(t146)
    t146c.id = "t146c"
    t146c.smirks = "[#6X4:1]-[#16X3+0:2]-[#7X3:3]-[#6X4:4]"
    torsion_handler.add_parameter(parameter=t146c, after=t146b.smirks)


    # t147
    t147 = torsion_handler.get_parameter({"id": "t147"})[0]
    assert t147.smirks == "[#6X3:1]-[#16X4,#16X3+0:2]-[#7X4,#7X3:3]-[#6X4:4]"
    t147.smirks = "[#6X3:1]-[#16X4:2]-[#7X4:3]-[#6X4:4]"

    t147a = copy.deepcopy(t147)
    t147a.id = "t147a"
    t147a.smirks = "[#6X3:1]-[#16X3+0:2]-[#7X4:3]-[#6X4:4]"
    torsion_handler.add_parameter(parameter=t147a, after=t147.smirks)

    t147b = copy.deepcopy(t147)
    t147b.id = "t147b"
    t147b.smirks = "[#6X3:1]-[#16X4:2]-[#7X3:3]-[#6X4:4]"
    torsion_handler.add_parameter(parameter=t147b, after=t147a.smirks)

    t147c = copy.deepcopy(t147)
    t147c.id = "t147c"
    t147c.smirks = "[#6X3:1]-[#16X3+0:2]-[#7X3:3]-[#6X4:4]"
    torsion_handler.add_parameter(parameter=t147c, after=t147b.smirks)

    # t148
    t148 = torsion_handler.get_parameter({"id": "t148"})[0]
    assert t148.smirks == "[#8X1:1]~[#16X4,#16X3+0:2]-[#7X4,#7X3:3]-[#1:4]"
    t148.smirks = "[#8X1:1]~[#16X4:2]-[#7X4:3]-[#1:4]"

    t148a = copy.deepcopy(t148)
    t148a.id = "t148a"
    t148a.smirks = "[#8X1:1]~[#16X3+0:2]-[#7X4:3]-[#1:4]"
    torsion_handler.add_parameter(parameter=t148a, after=t148.smirks)

    t148b = copy.deepcopy(t148)
    t148b.id = "t148b"
    t148b.smirks = "[#8X1:1]~[#16X4:2]-[#7X3:3]-[#1:4]"
    torsion_handler.add_parameter(parameter=t148b, after=t148a.smirks)

    t148c = copy.deepcopy(t148)
    t148c.id = "t148c"
    t148c.smirks = "[#8X1:1]~[#16X3+0:2]-[#7X3:3]-[#1:4]"
    torsion_handler.add_parameter(parameter=t148c, after=t148b.smirks)

    # t149
    t149 = torsion_handler.get_parameter({"id": "t149"})[0]
    assert t149.smirks == "[#8X1:1]~[#16X4,#16X3+0:2]-[#7X4,#7X3:3]-[#6X4:4]"
    t149.smirks = "[#8X1:1]~[#16X4:2]-[#7X4:3]-[#6X4:4]"

    t149a = copy.deepcopy(t149)
    t149a.id = "t149a"
    t149a.smirks = "[#8X1:1]~[#16X3+0:2]-[#7X4:3]-[#6X4:4]"
    torsion_handler.add_parameter(parameter=t149a, after=t149.smirks)

    t149b = copy.deepcopy(t149)
    t149b.id = "t149b"
    t149b.smirks = "[#8X1:1]~[#16X4:2]-[#7X3:3]-[#6X4:4]"
    torsion_handler.add_parameter(parameter=t149b, after=t149a.smirks)

    t149c = copy.deepcopy(t149)
    t149c.id = "t149c"
    t149c.smirks = "[#8X1:1]~[#16X3+0:2]-[#7X3:3]-[#6X4:4]"
    torsion_handler.add_parameter(parameter=t149c, after=t149b.smirks)

    # t150
    t150 = torsion_handler.get_parameter({"id": "t150"})[0]
    assert t150.smirks == "[#6X3:1]-[#16X4,#16X3+0:2]-[#7X3:3]-[#6X3:4]"
    t150.smirks = "[#6X3:1]-[#16X4:2]-[#7X3:3]-[#6X3:4]"

    t150a = copy.deepcopy(t150)
    t150a.id = "t150a"
    t150a.smirks = "[#6X3:1]-[#16X3+0:2]-[#7X3:3]-[#6X3:4]"
    torsion_handler.add_parameter(parameter=t150a, after=t150.smirks)

    # t151
    t151 = torsion_handler.get_parameter({"id": "t151"})[0]
    assert t151.smirks == "[#6X4:1]-[#16X4,#16X3+0:2]-[#7X3:3]-[#6X3:4]"
    t151.smirks = "[#6X4:1]-[#16X4:2]-[#7X3:3]-[#6X3:4]"

    t151a = copy.deepcopy(t151)
    t151a.id = "t151a"
    t151a.smirks = "[#6X4:1]-[#16X3+0:2]-[#7X3:3]-[#6X3:4]"
    torsion_handler.add_parameter(parameter=t151a, after=t151.smirks)

    # t152
    t152 = torsion_handler.get_parameter({"id": "t152"})[0]
    assert t152.smirks == "[#8X1:1]~[#16X4,#16X3+0:2]-[#7X3:3]-[#6X3:4]"
    t152.smirks = "[#8X1:1]~[#16X4:2]-[#7X3:3]-[#6X3:4]"

    t152a = copy.deepcopy(t152)
    t152a.id = "t152a"
    t152a.smirks = "[#8X1:1]~[#16X3+0:2]-[#7X3:3]-[#6X3:4]"
    torsion_handler.add_parameter(parameter=t152a, after=t152.smirks)

    # t153
    t153 = torsion_handler.get_parameter({"id": "t153"})[0]
    assert t153.smirks == "[#8X1:1]~[#16X4,#16X3+0:2]-[#7X3:3]-[#7X2:4]"
    t153.smirks = "[#8X1:1]~[#16X4:2]-[#7X3:3]-[#7X2:4]"

    t153a = copy.deepcopy(t153)
    t153a.id = "t153a"
    t153a.smirks = "[#8X1:1]~[#16X3+0:2]-[#7X3:3]-[#7X2:4]"
    torsion_handler.add_parameter(parameter=t153a, after=t153.smirks)

    # t154
    t154 = torsion_handler.get_parameter({"id": "t154"})[0]
    assert t154.smirks == "[*:1]~[#16X4,#16X3+0:2]=,:[#7X2:3]-,:[*:4]"
    t154.smirks = "[*:1]~[#16X4:2]=,:[#7X2:3]-,:[*:4]"

    t154a = copy.deepcopy(t154)
    t154a.id = "t154a"
    t154a.smirks = "[*:1]~[#16X3+0:2]=,:[#7X2:3]-,:[*:4]"
    torsion_handler.add_parameter(parameter=t154a, after=t154.smirks)

    # t155
    t155 = torsion_handler.get_parameter({"id": "t155"})[0]
    assert t155.smirks == "[#6X4:1]-[#16X4,#16X3+0:2]-[#7X2:3]~[#6X3:4]"
    t155.smirks = "[#6X4:1]-[#16X4:2]-[#7X2:3]~[#6X3:4]"

    t155a = copy.deepcopy(t155)
    t155a.id = "t155a"
    t155a.smirks = "[#6X4:1]-[#16X3+0:2]-[#7X2:3]~[#6X3:4]"
    torsion_handler.add_parameter(parameter=t155a, after=t155.smirks)

    # t156
    t156 = torsion_handler.get_parameter({"id": "t156"})[0]
    assert t156.smirks == "[#8X1:1]~[#16X4,#16X3+0:2]-[#7X2:3]~[#6X3:4]"
    t156.smirks = "[#8X1:1]~[#16X4:2]-[#7X2:3]~[#6X3:4]"

    t156a = copy.deepcopy(t156)
    t156a.id = "t156a"
    t156a.smirks = "[#8X1:1]~[#16X3+0:2]-[#7X2:3]~[#6X3:4]"
    torsion_handler.add_parameter(parameter=t156a, after=t156.smirks)

    # t157
    t157 = torsion_handler.get_parameter({"id": "t157"})[0]
    assert t157.smirks == "[*:1]~[#16X4,#16X3+0:2]-[#8X2:3]-[*:4]"
    t157.smirks = "[*:1]~[#16X4:2]-[#8X2:3]-[*:4]"

    t157a = copy.deepcopy(t157)
    t157a.id = "t157a"
    t157a.smirks = "[*:1]~[#16X3+0:2]-[#8X2:3]-[*:4]"
    torsion_handler.add_parameter(parameter=t157a, after=t157.smirks)

    # t158
    t158 = torsion_handler.get_parameter({"id": "t158"})[0]
    assert t158.smirks == "[*:1]-[#16X2,#16X3+1:2]-[#16X2,#16X3+1:3]-[*:4]"
    t158.smirks = "[*:1]-[#16X2:2]-[#16X2:3]-[*:4]"

    t158a = copy.deepcopy(t158)
    t158a.id = "t158a"
    t158a.smirks = "[*:1]-[#16X3+1:2]-[#16X2:3]-[*:4]"
    torsion_handler.add_parameter(parameter=t158a, after=t158.smirks)

    t158b = copy.deepcopy(t158)
    t158b.id = "t158b"
    t158b.smirks = "[*:1]-[#16X2:2]-[#16X3+1:3]-[*:4]"
    torsion_handler.add_parameter(parameter=t158b, after=t158a.smirks)

    t158c = copy.deepcopy(t158)
    t158c.id = "t158c"
    t158c.smirks = "[*:1]-[#16X3+1:2]-[#16X3+1:3]-[*:4]"
    torsion_handler.add_parameter(parameter=t158c, after=t158b.smirks)
    

def modify_cresset_torsions(torsion_handler):
    # cresset suggestion
    t47 = torsion_handler.get_parameter({"id": "t47"})[0]
    assert t47.smirks == "[*:1]~[#6X3:2]-[#6X3$(*=[#8,#16,#7]):3]~[*:4]"

    t47a = copy.deepcopy(t47)
    t47a.id = "t47a"
    t47a.smirks = "[#1:1]~[#6X3:2]-[#6X3$(*=[#8,#16,#7]):3]~[*:4]"
    torsion_handler.add_parameter(parameter=t47a, after=t47.smirks)

    t17 = torsion_handler.get_parameter({"id": "t17"})[0]
    assert t17.smirks == "[*:1]~[#6X3:2]-[#6X4:3]-[*:4]"

    t17a = copy.deepcopy(t17)
    t17a.id = "t17a"
    t17a.smirks = "[aR:1]~[#6X3aR:2]-[#6X4:3]-[*:4]"
    torsion_handler.add_parameter(parameter=t17a, after=t17.smirks)

    # biaryl torsions
    t43 = torsion_handler.get_parameter({"id": "t43"})[0]
    assert t43.smirks == "[*:1]~[#6X3:2]-[#6X3:3]~[*:4]"

    t43a = copy.deepcopy(t43)
    t43a.id = "t43a"
    t43a.smirks = "[*:1]~[#6X3aR:2]-[#6X3aR:3]~[*:4]"
    torsion_handler.add_parameter(parameter=t43a, after=t43.smirks)


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
    default="output/initial-force-field-v6.offxml",
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

    # angles
    angle_handler = force_field.get_parameter_handler("Angles")
    modify_r3_angles(angle_handler)

    split_msm_angles(angle_handler)

    bond_handler = force_field.get_parameter_handler("Bonds")
    split_msm_bonds(bond_handler)

    fix_problematic_angles(angle_handler)
    move_r4_angles_to_end(angle_handler)
    fix_problematic_bonds(bond_handler)


    # Write out file
    force_field.to_file(output_path, discard_cosmetic_attributes=True)
    print(f"Force field written to {output_path}")


if __name__ == "__main__":
    generate_force_field()
