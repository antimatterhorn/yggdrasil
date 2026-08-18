# Copyright (C) 2026  Cody Raskin

from PYB11Generator import *
from physics import *

@PYB11template("dim")
class FMMGravity(Physics):
    def pyinit(self,
               nodeList="NodeList*",
               constants="PhysicalConstants&",
               plummerLength="double",
               maxSourcesPerLeaf=("int", 16)):
        return

FMMGravity1d = PYB11TemplateClass(FMMGravity,
                              template_parameters = ("1"),
                              cppname = "FMMGravity<1>",
                              pyname = "FMMGravity1d",
                              docext = " (1D).")
FMMGravity2d = PYB11TemplateClass(FMMGravity,
                              template_parameters = ("2"),
                              cppname = "FMMGravity<2>",
                              pyname = "FMMGravity2d",
                              docext = " (2D).")
FMMGravity3d = PYB11TemplateClass(FMMGravity,
                              template_parameters = ("3"),
                              cppname = "FMMGravity<3>",
                              pyname = "FMMGravity3d",
                              docext = " (3D).")
