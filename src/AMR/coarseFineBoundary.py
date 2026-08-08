# Copyright (C) 2026  Cody Raskin

from PYB11Generator import *
from boundary import Boundary

@PYB11template("dim")
class CoarseFineBoundary(Boundary):
    def pyinit(self,
               coarseNodeList="NodeList*",
               ghostIds="std::vector<int>",
               coarseIds="std::vector<int>"):
        return

CoarseFineBoundary1d = PYB11TemplateClass(CoarseFineBoundary,
                              template_parameters = ("1"),
                              cppname = "AMR::CoarseFineBoundary<1>",
                              pyname = "CoarseFineBoundary1d",
                              docext = " (1D).")
CoarseFineBoundary2d = PYB11TemplateClass(CoarseFineBoundary,
                              template_parameters = ("2"),
                              cppname = "AMR::CoarseFineBoundary<2>",
                              pyname = "CoarseFineBoundary2d",
                              docext = " (2D).")
CoarseFineBoundary3d = PYB11TemplateClass(CoarseFineBoundary,
                              template_parameters = ("3"),
                              cppname = "AMR::CoarseFineBoundary<3>",
                              pyname = "CoarseFineBoundary3d",
                              docext = " (3D).")
