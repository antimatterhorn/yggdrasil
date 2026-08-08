# Copyright (C) 2026  Cody Raskin

from PYB11Generator import *
# Only Boundary itself, not boundary.py's Boundary1d/2d/3d, which would collide with Boundaries' own registration.
from boundary import Boundary

@PYB11template("dim")
class PatchNeighborBoundary(Boundary):
    def pyinit(self,
               neighborNodeList="NodeList*",
               myGhostIds="std::vector<int>",
               neighborInteriorIds="std::vector<int>"):
        return

PatchNeighborBoundary1d = PYB11TemplateClass(PatchNeighborBoundary,
                              template_parameters = ("1"),
                              cppname = "AMR::PatchNeighborBoundary<1>",
                              pyname = "PatchNeighborBoundary1d",
                              docext = " (1D).")
PatchNeighborBoundary2d = PYB11TemplateClass(PatchNeighborBoundary,
                              template_parameters = ("2"),
                              cppname = "AMR::PatchNeighborBoundary<2>",
                              pyname = "PatchNeighborBoundary2d",
                              docext = " (2D).")
PatchNeighborBoundary3d = PYB11TemplateClass(PatchNeighborBoundary,
                              template_parameters = ("3"),
                              cppname = "AMR::PatchNeighborBoundary<3>",
                              pyname = "PatchNeighborBoundary3d",
                              docext = " (3D).")
