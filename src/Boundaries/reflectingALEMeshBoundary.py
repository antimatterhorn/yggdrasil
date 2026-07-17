# Copyright (C) 2026  Cody Raskin

from PYB11Generator import *
from boundary import *

@PYB11template("dim")
class ReflectingALEMeshBoundary(Boundary):
    def pyinit(self, mesh="Mesh::ALEMesh<%(dim)s>*"):
        return
    def setFaces(self, ids="std::vector<size_t>"):
        return
    def setAllBoundaryFaces(self):
        return

ReflectingALEMeshBoundary2d = PYB11TemplateClass(ReflectingALEMeshBoundary,
                              template_parameters = ("2"),
                              cppname = "ReflectingALEMeshBoundary<2>",
                              pyname = "ReflectingALEMeshBoundary2d",
                              docext = " (2D).")
