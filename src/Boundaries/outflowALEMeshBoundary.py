# Copyright (C) 2026  Cody Raskin

from PYB11Generator import *
from boundary import *

@PYB11template("dim")
class OutflowALEMeshBoundary(Boundary):
    def pyinit(self, mesh="Mesh::ALEMesh<%(dim)s>*"):
        return
    def setFaces(self, ids="std::vector<size_t>"):
        return
    def setAllBoundaryFaces(self):
        return

OutflowALEMeshBoundary2d = PYB11TemplateClass(OutflowALEMeshBoundary,
                              template_parameters = ("2"),
                              cppname = "OutflowALEMeshBoundary<2>",
                              pyname = "OutflowALEMeshBoundary2d",
                              docext = " (2D).")
