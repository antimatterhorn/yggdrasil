# Copyright (C) 2026  Cody Raskin

from PYB11Generator import *
from boundary import *

@PYB11template("dim")
class DirichletALEMeshBoundary(Boundary):
    def pyinit(self, mesh="Mesh::ALEMesh<%(dim)s>*",
               rho="double", v="Lin::Vector<%(dim)s>", u="double"):
        return
    def setFaces(self, ids="std::vector<size_t>"):
        return
    def setAllBoundaryFaces(self):
        return

DirichletALEMeshBoundary2d = PYB11TemplateClass(DirichletALEMeshBoundary,
                              template_parameters = ("2"),
                              cppname = "DirichletALEMeshBoundary<2>",
                              pyname = "DirichletALEMeshBoundary2d",
                              docext = " (2D).")
