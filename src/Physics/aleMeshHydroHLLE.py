# Copyright (C) 2026  Cody Raskin

from PYB11Generator import *
from hydro import *

@PYB11template("dim")
class ALEMeshHydroHLLE(Hydro):
    def pyinit(self,
               nodeList="NodeList*",
               constants="PhysicalConstants&",
               eos="EquationOfState*",
               mesh="Mesh::ALEMesh<%(dim)s>*"):
        return

ALEMeshHydroHLLE2d = PYB11TemplateClass(ALEMeshHydroHLLE,
                              template_parameters = ("2"),
                              cppname = "ALEMeshHydroHLLE<2>",
                              pyname = "ALEMeshHydroHLLE2d",
                              docext = " (2D).")
