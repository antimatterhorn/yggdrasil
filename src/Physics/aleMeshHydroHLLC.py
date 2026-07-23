# Copyright (C) 2026  Cody Raskin

from PYB11Generator import *
from hydro import *

@PYB11template("dim")
class ALEMeshHydroHLLC(Hydro):
    def pyinit(self,
               nodeList="NodeList*",
               constants="PhysicalConstants&",
               eos="EquationOfState*",
               mesh="Mesh::ALEMesh<%(dim)s>*",
               nodeVelocities=("NodeList*", "nullptr")):
        "nodeVelocities: optional node-centered NodeList (one entry per mesh node) with a Vector 'velocity' field prescribing mesh motion. Omit (or pass None) for a static (Eulerian) mesh."
        return

ALEMeshHydroHLLC2d = PYB11TemplateClass(ALEMeshHydroHLLC,
                              template_parameters = ("2"),
                              cppname = "ALEMeshHydroHLLC<2>",
                              pyname = "ALEMeshHydroHLLC2d",
                              docext = " (2D).")
