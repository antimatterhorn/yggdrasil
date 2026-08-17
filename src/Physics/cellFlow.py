# Copyright (C) 2026  Cody Raskin

from PYB11Generator import *
from physics import *

@PYB11template("dim")
class CellFlowPhysics(Physics):
    def pyinit(self,
               nodeList="NodeList*",
               constants="PhysicalConstants&",
               numTypes="int",
               forceTable="std::vector<double>",
               radiusByType="std::vector<double>",
               baseRadius="double",
               repulsion="double",
               attraction="double",
               k="double",
               forceMultiplier="double",
               damping="double"):
        return

CellFlowPhysics1d = PYB11TemplateClass(CellFlowPhysics,
                              template_parameters = ("1"),
                              cppname = "CellFlowPhysics<1>",
                              pyname = "CellFlowPhysics1d",
                              docext = " (1D).")
CellFlowPhysics2d = PYB11TemplateClass(CellFlowPhysics,
                              template_parameters = ("2"),
                              cppname = "CellFlowPhysics<2>",
                              pyname = "CellFlowPhysics2d",
                              docext = " (2D).")
CellFlowPhysics3d = PYB11TemplateClass(CellFlowPhysics,
                              template_parameters = ("3"),
                              cppname = "CellFlowPhysics<3>",
                              pyname = "CellFlowPhysics3d",
                              docext = " (3D).")
