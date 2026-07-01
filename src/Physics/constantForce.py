# Copyright (C) 2026  Cody Raskin

from PYB11Generator import *
from physics import *

@PYB11template("dim")
class ConstantForce(Physics):
    def pyinit(self,
               nodeList="NodeList*",
               constants="PhysicalConstants&",
               forceVector="Lin::Vector<%(dim)s>&"):
        return

ConstantForce1d = PYB11TemplateClass(ConstantForce,
                              template_parameters = ("1"),
                              cppname = "ConstantForce<1>",
                              pyname = "ConstantForce1d",
                              docext = " (1D).")
ConstantForce2d = PYB11TemplateClass(ConstantForce,
                              template_parameters = ("2"),
                              cppname = "ConstantForce<2>",
                              pyname = "ConstantForce2d",
                              docext = " (2D).")
ConstantForce3d = PYB11TemplateClass(ConstantForce,
                              template_parameters = ("3"),
                              cppname = "ConstantForce<3>",
                              pyname = "ConstantForce3d",
                              docext = " (3D).") 