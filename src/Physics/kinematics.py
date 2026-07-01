from PYB11Generator import *
from physics import *

@PYB11template("dim")
class Kinematics(Physics):
    def pyinit(self,
               nodeList="NodeList*",
               constants="PhysicalConstants&"):
        return

Kinematics1d = PYB11TemplateClass(Kinematics,
                             template_parameters = ("1"),
                             cppname = "Kinematics<1>",
                             pyname = "Kinematics1d",
                             docext = " (1D).")
Kinematics2d = PYB11TemplateClass(Kinematics,
                             template_parameters = ("2"),
                             cppname = "Kinematics<2>",
                             pyname = "Kinematics2d",
                             docext = " (2D).")
Kinematics3d = PYB11TemplateClass(Kinematics,
                             template_parameters = ("3"),
                             cppname = "Kinematics<3>",
                             pyname = "Kinematics3d",
                             docext = " (3D).")
