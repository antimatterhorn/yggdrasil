from PYB11Generator import *
from physics import *

@PYB11template("dim")
class DEMConstantForce(Physics):
    def pyinit(self,
               nodeList="NodeList*",
               constants="PhysicalConstants&",
               bodyForce="Lin::Vector<%(dim)s>"):
        return

DEMConstantForce1d = PYB11TemplateClass(DEMConstantForce,
                              template_parameters = ("1"),
                              cppname = "DEMConstantForce<1>",
                              pyname = "DEMConstantForce1d",
                              docext = " (1D).")
DEMConstantForce2d = PYB11TemplateClass(DEMConstantForce,
                              template_parameters = ("2"),
                              cppname = "DEMConstantForce<2>",
                              pyname = "DEMConstantForce2d",
                              docext = " (2D).")
DEMConstantForce3d = PYB11TemplateClass(DEMConstantForce,
                              template_parameters = ("3"),
                              cppname = "DEMConstantForce<3>",
                              pyname = "DEMConstantForce3d",
                              docext = " (3D).")
