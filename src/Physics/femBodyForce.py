from PYB11Generator import *
from physics import *

@PYB11template("dim")
class ConstantBodyForce(Physics):
    def pyinit(self,
               nodeList="NodeList*",
               constants="PhysicalConstants&",
               g="Lin::Vector<%(dim)s>"):
        return

ConstantBodyForce1d = PYB11TemplateClass(ConstantBodyForce,
                                         template_parameters=("1"),
                                         cppname="ConstantBodyForce<1>",
                                         pyname="ConstantBodyForce1d",
                                         docext=" (1D).")
ConstantBodyForce2d = PYB11TemplateClass(ConstantBodyForce,
                                         template_parameters=("2"),
                                         cppname="ConstantBodyForce<2>",
                                         pyname="ConstantBodyForce2d",
                                         docext=" (2D).")
ConstantBodyForce3d = PYB11TemplateClass(ConstantBodyForce,
                                         template_parameters=("3"),
                                         cppname="ConstantBodyForce<3>",
                                         pyname="ConstantBodyForce3d",
                                         docext=" (3D).")
