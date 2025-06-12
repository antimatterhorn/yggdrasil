from PYB11Generator import *
from physics import *

@PYB11template("dim")
class PhaseCoupling(Physics):
    def pyinit(self,
               nodeList="NodeList*",
               constants="PhysicalConstants&",
               couplingConstant="double"):
        return
    def pyinit1(self,
               nodeList="NodeList*",
               constants="PhysicalConstants&",
               couplingConstant="double",
               searchRadius="double"):
        return

PhaseCoupling1d = PYB11TemplateClass(PhaseCoupling,
                              template_parameters = ("1"),
                              cppname = "PhaseCoupling<1>",
                              pyname = "PhaseCoupling1d",
                              docext = " (1D).")
PhaseCoupling2d = PYB11TemplateClass(PhaseCoupling,
                              template_parameters = ("2"),
                              cppname = "PhaseCoupling<2>",
                              pyname = "PhaseCoupling2d",
                              docext = " (2D).")
PhaseCoupling3d = PYB11TemplateClass(PhaseCoupling,
                              template_parameters = ("3"),
                              cppname = "PhaseCoupling<3>",
                              pyname = "PhaseCoupling3d",
                              docext = " (3D).") 