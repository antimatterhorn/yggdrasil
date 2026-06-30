from PYB11Generator import *
from physics import *

@PYB11template("dim")
class DEM(Physics):
    def pyinit(self,
               nodeList="NodeList*",
               constants="PhysicalConstants&",
               kn="double",
               kt="double",
               gammaN="double",
               gammaT="double",
               mu="double"):
        return

DEM1d = PYB11TemplateClass(DEM,
                              template_parameters = ("1"),
                              cppname = "DEM<1>",
                              pyname = "DEM1d",
                              docext = " (1D).")
DEM2d = PYB11TemplateClass(DEM,
                              template_parameters = ("2"),
                              cppname = "DEM<2>",
                              pyname = "DEM2d",
                              docext = " (2D).")
DEM3d = PYB11TemplateClass(DEM,
                              template_parameters = ("3"),
                              cppname = "DEM<3>",
                              pyname = "DEM3d",
                              docext = " (3D).")
