from PYB11Generator import *
from physics import *

@PYB11template("dim")
class ThermalConduction(Physics):
    def pyinit(self,
               nodeList="NodeList*",
               constants="PhysicalConstants&",
               eos="EquationOfState*",
               opac="OpacityModel*",
               grid="Mesh::Grid<%(dim)s>*"):
        return
    @PYB11cppname("getCell")
    def getCell2d(self,i="int",j="int",fieldName="std::string"):
        return

ThermalConduction1d = PYB11TemplateClass(ThermalConduction,
                              template_parameters = ("1"),
                              cppname = "ThermalConduction<1>",
                              pyname = "ThermalConduction1d",
                              docext = " (1D).")
ThermalConduction2d = PYB11TemplateClass(ThermalConduction,
                              template_parameters = ("2"),
                              cppname = "ThermalConduction<2>",
                              pyname = "ThermalConduction2d",
                              docext = " (2D).")
ThermalConduction3d = PYB11TemplateClass(ThermalConduction,
                              template_parameters = ("3"),
                              cppname = "ThermalConduction<3>",
                              pyname = "ThermalConduction3d",
                              docext = " (3D).") 