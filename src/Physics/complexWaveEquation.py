from PYB11Generator import *
from physics import *

@PYB11template("dim")
class ComplexWaveEquation(Physics):
    def pyinit(self,
               nodeList="NodeList*",
               constants="PhysicalConstants&",
               grid="Mesh::Grid<%(dim)s>*",
               C="double"):
        return
    @PYB11cppname("getCell")
    def getCell2d(self,i="int",j="int",fieldName="std::string"):
        return

ComplexWaveEquation1d = PYB11TemplateClass(ComplexWaveEquation,
                              template_parameters = ("1"),
                              cppname = "ComplexWaveEquation<1>",
                              pyname = "ComplexWaveEquation1d",
                              docext = " (1D).")
ComplexWaveEquation2d = PYB11TemplateClass(ComplexWaveEquation,
                              template_parameters = ("2"),
                              cppname = "ComplexWaveEquation<2>",
                              pyname = "ComplexWaveEquation2d",
                              docext = " (2D).")
ComplexWaveEquation3d = PYB11TemplateClass(ComplexWaveEquation,
                              template_parameters = ("3"),
                              cppname = "ComplexWaveEquation<3>",
                              pyname = "ComplexWaveEquation3d",
                              docext = " (3D).") 