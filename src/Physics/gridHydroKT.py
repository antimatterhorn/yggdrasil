from PYB11Generator import *
from hydro import *

@PYB11template("dim")
class GridHydroKT(Hydro):
    def pyinit(self,
               nodeList="NodeList*",
               constants="PhysicalConstants&",
               eos="EquationOfState*",
               grid="Mesh::Grid<%(dim)s>*"):
        return
    @PYB11cppname("getCell")
    def getCell2d(self,i="int",j="int",fieldName="std::string"):
        return


GridHydroKT1d = PYB11TemplateClass(GridHydroKT,
                              template_parameters = ("1"),
                              cppname = "GridHydroKT<1>",
                              pyname = "GridHydroKT1d",
                              docext = " (1D).")
GridHydroKT2d = PYB11TemplateClass(GridHydroKT,
                              template_parameters = ("2"),
                              cppname = "GridHydroKT<2>",
                              pyname = "GridHydroKT2d",
                              docext = " (2D).")
GridHydroKT3d = PYB11TemplateClass(GridHydroKT,
                              template_parameters = ("3"),
                              cppname = "GridHydroKT<3>",
                              pyname = "GridHydroKT3d",
                              docext = " (3D).") 