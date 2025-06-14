from PYB11Generator import *
from gridBoundary import *

@PYB11template("dim")
class OutflowGridBoundary(GridBoundary):
    def pyinit(self,grid="Mesh::Grid<%(dim)s>*"):
        return
    def pyinit1(self,grid="Mesh::Grid<%(dim)s>*",derivative="std::string"):
        return
    
OutflowGridBoundary1d = PYB11TemplateClass(OutflowGridBoundary,
                              template_parameters = ("1"),
                              cppname = "OutflowGridBoundary<1>",
                              pyname = "OutflowGridBoundary1d",
                              docext = " (1D).")
OutflowGridBoundary2d = PYB11TemplateClass(OutflowGridBoundary,
                              template_parameters = ("2"),
                              cppname = "OutflowGridBoundary<2>",
                              pyname = "OutflowGridBoundary2d",
                              docext = " (2D).")
OutflowGridBoundary3d = PYB11TemplateClass(OutflowGridBoundary,
                              template_parameters = ("3"),
                              cppname = "OutflowGridBoundary<3>",
                              pyname = "OutflowGridBoundary3d",
                              docext = " (3D).") 