from PYB11Generator import *
from boundary import *

@PYB11template("dim")
class GridBoundary(Boundary):
    def pyinit(self,grid="Mesh::Grid<%(dim)s>*"):
        return
    
GridBoundary1d = PYB11TemplateClass(GridBoundary,
                              template_parameters = ("1"),
                              cppname = "GridBoundary<1>",
                              pyname = "GridBoundary1d",
                              docext = " (1D).")
GridBoundary2d = PYB11TemplateClass(GridBoundary,
                              template_parameters = ("2"),
                              cppname = "GridBoundary<2>",
                              pyname = "GridBoundary2d",
                              docext = " (2D).")
GridBoundary3d = PYB11TemplateClass(GridBoundary,
                              template_parameters = ("3"),
                              cppname = "GridBoundary<3>",
                              pyname = "GridBoundary3d",
                              docext = " (3D).") 