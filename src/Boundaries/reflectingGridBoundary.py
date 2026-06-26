from PYB11Generator import *
from gridBoundary import *

@PYB11template("dim")
class ReflectingGridBoundary(GridBoundary):
    def pyinit(self,grid="Mesh::Grid<%(dim)s>*"):
        return
    
ReflectingGridBoundary1d = PYB11TemplateClass(ReflectingGridBoundary,
                              template_parameters = ("1"),
                              cppname = "ReflectingGridBoundary<1>",
                              pyname = "ReflectingGridBoundary1d",
                              docext = " (1D).")
ReflectingGridBoundary2d = PYB11TemplateClass(ReflectingGridBoundary,
                              template_parameters = ("2"),
                              cppname = "ReflectingGridBoundary<2>",
                              pyname = "ReflectingGridBoundary2d",
                              docext = " (2D).")
ReflectingGridBoundary3d = PYB11TemplateClass(ReflectingGridBoundary,
                              template_parameters = ("3"),
                              cppname = "ReflectingGridBoundary<3>",
                              pyname = "ReflectingGridBoundary3d",
                              docext = " (3D).") 