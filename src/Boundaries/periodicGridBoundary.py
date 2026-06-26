from PYB11Generator import *
from gridBoundary import *

@PYB11template("dim")
class PeriodicGridBoundary(GridBoundary):
    def pyinit(self,grid="Mesh::Grid<%(dim)s>*"):
        return
    
PeriodicGridBoundary1d = PYB11TemplateClass(PeriodicGridBoundary,
                              template_parameters = ("1"),
                              cppname = "PeriodicGridBoundary<1>",
                              pyname = "PeriodicGridBoundary1d",
                              docext = " (1D).")
PeriodicGridBoundary2d = PYB11TemplateClass(PeriodicGridBoundary,
                              template_parameters = ("2"),
                              cppname = "PeriodicGridBoundary<2>",
                              pyname = "PeriodicGridBoundary2d",
                              docext = " (2D).")
PeriodicGridBoundary3d = PYB11TemplateClass(PeriodicGridBoundary,
                              template_parameters = ("3"),
                              cppname = "PeriodicGridBoundary<3>",
                              pyname = "PeriodicGridBoundary3d",
                              docext = " (3D).") 