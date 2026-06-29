from PYB11Generator import *
from gridBoundary import *

@PYB11template("dim")
class ReflectingGridBoundary(GridBoundary):
    def pyinit(self, grid="Mesh::Grid<%(dim)s>*"):
        return
    def addBox(self, p1="Lin::Vector<%(dim)s>", p2="Lin::Vector<%(dim)s>"):
        return
    def removeBox(self, p1="Lin::Vector<%(dim)s>", p2="Lin::Vector<%(dim)s>"):
        return
    def addSphere(self, p="Lin::Vector<%(dim)s>", radius="double"):
        return
    def removeSphere(self, p="Lin::Vector<%(dim)s>", radius="double"):
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