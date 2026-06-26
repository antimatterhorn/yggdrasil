from PYB11Generator import *
from gridBoundary import *

@PYB11template("dim")
class DirichletGridBoundary(GridBoundary):
    def pyinit(self,grid="Mesh::Grid<%(dim)s>*"):
        return
    def addBox(self,p1="Lin::Vector<%(dim)s>",p2="Lin::Vector<%(dim)s>"):
        return
    def removeBox(self,p1="Lin::Vector<%(dim)s>",p2="Lin::Vector<%(dim)s>"):
        return
    def addSphere(self,p="Lin::Vector<%(dim)s>",radius="double"):
        return
    def removeSphere(self,p="Lin::Vector<%(dim)s>",radius="double"):
        return
    def addDomain(self):
        return
    def boundaryIds(self):
        return "std::vector<int>"
    
DirichletGridBoundary1d = PYB11TemplateClass(DirichletGridBoundary,
                              template_parameters = ("1"),
                              cppname = "DirichletGridBoundary<1>",
                              pyname = "DirichletGridBoundary1d",
                              docext = " (1D).")
DirichletGridBoundary2d = PYB11TemplateClass(DirichletGridBoundary,
                              template_parameters = ("2"),
                              cppname = "DirichletGridBoundary<2>",
                              pyname = "DirichletGridBoundary2d",
                              docext = " (2D).")
DirichletGridBoundary3d = PYB11TemplateClass(DirichletGridBoundary,
                              template_parameters = ("3"),
                              cppname = "DirichletGridBoundary<3>",
                              pyname = "DirichletGridBoundary3d",
                              docext = " (3D).") 