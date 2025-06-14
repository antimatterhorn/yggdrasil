from PYB11Generator import *
from boundaries import *

@PYB11template("dim")
class Constraint(Boundaries):
    def pyinit(self,nodeList="NodeList*", nodeIndicies="std::vector<int>"):
        return
    
Constraint1d = PYB11TemplateClass(Constraint,
                              template_parameters = ("1"),
                              cppname = "Constraint<1>",
                              pyname = "Constraint1d",
                              docext = " (1D).")
Constraint2d = PYB11TemplateClass(Constraint,
                              template_parameters = ("2"),
                              cppname = "Constraint<2>",
                              pyname = "Constraint2d",
                              docext = " (2D).")
Constraint3d = PYB11TemplateClass(Constraint,
                              template_parameters = ("3"),
                              cppname = "Constraint<3>",
                              pyname = "Constraint3d",
                              docext = " (3D).") 