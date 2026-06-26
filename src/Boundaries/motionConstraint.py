from PYB11Generator import *
from constraint import *

@PYB11template("dim")
class MotionConstraint(Constraint):
    def pyinit(self,nodeList="NodeList*", nodeIndicies="std::vector<int>", direction="Lin::Vector<%(dim)s>"):
        return
    
MotionConstraint1d = PYB11TemplateClass(MotionConstraint,
                              template_parameters = ("1"),
                              cppname = "MotionConstraint<1>",
                              pyname = "MotionConstraint1d",
                              docext = " (1D).")
MotionConstraint2d = PYB11TemplateClass(MotionConstraint,
                              template_parameters = ("2"),
                              cppname = "MotionConstraint<2>",
                              pyname = "MotionConstraint2d",
                              docext = " (2D).")
MotionConstraint3d = PYB11TemplateClass(MotionConstraint,
                              template_parameters = ("3"),
                              cppname = "MotionConstraint<3>",
                              pyname = "MotionConstraint3d",
                              docext = " (3D).") 