from PYB11Generator import *
from physics import *

@PYB11template("dim")
class TreeGravity(Physics):
    def pyinit(self,
               nodeList="NodeList*",
               constants="PhysicalConstants&",
               plummerLength="double"):
        return

TreeGravity1d = PYB11TemplateClass(TreeGravity,
                              template_parameters = ("1"),
                              cppname = "TreeGravity<1>",
                              pyname = "TreeGravity1d",
                              docext = " (1D).")
TreeGravity2d = PYB11TemplateClass(TreeGravity,
                              template_parameters = ("2"),
                              cppname = "TreeGravity<2>",
                              pyname = "TreeGravity2d",
                              docext = " (2D).")
TreeGravity3d = PYB11TemplateClass(TreeGravity,
                              template_parameters = ("3"),
                              cppname = "TreeGravity<3>",
                              pyname = "TreeGravity3d",
                              docext = " (3D).") 