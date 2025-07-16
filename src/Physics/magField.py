from PYB11Generator import *
from physics import *

@PYB11template("dim")
class MagField(Physics):
    def pyinit(self,
               nodeList="NodeList*",
               constants="PhysicalConstants&",
               magF="Lin::Vector<%(dim)s>&"):
        return
    def SetB(self, mF="Lin::Vector<%(dim)s>&"):
        return

# MagField1d = PYB11TemplateClass(MagField,
#                               template_parameters = ("1"),
#                               cppname = "MagField<1>",
#                               pyname = "MagField1d",
#                               docext = " (1D).")
# MagField2d = PYB11TemplateClass(MagField,
#                               template_parameters = ("2"),
#                               cppname = "MagField<2>",
#                               pyname = "MagField2d",
#                               docext = " (2D).")
MagField3d = PYB11TemplateClass(MagField,
                              template_parameters = ("3"),
                              cppname = "MagField<3>",
                              pyname = "MagField3d",
                              docext = " (3D).") 