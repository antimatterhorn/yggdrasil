from PYB11Generator import *

@PYB11template("dim")
class SpatialTree:
    def pyinit(self,
               pos="Field<Lin::Vector<%(dim)s>>*",
               mass="Field<double>*"):
        return

# SpatialTree1d = PYB11TemplateClass(SpatialTree,
#                               template_parameters = ("1"),
#                               cppname = "SpatialTree<1>",
#                               pyname = "SpatialTree1d",
#                               docext = " (1D).")
SpatialTree2d = PYB11TemplateClass(SpatialTree,
                              template_parameters = ("2"),
                              cppname = "SpatialTree<2>",
                              pyname = "SpatialTree2d",
                              docext = " (2D).")
SpatialTree3d = PYB11TemplateClass(SpatialTree,
                              template_parameters = ("3"),
                              cppname = "SpatialTree<3>",
                              pyname = "SpatialTree3d",
                              docext = " (3D).") 