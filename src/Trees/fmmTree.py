# Copyright (C) 2026  Cody Raskin

from PYB11Generator import *

@PYB11template("dim")
class FMMTree:
    def pyinit(self,
               pos="Field<Lin::Vector<%(dim)s>>*",
               mass="Field<double>*",
               maxSourcesPerLeaf=("int", 16),
               maxDepth=("int", 40)):
        return

FMMTree1d = PYB11TemplateClass(FMMTree,
                              template_parameters = ("1"),
                              cppname = "FMMTree<1>",
                              pyname = "FMMTree1d",
                              docext = " (1D).")
FMMTree2d = PYB11TemplateClass(FMMTree,
                              template_parameters = ("2"),
                              cppname = "FMMTree<2>",
                              pyname = "FMMTree2d",
                              docext = " (2D).")
FMMTree3d = PYB11TemplateClass(FMMTree,
                              template_parameters = ("3"),
                              cppname = "FMMTree<3>",
                              pyname = "FMMTree3d",
                              docext = " (3D).")
