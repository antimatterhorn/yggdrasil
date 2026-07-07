# Copyright (C) 2026  Cody Raskin

from PYB11Generator import *

@PYB11template("dim")
class FEMMaterial:
    @PYB11const
    def name(self):
        return "std::string"

FEMMaterial2d = PYB11TemplateClass(FEMMaterial,
                                   template_parameters=("2"),
                                   cppname="FEMMaterial<2>",
                                   pyname="FEMMaterial2d",
                                   docext=" (2D).")
FEMMaterial3d = PYB11TemplateClass(FEMMaterial,
                                   template_parameters=("3"),
                                   cppname="FEMMaterial<3>",
                                   pyname="FEMMaterial3d",
                                   docext=" (3D).")
