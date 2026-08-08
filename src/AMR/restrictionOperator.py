# Copyright (C) 2026  Cody Raskin

from PYB11Generator import *

@PYB11template("dim")
class RestrictionOperator:
    def pyinit(self,
               fineNodes="NodeList*",
               coarseNodes="NodeList*",
               fineGrid="Mesh::Grid<%(dim)s>*",
               coarseGrid="Mesh::Grid<%(dim)s>*",
               eos="EquationOfState*",
               fineIds="std::vector<int>",
               coarseIds="std::vector<int>"):
        return
    def apply(self):
        "Average the fine patch's interior cells onto the coarse cells they cover."
        return

RestrictionOperator1d = PYB11TemplateClass(RestrictionOperator,
                              template_parameters = ("1"),
                              cppname = "AMR::RestrictionOperator<1>",
                              pyname = "RestrictionOperator1d",
                              docext = " (1D).")
RestrictionOperator2d = PYB11TemplateClass(RestrictionOperator,
                              template_parameters = ("2"),
                              cppname = "AMR::RestrictionOperator<2>",
                              pyname = "RestrictionOperator2d",
                              docext = " (2D).")
RestrictionOperator3d = PYB11TemplateClass(RestrictionOperator,
                              template_parameters = ("3"),
                              cppname = "AMR::RestrictionOperator<3>",
                              pyname = "RestrictionOperator3d",
                              docext = " (3D).")
