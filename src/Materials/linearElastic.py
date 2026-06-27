from PYB11Generator import *
from material import *

@PYB11template("dim")
class IsotropicLinearElastic(FEMMaterial):
    def pyinit(self,
               E="double",
               nu="double",
               condition="PlaneCondition"):
        return

    @PYB11const
    def youngsModulus(self):
        return "double"

    @PYB11const
    def poissonsRatio(self):
        return "double"

IsotropicLinearElastic2d = PYB11TemplateClass(IsotropicLinearElastic,
                                              template_parameters=("2"),
                                              cppname="IsotropicLinearElastic<2>",
                                              pyname="IsotropicLinearElastic2d",
                                              docext=" (2D).")
IsotropicLinearElastic3d = PYB11TemplateClass(IsotropicLinearElastic,
                                              template_parameters=("3"),
                                              cppname="IsotropicLinearElastic<3>",
                                              pyname="IsotropicLinearElastic3d",
                                              docext=" (3D).")
