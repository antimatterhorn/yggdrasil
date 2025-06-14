from PYB11Generator import *

@PYB11template("dim")
class Boundary:
    def pyinit(self):
        return
    
Boundary1d = PYB11TemplateClass(Boundary,
                              template_parameters = ("1"),
                              cppname = "Boundary<1>",
                              pyname = "Boundary1d",
                              docext = " (1D).")
Boundary2d = PYB11TemplateClass(Boundary,
                              template_parameters = ("2"),
                              cppname = "Boundary<2>",
                              pyname = "Boundary2d",
                              docext = " (2D).")
Boundary3d = PYB11TemplateClass(Boundary,
                              template_parameters = ("3"),
                              cppname = "Boundary<3>",
                              pyname = "Boundary3d",
                              docext = " (3D).") 