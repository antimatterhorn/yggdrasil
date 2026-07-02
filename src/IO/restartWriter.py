# Copyright (C) 2026  Cody Raskin

from PYB11Generator import *

@PYB11template("dim")
class RestartWriter:
    def pyinit(self,
               nodeList="NodeList&",
               integrator="Integrator<%(dim)s>&"):
        return
    def write(self, fileName="std::string&"):
        return

RestartWriter1d = PYB11TemplateClass(RestartWriter,
                              template_parameters = ("1"),
                              cppname = "RestartWriter<1>",
                              pyname = "RestartWriter1d",
                              docext = " (1D).")
RestartWriter2d = PYB11TemplateClass(RestartWriter,
                              template_parameters = ("2"),
                              cppname = "RestartWriter<2>",
                              pyname = "RestartWriter2d",
                              docext = " (2D).")
RestartWriter3d = PYB11TemplateClass(RestartWriter,
                              template_parameters = ("3"),
                              cppname = "RestartWriter<3>",
                              pyname = "RestartWriter3d",
                              docext = " (3D).")
