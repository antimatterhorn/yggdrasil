# Copyright (C) 2026  Cody Raskin

from PYB11Generator import *

@PYB11template("dim")
class RestartReader:
    def pyinit(self,
               nodeList="NodeList&",
               integrator="Integrator<%(dim)s>&"):
        return
    def read(self, fileName="std::string&"):
        return

RestartReader1d = PYB11TemplateClass(RestartReader,
                              template_parameters = ("1"),
                              cppname = "RestartReader<1>",
                              pyname = "RestartReader1d",
                              docext = " (1D).")
RestartReader2d = PYB11TemplateClass(RestartReader,
                              template_parameters = ("2"),
                              cppname = "RestartReader<2>",
                              pyname = "RestartReader2d",
                              docext = " (2D).")
RestartReader3d = PYB11TemplateClass(RestartReader,
                              template_parameters = ("3"),
                              cppname = "RestartReader<3>",
                              pyname = "RestartReader3d",
                              docext = " (3D).")
