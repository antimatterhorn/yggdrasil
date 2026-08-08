# Copyright (C) 2026  Cody Raskin

from PYB11Generator import *

@PYB11template("dim")
class FluxObserver:
    "Abstract tap on a finite-volume solver's face fluxes; implemented by AMR::FluxRegister."

FluxObserver1d = PYB11TemplateClass(FluxObserver,
                              template_parameters = ("1"),
                              cppname = "FluxObserver<1>",
                              pyname = "FluxObserver1d",
                              docext = " (1D).")
FluxObserver2d = PYB11TemplateClass(FluxObserver,
                              template_parameters = ("2"),
                              cppname = "FluxObserver<2>",
                              pyname = "FluxObserver2d",
                              docext = " (2D).")
FluxObserver3d = PYB11TemplateClass(FluxObserver,
                              template_parameters = ("3"),
                              cppname = "FluxObserver<3>",
                              pyname = "FluxObserver3d",
                              docext = " (3D).")
