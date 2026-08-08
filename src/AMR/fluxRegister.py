# Copyright (C) 2026  Cody Raskin

from PYB11Generator import *
# Only FluxObserver itself, not its 1d/2d/3d instantiations, which would collide with Physics' own registration.
from fluxObserver import FluxObserver

@PYB11template("dim")
class FluxRegister(FluxObserver):
    def pyinit(self, numCells="int"):
        return
    def registerFace(self, cell="int", axis="int", plusSide="bool"):
        "Start accumulating the face on `plusSide` of `cell` along `axis`; returns its slot."
        return "int"
    def reset(self):
        return
    def slotCount(self):
        return "int"

@PYB11template("dim")
class Refluxer:
    def pyinit(self,
               coarseNodes="NodeList*",
               coarseGrid="Mesh::Grid<%(dim)s>*",
               eos="EquationOfState*",
               coarseRegister="AMR::FluxRegister<%(dim)s>*",
               fineRegister="AMR::FluxRegister<%(dim)s>*",
               coarseSlots="std::vector<int>",
               fineSlots="std::vector<int>"):
        return
    def apply(self, dt="double"):
        "Correct the coarse cells along the interface for one coarse step, then clear both registers."
        return

FluxRegister1d = PYB11TemplateClass(FluxRegister,
                              template_parameters = ("1"),
                              cppname = "AMR::FluxRegister<1>",
                              pyname = "FluxRegister1d",
                              docext = " (1D).")
FluxRegister2d = PYB11TemplateClass(FluxRegister,
                              template_parameters = ("2"),
                              cppname = "AMR::FluxRegister<2>",
                              pyname = "FluxRegister2d",
                              docext = " (2D).")
FluxRegister3d = PYB11TemplateClass(FluxRegister,
                              template_parameters = ("3"),
                              cppname = "AMR::FluxRegister<3>",
                              pyname = "FluxRegister3d",
                              docext = " (3D).")

Refluxer1d = PYB11TemplateClass(Refluxer,
                              template_parameters = ("1"),
                              cppname = "AMR::Refluxer<1>",
                              pyname = "Refluxer1d",
                              docext = " (1D).")
Refluxer2d = PYB11TemplateClass(Refluxer,
                              template_parameters = ("2"),
                              cppname = "AMR::Refluxer<2>",
                              pyname = "Refluxer2d",
                              docext = " (2D).")
Refluxer3d = PYB11TemplateClass(Refluxer,
                              template_parameters = ("3"),
                              cppname = "AMR::Refluxer<3>",
                              pyname = "Refluxer3d",
                              docext = " (3D).")
