from PYB11Generator import *

@PYB11template("dim")
class Integrator:
    def pyinit(self,
               packages="std::vector<Physics<%(dim)s>*>",
               dtmin="double"):
        return
    def pyinit1(self,
               packages="std::vector<Physics<%(dim)s>*>",
               dtmin="double", verbose="bool"):
        return
    def Step(self):
        return
    def Time(self):
        return
    def Cycle(self):
        return
    def getPackages(self):
        return "std::vector<Physics<%(dim)s>*>"
    
    dt = PYB11property("double", getter="Dt", doc="timestep")
    time = PYB11property("double", getter="Time", doc="The time.")
    cycle = PYB11property("double", getter="Cycle", doc="The cycle.")
    
Integrator1d = PYB11TemplateClass(Integrator,
                              template_parameters = ("1"),
                              cppname = "Integrator<1>",
                              pyname = "Integrator1d",
                              docext = " (1D).")
Integrator2d = PYB11TemplateClass(Integrator,
                              template_parameters = ("2"),
                              cppname = "Integrator<2>",
                              pyname = "Integrator2d",
                              docext = " (2D).")
Integrator3d = PYB11TemplateClass(Integrator,
                              template_parameters = ("3"),
                              cppname = "Integrator<3>",
                              pyname = "Integrator3d",
                              docext = " (3D).") 