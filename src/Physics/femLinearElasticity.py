from PYB11Generator import *
from physics import *

@PYB11template("dim")
class FEMLinearElasticity(Physics):
    def pyinit(self,
               nodeList="NodeList*",
               constants="PhysicalConstants&",
               mesh="Mesh::FEMesh<%(dim)s>*",
               material="FEMMaterial<%(dim)s>*",
               rho="double",
               bodyForce="Lin::Vector<%(dim)s>",
               alpha="double",
               beta="double"):
        return

FEMLinearElasticity2d = PYB11TemplateClass(FEMLinearElasticity,
                                           template_parameters=("2"),
                                           cppname="FEMLinearElasticity<2>",
                                           pyname="FEMLinearElasticity2d",
                                           docext=" (2D).")
FEMLinearElasticity3d = PYB11TemplateClass(FEMLinearElasticity,
                                           template_parameters=("3"),
                                           cppname="FEMLinearElasticity<3>",
                                           pyname="FEMLinearElasticity3d",
                                           docext=" (3D).")
