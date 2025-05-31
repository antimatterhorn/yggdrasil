from PYB11Generator import *

@PYB11template("dim")
class VoronoiMesh:
    def pyinit(self, seedPoints="std::vector<Lin::Vector<%(dim)s>>&"):
        return

    
VoronoiMesh2d = PYB11TemplateClass(VoronoiMesh,
                              template_parameters = ("2"),
                              cppname = "Mesh::VoronoiMesh<2>",
                              pyname = "VoronoiMesh2d",
                              docext = " (2D).")
VoronoiMesh3d = PYB11TemplateClass(VoronoiMesh,
                              template_parameters = ("3"),
                              cppname = "Mesh::VoronoiMesh<3>",
                              pyname = "VoronoiMesh3d",
                              docext = " (3D).")

# no 1D instantiation because that doesn't make sense anyway
