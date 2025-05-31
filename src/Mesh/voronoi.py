from PYB11Generator import *

@PYB11template("dim")
class VoronoiCell:
    def pyinit(self, generator = "Lin::Vector<%(dim)s>", vertices = "std::vector<Lin::Vector<%(dim)s>>"):
        return
    generator = PYB11property("Lin::Vector<%(dim)s>", getter="getGenerator")
    vertices = PYB11property("std::vector<Lin::Vector<%(dim)s>>&", getter="getVertices")

VoronoiCell2d = PYB11TemplateClass(VoronoiCell,
                              template_parameters = ("2"),
                              cppname = "Mesh::VoronoiCell<2>",
                              pyname = "VoronoiCell2d",
                              docext = " (2D).")
VoronoiCell3d = PYB11TemplateClass(VoronoiCell,
                              template_parameters = ("3"),
                              cppname = "Mesh::VoronoiCell<3>",
                              pyname = "VoronoiCell3d",
                              docext = " (3D).")

@PYB11template("dim")
class VoronoiMesh:

    def pyinit(self, seedPoints="Field<Lin::Vector<%(dim)s>>&"):
        return
    def getCells(self):
        return "const std::vector<typename Mesh::VoronoiCell<%(dim)s>>&"


# VoronoiMesh1d = PYB11TemplateClass(VoronoiMesh,
#                               template_parameters = ("1"),
#                               cppname = "Mesh::VoronoiMesh<1>",
#                               pyname = "VoronoiMesh1d",
#                               docext = " (1D).")    
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
