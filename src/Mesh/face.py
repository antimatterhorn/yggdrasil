# Copyright (C) 2026  Cody Raskin

from PYB11Generator import *

@PYB11template("dim")
class Face:
    def pyinit(self):
        return

    nodeIndices = PYB11property("std::vector<size_t>", getter="getNodeIndices", doc="Mesh node indices shared by the two adjacent cells.")
    normal      = PYB11property("Lin::Vector<%(dim)s>", getter="getNormal", doc="Unit face normal, pointing from leftCell toward rightCell.")
    area        = PYB11property("double", getter="getArea", doc="Face area (edge length in 2D).")
    centroid    = PYB11property("Lin::Vector<%(dim)s>", getter="getCentroid", doc="Face centroid.")
    leftCell    = PYB11property("size_t", getter="getLeftCell", doc="Index of the cell the normal points away from.")
    rightCell   = PYB11property("size_t", getter="getRightCell", doc="Index of the cell the normal points toward.")

Face2d = PYB11TemplateClass(Face,
                              template_parameters = ("2"),
                              cppname = "Mesh::Face<2>",
                              pyname = "Face2d",
                              docext = " (2D).")
