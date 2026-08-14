# Copyright (C) 2026  Cody Raskin

from PYB11Generator import *

# Expose the coordinate-system enum from C++ (namespace Mesh, not templated).
Geometry = PYB11enum(
    name="Geometry",
    values=[
        ("Cartesian",     "Mesh::Geometry::Cartesian"),
        ("CylindricalRZ", "Mesh::Geometry::CylindricalRZ"),
    ],
    namespace="Mesh"
)

@PYB11template("dim")
class Grid:
    def pyinit1(self,nx="int",sx="double"):
        return
    def pyinit2(self,nx="int",ny="int",sx="double",sy="double",
                geometry=("Mesh::Geometry","Mesh::Geometry::Cartesian")):
        return
    def pyinit3(self,nx="int",ny="int",nz="int",sx="double",sy="double",sz="double"):
        return
    @PYB11pyname("position")
    def getPosition(self,id="int"):
        return
    def setOrigin(self,origin="Lin::Vector<%(dim)s>"):
        return
    def size(self):
        return
    def logicalSize(self):
        return
    def getNeighboringCells(self,idx="int"):
        return "std::array<int, 2*%(dim)s>"
    def index(self,i="int",j=("int",0),k=("int",0)):
        return
    def indexToCoordinates(self,idx="int"):
        return "std::array<int, %(dim)s>"
    def assignPositions(self,nodeList="NodeList*"):
        return
    def lowMost(self,axis="int"):
        return "std::vector<int>"
    def highMost(self,axis="int"):
        return "std::vector<int>"
    def isGhost(self,idx="int"):
        return
    def isInterior(self,idx="int"):
        return
    def geometry(self):
        return "Mesh::Geometry"
    def cellVolume(self,idx="int"):
        return "double"
    def cellRadius(self,idx="int"):
        return "double"

    nx = PYB11property("int", getter="getnx", doc="The number of x coords.")
    ny = PYB11property("int", getter="getny", doc="The number of y coords.")
    nz = PYB11property("int", getter="getnz", doc="The number of z coords.")
    dx = PYB11property("double", getter="getdx", doc="The number of x coords.")
    dy = PYB11property("double", getter="getdy", doc="The number of y coords.")
    dz = PYB11property("double", getter="getdz", doc="The number of z coords.")

Grid1d = PYB11TemplateClass(Grid,
                              template_parameters = ("1"),
                              cppname = "Mesh::Grid<1>",
                              pyname = "Grid1d",
                              docext = " (1D).")
Grid2d = PYB11TemplateClass(Grid,
                              template_parameters = ("2"),
                              cppname = "Mesh::Grid<2>",
                              pyname = "Grid2d",
                              docext = " (2D).")
Grid3d = PYB11TemplateClass(Grid,
                              template_parameters = ("3"),
                              cppname = "Mesh::Grid<3>",
                              pyname = "Grid3d",
                              docext = " (3D).") 