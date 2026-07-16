# Copyright (C) 2026  Cody Raskin

from PYB11Generator import *

@PYB11template("dim")
class ALEMesh:
    def pyinit(self):
        return
    def addNode(self,position="Lin::Vector<%(dim)s>"):
        return
    def setNodePosition(self,nodeId="size_t",position="Lin::Vector<%(dim)s>"):
        return
    def addCell(self,nodeIndices="std::vector<size_t>"):
        return "size_t"
    def getNodes(self):
        return "std::vector<Lin::Vector<%(dim)s>>&"
    def getConnectivityMap(self):
        return "std::vector<std::vector<size_t>>&"
    def computeNeighbors(self):
        return
    def identifyBoundaryNodes(self):
        return
    def getNeighbors(self,nodeId="size_t"):
        return "std::vector<size_t>"
    def getBoundaryNodes(self):
        return "std::vector<size_t>&"
    def computeConnectivityMap(self):
        return
    def computeFaces(self):
        return
    def updateGeometry(self):
        return
    def getFaces(self):
        return "std::vector<Mesh::Face<%(dim)s>>&"
    def cellVolume(self,cellIndex="size_t"):
        return "double"
    def cellCentroid(self,cellIndex="size_t"):
        return "Lin::Vector<%(dim)s>"
    def writeVTK(self, filepath="std::string"):
        return "void"

ALEMesh2d = PYB11TemplateClass(ALEMesh,
                              template_parameters = ("2"),
                              cppname = "Mesh::ALEMesh<2>",
                              pyname = "ALEMesh2d",
                              docext = " (2D).")
