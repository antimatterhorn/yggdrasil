from yggdrasil import *
import matplotlib.pyplot as plt
from Mesh import *
from numpy import sqrt

if __name__ == "__main__":
    commandLine = CommandLineArguments(numNodes = 100,
                                       bmin = [-4,-4],
                                       bmax = [4,4],
                                       method = "fibonacci")
    
    assert method in ["random", "fibonacci", "glass", 
                        "constantDTheta", "poisson", 
                        "glassDisk", "poissonDisk", "lattice",
                        "cvt", "cvtDisk","hcp"]


    if method == "random":
        from RandomNodeGenerator import RandomNodeGenerator2d
        posF = RandomNodeGenerator2d(numNodes).positions
    elif method == "fibonacci":
        from FibonacciNodeGenerator import FibonacciDisk2d
        posF = FibonacciDisk2d(numNodes).positions
    elif method == "glass":
        from GlassNodeGenerator import GlassNodeGenerator2d
        posF = GlassNodeGenerator2d(numNodes).positions
        numNodes = len(posF) 
    elif method == "constantDTheta":
        from ConstantDThetaGenerator import ConstantDThetaDisk2d
        posF = ConstantDThetaDisk2d(numNodes).positions
        numNodes = len(posF) # needs to be fixed because constantDTheta may overfill
    elif method == "poissonDisk":
        from PoissonNodeGenerator import PoissonDisk2d
        posF = PoissonDisk2d(numNodes).positions
        numNodes = len(posF) 
    elif method == "poisson":
        from PoissonNodeGenerator import PoissonNodeGenerator2d
        posF = PoissonNodeGenerator2d(numNodes).positions
        numNodes = len(posF) 
    elif method == "cvt":
        from CVTNodeGenerator import CVTNodeGenerator2d
        posF = CVTNodeGenerator2d(numNodes).positions
        numNodes = len(posF)
    elif method == "hcp":
        from HCPNodeGenerator import HCPNodeGenerator2d
        posF = HCPNodeGenerator2d(int(sqrt(numNodes)),int(sqrt(numNodes))).positions
        numNodes = len(posF) 
    elif method == "cvtDisk":
        from CVTNodeGenerator import CVTDiskGenerator2d
        posF = CVTDiskGenerator2d(numNodes).positions
        numNodes = len(posF) 
    elif method == "glassDisk":
        from GlassNodeGenerator import GlassDisk2d
        posF = GlassDisk2d(numNodes).positions
        numNodes = len(posF)
    elif method == "lattice":
        from LatticeNodeGenerator import Lattice2d
        posF = Lattice2d(int(sqrt(numNodes)),int(sqrt(numNodes))).positions

    myNodeList = NodeList(numNodes)
    myNodeList.insertFieldVector2d("position")

    print("field names =",myNodeList.fieldNames)

    pos = myNodeList.getFieldVector2d("position")

    for i in range(numNodes):
        pos.setValue(i,Vector2d(posF[i][0], posF[i][1])*0.5*(bmax[1]-bmin[1]))

    print("Creating VoronoiMesh2d")
    vor = VoronoiMesh2d(pos)
    print(len(vor.getCells()),"cells")



    # Compute cell areas
    cells = vor.getCells()
    areas = [cell.area for cell in cells]
    mean_area = sum(areas) / len(areas)
    # Or to center at mean:

    femesh = vor.generateDualFEMesh()
    print("Created FEMesh")

    femesh.writeVTK("fibonacci.vtk")
    # print("FEMesh has", len(femesh.getNodes()), "nodes")
    # print("FEMesh has", len(femesh.getElements()), "elements")
