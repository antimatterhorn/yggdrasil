from yggdrasil import *
import matplotlib.pyplot as plt
from Mesh import VoronoiMesh2d
from numpy import sqrt

if __name__ == "__main__":
    commandLine = CommandLineArguments(numNodes = 100,
                                       bmin = [-4,-4],
                                       bmax = [4,4],
                                       method = "random")
    
    assert method in ["random", "fibonacci", "glass", 
                        "constantDTheta", "poisson", 
                        "glassDisk", "poissonDisk", "lattice",
                        "cvt", "cvtDisk"]


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

    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    import matplotlib.collections as mc

    from matplotlib import colormaps
    from matplotlib.colors import Normalize
    from matplotlib.patches import Polygon
    from matplotlib.colors import TwoSlopeNorm


    fig, ax = plt.subplots(figsize=(6, 6))

    # Compute cell areas
    cells = vor.getCells()
    areas = [cell.area for cell in cells]
    mean_area = sum(areas) / len(areas)
    # Or to center at mean:
    norm = TwoSlopeNorm(vmin=min(areas), vcenter=mean_area, vmax=max(areas))

    # norm = Normalize(vmin=min(areas), vmax=max(areas))
    cmap = colormaps["bwr"]

    for cell, area in zip(cells, areas):
        verts = [(v.x, v.y) for v in cell.vertices]
        if len(verts) >= 3:
            color = cmap(norm(area))
            poly = Polygon(verts, edgecolor='black', facecolor=color, alpha=0.7)
            ax.add_patch(poly)

        # Optionally mark the generator point
        gen = cell.generator
        ax.plot(gen.x, gen.y, 'ro', markersize=3)

    ax.set_xlim(bmin[0], bmax[0])
    ax.set_ylim(bmin[1], bmax[1])
    ax.set_aspect('equal')
    plt.title("Voronoi Mesh Colored by Area")

    # Add a colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # Required for matplotlib < 3.1
    plt.colorbar(sm, ax=ax, label="Cell Area")

    plt.show()
