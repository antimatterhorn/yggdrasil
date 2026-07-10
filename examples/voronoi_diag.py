from yggdrasil import *
import matplotlib.pyplot as plt
from Mesh import VoronoiMesh2d
from numpy import sqrt
from typing import Any, Callable

Func = Callable[[int], Any]

def RandomFunc(n):
    from RandomNodeGenerator import RandomNodeGenerator2d
    return RandomNodeGenerator2d(n)

def FibonacciFunc(n):
    from FibonacciNodeGenerator import FibonacciDisk2d
    return FibonacciDisk2d(n)

def GlassFunc(n):
    from GlassNodeGenerator import GlassNodeGenerator2d
    return GlassNodeGenerator2d(n)

def ConstantDThetaFunc(n):
    from ConstantDThetaGenerator import ConstantDThetaDisk2d
    return ConstantDThetaDisk2d(n)

def PoissonDiskFunc(n):
    from PoissonNodeGenerator import PoissonDisk2d
    return PoissonDisk2d(n)

def PoissonFunc(n):
    from PoissonNodeGenerator import PoissonNodeGenerator2d
    return PoissonNodeGenerator2d(n)

def CVTFunc(n):
    from CVTNodeGenerator import CVTNodeGenerator2d
    return CVTNodeGenerator2d(n)

def HCPFunc(n):
    from HCPNodeGenerator import HCPNodeGenerator2d
    m = int(sqrt(n))
    return HCPNodeGenerator2d(m, m)

def CVTDiskFunc(n):
    from CVTNodeGenerator import CVTDiskGenerator2d
    return CVTDiskGenerator2d(n)

def GlassDiskFunc(n):
    from GlassNodeGenerator import GlassDisk2d
    return GlassDisk2d(n)

def LatticeFunc(n):
    from LatticeNodeGenerator import Lattice2d
    m = int(sqrt(n))
    return Lattice2d(m, m)

def TanhFunc(n):
    from TanhNodeGenerator import TanhNodeGenerator2d
    m = int(sqrt(n))
    return TanhNodeGenerator2d(m, m, beta_x=3.0, beta_y=3.0)

getPoints: dict[str, Func] = {
    "random": RandomFunc,
    "fibonacci": FibonacciFunc,
    "glass": GlassFunc,
    "constantDTheta": ConstantDThetaFunc,
    "poissonDisk": PoissonDiskFunc,
    "poisson": PoissonFunc,
    "cvt": CVTFunc,
    "hcp": HCPFunc,
    "cvtDisk": CVTDiskFunc,
    "glassDisk": GlassDiskFunc,
    "lattice": LatticeFunc,
    "tanh": TanhFunc,
}

if __name__ == "__main__":
    commandLine = CommandLineArguments(numNodes = 100,
                                       bmin = [-4,-4],
                                       bmax = [4,4],
                                       method = "fibonacci")

    assert method in getPoints, "unknown method '%s', must be one of %s" % (method, sorted(getPoints.keys()))

    posF = getPoints[method](numNodes).positions
    numNodes = len(posF)  # generators aren't guaranteed to return exactly the requested count

    myNodeList = NodeList(numNodes)
    myNodeList.insertFieldVector2d("position")

    print("field names =",myNodeList.fieldNames)

    pos = myNodeList.position

    for i in range(numNodes):
        pos[i] = Vector2d(posF[i][0], posF[i][1])*0.5*(bmax[1]-bmin[1])

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
