import random
from yggdrasil import *
import matplotlib.pyplot as plt

if __name__ == "__main__":
    commandLine = CommandLineArguments(numNodes = 100,
                                       bmin = [-4,-4],
                                       bmax = [4,4],
                                       method = "random")
    
    assert method in ["random", "fibonacci"]

    myNodeList = NodeList(numNodes)
    myNodeList.insertFieldVector2d("position")

    print("field names =",myNodeList.fieldNames)

    pos = myNodeList.getFieldVector2d("position")
    if method == "random":
        for i in range(numNodes):
            pos.setValue(i,Vector2d(random.uniform(bmin[0],bmax[0]), random.uniform(bmin[1],bmax[1])))
    else:
        from FibonacciNodeGenerator import FibonacciDisk2d
        posF = FibonacciDisk2d(numNodes).positions
        for i in range(numNodes):
            pos.setValue(i,Vector2d(posF[i][0], posF[i][1])*0.5*(bmax[1]-bmin[1]))

    print("Creating VoronoiMesh2d")
    vor = VoronoiMesh2d(pos)
    print(len(vor.getCells()),"cells")

    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    import matplotlib.collections as mc

    fig, ax = plt.subplots(figsize=(6, 6))

    for cell in vor.getCells():
        verts = [(v.x, v.y) for v in cell.vertices]
        if len(verts) >= 3:
            poly = Polygon(verts, edgecolor='black', facecolor='lightblue', alpha=0.4)
            ax.add_patch(poly)

        # Optionally mark the generator point
        gen = cell.generator
        ax.plot(gen.x, gen.y, 'ro', markersize=3)

    ax.set_xlim(bmin[0], bmax[0])
    ax.set_ylim(bmin[1], bmax[1])
    ax.set_aspect('equal')
    plt.title("Voronoi Mesh")
    plt.show()