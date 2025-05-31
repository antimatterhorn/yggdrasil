import random
from yggdrasil import *
import matplotlib.pyplot as plt

if __name__ == "__main__":
    commandLine = CommandLineArguments(numNodes = 5,
                                       bmin = [-1,-1],
                                       bmax = [1,1])
    
    myNodeList = NodeList(numNodes)
    myNodeList.insertFieldVector2d("position")

    print("field names =",myNodeList.fieldNames)

    pos = myNodeList.getFieldVector2d("position")
    for i in range(numNodes):
        pos.setValue(i,Vector2d(random.random() * (bmax[0]-bmin[0]) + bmin[0],
                  random.random() * (bmax[1]-bmin[1]) + bmin[1]))

    print("Creating VoronoiMesh2d")
    vor = VoronoiMesh2d(pos)
    print(len(vor.getCells()),"cells")

    for cell in vor.getCells():
        gen = np.array([cell.generator.x, cell.generator.y])
        verts = np.array([[v.x, v.y] for v in cell.vertices])
        centroid = np.mean(verts, axis=0)
        if np.linalg.norm(gen - centroid) > 1.0:
            print(f"WARNING: generator {gen} far from cell centroid {centroid}")

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