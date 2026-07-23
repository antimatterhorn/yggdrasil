from yggdrasil import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon
from scipy.spatial import Delaunay
from CircleQuadPolyGenerator import CircleQuadPolyGenerator2d

if __name__ == "__main__":
    commandLine = CommandLineArguments(nx = 20, ny = 20, k=1.0,
                                       angleTol = 25.0)

    gen = CircleQuadPolyGenerator2d(nx, ny, k=k)
    positions = gen.positions
    quads = gen.cells

    points = np.array(positions)

    fig, ax = plt.subplots(figsize=(7, 7))

    for q in quads:
        poly = MplPolygon(points[q], closed=True, fill=True,
                        facecolor="navajowhite", edgecolor="darkorange", linewidth=1.2)
        ax.add_patch(poly)

    ax.plot(points[:, 0], points[:, 1], "o", color="black", markersize=2)
    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(-1.1, 1.1)
    ax.set_aspect("equal")
    #ax.set_title(f"Quad-dominant mesh: {len(quads)} quads, {len(leftover_triangles)} triangles")
    plt.tight_layout()
    plt.show()