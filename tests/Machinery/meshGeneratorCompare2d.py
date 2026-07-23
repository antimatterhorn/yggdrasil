from yggdrasil import *
import numpy as np
from Utilities import MeshViz
from CircleQuadPolyGenerator import CircleQuadPolyGenerator2d

if __name__ == "__main__":
    commandLine = CommandLineArguments(nx = 20, ny = 20, k=1.0,
                                       angleTol = 25.0)

    gen = CircleQuadPolyGenerator2d(nx, ny, k=k)
    mesh_viz = MeshViz(gen)
    mesh_viz.plot()