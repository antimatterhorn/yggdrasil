from yggdrasil import *
from Mesh import Grid2d
from Calculators import Mandelbrot

import numpy as np
import matplotlib.pyplot as plt

commandLine = CommandLineArguments(ne = 200)

nx = ny = ne
dx = dy = 2./(nx//2)

numNodes = nx*ny

nodeList = NodeList(numNodes)
myGrid = Grid2d(nx,ny,dx,dy)

origin = Vector2d(2,2)
myGrid.setOrigin(origin)

mand   = Mandelbrot(nodeList,myGrid)
mand.compute()

mp     = nodeList.getFieldDouble("mandelbrot")

# Convert to 2D array for plotting
data = np.array([mp[i] for i in range(numNodes)])
data2D = data.reshape((ny, nx))

plt.imshow(data2D, extent=[-2, 2, -2, 2], origin='lower', cmap='inferno')
plt.title("Smoothed Mandelbrot Escape Time")
plt.colorbar(label="Smoothed Iteration Count")
plt.tight_layout()
plt.show()