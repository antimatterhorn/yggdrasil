from yggdrasil import *
from Mesh import Grid2d
from Calculators import Mandelbrot

import numpy as np
import matplotlib.pyplot as plt

# Command-line arguments
commandLine = CommandLineArguments(ne = 500,
                                   xmax = -1.105945,
                                   xmin = -1.106815,
                                   ymin = 0.240362,
                                   ymax = 0.241187)

# Extract bounds and resolution
nx = ny = ne

# Compute spacing
dx = (xmax - xmin) / nx
dy = (ymax - ymin) / ny
numNodes = nx * ny

# Create grid
print("Creating Grid...")
myGrid = Grid2d(nx, ny, dx, dy)
nodeList = NodeList(myGrid.size())

# Set grid origin to upper-right corner (consistent with your system)
# so that coordinates descend as i,j increase
origin = Vector2d(-xmin,-ymin)
myGrid.setOrigin(origin)

print(xmin,ymin,xmax,ymax)

# Compute Mandelbrot set
print("Computing Mandelbrot...")
mand = Mandelbrot(nodeList, myGrid)
mand.compute()

# Extract and reshape results
mp = nodeList.mandelbrot
print("Generating Plot...")
data2D = np.zeros((ny, nx))  # NOTE: row = y, col = x
for j in range(ny):
    for i in range(nx):
        data2D[j, i] = mp[myGrid.index(i, j, 0)]

# Plot using user-specified domain extents
plt.imshow(data2D, extent=[xmin, xmax, ymin, ymax], origin='lower', cmap='inferno')
plt.title("Smoothed Mandelbrot Escape Time")
plt.colorbar(label="Smoothed Iteration Count")
plt.tight_layout()
plt.show()
