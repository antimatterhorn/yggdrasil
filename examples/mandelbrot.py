from yggdrasil import *
from Mesh import Grid2d
from Calculators import Mandelbrot

import numpy as np
import matplotlib.pyplot as plt

# Command-line arguments
commandLine = CommandLineArguments(ne = 1000,
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
nodeList = NodeList(numNodes)
myGrid = Grid2d(nx, ny, dx, dy)

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
mp = nodeList.getFieldDouble("mandelbrot")
print("Generating Plot...")
data = np.array([mp[i] for i in range(numNodes)])
data2D = data.reshape((ny, nx))  # NOTE: row = y, col = x

# Plot using user-specified domain extents
plt.imshow(data2D, extent=[xmin, xmax, ymin, ymax], origin='lower', cmap='inferno')
plt.title("Smoothed Mandelbrot Escape Time")
plt.colorbar(label="Smoothed Iteration Count")
plt.tight_layout()
plt.show()
