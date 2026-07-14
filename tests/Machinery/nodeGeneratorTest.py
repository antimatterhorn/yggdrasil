

import numpy as np
import matplotlib.pyplot as plt
from yggdrasil import *
from typing import Any,Callable

commandLine = CommandLineArguments(method = "Fibonacci",
                                   n = 500)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

def RPRFunc(n):
    from RPRPSNodeGenerator import RecursivePrimitiveRefinementSurface3d
    return RecursivePrimitiveRefinementSurface3d(n)

def PSFunc(n):
    from RPRPSNodeGenerator import ParameterizedSpiralSurface3d
    return ParameterizedSpiralSurface3d(n)

def FibonacciFunc(n):
    from FibonacciNodeGenerator import FibonacciSurface3d
    return FibonacciSurface3d(n)

def SEAGenFunc(n):
    from SEANodeGenerator import SEAGenSurface3d
    return SEAGenSurface3d(n)

Func = Callable[[int], Any]

getPoints: dict[str,Func] ={
    "RPR": RPRFunc,
    "PS": PSFunc,
    "Fibonacci": FibonacciFunc,
    "SEAGen": SEAGenFunc
}

points = getPoints[method](n).positions



from scipy.spatial import distance_matrix

points = np.array(points)

# Compute pairwise distances on the sphere
dists = distance_matrix(points, points)

# Define a Gaussian-like weight function
sigma = 0.1  # Tune this for best results
weights = np.exp(-dists**2 / (2 * sigma**2))

# Normalize so that the sum of weights approximates the total area (4π for a unit sphere)
local_density = weights.sum(axis=1)
As = (4 * np.pi) / local_density    # Approximate area per point
As = As/As.mean()                   # Normalize areas

print("Standard Deviation from Equal Area: %3.4f"%(np.std(As)))

# Extract coordinates
xs = points[:, 0]
ys = points[:, 1]
zs = points[:, 2]

# 3D Scatter plot
ax.scatter(xs, ys, zs, c=As, cmap="viridis")
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
ax.set_title(method)

# Create a separate 2D scatter plot of lattitude vs. area
# This is with knowledge that the worst portions of any tiling method are at the poles
fig2 = plt.figure()
ax2 = fig2.add_subplot()
ax2.scatter(zs, As, c=As, cmap="viridis")
ax2.set_xlabel("Z Coordinate")
ax2.set_ylabel("Area per Point (As)")
ax2.set_title("Z Coordinate vs Area")

plt.show()