import sys
import numpy as np
import matplotlib.pyplot as plt
from yggdrasil import *
from scipy.spatial import distance_matrix
from progressBar import ProgressBar

def stDev(method):
    ns = []
    sd = []
    for n in range(10,500,10):
        p = method(n).positions
        p = np.array(p)

        # Compute pairwise distances on the sphere
        dists = distance_matrix(p, p)

        # Define a Gaussian-like weight function
        sigma = 0.1  # Tune this for best results
        weights = np.exp(-dists**2 / (2 * sigma**2))

        # Normalize so that the sum of weights approximates the total area (4π for a unit sphere)
        local_density = weights.sum(axis=1)
        As = (4 * np.pi) / local_density    # Approximate area per point
        As = As/As.mean()                   # Normalize areas

        ns.append(len(p))
        sd.append(np.std(As))
        ProgressBar((n+10)/500,method.__name__)
    print("\n")
    return ns,sd

fig = plt.figure()
ax = fig.add_subplot()


from FibonacciNodeGenerator import FibonacciDisk2d
method = FibonacciDisk2d

xs,ys = stDev(method)
scatter3 = ax.plot(xs, ys, c="blue", label="Fibonacci")

from GlassNodeGenerator import GlassNodeGenerator2d
method  = GlassNodeGenerator2d

xs,ys = stDev(method)
scatter4 = ax.plot(xs, ys, c="black", label="Glass")

from ConstantDThetaGenerator import ConstantDThetaDisk2d
method  = ConstantDThetaDisk2d

xs,ys = stDev(method)
scatter4 = ax.plot(xs, ys, c="green", label="Constant DTheta")

from PoissonNodeGenerator import PoissonDisk2d
method   = PoissonDisk2d

xs,ys = stDev(method)
scatter4 = ax.plot(xs, ys, c="red", label="Poisson")

ax.set_xlabel("number of points")
ax.set_ylabel("standard deviation")

# Add legend
ax.legend()

plt.show()