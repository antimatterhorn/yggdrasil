from yggdrasil import *
from Mesh import Grid2d
from Calculators import Mandelbrot
from Animation import AnimateGrid2d

import numpy as np

class AnimatedMandelbrotZoom:
    def __init__(self, nx, ny, xmin, xmax, ymin, ymax, zoom_rate=0.98):
        self.nx = nx
        self.ny = ny
        self.zoom_rate = zoom_rate
        self.center_x = (xmin + xmax) / 2
        self.center_y = (ymin + ymax) / 2
        self.width = xmax - xmin
        self.height = ymax - ymin
        self.num_nodes = nx * ny
        self.t = 0

        self.nodeList = NodeList(self.num_nodes)
        self.grid = Grid2d(nx, ny, self.width / nx, self.height / ny)
        origin = Vector2d(-xmin, -ymin)
        self.grid.setOrigin(origin)
        self.compute()

    def compute(self):
        mand = Mandelbrot(self.nodeList, self.grid)
        mand.compute()
        self.values = [self.nodeList.mandelbrot[i] for i in range(self.num_nodes)]

    def module_stepper(self):
        self.t += 1
        self.width *= self.zoom_rate
        self.height *= self.zoom_rate
        xmin = self.center_x - self.width / 2
        ymin = self.center_y - self.height / 2
        dx = self.width / self.nx
        dy = self.height / self.ny

        self.grid = Grid2d(self.nx, self.ny, dx, dy)
        origin = Vector2d(-xmin, -ymin)
        self.grid.setOrigin(origin)
        self.nodeList = NodeList(self.num_nodes)
        self.compute()

    def module_call(self, i, j, fieldName):
        idx = j * self.nx + i
        return self.values[idx]

    def module_title(self):
        return f"Zoom Step {self.t}"

    def __call__(self, i, j):
        return self.module_call(i, j, None)

# Initial bounds
xmin, xmax = -1.106815, -1.105945
ymin, ymax = 0.240362, 0.241187
ne = 300

sim = AnimatedMandelbrotZoom(ne, ne, xmin, xmax, ymin, ymax)
AnimateGrid2d((ne, ne), sim, frames=150, interval=50, cmap='inferno', extremis=(0, 100))
