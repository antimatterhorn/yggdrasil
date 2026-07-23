# Copyright (C) 2026  Cody Raskin

from LinearAlgebra import *
import numpy as np
from LatticeNodeGenerator import Lattice2d

class CircleQuadPolyGenerator2d:
    def __init__(self, nx,ny,k=1.0):
        if nx %2 == 0:
            nx += 1
        if ny %2 == 0:
            ny += 1
        self.positions = []
        for i in range(nx):
            for j in range(ny):
                x = -1+ 2*i / float(nx - 1)
                y = -1+ 2*j / float(ny - 1)
                self.positions.append((x,y))

        for i in range(len(self.positions)):
            x,y = self.positions[i]
            r_euc = np.sqrt(x**2+y**2)
            r_inf = max(abs(x), abs(y))             # square radius 
            t = r_inf**k    # k > 1: more square, k < 1: more circular
            if r_euc > 0:
                scale = (1 - t) + t * (r_inf / r_euc)
                x_new = x * scale
                y_new = y * scale
            else:
                x_new, y_new = 0.0, 0.0
            self.positions[i] = (x_new, y_new)


        # Create a quad mesh from the points.
        self.cells = []
        def idx(i, j):
            return i * ny + j

        for i in range(nx - 1):
            for j in range(ny - 1):
                self.cells.append([
                    idx(i, j),
                    idx(i+1, j),
                    idx(i+1, j+1),
                    idx(i, j+1)
                ])

        self.npoints = len(self.positions)
        self.ncells = len(self.cells)