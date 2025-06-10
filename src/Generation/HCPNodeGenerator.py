from yggdrasil import *

class HCPNodeGenerator2d:
    def __init__(self, nx,ny):
        self.positions = []

        # Step 1: Generate initial seed using a Lattice
        from LatticeNodeGenerator import Lattice2d
        self.positions = Lattice2d(nx,ny).positions

        # Step 2: Alternate shifting rows left or right by dx/2
        dx = 1.0/nx
        for i in range(len(self.positions)):
            j = i//nx
            ddx = dx/2.0*((-1)**j)
            self.positions[i] = (self.positions[i][0]+ddx,self.positions[i][1])