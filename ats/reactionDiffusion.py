import sys
sys.path.append("../build/src/Utilities")

from yggdrasil import *
from Animation import *
import numpy as np
from Mesh import Grid2d
from Physics import ReactionDiffusion
from Boundaries import PeriodicGridBoundary2d

import random

def initialize_reaction_diffusion(rd: ReactionDiffusion, nodeList, grid):
    # Seed the RNG for reproducibility
    random.seed(42)
    
    c1 = nodeList.getFieldDouble("c1")
    c2 = nodeList.getFieldDouble("c2")
    c3 = nodeList.getFieldDouble("c3")

    sx = grid.nx
    sy = grid.ny

    for j in range(sy):
        for i in range(sx):
            idx = grid.index(i, j, 0)
            choice = random.randint(0, 2)
            if choice == 0:
                c1.setValue(idx, 1.0)
            elif choice == 1:
                c2.setValue(idx, 1.0)
            else:
                c3.setValue(idx, 1.0)

def run():

    animate = False
    nx = 100
    ny = 100
    A = 0.8
    D = 0.4
    cycles = 100

    constants = MKS()
    numNodes = nx*ny
    myGrid = Grid2d(nx, ny,1,1)
    myNodes = NodeList(numNodes)
    rps = ReactionDiffusion(myNodes,constants,grid=myGrid,A=A,D=D)

    pm = PeriodicGridBoundary2d(grid=myGrid)
    rps.addBoundary(pm)

    initialize_reaction_diffusion(rps,myNodes,myGrid)

    integrator = Integrator2d(packages=[rps],dtmin=0.05,verbose=False)

    controller = Controller(integrator=integrator,
                            statStep=200,
                            periodicWork=[])

    # ------------------------------------------------
    # Step
    # ------------------------------------------------
    if(animate):
        title = MakeTitle(controller,"time","time")

        bounds = (nx,ny)
        update_method = AnimationUpdateMethod2d(call=rps.getCell,
                                                stepper=controller.Step,
                                                title=title)
        AnimateGrid2d(bounds,update_method,threeColors=True)
    else:
        out = []
        controller.Step(cycles)
        out.append(integrator.time)
        rgb_grid = []
        max_values = []
        for j in range(ny):
            maxi = 0
            for i in range(nx):
                k = j*ny + i
                rgb_grid.append(rps.getCell(i % nx, j % ny, "maxphi"))
                if i == nx // 2:
                    maxi = float(rgb_grid[k][0])
            max_values.append(maxi)
        return max_values

if __name__ == "__main__":
    print(run())
