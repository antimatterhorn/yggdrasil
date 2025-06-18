import sys
sys.path.append("../build/src/Utilities")

from yggdrasil import *
from Animation import *
from Mesh import Grid2d
from Physics import ReactionDiffusion
from Boundaries import PeriodicGridBoundary2d

import random

def initialize_reaction_diffusion(rd: ReactionDiffusion, nodeList, grid):
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

if __name__ == "__main__":

    commandLine = CommandLineArguments(animate = True,
                                        nx = 100,
                                        ny = 100,
                                        A = 0.8,
                                        D = 0.4)


    constants = MKS()
    numNodes = nx*ny
    myGrid = Grid2d(nx, ny,1,1)
    myNodes = NodeList(numNodes)
    print("%i x %i grid"%(nx,ny))
    rps = ReactionDiffusion(myNodes,constants,grid=myGrid,A=A,D=D)

    pm = PeriodicGridBoundary2d(grid=myGrid)
    rps.addBoundary(pm)

    initialize_reaction_diffusion(rps,myNodes,myGrid)

    integrator = Integrator2d(packages=[rps],dtmin=0.05,verbose=False)

    print(integrator)
    print("numNodes =",myNodes.numNodes)
    print("field names =",myNodes.fieldNames)

    controller = Controller(integrator=integrator,
                            statStep=10,
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
        controller.Step(5000)
