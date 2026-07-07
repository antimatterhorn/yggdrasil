import numpy as np
from yggdrasil import *


def run():
    siloDump = False
    nx = 100
    ny = 50
    dx = 0.02
    dy = 0.02
    dtmin = 1.0e-7
    tstop = 0.3

    myGrid = Grid2d(nx,ny,dx,dy)

    myNodeList = NodeList(nx*ny)

    constants = MKS()
    eos = IdealGasEOS(1.4,constants)

    hydro = GridHydroKT2d(myNodeList,constants,eos,myGrid)

    box = PeriodicGridBoundary2d(grid=myGrid)
    hydro.addBoundary(box)

    integrator = RungeKutta4Integrator2d([hydro],dtmin=dtmin,verbose=False)

    density  = myNodeList.getFieldDouble("density")
    energy   = myNodeList.getFieldDouble("specificInternalEnergy")
    velocity = myNodeList.getFieldVector2d("velocity")
    position = myNodeList.getFieldVector2d("position")

    p0 = 2.5
    gamma = eos.gamma

    for j in range(ny):
        for i in range(nx):
            idx = myGrid.index(i,j,0)
            pos = position[idx]
            x = pos.x

            if j < ny // 4:
                rho = 1.0
                vx = -0.5
            elif j < 3 * ny // 4:
                rho = 2.0
                vx = 0.5
            else:
                rho = 1.0
                vx = -0.5

            vy = 0.1 * np.sin(4 * np.pi * x)

            velocity.setValue(idx, Vector2d(vx, vy))
            density.setValue(idx, rho)
            energy.setValue(idx, p0 / ((gamma - 1.0) * rho))

    periodicWork = []

    if siloDump:
        meshWriter = SiloDump(baseName="KH",
                                nodeList=myNodeList,
                                fieldNames=["density","specificInternalEnergy","pressure","velocity"],
                                dumpCycle=50)
        periodicWork += [meshWriter]

    controller = Controller(integrator=integrator,periodicWork=periodicWork,statStep=10000,tstop=tstop)

    controller.Step(1000000)

    # vertical lineout of density nearest x=0.5
    i_lineout = int(round(0.5/dx - 0.5))

    ys = []
    ys.append(integrator.time)
    for j in range(ny):
        idx = myGrid.index(i_lineout, j, 0)
        ys.append(density[idx])
    return ys

if __name__ == "__main__":
    print(run())
