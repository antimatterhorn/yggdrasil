from yggdrasil import *


def run():
    siloDump = False
    cycles = 100
    nx = 100
    ny = 20
    dx = 1
    dy = 1
    dtmin = 0.001
    myGrid = Grid2d(nx,ny,dx,dy)
    
    myNodeList = NodeList(nx*ny)

    constants = MKS()
    eos = IdealGasEOS(1.4,constants)

    hydro = GridHydroHLLE2d(myNodeList,constants,eos,myGrid)

    box = ReflectingGridBoundary2d(grid=myGrid)
    hydro.addBoundary(box)

    integrator = RungeKutta4Integrator2d([hydro],dtmin=dtmin,verbose=False)

    density = myNodeList.getFieldDouble("density")
    energy  = myNodeList.getFieldDouble("specificInternalEnergy")

    for j in range(ny):
        for i in range(nx):
            idx = myGrid.index(i,j,0)
            if i < nx // 2:
                density.setValue(idx, 1.0)
                energy.setValue(idx, 2.5)   # high pressure side
            else:
                density.setValue(idx, 0.125)
                energy.setValue(idx, 2.0)   # low pressure side

    periodicWork = []
    
    if siloDump:
        meshWriter = SiloDump(baseName="HLL",
                                nodeList=myNodeList,
                                fieldNames=["density","specificInternalEnergy","pressure","velocity"],
                                dumpCycle=50)
        periodicWork += [meshWriter]

    controller = Controller(integrator=integrator,periodicWork=periodicWork,statStep=500)

    controller.Step(cycles)

    ys = []
    ys.append(integrator.time)
    position = myNodeList.getFieldVector2d("position")
    for i in range(nx*ny):
        if position[i].y == ((ny/2.0)+(dy/2.0)):
            ys.append(density[i])
    return ys

if __name__ == "__main__":
    print(run())