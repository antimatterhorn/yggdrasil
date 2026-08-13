import os
import numpy as np
from yggdrasil import *
from Animation import *
from Mesh import Grid2d
from Physics import GridHydroKT2d
from EOS import IdealGasEOS
from Boundaries import PeriodicGridBoundary2d
from Utilities import SiloDump
from IO import RestartWriter

if __name__ == "__main__":
    commandLine = CommandLineArguments(animate = False,
                                        siloDump = True,
                                        statStep = 10,
                                        dumpCycle= 50,
                                        cycles=9999,
                                        tstop = 2.0,
                                        res = 200,
                                        dtmin = 0.1e-7,
                                        intVerbose = False,
                                        restartCycle = 250,
                                        rootName = "kelvin-helmholtz",
                                        vizDir = "viz",
                                        restartDir = "restart",
                                        restoreCycle = None)

    nx = res
    ny = res // 2
    
    dx = dy = 1.0/res

    dumpDir = rootName + "/" + "nx=%d_ny=%d" % (nx, ny)
    
    vizDir = dumpDir + "/" + vizDir
    restartDir = dumpDir + "/" + restartDir
    if siloDump: 
        os.makedirs(vizDir, exist_ok=True)
    if restartCycle > 0:
        os.makedirs(restartDir, exist_ok=True)

    myGrid = Grid2d(nx,ny,dx,dy)
    print("grid size:",myGrid.size())
    
    myNodeList = NodeList(nx*ny)
    print("numNodes =",myNodeList.numNodes)
    print("field names =",myNodeList.fieldNames)

    constants = MKS()
    eos = IdealGasEOS(1.4,constants)
    print(eos,"gamma =",eos.gamma)

    hydro = GridHydroKT2d(myNodeList,constants,eos,myGrid) 
    #hydro = GridHydroHLLE2d(myNodeList,constants,eos,myGrid) 
    print("numNodes =",myNodeList.numNodes)
    print("field names =",myNodeList.fieldNames)

    box = PeriodicGridBoundary2d(grid=myGrid)
    hydro.addBoundary(box)

    integrator = RungeKutta4Integrator2d([hydro],dtmin=dtmin,verbose=intVerbose)

    density = myNodeList.getFieldDouble("density")
    energy  = myNodeList.getFieldDouble("specificInternalEnergy")
    velocity = myNodeList.getFieldVector2d("velocity")
    position = myNodeList.getFieldVector2d("position")

    p0 = 2.5
    gamma = eos.gamma

    numSeeds = 2
    Lx = nx * dx

    for j in range(ny):
        for i in range(nx):
            idx = myGrid.index(i, j, 0)
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

            vy = 0.1 * np.sin(2 * np.pi * numSeeds * x / Lx)

            velocity.setValue(idx, Vector2d(vx, vy))
            density.setValue(idx, rho)
            energy.setValue(idx, p0 / ((gamma - 1.0) * rho))

    restoreIfAvailable(myNodeList, integrator, restartDir, rootName, restoreCycle)

    periodicWork = []

    if siloDump:
        meshWriter = SiloDump(baseName=os.path.join(vizDir, rootName),
                                nodeList=myNodeList,
                                fieldNames=["density","specificInternalEnergy","pressure","velocity"],
                                dumpCycle=dumpCycle,
                                grid=myGrid)
        periodicWork += [meshWriter]

    if restartCycle > 0:
        restartWriter = RestartWriter(myNodeList, integrator)
        def dropRestart(cycle, time, dt):
            restartWriter.write(restartFileName(restartDir, rootName, cycle))
        dropRestart.cycle = restartCycle
        periodicWork += [dropRestart]

    controller = Controller(integrator=integrator,periodicWork=periodicWork,statStep=statStep,tstop=tstop)

    if(animate):
        title = MakeTitle(controller,"time","time")

        bounds = (nx,ny)
        update_method = AnimationUpdateMethod2d(call=hydro.getCell2d,
                                                stepper=controller.Step,
                                                title=title,
                                                fieldName="density")
        AnimateGrid2d(bounds,update_method,extremis=[0,5],frames=cycles,cmap="plasma")
    else:
        controller.Step(9999999) # run until tstop
