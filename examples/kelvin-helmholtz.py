import os
import numpy as np
from yggdrasil import *
from Animation import *
from Mesh import Grid2d
from Physics import GridHydroKT2d
from EOS import IdealGasEOS
from Boundaries import PeriodicGridBoundary2d
from Utilities import SiloDump
from IO import RestartWriter2d, RestartReader2d

def restartFileName(restartDir, rootName, cycle):
    return os.path.join(restartDir, "%s-restart-cycle%d.ygr" % (rootName, cycle))

if __name__ == "__main__":
    commandLine = CommandLineArguments(animate = False,
                                        siloDump = True,
                                        dumpCycle= 50,
                                        cycles = 3000,
                                        nx = 200,
                                        ny = 100,
                                        dx = 0.01,
                                        dy = 0.01,
                                        dtmin = 0.1e-7,
                                        intVerbose = False,
                                        restartCycle = 250,
                                        rootName = "kelvin-helmholtz",
                                        vizDir = "KH/viz",
                                        restartDir = "KH/restart",
                                        restoreCycle = None)

    os.makedirs(vizDir, exist_ok=True)
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

    for j in range(ny):
        for i in range(nx):
            idx = myGrid.index(i, j, 0)
            pos = position[idx]
            x = pos.x
            y = pos.y

            y0 = j * dy

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

    # Pick up where a previous run left off. --restoreCycle=N picks up that
    # exact checkpoint; otherwise (the default, restoreCycle=None) we pick up
    # the highest-cycle checkpoint found in restartDir, if any. Construction
    # above is deterministic, so it's safe to always build the fresh IC
    # first and then overwrite it.
    if restoreCycle is not None:
        restoreFrom = restartFileName(restartDir, rootName, int(restoreCycle))
        if not os.path.exists(restoreFrom):
            print("--restoreCycle=%s requested but %s does not exist" %
                                     (restoreCycle, restoreFrom))
            restorFrom = None
    else:
        restoreFrom = findLatestRestart(restartDir, rootName)

    if restoreFrom is not None:
        RestartReader2d(myNodeList, integrator).read(restoreFrom)
        print("Restored from %s at cycle %d, time %g" %
              (restoreFrom, integrator.Cycle(), integrator.Time()))

    periodicWork = []

    if siloDump:
        meshWriter = SiloDump(baseName=os.path.join(vizDir, rootName),
                                nodeList=myNodeList,
                                fieldNames=["density","specificInternalEnergy","pressure","velocity"],
                                dumpCycle=dumpCycle)
        periodicWork += [meshWriter]

    # Drop a new checkpoint (named for the cycle it was taken at) every
    # restartCycle steps, so old checkpoints are never overwritten and a run
    # can be resumed from any of them via --restoreCycle.
    restartWriter = RestartWriter2d(myNodeList, integrator)
    def dropRestart(cycle, time, dt):
        restartWriter.write(restartFileName(restartDir, rootName, cycle))
    dropRestart.cycle = restartCycle
    periodicWork += [dropRestart]

    controller = Controller(integrator=integrator,periodicWork=periodicWork,statStep=1)

    if(animate):
        title = MakeTitle(controller,"time","time")

        bounds = (nx,ny)
        update_method = AnimationUpdateMethod2d(call=hydro.getCell2d,
                                                stepper=controller.Step,
                                                title=title,
                                                fieldName="density")
        AnimateGrid2d(bounds,update_method,extremis=[0,5],frames=cycles,cmap="plasma")
    else:
        controller.Step(cycles)
