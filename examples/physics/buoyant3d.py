import os
from yggdrasil import *
from Animation import *
from Mesh import Grid3d
from Physics import GridHydroKT3d,ConstantForce3d
from EOS import IdealGasEOS
from Boundaries import ReflectingGridBoundary3d
from Utilities import SiloDump
from math import log
from IO import RestartWriter

if __name__ == "__main__":
    commandLine = CommandLineArguments(siloDump = True,
                                        dumpCycle = 100,
                                        cycles = 6000,
                                        nx = 20,
                                        ny = 20,
                                        nz = 80,
                                        dx = 1,
                                        dy = 1,
                                        dz = 1,
                                        dtmin = 0.1e-7,
                                        g  = -10,
                                        rho1 = 1,
                                        rho2 = 5,
                                        intVerbose = False,
                                        restartCycle = 50,
                                        rootName = "buoy",
                                        vizDir = "buoy/viz",
                                        restartDir = "buoy/restart",
                                        restoreCycle = None)

    os.makedirs(vizDir, exist_ok=True)
    os.makedirs(restartDir, exist_ok=True)

    myGrid = Grid3d(nx,ny,nz,dx,dy,dz)
    print("grid size:",myGrid.size())

    myNodeList = NodeList(myGrid.size())
    print("numNodes =",myNodeList.numNodes)
    print("field names =",myNodeList.fieldNames)

    constants = CGS()
    eos = IdealGasEOS(1.4,constants)
    print(eos,"gamma =",eos.gamma)

    hydro = GridHydroKT3d(myNodeList,constants,eos,myGrid) 
    #hydro = GridHydroHLLE2d(myNodeList,constants,eos,myGrid) 
    print("numNodes =",myNodeList.numNodes)
    print("field names =",myNodeList.fieldNames)

    box = ReflectingGridBoundary3d(grid=myGrid)
    hydro.addBoundary(box)

    gravityVector = Vector3d(0.,0.,g)
    gravity  = ConstantForce3d(myNodeList,constants,gravityVector)

    integrator = RungeKutta4Integrator3d([hydro,gravity],dtmin=dtmin,verbose=intVerbose)

    density = myNodeList.density
    energy  = myNodeList.specificInternalEnergy
    velocity = myNodeList.velocity
    position = myNodeList.position

    p0 = 2.5
    gamma = eos.gamma

    loc = [int(nx//2),int(ny//2),int(nz//10)]

    p0 = 2.5
    p_top = p0

    p = np.zeros((nx, ny, nz))
    rho = np.zeros((nx, ny, nz))

    for i in range(nx):
        for j in range(ny):
            for k in reversed(range(nz)):
                idx = myGrid.index(i, j, k)
                pos = position[idx]
                x = pos.x
                y = pos.y
                z = pos.z

                r = np.sqrt((i-loc[0])**2 + (j-loc[1])**2 + (z-loc[2])**2)

                if r > (nx//20):
                    rhoij = rho2
                else:
                    rhoij = rho1

                rho[i, j, k] = rhoij

                if k == nz - 1:
                    p[i, j, k] = p_top
                else:
                    p[i, j, k] = p[i, j, k+1] - rho[i, j, k+1] * g * dz

                density.setValue(idx, rhoij)
                u = p[i, j, k] / ((gamma - 1.0) * rho[i, j, k])
                energy.setValue(idx, u)

    restoreIfAvailable(myNodeList, integrator, restartDir, rootName, restoreCycle)

    periodicWork = []

    if siloDump:
        meshWriter = SiloDump(baseName=os.path.join(vizDir, rootName),
                                nodeList=myNodeList,
                                fieldNames=["density","specificInternalEnergy","pressure","velocity"],
                                dumpCycle=dumpCycle,
                                grid=myGrid)
        periodicWork += [meshWriter]

    restartWriter = RestartWriter(myNodeList, integrator)
    def dropRestart(cycle, time, dt):
        restartWriter.write(restartFileName(restartDir, rootName, cycle))
    dropRestart.cycle = restartCycle
    periodicWork += [dropRestart]

    controller = Controller(integrator=integrator,periodicWork=periodicWork,statStep=1)

    controller.Step(cycles)
