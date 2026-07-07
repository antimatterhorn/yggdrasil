from yggdrasil import *
from Animation import *
from Mesh import Grid2d
from Physics import GridHydroKT2d,GridHydroHLLC2d,GridHydroHLLE2d
from EOS import IdealGasEOS
from Boundaries import PeriodicGridBoundary2d
from Utilities import SiloDump
from math import log

if __name__ == "__main__":
    commandLine = CommandLineArguments(animate = True,
                                        siloDump = False,
                                        cycles = 3000,
                                        nx = 100,
                                        ny = 100,
                                        dx = 1,
                                        dy = 1,
                                        dtmin = 0.1e-7,
                                        intVerbose = False)

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

    loc = [int(nx//2),int(ny//2)]

    a = 2.0/log(2)

    for j in range(ny):
        for i in range(nx):
            idx = myGrid.index(i, j, 0)
            pos = position[idx]
            x = pos.x
            y = pos.y

            r = np.sqrt((i-loc[0])**2 + (j-loc[1])**2)

            if r <= 2.0:
                rho = 10
                vx = 5
            else:
                rho = 20*np.exp(-r/a)
                if rho < 1: 
                    rho = 1
                    vx = 0
                    vy = 0
                else:
                    vx = 5
                    vy = 5

            velocity.setValue(idx, Vector2d(vx, vy))
            density.setValue(idx, rho)
            energy.setValue(idx, p0 / ((gamma - 1.0) * rho))

    periodicWork = []

    if siloDump:
        meshWriter = SiloDump(baseName="HLL",
                                nodeList=myNodeList,
                                fieldNames=["density","specificInternalEnergy","pressure","velocity"],
                                dumpCycle=50)
        periodicWork += [meshWriter]

    controller = Controller(integrator=integrator,periodicWork=periodicWork,statStep=20)

    if(animate):
        title = MakeTitle(controller,"time","time")

        bounds = (nx,ny)
        # call = lambda i, j, fieldName: hydro.getCell2dComponent(i, j, 0, "velocity"),
        update_method = AnimationUpdateMethod2d(call = hydro.getCell2d,
                                                stepper=controller.Step,
                                                title=title,
                                                fieldName="density")
        AnimateGrid2d(bounds,update_method,extremis=[0,5],frames=cycles,cmap="plasma")
    else:
        controller.Step(cycles)
