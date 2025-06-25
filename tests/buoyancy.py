from yggdrasil import *
from Animation import *
from Mesh import Grid2d
from Physics import GridHydroKT2d,GridHydroHLLC2d,GridHydroHLLE2d,ConstantGridAccel2d
from EOS import IdealGasEOS
from Boundaries import ReflectingGridBoundary2d
from Utilities import SiloDump
from math import log

if __name__ == "__main__":
    commandLine = CommandLineArguments(animate = True,
                                        siloDump = False,
                                        dumpCycle = 100,
                                        cycles = 6000,
                                        nx = 100,
                                        ny = 100,
                                        dx = 1,
                                        dy = 1,
                                        dtmin = 0.1e-7,
                                        g  = -10,
                                        rho1 = 1,
                                        rho2 = 5,
                                        intVerbose = False)

    myGrid = Grid2d(nx,ny,dx,dy)
    print("grid size:",myGrid.size())
    
    myNodeList = NodeList(nx*ny)
    print("numNodes =",myNodeList.numNodes)
    print("field names =",myNodeList.fieldNames)

    constants = CGS()
    eos = IdealGasEOS(1.4,constants)
    print(eos,"gamma =",eos.gamma)

    hydro = GridHydroKT2d(myNodeList,constants,eos,myGrid) 
    #hydro = GridHydroHLLE2d(myNodeList,constants,eos,myGrid) 
    print("numNodes =",myNodeList.numNodes)
    print("field names =",myNodeList.fieldNames)

    box = ReflectingGridBoundary2d(grid=myGrid)
    hydro.addBoundary(box)

    gravityVector = Vector2d(0.,g)
    gravity  = ConstantGridAccel2d(myNodeList,constants,gravityVector)

    integrator = RungeKutta4Integrator2d([hydro,gravity],dtmin=dtmin,verbose=intVerbose)

    density = myNodeList.getFieldDouble("density")
    energy  = myNodeList.getFieldDouble("specificInternalEnergy")
    velocity = myNodeList.getFieldVector2d("velocity")
    position = myNodeList.getFieldVector2d("position")

    p0 = 2.5
    gamma = eos.gamma

    loc = [int(nx//2),int(ny//10)]

    p0 = 2.5
    p_top = p0

    p = np.zeros((nx, ny))
    rho = np.zeros((nx, ny))

    for i in range(nx):
        for j in reversed(range(ny)):
            idx = myGrid.index(i, j, 0)
            pos = position[idx]
            x = pos.x
            y = pos.y

            r = np.sqrt((i-loc[0])**2 + (j-loc[1])**2)

            if r > (nx//20):
                rhoij = rho2
            else:
                rhoij = rho1

            rho[i, j] = rhoij

            if j == ny - 1:
                p[i, j] = p_top
            else:
                p[i, j] = p[i, j+1] - rho[i, j+1] * g * dy

            density.setValue(idx, rhoij)
            u = p[i, j] / ((gamma - 1.0) * rho[i, j])
            energy.setValue(idx, u)

    periodicWork = []

    if siloDump:
        meshWriter = SiloDump(baseName="buoyancy",
                                nodeList=myNodeList,
                                fieldNames=["density","specificInternalEnergy","pressure","velocity"],
                                dumpCycle=dumpCycle)
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
        AnimateGrid2d(bounds,update_method,extremis=[0,1.1*rho2],frames=cycles,cmap="YlOrRd")
    else:
        controller.Step(cycles)
