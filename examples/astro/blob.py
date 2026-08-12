import os
from yggdrasil import *
from Animation import *
from Mesh import Grid2d
from Physics import GridHydroKT2d,GridHydroHLLC2d,GridHydroHLLE2d
from EOS import IdealGasEOS
from Boundaries import PeriodicGridBoundary2d, OutflowGridBoundary2d
from Utilities import SiloDump
from math import log

if __name__ == "__main__":
    commandLine = CommandLineArguments(animate = True,
                                        siloDump = False,
                                        cycles = 3000,
                                        nx = 200,
                                        ny = 100,
                                        dx = 1,
                                        dy = 1,
                                        dtmin = 0.1e-7,
                                        intVerbose = False,
                                        dumpCycle= 50,
                                        rootName = "blob",
                                        vizDir = "viz")
    vizDir = rootName + "/" + vizDir
    if siloDump:
        os.makedirs(vizDir, exist_ok=True)

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

    pb = PeriodicGridBoundary2d(grid=myGrid)
    hydro.addBoundary(pb)

    integrator = RungeKutta4Integrator2d([hydro],dtmin=dtmin,verbose=intVerbose)

    density = myNodeList.getFieldDouble("density")
    energy  = myNodeList.getFieldDouble("specificInternalEnergy")
    velocity = myNodeList.getFieldVector2d("velocity")
    position = myNodeList.getFieldVector2d("position")

    p0 = 2.5
    gamma = eos.gamma

    blobloc = [int(3*nx//4),int(ny//2)]
    shockfr = nx-2

    a = 2.0/log(2)

    for j in range(ny):
        for i in range(nx):
            idx = myGrid.index(i, j, 0)
            pos = position[idx]
            x = pos.x
            y = pos.y

            r = np.sqrt((i-blobloc[0])**2 + (j-blobloc[1])**2)

            vx = vy = 0

            if r <= 2.0:
                rho = 10
            else:
                rho = max(20*np.exp(-r/a),1)

            if x > shockfr:
                rho = 4
                vx  = -10

            velocity.setValue(idx, Vector2d(vx, vy))
            density.setValue(idx, rho)
            energy.setValue(idx, p0 / ((gamma - 1.0) * rho))

    periodicWork = []

    if siloDump:
        meshWriter = SiloDump(baseName=os.path.join(vizDir, rootName),
                                nodeList=myNodeList,
                                fieldNames=["density","specificInternalEnergy","pressure","velocity"],
                                dumpCycle=dumpCycle,
                                grid=myGrid)
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
        AnimateGrid2d(bounds,update_method,extremis=[0,5],frames=cycles,cmap="plasma",lineout_axis='y', lineout_index=ny // 2)
    else:
        controller.Step(cycles)
