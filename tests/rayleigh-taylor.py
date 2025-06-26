from yggdrasil import *
import matplotlib.pyplot as plt
from Animation import *
from Mesh import Grid2d
from Physics import GridHydroKT2d, GridHydroHLLE2d, ConstantGridAccel2d
from EOS import IdealGasEOS
from Boundaries import ReflectingGridBoundary2d,DirichletGridBoundary2d

from matplotlib.colors import LinearSegmentedColormap
colors = [(1,0,0), (0, 0, 0), (0,0,1)]  # Red -> Black -> Blue
cmap = LinearSegmentedColormap.from_list('rbbl', colors, N=256)

if __name__ == "__main__":
    commandLine = CommandLineArguments(animate = True,
                                        siloDump = False,
                                        cycles = 3000,
                                        nx = 100,
                                        ny = 100,
                                        dx = 1,
                                        dy = 1,
                                        g  = -980,
                                        dtmin = 1e-6,
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
    p_top = p0
    gamma = eos.gamma

    p = np.zeros((nx, ny))
    rho = np.zeros((nx, ny))

    # Set up the interface shape
    interface_y = np.zeros(nx)
    for i in range(nx):
        x = (i + 0.5) * dx
        interface_y[i] = ny / 2 + 4.0 * np.sin(9.0 * np.pi * x / nx + np.pi)

    for i in range(nx):
        for j in reversed(range(ny)):
            y = (j + 0.5) * dy
            idx = myGrid.index(i, j, 0)

            if y >= interface_y[i]:
                rho_ij = 5.0
            else:
                rho_ij = 1.0

            rho[i, j] = rho_ij

            if j == ny - 1:
                p[i, j] = p_top
            else:
                p[i, j] = p[i, j+1] - rho[i, j+1] * g * dy  # integrate dp = -ρg dy

            # Now set fields
            density.setValue(idx, rho[i, j])
            u = p[i, j] / ((gamma - 1.0) * rho[i, j])
            energy.setValue(idx, u)

    periodicWork = []

    if siloDump:
        meshWriter = SiloDump(baseName="Rayleigh-Taylor",
                                nodeList=myNodeList,
                                fieldNames=["density","specificInternalEnergy","pressure","velocity","acceleration"],
                                dumpCycle=50)
        periodicWork += [meshWriter]

    controller = Controller(integrator=integrator,periodicWork=periodicWork,statStep=50)

    if(animate):
        title = MakeTitle(controller,"time","time")

        bounds = (nx,ny)
        update_method = AnimationUpdateMethod2d(call=hydro.getCell2d,
                                                stepper=controller.Step,
                                                title=title,
                                                fieldName="density")
        AnimateGrid2d(bounds,update_method,extremis=[1,3],frames=cycles,cmap=cmap)
    else:
        controller.Step(cycles)
