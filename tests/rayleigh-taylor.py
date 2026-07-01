from yggdrasil import *
import matplotlib.pyplot as plt
from Animation import *
from Mesh import Grid2d
from Physics import GridHydroKT2d, GridHydroHLLE2d, ConstantGravity2d
from EOS import IdealGasEOS
from Boundaries import ReflectingGridBoundary2d,DirichletGridBoundary2d
from Utilities import SiloDump

if __name__ == "__main__":
    commandLine = CommandLineArguments(animate = False,
                                        siloDump = True,
                                        dumpCycle = 50,
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
    gravity  = ConstantGravity2d(myNodeList,constants,gravityVector)

    integrator = RungeKutta4Integrator2d([hydro,gravity],dtmin=dtmin,verbose=intVerbose)

    density = myNodeList.getFieldDouble("density")
    energy  = myNodeList.getFieldDouble("specificInternalEnergy")
    velocity = myNodeList.getFieldVector2d("velocity")
    position = myNodeList.getFieldVector2d("position")

    # p0 must be large enough that cs >> v_RT everywhere.
    # v_RT ~ sqrt(At * |g| * H) ~ sqrt(2/3 * 980 * 100) ~ 255 cm/s.
    # With rho_heavy=5 at the top, cs_top = sqrt(gamma*p0/rho_heavy).
    # p0 = 1e6 gives cs_top ~ 529 cm/s (Ma ~ 0.5, marginal).
    # p0 = 1e7 gives cs_top ~ 1673 cm/s (Ma ~ 0.15, safely subsonic).
    p0    = 1.0e7
    gamma = eos.gamma

    rho = np.zeros((nx, ny))

    # Interface shape: sinusoidal perturbation centered at mid-domain
    interface_y = np.zeros(nx)
    for i in range(nx):
        x = (i + 0.5) * dx
        interface_y[i] = ny / 2 + 1.0 * np.sin(3.0 * np.pi * x / nx + np.pi)

    # First pass: assign density from the interface (heavy above, light below)
    for i in range(nx):
        for j in range(ny):
            y = (j + 0.5) * dy
            rho[i, j] = 5.0 if y >= interface_y[i] else 1.0

    # Build a HORIZONTALLY UNIFORM 1D hydrostatic pressure by integrating
    # the row-averaged density from the top downward.  Using a column-local
    # density would create spurious horizontal pressure gradients in the IC
    # that drive immediate horizontal flow before the instability can develop.
    p_1d = np.zeros(ny)
    p_1d[ny - 1] = p0
    for j in range(ny - 2, -1, -1):
        rho_row = np.mean(rho[:, j + 1])
        p_1d[j] = p_1d[j + 1] - rho_row * g * dy   # g < 0, so pressure increases downward

    # Second pass: set density and energy fields
    for i in range(nx):
        for j in range(ny):
            idx = myGrid.index(i, j, 0)
            density.setValue(idx, rho[i, j])
            energy.setValue(idx, p_1d[j] / ((gamma - 1.0) * rho[i, j]))

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
        AnimateGrid2d(bounds,update_method,extremis=[1,5],frames=cycles,cmap='RdBu_r')
    else:
        controller.Step(cycles)
