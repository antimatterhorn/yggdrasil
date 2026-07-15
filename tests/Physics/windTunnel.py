import os
from yggdrasil import *
from Animation import *
from Mesh import Grid2d
from Physics import GridHydroKT2d, GridHydroHLLC2d
from EOS import IdealGasEOS
from Boundaries import DirichletGridBoundary2d, ReflectingGridBoundary2d, OutflowGridBoundary2d
from Utilities import SiloDump

if __name__ == "__main__":
    commandLine = CommandLineArguments(animate = True,
                                        siloDump = False,
                                        cycles = 3000,
                                        nx = 260,
                                        ny = 160,
                                        dx = 1.0,
                                        dy = 1.0,
                                        dtmin = 1e-7,
                                        intVerbose = False,
                                        dumpCycle = 25,
                                        rootName = "windTunnel",
                                        vizDir = "viz")

    vizDir = rootName + "/" + vizDir
    os.makedirs(vizDir, exist_ok=True)

    myGrid = Grid2d(nx, ny, dx, dy)
    print("grid size:", myGrid.size())

    myNodeList = NodeList(nx * ny)

    constants = MKS()
    eos = IdealGasEOS(1.4, constants)
    gamma = eos.gamma

    # GridHydroKT2d: the Kurganov-Tadmor central scheme's extra numerical
    # dissipation is what keeps the immediate wake of the obstacle stable.
    # GridHydroHLLC2d's exact-contact-preserving second-order reconstruction
    # is fragile right at a masked (non-body-fitted) obstacle boundary --
    # it develops a runaway velocity spike in the trailing-edge wake cells
    # regardless of Mach number (confirmed from Mach ~0.08 to ~0.3), because
    # the obstacle is represented as flat masked cells rather than a proper
    # body-fitted or ghost-fluid mirrored boundary. Swap this line if you
    # want to see it yourself.
    hydro = GridHydroKT2d(myNodeList, constants, eos, myGrid)
    # hydro = GridHydroHLLC2d(myNodeList, constants, eos, myGrid)

    # ------------------------------------------------------------------
    # Freestream conditions (subsonic, Mach ~0.50)
    # ------------------------------------------------------------------
    rhoInf = 1.0
    pInf   = 1.0
    vxInf  = 0.6
    uInf   = pInf / ((gamma - 1.0) * rhoInf)

    density  = myNodeList.getFieldDouble("density")
    energy   = myNodeList.getFieldDouble("specificInternalEnergy")
    velocity = myNodeList.getFieldVector2d("velocity")

    obstacleCenter = Vector2d(nx * dx * 0.3, ny * dy * 0.5)
    obstacleRadius = 0.05 * ny * dy

    for j in range(ny):
        for i in range(nx):
            idx = myGrid.index(i, j, 0)
            density.setValue(idx, rhoInf)
            energy.setValue(idx, uInf)

            # Start at rest inside the obstacle footprint so the boundary's
            # interior fill (velocity pinned to zero, see ReflectingGridBoundary)
            # isn't an instantaneous full-magnitude discontinuity on step one --
            # that "materializing a solid disk into a moving fluid" transient
            # is itself enough to destabilize the flow for a few cycles.
            px, py = i * dx + 0.5 * dx, j * dy + 0.5 * dy
            inObstacle = (px - obstacleCenter.x) ** 2 + (py - obstacleCenter.y) ** 2 <= obstacleRadius ** 2
            velocity.setValue(idx, Vector2d(0.0, 0.0) if inObstacle else Vector2d(vxInf, 0.0))

    # ------------------------------------------------------------------
    # Boundaries:
    #   - left face:              Dirichlet, pinned to freestream (inflow)
    #   - right/top/bottom faces: Outflow (linear extrapolation) -- an open
    #                             field rather than a confined channel, so
    #                             flow can freely leave through the sides
    #                             instead of reflecting off them
    #   - circular obstacle:      Reflecting interior obstacle, on its own
    #                             boundary object with face mode disabled
    #                             (setFaces([])) so it contributes only the
    #                             interior fill, not a second face treatment
    #                             on top of Outflow's. GridHydroBase excludes
    #                             obstacle cells from its own flux loop (see
    #                             ZeroTimeInitialize) and leaves them for
    #                             this boundary's Apply to fill each step.
    # ------------------------------------------------------------------
    inflow = DirichletGridBoundary2d(grid=myGrid)
    inflow.addBox(Vector2d(-10.0, -10.0), Vector2d(dx, ny * dy + 10.0))
    hydro.addBoundary(inflow)

    outflow = OutflowGridBoundary2d(grid=myGrid)
    outflow.setFaces(["right", "top", "bottom"])
    hydro.addBoundary(outflow)

    obstacle = ReflectingGridBoundary2d(grid=myGrid)
    obstacle.setFaces([])
    obstacle.addSphere(obstacleCenter, obstacleRadius)
    hydro.addBoundary(obstacle)

    integrator = RungeKutta4Integrator2d([hydro], dtmin=dtmin, verbose=intVerbose)

    periodicWork = []
    if siloDump:
        meshWriter = SiloDump(baseName=os.path.join(vizDir, rootName),
                               nodeList=myNodeList,
                               fieldNames=["density", "specificInternalEnergy", "pressure", "velocity"],
                               dumpCycle=dumpCycle,
                               grid=myGrid)
        periodicWork += [meshWriter]

    controller = Controller(integrator=integrator, periodicWork=periodicWork, statStep=20)

    if animate:
        title = MakeTitle(controller, "time", "time")
        bounds = (nx, ny)
        update_method = AnimationUpdateMethod2d(call=hydro.getCell2d,
                                                 stepper=controller.Step,
                                                 title=title,
                                                 fieldName="density")
        AnimateGrid2d(bounds, update_method, extremis=[0, 2], frames=cycles, cmap="plasma")
    else:
        controller.Step(cycles)
