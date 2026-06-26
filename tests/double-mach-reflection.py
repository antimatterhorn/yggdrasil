import numpy as np
from math import sqrt
from yggdrasil import *
from Animation import *
from Mesh import Grid2d
from Physics import GridHydroKT2d
from EOS import IdealGasEOS
from Boundaries import ReflectingGridBoundary2d, OutflowGridBoundary2d
from Utilities import SiloDump

# Woodward-Colella double Mach reflection benchmark (JCP 1984).
# A Mach 10 shock inclined at 60° to the horizontal strikes a reflecting wedge.

if __name__ == "__main__":
    commandLine = CommandLineArguments(
        animate    = False,
        siloDump   = True,
        dumpCycle  = 50,
        cycles     = 3000,
        nx         = 400,
        ny         = 100,
        dx         = 0.01,
        dy         = 0.01,
        dtmin      = 1e-5,
        intVerbose = False)

    myGrid     = Grid2d(nx, ny, dx, dy)
    myNodeList = NodeList(nx * ny)

    constants = MKS()
    eos       = IdealGasEOS(1.4, constants)
    gamma     = eos.gamma

    hydro = GridHydroKT2d(myNodeList, constants, eos, myGrid)

    # Bottom wall only: reflecting (models the 30-degree wedge surface).
    wall = ReflectingGridBoundary2d(grid=myGrid)
    wall.setFaces(["bottom"])
    hydro.addBoundary(wall)

    # Right wall only: outflow (pre-shock gas exits without reflection).
    outflow = OutflowGridBoundary2d(grid=myGrid)
    outflow.setFaces(["right"])
    hydro.addBoundary(outflow)

    integrator = RungeKutta4Integrator2d([hydro], dtmin=dtmin, verbose=intVerbose)

    density  = myNodeList.getFieldDouble("density")
    energy   = myNodeList.getFieldDouble("specificInternalEnergy")
    velocity = myNodeList.getFieldVector2d("velocity")

    # Pre-shock state.  rho1=1.4 is chosen so that c1 = sqrt(gamma*p1/rho1) = 1
    # exactly, making the shock Mach number = 10 = shock speed.
    rho1 = 1.4
    p1   = 1.0
    u1   = p1 / ((gamma - 1.0) * rho1)

    # Post-shock state from Rankine-Hugoniot at Mach 10, gamma=1.4.
    # The shock normal is n = (sqrt(3)/2, -1/2), so the post-shock velocity
    # is 8.25 in the normal direction, decomposed into lab-frame (vx2, vy2).
    rho2 = 8.0
    p2   = 116.5
    u2   = p2 / ((gamma - 1.0) * rho2)
    vx2  = 4.125 * sqrt(3.0)   # 8.25 * cos(30 deg)
    vy2  = -4.125               # -8.25 * sin(30 deg)

    inv_sqrt3 = 1.0 / sqrt(3.0)
    x0        = 1.0 / 6.0      # x-position of the shock foot at t=0

    # Set initial conditions: post-shock left of x = x0 + y/sqrt(3), pre-shock right.
    for j in range(ny):
        for i in range(nx):
            idx = myGrid.index(i, j, 0)
            x   = (i + 0.5) * dx
            y   = (j + 0.5) * dy
            if x < x0 + y * inv_sqrt3:
                density.setValue(idx, rho2)
                velocity.setValue(idx, Vector2d(vx2, vy2))
                energy.setValue(idx, u2)
            else:
                density.setValue(idx, rho1)
                velocity.setValue(idx, Vector2d(0.0, 0.0))
                energy.setValue(idx, u1)

    # Speed at which the shock sweeps along x at any fixed y.
    # From the shock plane equation n*r = d0 + vs*t with vs=10, n_x=sqrt(3)/2:
    # dx/dt|_y = vs / n_x = 10 / (sqrt(3)/2) = 20/sqrt(3).
    vs_x = 20.0 * inv_sqrt3

    class DMRBoundary:
        """
        Handles the three time-varying inflow boundaries each cycle.
        The right wall is handled statically by OutflowGridBoundary.
        """
        def __init__(self):
            self.cycle = 1
            self._rho  = myNodeList.getFieldDouble("density")
            self._u    = myNodeList.getFieldDouble("specificInternalEnergy")
            self._v    = myNodeList.getFieldVector2d("velocity")

        def _post(self, idx):
            self._rho.setValue(idx, rho2)
            self._v.setValue(idx, Vector2d(vx2, vy2))
            self._u.setValue(idx, u2)

        def _pre(self, idx):
            self._rho.setValue(idx, rho1)
            self._v.setValue(idx, Vector2d(0.0, 0.0))
            self._u.setValue(idx, u1)

        def __call__(self, cycle, time, dt):
            x_foot = x0 + vs_x * time
            x_top  = x0 + ny * dy * inv_sqrt3 + vs_x * time

            # Left wall: steady post-shock inflow.
            for j in range(ny):
                self._post(myGrid.index(0, j, 0))

            # Bottom wall: post-shock inflow left of the moving shock foot;
            # ReflectingGridBoundary handles the wedge to the right.
            for i in range(nx):
                if (i + 0.5) * dx < x_foot:
                    self._post(myGrid.index(i, 0, 0))

            # Top wall: post-shock left of the shock, pre-shock right.
            for i in range(nx):
                idx = myGrid.index(i, ny - 1, 0)
                if (i + 0.5) * dx < x_top:
                    self._post(idx)
                else:
                    self._pre(idx)

    dmr_bc = DMRBoundary()

    periodicWork = [dmr_bc]

    if siloDump:
        meshWriter = SiloDump(
            baseName   = "DMR",
            nodeList   = myNodeList,
            fieldNames = ["density", "specificInternalEnergy", "pressure", "velocity"],
            dumpCycle  = dumpCycle)
        periodicWork += [meshWriter]

    controller = Controller(integrator=integrator, periodicWork=periodicWork, statStep=10)

    if animate:
        title = MakeTitle(controller, "time", "time")
        bounds = (nx, ny)
        update_method = AnimationUpdateMethod2d(
            call      = hydro.getCell2d,
            stepper   = controller.Step,
            title     = title,
            fieldName = "density")
        AnimateGrid2d(bounds, update_method, extremis=[1, 22], frames=cycles, cmap="plasma")
    else:
        controller.Step(cycles)
