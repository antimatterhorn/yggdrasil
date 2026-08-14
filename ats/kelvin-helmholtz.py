import numpy as np
from yggdrasil import *

# Kelvin-Helmholtz shear instability. There is no simple closed-form solution,
# but the defining behavior is that a seeded perturbation *grows*: the initially
# flat, x-uniform density interfaces (density depends only on y at t=0) roll up,
# so the spread of density along x at the interface rows climbs from zero to a
# clear finite amplitude. This checks that growth (and that it stays bounded),
# rather than snapshotting a particular late-time state.


def _interface_spread(nodes, grid, nx, ny):
    density = nodes.getFieldDouble("density")
    spread = 0.0
    for j in (ny // 4, 3 * ny // 4):                 # the two shear interfaces
        row = [density[grid.index(i, j, 0)] for i in range(nx)]
        spread += float(np.std(row))
    return 0.5 * spread


def run():
    nx, ny, dx, dy, tstop = 100, 50, 0.02, 0.02, 0.3
    grid = Grid2d(nx, ny, dx, dy)
    nodes = NodeList(grid.size())
    constants = MKS()
    eos = IdealGasEOS(1.4, constants)
    hydro = GridHydroKT2d(nodes, constants, eos, grid)
    boundary = PeriodicGridBoundary2d(grid=grid)
    hydro.addBoundary(boundary)
    integrator = RungeKutta4Integrator2d([hydro], dtmin=1.0e-7, verbose=False)

    density = nodes.getFieldDouble("density")
    energy = nodes.getFieldDouble("specificInternalEnergy")
    velocity = nodes.getFieldVector2d("velocity")
    position = nodes.getFieldVector2d("position")
    p0, gamma = 2.5, eos.gamma
    for j in range(ny):
        for i in range(nx):
            idx = grid.index(i, j, 0)
            x = position[idx].x
            if j < ny // 4 or j >= 3 * ny // 4:
                rho, vx = 1.0, -0.5
            else:
                rho, vx = 2.0, 0.5
            velocity[idx] = Vector2d(vx, 0.1 * np.sin(4 * np.pi * x))
            density[idx] = rho
            energy[idx] = p0 / ((gamma - 1.0) * rho)

    spread0 = _interface_spread(nodes, grid, nx, ny)     # ~0: interfaces flat

    Controller(integrator=integrator, periodicWork=[], statStep=100000, tstop=tstop).Step(1000000)
    spread1 = _interface_spread(nodes, grid, nx, ny)

    return {"mode": "analytic", "checks": [
        ("interfaces start flat", bool(spread0 < 1e-9), f"initial x-spread={spread0:.2e}"),
        ("instability grows (interface rolls up)", bool(spread1 > 0.05),
         f"interface x-spread {spread0:.2e} -> {spread1:.3f}"),
        ("stays bounded", bool(np.isfinite(spread1) and spread1 < 10.0),
         f"final spread={spread1:.3f} finite"),
    ]}


if __name__ == "__main__":
    print(run())
