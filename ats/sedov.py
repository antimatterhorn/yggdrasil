from yggdrasil import *
import numpy as np

# Sedov-Taylor point blast on a 2D Cartesian grid (cylindrically symmetric,
# d=2), verified against the self-similar solution: the shock radius follows
# R(t) ~ t^(2/(d+2)) = t^(1/2). The exponent is coefficient-free (a planar
# solver would give 2/3, a spherical one 2/5), so it is the robust analytical
# check. Total energy is conserved by the reflecting-walled FV scheme.

GAMMA = 1.4
SEDOV_EXPONENT = 0.5


def _measure(grid, density, velocity, energy, nx, ny, x0, y0, amb):
    r_shock = rho_peak = E = 0.0
    for j in range(ny):
        for i in range(nx):
            idx = grid.index(i, j, 0)
            rho = density[idx]
            if rho > 1.5 * amb:
                r_shock = max(r_shock, ((i - x0) ** 2 + (j - y0) ** 2) ** 0.5)
            rho_peak = max(rho_peak, rho)
            v = velocity[idx]
            E += rho * (energy[idx] + 0.5 * (v.x * v.x + v.y * v.y)) * grid.cellVolume(idx)
    return r_shock, rho_peak, E


def run():
    nx, ny = 100, 100
    cyclesA, cyclesB = 300, 600
    grid = Grid2d(nx, ny, 1.0, 1.0)
    nodes = NodeList(nx * ny)
    constants = MKS()
    eos = IdealGasEOS(GAMMA, constants)
    hydro = GridHydroHLLE2d(nodes, constants, eos, grid)
    boundary = ReflectingGridBoundary2d(grid=grid)
    hydro.addBoundary(boundary)
    integrator = RungeKutta4Integrator2d([hydro], dtmin=1e-7, verbose=False)

    density = nodes.getFieldDouble("density")
    energy = nodes.getFieldDouble("specificInternalEnergy")
    velocity = nodes.getFieldVector2d("velocity")
    x0, y0, amb = nx // 2, ny // 2, 1.0
    for j in range(ny):
        for i in range(nx):
            idx = grid.index(i, j, 0)
            r2 = (i - x0) ** 2 + (j - y0) ** 2
            energy[idx] = 1000.0 * np.exp(-r2 / (2.0 * 1.5 ** 2)) + 1e-3
            density[idx] = amb

    controller = Controller(integrator=integrator, periodicWork=[], statStep=100000)
    E0 = _measure(grid, density, velocity, energy, nx, ny, x0, y0, amb)[2]
    controller.Step(cyclesA)
    t1 = integrator.time
    r1, _, _ = _measure(grid, density, velocity, energy, nx, ny, x0, y0, amb)
    controller.Step(cyclesB - cyclesA)
    t2 = integrator.time
    r2, peak, E2 = _measure(grid, density, velocity, energy, nx, ny, x0, y0, amb)

    exponent = float(np.log(r2 / r1) / np.log(t2 / t1))
    dE = abs(E2 - E0) / E0

    return {"mode": "analytic", "checks": [
        ("self-similar radius exponent", bool(abs(exponent - SEDOV_EXPONENT) < 0.07),
         f"{exponent:.3f} vs exact {SEDOV_EXPONENT:.3f}"),
        ("shock present, interior", bool(peak > 1.5 * amb and r2 < 0.45 * nx),
         f"peak rho={peak:.2f}, R={r2:.1f} < {0.45*nx:.0f}"),
        ("energy conserved", bool(dE < 2.0e-2), f"dE/E={dE:.2e} < 2e-2"),
    ]}


if __name__ == "__main__":
    print(run())
