from yggdrasil import *
import numpy as np
from Animation import *
from Physics import GridHydroHLLE2d
from Mesh import Grid2d
from EOS import IdealGasEOS
from Boundaries import ReflectingGridBoundary2d

# Sedov-Taylor point blast on a 2D Cartesian grid (which is cylindrically
# symmetric, d = 2), verified against the self-similar solution. A concentrated
# energy deposit drives a blast whose shock radius must follow the Sedov law
#   R(t) proportional to t^(2/(d+2)) = t^(1/2)   for d = 2.
# The exponent is coefficient-free, so it is the robust analytical check (a
# planar solver would give 2/3, a spherical one 2/5). Total energy, conserved
# by the reflecting-walled finite-volume scheme, is checked as an invariant.
# (The peak density stays below the ideal (g+1)/(g-1) jump because the thin
# shell is under-resolved on a coarse grid, so it is not asserted tightly.)

GAMMA = 1.4
D_GEOM = 2                                  # cylindrical symmetry of a 2D blast
SEDOV_EXPONENT = 2.0 / (D_GEOM + 2.0)       # = 0.5


def measure(grid, density, velocity, energy, nx, ny, x0, y0, rho_ambient):
    """Shock radius (outermost cell above 1.5*ambient), peak density, and total
    energy sum(rho*(u+0.5 v^2)*cellVolume)."""
    r_shock = 0.0
    rho_peak = 0.0
    E_tot = 0.0
    for j in range(ny):
        for i in range(nx):
            idx = grid.index(i, j, 0)
            rho = density[idx]
            if rho > 1.5 * rho_ambient:
                r_shock = max(r_shock, ((i - x0) ** 2 + (j - y0) ** 2) ** 0.5)
            rho_peak = max(rho_peak, rho)
            v = velocity[idx]
            E_tot += rho * (energy[idx] + 0.5 * (v.x * v.x + v.y * v.y)) * grid.cellVolume(idx)
    return r_shock, rho_peak, E_tot


if __name__ == "__main__":
    commandLine = CommandLineArguments(animate = False,
                                       cyclesA = 300,      # first measurement
                                       cyclesB = 600,      # second measurement
                                       nx = 100,
                                       ny = 100,
                                       dx = 1.0,
                                       dy = 1.0,
                                       dtmin = 1e-7)

    myGrid = Grid2d(nx, ny, dx, dy)
    myNodeList = NodeList(nx * ny)
    constants = MKS()
    eos = IdealGasEOS(GAMMA, constants)
    hydro = GridHydroHLLE2d(myNodeList, constants, eos, myGrid)
    hydro.addBoundary(ReflectingGridBoundary2d(grid=myGrid))
    integrator = RungeKutta4Integrator2d([hydro], dtmin=dtmin, verbose=False)

    density = myNodeList.getFieldDouble("density")
    energy = myNodeList.getFieldDouble("specificInternalEnergy")
    velocity = myNodeList.getFieldVector2d("velocity")

    x0, y0 = nx // 2, ny // 2
    sigma = 1.5
    e0 = 1000.0
    rho_ambient = 1.0
    for j in range(ny):
        for i in range(nx):
            idx = myGrid.index(i, j, 0)
            r2 = (i - x0) ** 2 + (j - y0) ** 2
            energy.setValue(idx, e0 * np.exp(-r2 / (2.0 * sigma ** 2)) + 1e-3)
            density.setValue(idx, rho_ambient)

    controller = Controller(integrator=integrator, periodicWork=[], statStep=200)

    if animate:
        title = MakeTitle(controller, "time", "time")
        update_method = AnimationUpdateMethod2d(call=hydro.getCell2d, stepper=controller.Step,
                                                title=title, fieldName="density")
        AnimateGrid2d((nx, ny), update_method, extremis=[0, 4], frames=cyclesB, cmap="plasma")
    else:
        E0 = measure(myGrid, density, velocity, energy, nx, ny, x0, y0, rho_ambient)[2]

        controller.Step(cyclesA)
        t1 = controller.time
        r1, peak1, _ = measure(myGrid, density, velocity, energy, nx, ny, x0, y0, rho_ambient)

        controller.Step(cyclesB - cyclesA)
        t2 = controller.time
        r2, peak2, E2 = measure(myGrid, density, velocity, energy, nx, ny, x0, y0, rho_ambient)

        exponent = np.log(r2 / r1) / np.log(t2 / t1)
        dE = abs(E2 - E0) / E0

        print(f"Sedov: R({t1:.3f})={r1:.2f}, R({t2:.3f})={r2:.2f} -> exponent={exponent:.3f} "
              f"(exact {SEDOV_EXPONENT:.3f}); peak rho={peak2:.2f}; dE/E={dE:.2e}")
        assert r2 < 0.45 * nx, "shock reached the domain boundary; measurement invalid"
        assert peak2 > 1.5 * rho_ambient, "no shock present"
        assert abs(exponent - SEDOV_EXPONENT) < 0.07, \
            f"Sedov radius scaling off: {exponent:.3f} vs {SEDOV_EXPONENT:.3f}"
        assert dE < 2.0e-2, f"energy not conserved: dE/E={dE:.2e}"
        print("[ok] Sedov blast follows the self-similar R ~ t^(1/2) law and conserves energy.")
