from yggdrasil import *
import numpy as np
from Animation import *
from Physics import GridHydroKT2d
from Mesh import Grid2d
from EOS import IdealGasEOS
from Boundaries import ReflectingGridBoundary2d

# Sod shock tube, verified against the exact Riemann solution. The classic
# initial state (rho,u,p) = (1,0,1) | (0.125,0,0.1) with gamma=1.4 develops a
# left rarefaction, a contact, and a right shock; the exact self-similar
# solution is sampled at the achieved time and compared to the numerical
# profile (L1 error in density and velocity).

GAMMA = 1.4


# --- Exact Riemann solver (Toro, Ch. 4) -----------------------------------
def _fK(p, rhoK, pK, cK, g):
    if p > pK:                                   # shock branch
        A = 2.0 / ((g + 1.0) * rhoK)
        B = (g - 1.0) / (g + 1.0) * pK
        return (p - pK) * np.sqrt(A / (p + B))
    return (2.0 * cK / (g - 1.0)) * ((p / pK) ** ((g - 1.0) / (2.0 * g)) - 1.0)


def _fKp(p, rhoK, pK, cK, g):
    if p > pK:
        A = 2.0 / ((g + 1.0) * rhoK)
        B = (g - 1.0) / (g + 1.0) * pK
        return np.sqrt(A / (B + p)) * (1.0 - (p - pK) / (2.0 * (B + p)))
    return (1.0 / (rhoK * cK)) * (p / pK) ** (-(g + 1.0) / (2.0 * g))


def exact_riemann(rhoL, uL, pL, rhoR, uR, pR, g):
    """Return a function S -> (rho, u, p) giving the exact solution at self-
    similar coordinate S = (x - x0)/t."""
    cL = np.sqrt(g * pL / rhoL)
    cR = np.sqrt(g * pR / rhoR)

    def f(p):
        return _fK(p, rhoL, pL, cL, g) + _fK(p, rhoR, pR, cR, g) + (uR - uL)

    def fp(p):
        return _fKp(p, rhoL, pL, cL, g) + _fKp(p, rhoR, pR, cR, g)

    p = max(1e-8, 0.5 * (pL + pR))               # Newton for star pressure
    for _ in range(100):
        pn = p - f(p) / fp(p)
        if pn <= 0.0:
            pn = 1e-8
        if abs(pn - p) < 1e-12 * max(1.0, p):
            p = pn
            break
        p = pn
    pstar = p
    ustar = 0.5 * (uL + uR) + 0.5 * (_fK(pstar, rhoR, pR, cR, g)
                                     - _fK(pstar, rhoL, pL, cL, g))

    def sample(S):
        if S <= ustar:                           # left of contact
            if pstar > pL:                       # left shock
                SL = uL - cL * np.sqrt((g + 1) / (2 * g) * pstar / pL
                                       + (g - 1) / (2 * g))
                if S <= SL:
                    return rhoL, uL, pL
                rho = rhoL * (pstar / pL + (g - 1) / (g + 1)) \
                    / ((g - 1) / (g + 1) * pstar / pL + 1.0)
                return rho, ustar, pstar
            cs = cL * (pstar / pL) ** ((g - 1) / (2 * g))   # left rarefaction
            if S <= uL - cL:
                return rhoL, uL, pL
            if S >= ustar - cs:
                return rhoL * (pstar / pL) ** (1.0 / g), ustar, pstar
            u = 2 / (g + 1) * (cL + (g - 1) / 2 * uL + S)
            c = 2 / (g + 1) * (cL + (g - 1) / 2 * (uL - S))
            rho = rhoL * (c / cL) ** (2 / (g - 1))
            return rho, u, pL * (c / cL) ** (2 * g / (g - 1))
        else:                                    # right of contact
            if pstar > pR:                       # right shock
                SR = uR + cR * np.sqrt((g + 1) / (2 * g) * pstar / pR
                                       + (g - 1) / (2 * g))
                if S >= SR:
                    return rhoR, uR, pR
                rho = rhoR * (pstar / pR + (g - 1) / (g + 1)) \
                    / ((g - 1) / (g + 1) * pstar / pR + 1.0)
                return rho, ustar, pstar
            cs = cR * (pstar / pR) ** ((g - 1) / (2 * g))   # right rarefaction
            if S >= uR + cR:
                return rhoR, uR, pR
            if S <= ustar + cs:
                return rhoR * (pstar / pR) ** (1.0 / g), ustar, pstar
            u = 2 / (g + 1) * (-cR + (g - 1) / 2 * uR + S)
            c = 2 / (g + 1) * (cR - (g - 1) / 2 * (uR - S))
            rho = rhoR * (c / cR) ** (2 / (g - 1))
            return rho, u, pR * (c / cR) ** (2 * g / (g - 1))

    return sample


if __name__ == "__main__":
    commandLine = CommandLineArguments(animate = False,
                                       cycles = 20000,
                                       nx = 200,
                                       ny = 10,
                                       dx = 0.005,
                                       dy = 0.005,
                                       tstop = 0.2,
                                       dtmin = 1e-7)

    myGrid = Grid2d(nx, ny, dx, dy)
    myNodeList = NodeList(nx * ny)
    constants = MKS()
    eos = IdealGasEOS(GAMMA, constants)
    hydro = GridHydroKT2d(myNodeList, constants, eos, myGrid)
    hydro.addBoundary(ReflectingGridBoundary2d(grid=myGrid))
    integrator = RungeKutta4Integrator2d([hydro], dtmin=dtmin, verbose=False)

    # Classic Sod state: (rho,p) = (1,1) left of x0, (0.125,0.1) right.
    x0 = 0.5 * nx * dx
    rhoL, pL, rhoR, pR = 1.0, 1.0, 0.125, 0.1
    density = myNodeList.getFieldDouble("density")
    energy = myNodeList.getFieldDouble("specificInternalEnergy")
    for j in range(ny):
        for i in range(nx):
            idx = myGrid.index(i, j, 0)
            if (i + 0.5) * dx < x0:
                density.setValue(idx, rhoL); energy.setValue(idx, pL / ((GAMMA - 1) * rhoL))
            else:
                density.setValue(idx, rhoR); energy.setValue(idx, pR / ((GAMMA - 1) * rhoR))

    controller = Controller(integrator=integrator, periodicWork=[], statStep=200, tstop=tstop)

    if animate:
        title = MakeTitle(controller, "time", "time")
        update_method = AnimationUpdateMethod2d(call=hydro.getCell2d, stepper=controller.Step,
                                                title=title, fieldName="density")
        AnimateGrid2d((nx, ny), update_method, extremis=[0, 1.1], frames=cycles, cmap="plasma")
    else:
        controller.Step(cycles)
        t = controller.time
        sample = exact_riemann(rhoL, 0.0, pL, rhoR, 0.0, pR, GAMMA)

        pressure = myNodeList.getFieldDouble("pressure")
        velocity = myNodeList.getFieldVector2d("velocity")
        jmid = ny // 2
        err_rho = err_u = err_p = 0.0
        n = 0
        for i in range(nx):
            idx = myGrid.index(i, jmid, 0)
            x = (i + 0.5) * dx
            rho_e, u_e, p_e = sample((x - x0) / t)
            err_rho += abs(density[idx] - rho_e)
            err_u += abs(velocity[idx].x - u_e)
            err_p += abs(pressure[idx] - p_e)
            n += 1
        err_rho /= n; err_u /= n; err_p /= n

        print(f"Sod at t={t:.4f}: L1 errors  rho={err_rho:.4e}  u={err_u:.4e}  p={err_p:.4e}")
        assert err_rho < 2.0e-2, f"density L1 error too large: {err_rho:.4e}"
        assert err_u < 2.0e-2, f"velocity L1 error too large: {err_u:.4e}"
        assert err_p < 2.0e-2, f"pressure L1 error too large: {err_p:.4e}"
        print("[ok] Sod shock tube matches the exact Riemann solution.")
