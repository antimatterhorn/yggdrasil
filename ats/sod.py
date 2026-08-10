from yggdrasil import *
import numpy as np
from Physics import GridHydroKT2d, GridHydroHLLC2d, GridHydroHLLE2d

# Sod shock tube verified against the exact Riemann solution, for all three grid
# solvers (KT, HLLC, HLLE). Replaces the old snapshot-only hllcsod/hllesod tests.
# Each solver's density/velocity/pressure lineout is compared to the exact
# self-similar solution sampled at the achieved time (L1 error).

GAMMA = 1.4
RHOL, PL, RHOR, PR = 1.0, 1.0, 0.125, 0.1


def _fK(p, rhoK, pK, cK, g):
    if p > pK:
        A = 2.0 / ((g + 1) * rhoK); B = (g - 1) / (g + 1) * pK
        return (p - pK) * np.sqrt(A / (p + B))
    return (2 * cK / (g - 1)) * ((p / pK) ** ((g - 1) / (2 * g)) - 1.0)


def _fKp(p, rhoK, pK, cK, g):
    if p > pK:
        A = 2.0 / ((g + 1) * rhoK); B = (g - 1) / (g + 1) * pK
        return np.sqrt(A / (B + p)) * (1 - (p - pK) / (2 * (B + p)))
    return (1.0 / (rhoK * cK)) * (p / pK) ** (-(g + 1) / (2 * g))


def exact_riemann(rhoL, uL, pL, rhoR, uR, pR, g):
    cL = np.sqrt(g * pL / rhoL); cR = np.sqrt(g * pR / rhoR)
    f = lambda p: _fK(p, rhoL, pL, cL, g) + _fK(p, rhoR, pR, cR, g) + (uR - uL)
    fp = lambda p: _fKp(p, rhoL, pL, cL, g) + _fKp(p, rhoR, pR, cR, g)
    p = max(1e-8, 0.5 * (pL + pR))
    for _ in range(100):
        pn = p - f(p) / fp(p)
        if pn <= 0: pn = 1e-8
        if abs(pn - p) < 1e-12 * max(1.0, p): p = pn; break
        p = pn
    pstar = p
    ustar = 0.5 * (uL + uR) + 0.5 * (_fK(pstar, rhoR, pR, cR, g) - _fK(pstar, rhoL, pL, cL, g))

    def sample(S):
        if S <= ustar:
            if pstar > pL:
                SL = uL - cL * np.sqrt((g + 1) / (2 * g) * pstar / pL + (g - 1) / (2 * g))
                if S <= SL: return rhoL, uL, pL
                return rhoL * (pstar / pL + (g - 1) / (g + 1)) / ((g - 1) / (g + 1) * pstar / pL + 1), ustar, pstar
            cs = cL * (pstar / pL) ** ((g - 1) / (2 * g))
            if S <= uL - cL: return rhoL, uL, pL
            if S >= ustar - cs: return rhoL * (pstar / pL) ** (1 / g), ustar, pstar
            u = 2 / (g + 1) * (cL + (g - 1) / 2 * uL + S)
            c = 2 / (g + 1) * (cL + (g - 1) / 2 * (uL - S))
            return rhoL * (c / cL) ** (2 / (g - 1)), u, pL * (c / cL) ** (2 * g / (g - 1))
        else:
            if pstar > pR:
                SR = uR + cR * np.sqrt((g + 1) / (2 * g) * pstar / pR + (g - 1) / (2 * g))
                if S >= SR: return rhoR, uR, pR
                return rhoR * (pstar / pR + (g - 1) / (g + 1)) / ((g - 1) / (g + 1) * pstar / pR + 1), ustar, pstar
            cs = cR * (pstar / pR) ** ((g - 1) / (2 * g))
            if S >= uR + cR: return rhoR, uR, pR
            if S <= ustar + cs: return rhoR * (pstar / pR) ** (1 / g), ustar, pstar
            u = 2 / (g + 1) * (-cR + (g - 1) / 2 * uR + S)
            c = 2 / (g + 1) * (cR - (g - 1) / 2 * (uR - S))
            return rhoR * (c / cR) ** (2 / (g - 1)), u, pR * (c / cR) ** (2 * g / (g - 1))
    return sample


def _solve(Solver, nx=200, ny=4, dx=0.005, tstop=0.2):
    grid = Grid2d(nx, ny, dx, dx)
    nodes = NodeList(nx * ny)
    constants = MKS()
    eos = IdealGasEOS(GAMMA, constants)
    hydro = Solver(nodes, constants, eos, grid)
    # Bound to a name, not passed inline: the package stores a non-owning pointer, so a
    # temporary would be collected on return and leave it dangling.
    boundary = ReflectingGridBoundary2d(grid=grid)
    hydro.addBoundary(boundary)
    integrator = RungeKutta4Integrator2d([hydro], dtmin=1e-7, verbose=False)

    x0 = 0.5 * nx * dx
    density = nodes.getFieldDouble("density")
    energy = nodes.getFieldDouble("specificInternalEnergy")
    for j in range(ny):
        for i in range(nx):
            idx = grid.index(i, j, 0)
            if (i + 0.5) * dx < x0:
                density[idx] = RHOL; energy[idx] = PL / ((GAMMA - 1) * RHOL)
            else:
                density[idx] = RHOR; energy[idx] = PR / ((GAMMA - 1) * RHOR)

    Controller(integrator=integrator, periodicWork=[], statStep=100000, tstop=tstop).Step(100000)
    t = integrator.time
    sample = exact_riemann(RHOL, 0.0, PL, RHOR, 0.0, PR, GAMMA)
    pressure = nodes.getFieldDouble("pressure")
    velocity = nodes.getFieldVector2d("velocity")
    jmid = ny // 2
    e_rho = e_u = e_p = 0.0
    for i in range(nx):
        idx = grid.index(i, jmid, 0)
        rho_e, u_e, p_e = sample(((i + 0.5) * dx - x0) / t)
        e_rho += abs(density[idx] - rho_e)
        e_u += abs(velocity[idx].x - u_e)
        e_p += abs(pressure[idx] - p_e)
    return e_rho / nx, e_u / nx, e_p / nx


def run():
    checks = []
    for label, Solver, tol in (("KT", GridHydroKT2d, 2.0e-2),
                               ("HLLC", GridHydroHLLC2d, 3.0e-2),
                               ("HLLE", GridHydroHLLE2d, 3.0e-2)):
        e_rho, e_u, e_p = _solve(Solver)
        ok = bool(e_rho < tol and e_u < tol and e_p < tol)
        checks.append((f"{label} vs exact Riemann",
                       ok, f"L1 rho={e_rho:.2e} u={e_u:.2e} p={e_p:.2e} (tol {tol:.0e})"))
    return {"mode": "analytic", "checks": checks}


if __name__ == "__main__":
    print(run())
