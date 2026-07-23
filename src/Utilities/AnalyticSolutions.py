import math
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar

# Self-similar Sedov-Taylor point-blast solution for a uniform ambient density
# (omega=0), planar/cylindrical/spherical geometry, general gamma. This is the
# "standard" branch of the closed-form parametric solution (Kamm & Timmes,
# LA-UR-07-2849), ported from LANL's peer-reviewed reference implementation
# (ExactPack, exactpack/solvers/sedov/sedov.py) with omega=0 substituted
# through -- the singular/vacuum branches are dropped since they never trigger
# for omega=0, gamma>1 (verified: v2 < vstar always holds in that regime).
# Verified against the well-known reference value alpha=0.851072 for the
# gamma=1.4 spherical case (the energy constant conventionally used to set
# eblast so a unit-density blast's shock reaches R=1 at t=1).


class SedovSolution:
    """geometry: 1=planar, 2=cylindrical, 3=spherical. eblast must match that
    dimensionality (e.g. energy per unit length for geometry=2). rho0: uniform
    ambient density."""

    def __init__(self, gamma, geometry, eblast, rho0):
        assert geometry in (1, 2, 3)
        assert gamma > 1.0
        self.gamma = gamma
        self.geometry = geometry
        self.eblast = eblast
        self.rho0 = rho0

        gamm1 = gamma - 1.0
        gamp1 = gamma + 1.0
        gpogm = gamp1 / gamm1
        xg2 = geometry + 2.0
        denom2 = 2.0 * gamm1 + geometry
        denom3 = geometry * (2.0 - gamma)

        v2 = 4.0 / (xg2 * gamp1)
        vstar = 2.0 / (gamm1 * geometry + 2.0)
        assert v2 < vstar - 1e-4, "outside the 'standard' Sedov solution regime"
        assert abs(denom2) > 1e-4 and abs(denom3) > 1e-4, \
            "gamma/geometry combination hits a special-singularity branch (unsupported)"

        self.gamm1, self.gamp1, self.gpogm, self.xg2 = gamm1, gamp1, gpogm, xg2

        a0 = 2.0 / xg2
        a2 = -gamm1 / denom2
        a1 = xg2 * gamma / (2.0 + geometry * gamm1) * \
            ((2.0 * geometry * (2.0 - gamma)) / (gamma * xg2 * xg2) - a2)
        a3 = geometry / denom2
        a4 = xg2 * geometry * a1 / denom3
        a5 = -2.0 * geometry / denom3
        self.a0, self.a1, self.a2, self.a3, self.a4, self.a5 = a0, a1, a2, a3, a4, a5

        self.a_val = 0.25 * xg2 * gamp1
        self.b_val = gpogm
        self.c_val = 0.5 * xg2 * gamma
        self.d_val = (xg2 * gamp1) / (xg2 * gamp1 - 2.0 * (2.0 + geometry * gamm1))
        self.e_val = 0.5 * (2.0 + geometry * gamm1)

        self.v0 = 2.0 / (xg2 * gamma)          # similarity variable V at the origin
        self.v2 = v2                           # similarity variable V at the shock

        eval1 = quad(self._efun1, self.v0, self.v2, epsabs=1e-12)[0]
        eval2 = quad(self._efun2, self.v0, self.v2, epsabs=1e-12)[0]
        if geometry == 1:
            self.alpha = 0.5 * eval1 + eval2 / gamm1
        else:
            self.alpha = (geometry - 1.0) * math.pi * (eval1 + 2.0 * eval2 / gamm1)

    def _sedovFunctions(self, v):
        """lambda=r/R(t), and the density/velocity/pressure ratios (relative
        to their post-shock values) as functions of the similarity variable v."""
        a0, a1, a2, a3, a4, a5 = self.a0, self.a1, self.a2, self.a3, self.a4, self.a5
        x1 = self.a_val * v
        dx1dv = self.a_val
        cbag = max(1e-30, self.c_val * v - 1.0)
        x2 = self.b_val * cbag
        dx2dv = self.b_val * self.c_val
        ebag = 1.0 - self.e_val * v
        x3 = self.d_val * ebag
        dx3dv = -self.d_val * self.e_val
        x4 = max(self.b_val * (1.0 - 0.5 * self.xg2 * v), 1e-12)
        dx4dv = -self.b_val * 0.5 * self.xg2

        lam = x1 ** (-a0) * x2 ** (-a2) * x3 ** (-a1)
        dlamdv = -(a0 * dx1dv / x1 + a2 * dx2dv / x2 + a1 * dx3dv / x3) * lam
        f = x1 * lam
        g = x2 ** a3 * x3 ** a4 * x4 ** a5
        h = x1 ** (a0 * self.geometry) * x3 ** (a4 - 2.0 * a1) * x4 ** (1.0 + a5)
        return lam, dlamdv, f, g, h

    def _efun1(self, v):
        lam, dlamdv, f, g, h = self._sedovFunctions(v)
        return dlamdv * lam ** (self.geometry + 1.0) * self.gpogm * g * v * v

    def _efun2(self, v):
        lam, dlamdv, f, g, h = self._sedovFunctions(v)
        z = 8.0 / (self.xg2 * self.xg2 * self.gamp1)
        return dlamdv * lam ** (self.geometry - 1.0) * h * z

    def shockRadius(self, t):
        return (self.eblast / (self.alpha * self.rho0)) ** (1.0 / self.xg2) * t ** (2.0 / self.xg2)

    def shockJump(self, t):
        """Exact post-shock (density, velocity, pressure) at time t -- the
        Rankine-Hugoniot values, independent of any radial sampling. Use these
        (not the max of a sampled/binned profile array, which can miss the
        narrow near-shock spike entirely on a coarse or unstructured sampling)
        to normalize profile errors."""
        r2 = self.shockRadius(t)
        us = (2.0 / self.xg2) * r2 / t
        u2 = 2.0 * us / self.gamp1
        rho2 = self.gpogm * self.rho0
        p2 = 2.0 * self.rho0 * us * us / self.gamp1
        return rho2, u2, p2

    def profile(self, t, radii):
        """Density, velocity, and pressure arrays at time t for an array of
        radii >= 0 (ambient values beyond the shock)."""
        radii = np.asarray(radii, dtype=float)
        r2 = self.shockRadius(t)
        rho2, u2, p2 = self.shockJump(t)

        density = np.full_like(radii, self.rho0)
        velocity = np.zeros_like(radii)
        pressure = np.zeros_like(radii)

        for i, r in enumerate(radii):
            if r >= r2 or r2 <= 0.0:
                continue
            lamWant = r / r2

            def resid(v):
                return (self._sedovFunctions(v)[0] - lamWant) ** 2

            v = minimize_scalar(resid, bounds=(self.v0, self.v2), method="bounded",
                                options={"xatol": 1e-13}).x
            _, _, f, g, h = self._sedovFunctions(v)
            density[i] = rho2 * g
            velocity[i] = u2 * f
            pressure[i] = p2 * h

        return density, velocity, pressure


# Exact 1D Riemann solution (Toro, "Riemann Solvers and Numerical Methods for
# Fluid Dynamics", Ch. 4) for a two-state shock-tube initial condition, general
# gamma. Star-region pressure/velocity found by Newton iteration on the
# standard f_L(p) + f_R(p) + (uR-uL) = 0 relation; each state's own branch
# (shock vs. rarefaction) is picked by comparing p* to that side's initial
# pressure.


class SodSolution:
    """Two-state shock-tube initial condition (rhoL,uL,pL) | (rhoR,uR,pR)."""

    def __init__(self, gamma, rhoL, uL, pL, rhoR, uR, pR):
        self.gamma = gamma
        self.rhoL, self.uL, self.pL = rhoL, uL, pL
        self.rhoR, self.uR, self.pR = rhoR, uR, pR
        self.cL = math.sqrt(gamma * pL / rhoL)
        self.cR = math.sqrt(gamma * pR / rhoR)

        p = max(1e-8, 0.5 * (pL + pR))            # Newton solve for star pressure
        for _ in range(100):
            fval = self._fK(p, rhoL, pL, self.cL) + self._fK(p, rhoR, pR, self.cR) + (uR - uL)
            fderiv = self._fKp(p, rhoL, pL, self.cL) + self._fKp(p, rhoR, pR, self.cR)
            pn = p - fval / fderiv
            if pn <= 0.0:
                pn = 1e-8
            if abs(pn - p) < 1e-12 * max(1.0, p):
                p = pn
                break
            p = pn
        self.pstar = p
        self.ustar = 0.5 * (uL + uR) + 0.5 * (self._fK(p, rhoR, pR, self.cR)
                                              - self._fK(p, rhoL, pL, self.cL))

    def _fK(self, p, rhoK, pK, cK):
        g = self.gamma
        if p > pK:                                   # shock branch
            A = 2.0 / ((g + 1.0) * rhoK)
            B = (g - 1.0) / (g + 1.0) * pK
            return (p - pK) * math.sqrt(A / (p + B))
        return (2.0 * cK / (g - 1.0)) * ((p / pK) ** ((g - 1.0) / (2.0 * g)) - 1.0)

    def _fKp(self, p, rhoK, pK, cK):
        g = self.gamma
        if p > pK:
            A = 2.0 / ((g + 1.0) * rhoK)
            B = (g - 1.0) / (g + 1.0) * pK
            return math.sqrt(A / (B + p)) * (1.0 - (p - pK) / (2.0 * (B + p)))
        return (1.0 / (rhoK * cK)) * (p / pK) ** (-(g + 1.0) / (2.0 * g))

    def sample(self, S):
        """(rho, u, p) at self-similar coordinate S = (x-x0)/t."""
        g = self.gamma
        rhoL, uL, pL, cL = self.rhoL, self.uL, self.pL, self.cL
        rhoR, uR, pR, cR = self.rhoR, self.uR, self.pR, self.cR
        pstar, ustar = self.pstar, self.ustar

        if S <= ustar:                                # left of contact
            if pstar > pL:                             # left shock
                SL = uL - cL * math.sqrt((g + 1) / (2 * g) * pstar / pL + (g - 1) / (2 * g))
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
        else:                                          # right of contact
            if pstar > pR:                             # right shock
                SR = uR + cR * math.sqrt((g + 1) / (2 * g) * pstar / pR + (g - 1) / (2 * g))
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

    def profile(self, t, x, x0=0.0):
        """Density, velocity, and pressure arrays at time t for an array of
        absolute positions x (x0: initial diaphragm location)."""
        S = (np.asarray(x, dtype=float) - x0) / t
        out = np.array([self.sample(s) for s in S])
        return out[:, 0], out[:, 1], out[:, 2]


# Exact Noh problem (converging cold gas stagnating into a shock), general
# gamma and geometry index (1=planar, 2=cylindrical, 3=spherical). Derived
# from first principles (pre-shock density compresses purely kinematically
# from geometric convergence -- rho ~ (1 + v0 t/r)^(geometry-1) -- then the
# standard normal-shock jump is applied to that *local* pre-shock density at
# the moment the shock passes each radius; since post-shock velocity is
# exactly zero, shocked material never moves again, so the post-shock state is
# uniform and time-independent once set) and cross-checked against LANL's
# ExactPack (exactpack/solvers/noh/noh1.py) -- matches exactly, including the
# well-known rho2/rho0 = ((gamma+1)/(gamma-1))^geometry compression (e.g. 64x
# for the classic spherical gamma=5/3 case).


class NohSolution:
    """geometry: 1=planar, 2=cylindrical, 3=spherical. v0: inbound speed
    magnitude (positive -- the gas moves at -v0, toward the origin). rho0:
    uniform ambient density."""

    def __init__(self, gamma, geometry, v0, rho0):
        assert geometry in (1, 2, 3)
        assert gamma > 1.0
        assert v0 > 0.0
        self.gamma = gamma
        self.geometry = geometry
        self.v0 = v0
        self.rho0 = rho0

        self.D = v0 * (gamma - 1.0) / 2.0                             # shock speed
        self.rho2 = rho0 * ((gamma + 1.0) / (gamma - 1.0)) ** geometry  # post-shock density
        self.p2 = (gamma - 1.0) * self.rho2 * v0 * v0 / 2.0             # post-shock pressure

    def shockRadius(self, t):
        return self.D * t

    def profile(self, t, radii):
        """Density, radial velocity (negative = inbound ahead of the shock,
        zero/stagnant behind it), and pressure arrays at time t for an array
        of radii >= 0."""
        radii = np.asarray(radii, dtype=float)
        shocked = radii < self.shockRadius(t)

        preshock = self.rho0 * (1.0 + self.v0 * t / np.maximum(radii, 1e-300)) ** (self.geometry - 1.0)
        density = np.where(shocked, self.rho2, preshock)
        velocity = np.where(shocked, 0.0, -self.v0)
        pressure = np.where(shocked, self.p2, 0.0)
        return density, velocity, pressure
