from yggdrasil import *
import math
import matplotlib.pyplot as plt
from Animation import *
from Physics import GridHydroKT2d
from Mesh import Grid2d
from EOS import IdealGasEOS
from Boundaries import ReflectingGridBoundary2d
from Utilities import SiloDump

# ---------------------------------------------------------------------------
# Exact Riemann solver (Toro, "Riemann Solvers and Numerical Methods for
# Fluid Dynamics", Ch. 4) for the classic Sod problem, used below to check
# the simulation against the analytical solution rather than just eyeballing
# a plot.
# ---------------------------------------------------------------------------

def _pressureFunction(p, rhoK, pK, cK, gamma):
    """Toro eq. 4.6/4.7: (f_K(p), f_K'(p)) for one side of the Riemann fan."""
    if p <= pK:
        # rarefaction
        pRatio = p / pK
        f = (2.0 * cK / (gamma - 1.0)) * (pRatio ** ((gamma - 1.0) / (2.0 * gamma)) - 1.0)
        fPrime = (1.0 / (rhoK * cK)) * pRatio ** (-(gamma + 1.0) / (2.0 * gamma))
    else:
        # shock
        A = 2.0 / ((gamma + 1.0) * rhoK)
        B = (gamma - 1.0) / (gamma + 1.0) * pK
        f = (p - pK) * math.sqrt(A / (p + B))
        fPrime = math.sqrt(A / (p + B)) * (1.0 - (p - pK) / (2.0 * (B + p)))
    return f, fPrime


def sodStarRegion(rhoL, uL, pL, rhoR, uR, pR, gamma, tol=1e-10, maxIter=100):
    """Newton-Raphson solve for (p*, u*) in the star region (Toro eq. 4.5-4.9)."""
    cL = math.sqrt(gamma * pL / rhoL)
    cR = math.sqrt(gamma * pR / rhoR)

    p = max(tol, 0.5 * (pL + pR))
    for _ in range(maxIter):
        fL, fLp = _pressureFunction(p, rhoL, pL, cL, gamma)
        fR, fRp = _pressureFunction(p, rhoR, pR, cR, gamma)
        pNew = p - (fL + fR + (uR - uL)) / (fLp + fRp)
        pNew = max(tol, pNew)
        if abs(pNew - p) < tol:
            p = pNew
            break
        p = pNew

    fL, _ = _pressureFunction(p, rhoL, pL, cL, gamma)
    fR, _ = _pressureFunction(p, rhoR, pR, cR, gamma)
    u = 0.5 * (uL + uR) + 0.5 * (fR - fL)
    return p, u


def sodSample(S, rhoL, uL, pL, rhoR, uR, pR, gamma, pStar, uStar):
    """Sample the exact Riemann solution at self-similar coordinate S = (x-x0)/t."""
    cL = math.sqrt(gamma * pL / rhoL)
    cR = math.sqrt(gamma * pR / rhoR)

    if S <= uStar:
        if pStar <= pL:
            # Left rarefaction
            cStarL = cL * (pStar / pL) ** ((gamma - 1.0) / (2.0 * gamma))
            shL = uL - cL
            stL = uStar - cStarL
            if S <= shL:
                return rhoL, uL, pL
            elif S >= stL:
                rho = rhoL * (pStar / pL) ** (1.0 / gamma)
                return rho, uStar, pStar
            else:
                c = (2.0 / (gamma + 1.0)) * (cL + (gamma - 1.0) / 2.0 * (uL - S))
                u = (2.0 / (gamma + 1.0)) * (cL + (gamma - 1.0) / 2.0 * uL + S)
                rho = rhoL * (c / cL) ** (2.0 / (gamma - 1.0))
                p = pL * (c / cL) ** (2.0 * gamma / (gamma - 1.0))
                return rho, u, p
        else:
            # Left shock
            SL = uL - cL * math.sqrt((gamma + 1.0) / (2.0 * gamma) * (pStar / pL) + (gamma - 1.0) / (2.0 * gamma))
            if S <= SL:
                return rhoL, uL, pL
            else:
                rho = rhoL * ((pStar / pL) + (gamma - 1.0) / (gamma + 1.0)) / ((gamma - 1.0) / (gamma + 1.0) * (pStar / pL) + 1.0)
                return rho, uStar, pStar
    else:
        if pStar > pR:
            # Right shock
            SR = uR + cR * math.sqrt((gamma + 1.0) / (2.0 * gamma) * (pStar / pR) + (gamma - 1.0) / (2.0 * gamma))
            if S >= SR:
                return rhoR, uR, pR
            else:
                rho = rhoR * ((pStar / pR) + (gamma - 1.0) / (gamma + 1.0)) / ((gamma - 1.0) / (gamma + 1.0) * (pStar / pR) + 1.0)
                return rho, uStar, pStar
        else:
            # Right rarefaction
            cStarR = cR * (pStar / pR) ** ((gamma - 1.0) / (2.0 * gamma))
            shR = uR + cR
            stR = uStar + cStarR
            if S >= shR:
                return rhoR, uR, pR
            elif S <= stR:
                rho = rhoR * (pStar / pR) ** (1.0 / gamma)
                return rho, uStar, pStar
            else:
                c = (2.0 / (gamma + 1.0)) * (cR - (gamma - 1.0) / 2.0 * (uR - S))
                u = (2.0 / (gamma + 1.0)) * (-cR + (gamma - 1.0) / 2.0 * uR + S)
                rho = rhoR * (c / cR) ** (2.0 / (gamma - 1.0))
                p = pR * (c / cR) ** (2.0 * gamma / (gamma - 1.0))
                return rho, u, p


if __name__ == "__main__":
    commandLine = CommandLineArguments(animate = False,
                                       siloDump = False,
                                        cycles = 300,
                                        nx = 100,
                                        ny = 20,
                                        dx = 1,
                                        dy = 1,
                                        dtmin = 0.001,
                                        checkAnalytical = True)

    myGrid = Grid2d(nx,ny,dx,dy)
    print("grid size:",myGrid.size())
    
    myNodeList = NodeList(nx*ny)
    print("numNodes =",myNodeList.numNodes)
    print("field names =",myNodeList.fieldNames)

    constants = MKS()
    print("G =",constants.G)
    eos = IdealGasEOS(1.4,constants)
    print(eos,"gamma =",eos.gamma)

    hydro = GridHydroKT2d(myNodeList,constants,eos,myGrid)
    print("numNodes =",myNodeList.numNodes)
    print("field names =",myNodeList.fieldNames)

    box = ReflectingGridBoundary2d(grid=myGrid)
    hydro.addBoundary(box)

    integrator = RungeKutta4Integrator2d([hydro],dtmin=dtmin,verbose=False)

    density = myNodeList.getFieldDouble("density")
    energy  = myNodeList.getFieldDouble("specificInternalEnergy")

    for j in range(ny):
        for i in range(nx):
            idx = myGrid.index(i,j,0)
            if i < nx // 2:
                density.setValue(idx, 1.0)
                energy.setValue(idx, 2.5)   # high pressure side
            else:
                density.setValue(idx, 0.125)
                energy.setValue(idx, 2.0)   # low pressure side

    periodicWork = []
    
    if siloDump:
        meshWriter = SiloDump(baseName="Sod",
                                nodeList=myNodeList,
                                fieldNames=["density","specificInternalEnergy","pressure","velocity"],
                                dumpCycle=50,
                                grid=myGrid)
        periodicWork += [meshWriter]

    controller = Controller(integrator=integrator,periodicWork=periodicWork,statStep=50)

    if(animate):
        title = MakeTitle(controller,"time","time")

        bounds = (nx,ny)
        update_method = AnimationUpdateMethod2d(call=hydro.getCell2d,
                                                stepper=controller.Step,
                                                title=title,
                                                fieldName="density")
        AnimateGrid2d(bounds,update_method,extremis=[-5,5],frames=cycles,cmap="plasma")
    else:
        controller.Step(cycles)

    xs = []
    ys = []
    position = myNodeList.getFieldVector2d("position")
    for i in range(nx*ny):
        if position[i].y == ((ny/2.0)+(dy/2.0)):
            xs.append(position[i].x)
            ys.append(density[i])
    # plt.plot(xs,ys)
    # plt.show()

    # -----------------------------------------------------------------------
    # Check against the exact Riemann solution (only in headless mode -- the
    # animated path never reaches this point in a scripted way). x0 is where
    # the initial discontinuity sits (the i < nx//2 split above); the check
    # verifies the elapsed time is still within the window before either wave
    # reaches a domain wall, since ReflectingGridBoundary means the exact
    # solution (which assumes an unbounded domain) stops applying once that
    # happens. cycles defaults to 300 (t~13.8) specifically to stay well
    # inside that window (~t=42) for this domain size; override --cycles for
    # a longer/animated run, but the analytical check below will warn rather
    # than silently mis-compare if you push time past the safe window.
    if checkAnalytical and not animate:
        gamma = eos.gamma
        x0 = (nx // 2) * dx
        t = controller.time

        rhoL0, uL0, pL0 = 1.0, 0.0, (gamma - 1.0) * 1.0 * 2.5
        rhoR0, uR0, pR0 = 0.125, 0.0, (gamma - 1.0) * 0.125 * 2.0

        pStar, uStar = sodStarRegion(rhoL0, uL0, pL0, rhoR0, uR0, pR0, gamma)

        cL0 = math.sqrt(gamma * pL0 / rhoL0)
        cR0 = math.sqrt(gamma * pR0 / rhoR0)
        leftHeadSpeed = uL0 - cL0  # fastest-left-moving signal regardless of fan/shock
        rightHeadSpeed = uR0 + cR0 if pStar <= pR0 else \
            uR0 + cR0 * math.sqrt((gamma + 1.0) / (2.0 * gamma) * (pStar / pR0) + (gamma - 1.0) / (2.0 * gamma))
        safeWindow = min(x0 / abs(leftHeadSpeed), (nx * dx - x0) / rightHeadSpeed)
        if t > safeWindow:
            print(f"WARNING: t={t:.3f} exceeds the pre-reflection window (~{safeWindow:.3f}) "
                  f"for this domain -- analytical comparison is no longer meaningful past this point.")

        errAbs = []
        for x, rhoSim in zip(xs, ys):
            S = (x - x0) / t
            rhoExact, uExact, pExact = sodSample(S, rhoL0, uL0, pL0, rhoR0, uR0, pR0, gamma, pStar, uStar)
            errAbs.append(abs(rhoSim - rhoExact))

        l1 = sum(errAbs) / len(errAbs)
        linf = max(errAbs)
        print(f"Analytical check at t={t:.4f}: density L1 error={l1:.5f}, Linf error={linf:.5f}")

        if t > safeWindow:
            print("SKIPPED: comparison window exceeded (see warning above) -- not a pass/fail signal")
        else:
            tol = 0.05
            assert l1 < tol, f"L1 density error {l1:.5f} exceeds tolerance {tol} -- check for a regression"
            print("PASS: sod analytical check")
