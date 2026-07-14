from yggdrasil import *
from EOS import IdealGasEOS
from Opac import KramersOpacity
from Physics import StellarEvolution
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Builds a toy 1D stellar structure (roughly Sun-sized/massed) with
# StellarEvolution, evolves it for a number of cycles under toy pp-chain
# burning, checks a few physical invariants, and (optionally) plots the
# radial structure and the surface luminosity history.
# ---------------------------------------------------------------------------

_results = []

def check(label, condition):
    status = "PASS" if condition else "FAIL"
    print(f"    [{status}] {label}")
    _results.append(condition)
    return condition

if __name__ == "__main__":
    commandLine = CommandLineArguments(numZones = 60,
                                        totalMass = 2.0e30,          # ~1 solar mass, kg
                                        radius = 7.0e8,               # ~1 solar radius, m
                                        centralTemperature = 1.5e7,   # K
                                        gamma = 5.0 / 3.0,
                                        kappa0 = 1.5e21,               # Kramers coefficient
                                        kappaES = 0.12,                # electron-scattering floor, m^2/kg
                                        cycles = 200,
                                        dtmin = 1.0,
                                        statStep = 40,
                                        plot = True)

    constants = MKS()
    print("G =", constants.G)

    nodeList = NodeList(numZones)
    eos = IdealGasEOS(gamma, constants)
    opac = KramersOpacity(kappa0, kappaES, constants)
    print(eos, "gamma =", eos.gamma)

    star = StellarEvolution(nodeList, constants, eos, opac,
                             totalMass=totalMass,
                             radius=radius,
                             centralTemperature=centralTemperature,
                             gamma=gamma)

    radiusField      = nodeList.getFieldDouble("radius")
    massField        = nodeList.getFieldDouble("mass")
    densityField     = nodeList.getFieldDouble("density")
    pressureField    = nodeList.getFieldDouble("pressure")
    temperatureField = nodeList.getFieldDouble("temperature")
    luminosityField  = nodeList.getFieldDouble("luminosity")
    energyField      = nodeList.getFieldDouble("specificInternalEnergy")

    print(f"\nInitial hydrostatic model ({numZones} zones):")
    print(f"  central density     = {densityField[0]:.6e} kg/m^3")
    print(f"  central pressure    = {pressureField[0]:.6e} Pa")
    print(f"  central temperature = {temperatureField[0]:.6e} K")
    print(f"  surface radius      = {radiusField[numZones-1]:.6e} m  (target {radius:.3e})")
    print(f"  surface mass        = {massField[numZones-1]:.6e} kg  (target {totalMass:.3e})")

    massError = abs(massField[numZones-1] - totalMass) / totalMass
    check("hydrostatic shooting converged to target mass (< 1%)", massError < 1e-2)
    check("density decreases monotonically outward at t=0",
          all(densityField[i] >= densityField[i+1] for i in range(numZones-1)))
    check("pressure decreases monotonically outward at t=0",
          all(pressureField[i] >= pressureField[i+1] for i in range(numZones-1)))
    check("temperature decreases monotonically outward at t=0",
          all(temperatureField[i] >= temperatureField[i+1] for i in range(numZones-1)))

    # Track surface luminosity and central conditions as the star burns
    history = {"time": [], "luminosity": [], "centralTemperature": []}

    class recordHistory:
        cycle = 1
        def __call__(self, cycle, time, dt):
            history["time"].append(time)
            history["luminosity"].append(luminosityField[numZones-1])
            history["centralTemperature"].append(temperatureField[0])

    integrator = RungeKutta2Integrator1d([star], dtmin=dtmin)
    controller = Controller(integrator=integrator,
                             periodicWork=[recordHistory()],
                             statStep=statStep)
    controller.Step(cycles)

    print(f"\nAfter {controller.cycle} cycles (t = {integrator.time:.6e} s):")
    print(f"  central temperature = {temperatureField[0]:.6e} K")
    print(f"  surface luminosity  = {luminosityField[numZones-1]:.6e} (toy units)")

    allFinite = all(
        f[i] == f[i] and abs(f[i]) != float("inf")
        for f in (densityField, pressureField, temperatureField, energyField)
        for i in range(numZones)
    )
    check("all zones finite (no NaN/Inf) after evolution", allFinite)
    check("specific internal energy stays positive", all(energyField[i] > 0 for i in range(numZones)))
    check("luminosity is non-negative and non-decreasing outward",
          all(luminosityField[i] <= luminosityField[i+1] and luminosityField[i] >= 0.0
              for i in range(numZones-1)))

    n_pass, n_total = sum(_results), len(_results)
    print()
    print(f"ALL {n_total} CHECKS PASSED" if n_pass == n_total
          else f"{n_total - n_pass} / {n_total} CHECKS FAILED")

    if plot:
        r = [radiusField[i] for i in range(numZones)]

        fig, axes = plt.subplots(2, 2, figsize=(10, 8))
        axes[0,0].plot(r, [densityField[i] for i in range(numZones)])
        axes[0,0].set_xlabel("radius (m)"); axes[0,0].set_ylabel("density (kg/m^3)")

        axes[0,1].plot(r, [pressureField[i] for i in range(numZones)])
        axes[0,1].set_xlabel("radius (m)"); axes[0,1].set_ylabel("pressure (Pa)")

        axes[1,0].plot(r, [temperatureField[i] for i in range(numZones)])
        axes[1,0].set_xlabel("radius (m)"); axes[1,0].set_ylabel("temperature (K)")

        axes[1,1].plot(r, [luminosityField[i] for i in range(numZones)])
        axes[1,1].set_xlabel("radius (m)"); axes[1,1].set_ylabel("luminosity (toy units)")

        fig.suptitle(f"StellarEvolution radial profile after {controller.cycle} cycles")
        fig.tight_layout()

        plt.figure()
        plt.plot(history["time"], history["luminosity"])
        plt.xlabel("time (s)")
        plt.ylabel("surface luminosity (toy units)")
        plt.title("Surface luminosity vs. time")

        plt.show()

    sys.exit(0 if n_pass == n_total else 1)
