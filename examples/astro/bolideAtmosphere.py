import os
import math
from yggdrasil import *
from Animation import *
from Mesh import Grid2d, Geometry
from Physics import GridHydroKT2d, ConstantForce2d
from EOS import IdealGasEOS
from Boundaries import ReflectingGridBoundary2d, OutflowGridBoundary2d
from Utilities import SiloDump

# Axisymmetric (r,z) Earth atmosphere in hydrostatic equilibrium, perturbed by
# a bolide's energy trail read from an external (time, altitude, energy) text
# file and sourced into the r=0 column as the simulation clock crosses each
# entry's timestamp. KT hydro on a CylindricalRZ Grid, constant downward
# gravity, isothermal ideal-gas atmosphere out to ~10km.


class BolideEnergySource:
    """Periodic-work callable: at construction, reads a (time, altitude,
    energy) text file describing a bolide's deposited-energy trail; on each
    call, deposits any event whose timestamp the simulation clock has reached
    into the atmospheric column within columnRadius of the r=0 axis and
    layerThickness of the event's altitude."""

    def __init__(self, nodeList, grid, filename, columnRadius, layerThickness, sourceCycle=1):
        self.density = nodeList.getFieldDouble("density")
        self.energy = nodeList.getFieldDouble("specificInternalEnergy")
        self.position = nodeList.getFieldVector2d("position")
        self.grid = grid
        self.cycle = sourceCycle
        self.columnRadius = columnRadius
        self.layerThickness = layerThickness
        self.events = self._readEvents(filename)
        self.nextEvent = 0

    @staticmethod
    def _readEvents(filename):
        events = []
        with open(filename) as f:
            for line in f:
                fields = line.strip().split(",")
                if len(fields) < 3:
                    continue
                try:
                    t, altitude, e = float(fields[0]), float(fields[1]), float(fields[2])
                except ValueError:
                    continue   # header/comment row
                events.append((t, altitude, e))
        events.sort(key=lambda ev: ev[0])
        return events

    def _columnCellIndices(self, altitude):
        halfThickness = 0.5 * self.layerThickness
        cells = []
        for idx in range(self.grid.nx * self.grid.ny):
            pos = self.position[idx]
            if pos.x <= self.columnRadius and abs(pos.y - altitude) <= halfThickness:
                cells.append(idx)
        return cells

    def __call__(self, cycle, time, dt):
        while self.nextEvent < len(self.events) and self.events[self.nextEvent][0] <= time:
            eventTime, altitude, totalEnergy = self.events[self.nextEvent]
            self.nextEvent += 1
            cells = self._columnCellIndices(altitude)
            if not cells:
                continue
            # grid.cellVolume() is "per radian" (the 2*pi revolve factor
            # cancels inside the hydro's flux divergence) -- restore it here
            # to get the column's true mass, since totalEnergy is a real
            # physical Joule quantity for the fully revolved event.
            trueMass = 2.0 * math.pi * sum(self.density[idx] * self.grid.cellVolume(idx) for idx in cells)
            if trueMass <= 0.0:
                continue
            deltaSie = totalEnergy / trueMass
            for idx in cells:
                self.energy.setValue(idx, self.energy[idx] + deltaSie)
            print("bolide source: t=%.3f (event t=%.3f) deposited %.3e J over %d cells at altitude %.1fm "
                  "(column mass=%.3e kg, du=%.3e J/kg)"
                  % (time, eventTime, totalEnergy, len(cells), altitude, trueMass, deltaSie))


if __name__ == "__main__":
    commandLine = CommandLineArguments(animate = False,
                                       siloDump = True,
                                       dumpCycle = 25,
                                       cycles = 20000,
                                       tstop = 5.0,
                                       nr = 100,
                                       nz = 100,
                                       dr = 100.0,
                                       dz = 100.0,
                                       dtmin = 1e-6,
                                       intVerbose = False,
                                       statStep = 25,
                                       rho0 = 1.225,        # sea-level air density, kg/m^3
                                       p0 = 101325.0,        # sea-level pressure, Pa
                                       gamma = 1.4,
                                       g = 9.81,             # m/s^2
                                       bolideFile = "bolide.txt",
                                       columnRadius = 300.0,  # m, radial extent of the r=0 source column
                                       layerThickness = 300.0,  # m, altitude window each event deposits into
                                       rootName = "bolideAtmosphere",
                                       vizDir = "bolideAtmosphere/viz")

    myGrid = Grid2d(nr, nz, dr, dz, Geometry.CylindricalRZ)
    print("grid size:", myGrid.size(), " domain: %gm x %gm (r x z)" % (nr * dr, nz * dz))

    myNodeList = NodeList(nr * nz)
    constants = MKS()
    eos = IdealGasEOS(gamma, constants)
    print(eos, "gamma =", eos.gamma)

    hydro = GridHydroKT2d(myNodeList, constants, eos, myGrid)

    # r=0 axis boundary is installed automatically for RZ grids. Ground (z=0)
    # is a solid reflecting floor; the outer radial and top edges are open so
    # the blast can vent rather than reflect back into the domain.
    ground = ReflectingGridBoundary2d(grid=myGrid)
    ground.setFaces(["bottom"])
    hydro.addBoundary(ground)

    sky = OutflowGridBoundary2d(grid=myGrid)
    sky.setFaces(["right", "top"])
    hydro.addBoundary(sky)

    gravityVector = Vector2d(0.0, -g)
    gravity = ConstantForce2d(myNodeList, constants, gravityVector)

    integrator = RungeKutta4Integrator2d([hydro, gravity], dtmin=dtmin, verbose=intVerbose)

    # Isothermal hydrostatic atmosphere: rho(z) = rho0*exp(-z/H), p(z) =
    # p0*exp(-z/H), with scale height H = p0/(rho0*g). Constant temperature
    # keeps specific internal energy p/((gamma-1)*rho) uniform with altitude,
    # so dp/dz = -rho*g is satisfied exactly, not by discrete integration.
    # GridHydroBase's gravity coupling isn't well-balanced, so this discrete
    # IC holds only approximately -- expect a modest ambient drift (a few
    # tens of m/s) over the run, well below the bolide's own signal.
    H = p0 / (rho0 * g)
    sie0 = p0 / ((gamma - 1.0) * rho0)
    print("scale height H = %.1fm, sea-level specific internal energy = %.3e J/kg" % (H, sie0))

    density = myNodeList.getFieldDouble("density")
    energy = myNodeList.getFieldDouble("specificInternalEnergy")
    position = myNodeList.getFieldVector2d("position")
    for i in range(nr):
        for j in range(nz):
            idx = myGrid.index(i, j, 0)
            z = position[idx].y
            density.setValue(idx, rho0 * math.exp(-z / H))
            energy.setValue(idx, sie0)

    periodicWork = []

    bolideSource = BolideEnergySource(myNodeList, myGrid, bolideFile,
                                      columnRadius=columnRadius,
                                      layerThickness=layerThickness,
                                      sourceCycle=1)
    periodicWork.append(bolideSource)

    if siloDump:
        os.makedirs(vizDir, exist_ok=True)
        meshWriter = SiloDump(baseName=os.path.join(vizDir, rootName),
                              nodeList=myNodeList,
                              fieldNames=["density", "specificInternalEnergy",
                                          "pressure", "velocity"],
                              dumpCycle=dumpCycle,
                              grid=myGrid)
        periodicWork.append(meshWriter)

    controller = Controller(integrator=integrator, periodicWork=periodicWork,
                            statStep=statStep, tstop=tstop)

    if animate:
        title = MakeTitle(controller, "time", "time")
        bounds = (nr, nz)
        update_method = AnimationUpdateMethod2d(call=hydro.getCell2d,
                                                stepper=controller.Step,
                                                title=title,
                                                fieldName="density")
        AnimateGrid2d(bounds, update_method, extremis=[0, 1.3 * rho0], frames=cycles, cmap="plasma")
    else:
        controller.Step(cycles)
