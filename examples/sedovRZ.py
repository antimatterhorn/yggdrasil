import os
from yggdrasil import *
from Animation import *
from Mesh import Grid2d, Geometry
from Physics import GridHydroHLLE2d
from EOS import IdealGasEOS
from Boundaries import ReflectingGridBoundary2d
from Utilities import SiloDump

# Sedov point blast in axisymmetric (r,z) coordinates. A pot of internal energy
# in the (r=0, z=0) corner drives a shock that expands as a hemisphere about the
# origin (the reflecting r=0 axis and z=0 base mirror the quarter-plane into a
# full sphere on revolution). Demonstrates the CylindricalRZ Grid, the p/r source
# term, and the auto-installed r=0 axis boundary, with live animation and Silo
# output for VisIt.

if __name__ == "__main__":
    commandLine = CommandLineArguments(animate = True,
                                       siloDump = False,
                                       dumpCycle = 25,
                                       cycles = 1500,
                                       nr = 120,
                                       nz = 120,
                                       dr = 0.01,
                                       dz = 0.01,
                                       dtmin = 1e-6,
                                       intVerbose = False,
                                       E_blob = 400.0,
                                       blob = 2,
                                       rootName = "sedovRZ",
                                       vizDir = "sedovRZ/viz")

    myGrid = Grid2d(nr, nz, dr, dz, Geometry.CylindricalRZ)
    print("grid size:", myGrid.size())

    myNodeList = NodeList(nr * nz)
    constants = MKS()
    eos = IdealGasEOS(1.4, constants)
    print(eos, "gamma =", eos.gamma)

    hydro = GridHydroHLLE2d(myNodeList, constants, eos, myGrid)

    # r=0 axis boundary is installed automatically for RZ grids; only the outer
    # walls are set here.
    box = ReflectingGridBoundary2d(grid=myGrid)
    box.setFaces(["right", "top", "bottom"])
    hydro.addBoundary(box)

    integrator = RungeKutta4Integrator2d([hydro], dtmin=dtmin, verbose=intVerbose)

    density = myNodeList.getFieldDouble("density")
    energy = myNodeList.getFieldDouble("specificInternalEnergy")
    for i in range(nr * nz):
        density.setValue(i, 1.0)
        energy.setValue(i, 1e-3)               # cold ambient

    # Deposit the blast energy in a small blob at the axis/base corner.
    for j in range(blob):
        for i in range(blob):
            energy.setValue(myGrid.index(i, j, 0), E_blob)

    periodicWork = []

    if siloDump:
        os.makedirs(vizDir, exist_ok=True)
        meshWriter = SiloDump(baseName=os.path.join(vizDir, rootName),
                              nodeList=myNodeList,
                              fieldNames=["density", "specificInternalEnergy",
                                          "pressure", "velocity"],
                              dumpCycle=dumpCycle,
                              grid=myGrid)
        periodicWork += [meshWriter]

    controller = Controller(integrator=integrator, periodicWork=periodicWork, statStep=25)

    if animate:
        title = MakeTitle(controller, "time", "time")
        bounds = (nr, nz)
        update_method = AnimationUpdateMethod2d(call=hydro.getCell2d,
                                                stepper=controller.Step,
                                                title=title,
                                                fieldName="density")
        AnimateGrid2d(bounds, update_method, extremis=[0, 4], frames=cycles, cmap="plasma")
    else:
        controller.Step(cycles)
