from yggdrasil import *
from Animation import *
from math import cos
from Mesh import Grid2d
from Physics import WaveEquation2d
from Boundaries import OutflowGridBoundary2d, DirichletGridBoundary2d, ReflectingGridBoundary2d
from Utilities import SiloDump


class oscillate:
    def __init__(self, nodeList, grid, nx, ny, workCycle=1):
        self.phi = nodeList.phi
        self.grid = grid
        self.cycle = workCycle
        # top source sits inside the horn, bottom source is in free space
        self.x_top = nx // 8
        self.y_top = 3 * ny // 4
        self.x_bot = nx // 8
        self.y_bot = ny // 4

    def __call__(self, cycle, time, dt):
        idx = self.grid.index(self.x_top, self.y_top, 0)
        self.phi.setValue(idx, 5.0 * cos(time))
        idx = self.grid.index(self.x_bot, self.y_bot, 0)
        self.phi.setValue(idx, 5.0 * cos(time))


if __name__ == "__main__":
    commandLine = CommandLineArguments(
        animate = True,
        cycles  = 20000,
        nx      = 200,
        ny      = 150,
    )

    constants = MKS()
    myNodeList = NodeList(nx * ny)
    grid = Grid2d(nx, ny, 1, 1)

    waveEqn = WaveEquation2d(nodeList=myNodeList, constants=constants,
                              grid=grid, C=1.0)

    outflow = OutflowGridBoundary2d(grid=grid)
    divWall = DirichletGridBoundary2d(grid=grid)   # hard barrier between halves
    cone    = ReflectingGridBoundary2d(grid=grid)  # Neumann (rigid) cone walls

    # ── dividing wall ─────────────────────────────────────────────────────────
    # Horizontal Dirichlet strip at y ≈ ny//2, full width.
    # Grid cell (i,j) has its centre at position (i+0.5, j+0.5) for dx=dy=1,
    # so addBox(p1, p2) covers j where p1.y ≤ j+0.5 ≤ p2.y.
    wall_y   = ny // 2      # = 75
    wall_bot = wall_y - 1   # p1.y: covers j ≥ 74
    wall_top = wall_y + 1   # p2.y: covers j ≤ 75  → two-cell-thick wall (j=74,75)
    divWall.addBox(Vector2d(0, wall_bot), Vector2d(nx, wall_top))

    # ── cone walls in top half (Neumann / reflecting) ─────────────────────────
    # An aperture wall at x=x_wall has a hole from y=y_aperture_lo to
    # y=y_aperture_hi.  From each hole edge the wall fans outward diagonally,
    # reaching the full top-half extent at x=nx.
    #
    #  top_bound ──────────────────────────────────────────
    #            |            |RRRR\                      |
    #            |   source   |RRRRR\  upper cone wall    |
    #            |    (•)     |RRRRRR\ (reflecting)       |
    #   aperture |    x=25    |hole   \                   |
    #   y=120 ── |............|.......>\  (cone interior) |
    #            |            |        >                  |
    #   y=104 ── |............|.......>/                  |
    #            |            |RRRRRR/                    |
    #            |            |RRRRR/  lower cone wall    |
    #            |            |RRRR/                      |
    #  first_free ─────────────────────────────────────────
    #            x=0         x=35                       x=199

    y_src_top = 3 * ny // 4      # = 112
    gap       = 8                # half-aperture in cells
    x_wall    = nx // 8 + 10     # = 35 — aperture wall x-position
    first_free = wall_top + 1    # = 77 — first free row above dividing wall
    top_bound  = ny - 1          # = 149

    y_aperture_lo = y_src_top - gap   # = 104
    y_aperture_hi = y_src_top + gap   # = 120

    horiz_span  = float(nx - x_wall)                 # = 165
    upper_coeff = float(ny - y_aperture_hi)           # = 30
    lower_coeff = float(y_aperture_lo - first_free)   # = 27

    # addBox(Vector2d(a, b), Vector2d(c, d)) covers cells i in [a, c) and j=b
    # when d = b+1.  x_right is the last inclusive column index, so we pass
    # x_right+1 as the exclusive upper bound, and y+1 to span exactly row y.

    # upper cone wall
    for y in range(y_aperture_hi + 1, top_bound + 1):
        x_right = int(x_wall + (y - y_aperture_hi) * horiz_span / upper_coeff - 1e-9)
        if x_right >= x_wall:
            cone.addBox(Vector2d(x_wall, y), Vector2d(x_right + 1, y + 1))

    # lower cone wall
    for y in range(first_free, y_aperture_lo):
        x_right = int(x_wall + (y_aperture_lo - y) * horiz_span / lower_coeff - 1e-9)
        if x_right >= x_wall:
            cone.addBox(Vector2d(x_wall, y), Vector2d(x_right + 1, y + 1))

    waveEqn.addBoundary(outflow)
    waveEqn.addBoundary(divWall)
    waveEqn.addBoundary(cone)

    integrator = RungeKutta4Integrator2d(packages=[waveEqn],
                                          dtmin=0.05, verbose=False)

    osc = oscillate(nodeList=myNodeList, grid=grid, nx=nx, ny=ny, workCycle=1)
    periodicWork = [osc]
    if not animate:
        dump = SiloDump("hornWave", myNodeList,
                        fieldNames=["phi", "xi"], dumpCycle=50)
        periodicWork.append(dump)

    controller = Controller(integrator=integrator,
                            statStep=100,
                            periodicWork=periodicWork)

    if animate:
        title = MakeTitle(controller, "time", "time")
        bounds = (nx, ny)
        update_method = AnimationUpdateMethod2d(
            call=waveEqn.getCell2d,
            stepper=controller.Step,
            title=title,
            fieldName="maxphi",
        )
        AnimateGrid2d(bounds, update_method,
                      extremis=[0, 0.1], frames=cycles, cmap='plasma')
    else:
        controller.Step(cycles)
