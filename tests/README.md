# Tests

Split into two subdirectories:

- **[Machinery/](Machinery/)** — lightweight checks of core infrastructure: mesh/tree/integrator/database classes, boundary conditions, node generators, I/O writers, restart. Most of these validate one specific class's behavior (often via a small physics scenario used only as a vehicle to exercise it).
- **[Physics/](Physics/)** — full physics-package acceptance tests: set up a simulation, run it, and inspect, plot, or animate the result.

Most scripts are standalone and runnable directly with `python <script>.py` from inside their subdirectory. Each subdirectory has its own `yggdrasil.py`, which is a symlink to [yggdrasil.py](yggdrasil.py) below — not a copy — so there is exactly one real bootstrap module; `from yggdrasil import *` resolves the same way regardless of which subdirectory a script runs from.

Two scripts don't use the Yggdrasil stack at all and stay at this top level: [gcodeCNC.py](gcodeCNC.py) and [poisson.py](poisson.py). `CantinaBand.wav` (audio fixture for `Physics/mic_test.py`) and `KH/` (default viz/restart output dirs for `Physics/kelvin-helmholtz.py`) moved alongside the scripts that use them.

## Machinery/

| Script | Description |
|---|---|
| [boundary_test.py](Machinery/boundary_test.py) | Pass/fail regression suite for grid boundary conditions — reflecting, outflow, periodic, and combined multi-face setups on `GridHydroKT2d`, checked via ghost-cell velocity/density values after one RK4 step; plus non-cubic 3D (`GridHydroKT3d`) coverage for Periodic/Outflow/Reflecting guarding against the axis-mislabeling bugs fixed 2026-07-02 |
| [centroid_test.py](Machinery/centroid_test.py) | Builds a `VoronoiMesh2d` from random seed points and prints each cell's centroid, catching any exceptions |
| [constraint_test.py](Machinery/constraint_test.py) | HCP-lattice ball drop under `ConstantForce2d`/`Kinetics2d` with a `MotionConstraint2d` pinning a subset of nodes, animated as a scatter plot |
| [dataBase_test.py](Machinery/dataBase_test.py) | Minimal smoke test that constructs a `WaveEquation2d` on a `NodeList`/`Grid2d` and confirms the `phi` field is reachable both by attribute and by `getFieldDouble` |
| [eosView.py](Machinery/eosView.py) | Evaluates a selected EOS (ideal gas, Helmholtz, Mie–Gruneisen, Tillotson, polytropic, isothermal) over a log-spaced density/energy grid and 3D-surface-plots the resulting pressure |
| [grid_test.py](Machinery/grid_test.py) | Prints `Grid2d` indexing, boundary flags, edge-node lists, and cell positions, and exercises `DirichletGridBoundary2d` domain/neighbor queries |
| [integrator_test.py](Machinery/integrator_test.py) | Compares Euler, RK2, RK4, and Crank–Nicolson integrators against an analytic ODE solution using `ImplicitPhysics2d`, plotting solutions and residuals |
| [nodeGeneratorCompare2d.py](Machinery/nodeGeneratorCompare2d.py) | Compares point-distribution uniformity (area-weight standard deviation) across 2D disk node generators — Fibonacci, constant-dtheta, Poisson, CVT |
| [nodeGeneratorCompare3d.py](Machinery/nodeGeneratorCompare3d.py) | Compares point-distribution uniformity across 3D spherical-surface node generators — RPR, parameterized spiral, Fibonacci, SEAGen |
| [nodeGeneratorTest.py](Machinery/nodeGeneratorTest.py) | Generates points on a spherical surface with a selectable generator and 3D-scatter-plots them colored by local area deviation from uniform |
| [outflowTest.py](Machinery/outflowTest.py) | `WaveEquation2d` on a grid with `OutflowGridBoundary2d`, driven by a localized pulse, animated to check non-reflecting boundary behavior |
| [restart_test.py](Machinery/restart_test.py) | Pass/fail regression suite for `RestartWriter2d`/`RestartReader2d` — checks that checkpointing a `GridHydroKT2d` run partway through and resuming from a fresh construction reproduces a continuous run's cycle/time/fields exactly, and that restoring into a mismatched-size `NodeList` raises rather than silently corrupting state |
| [silo_test.py](Machinery/silo_test.py) | Minimal `SiloMeshWriter2d` smoke test writing a `WaveEquation2d` field to a `.silo` file (contains a pre-existing `nodeList`/`myNodeList` naming bug) |
| [tree_test.py](Machinery/tree_test.py) | Builds a `KDTree3d` over a `Grid3d`'s node positions and looks up nearest neighbors of a central point |
| [voronoi_test.py](Machinery/voronoi_test.py) | Builds a `VoronoiMesh2d` from a selectable 2D node generator, computes cell areas, and generates the dual `FEMesh` |
| [vtk_test.py](Machinery/vtk_test.py) | Minimal `VTKMeshWriter2d` smoke test writing a `WaveEquation2d` field to a `.vtk` file |

## Physics/

| Script | Description |
|---|---|
| [bowshock.py](Physics/bowshock.py) | 2D grid-hydro simulation of a dense, fast-moving blob plunging into an exponential ambient density field (bow shock formation), animated with `GridHydroKT2d` and periodic boundaries |
| [buoyancy.py](Physics/buoyancy.py) | 2D grid-hydro Rayleigh–Taylor-style buoyancy setup: a dense circular region embedded in a hydrostatic two-density atmosphere under `ConstantForce2d` gravity |
| [cncTest.py](Physics/cncTest.py) | Exercises `CNCPathPhysics2d` by programming a rectangular toolpath, then samples and plots the path kinematically in pure Python for comparison |
| [conduction.py](Physics/conduction.py) | `ThermalConduction2d` diffusion of a hot circular region on a reflecting grid, using `IdealGasEOS` and `ConstantOpacity`, animated as a temperature field |
| [double-mach-reflection.py](Physics/double-mach-reflection.py) | Woodward–Colella double Mach reflection benchmark: a Mach-10 oblique shock striking a reflecting wedge, with a custom time-varying inflow boundary class and Silo dumps |
| [kelvin-helmholtz.py](Physics/kelvin-helmholtz.py) | Kelvin–Helmholtz shear-layer instability: three density/velocity bands with a sinusoidal velocity perturbation, evolved with `GridHydroKT2d` and periodic boundaries, dumped to Silo under `--vizDir` (default `KH/viz`). Also drops a new `RestartWriter2d` checkpoint (named `<rootName>-restart-cycle<N>.ygr`) into `--restartDir` (default `KH/restart`) every `--restartCycle` steps, and on startup picks up the highest-numbered checkpoint found there automatically -- or an exact one via `--restoreCycle=N` -- so re-running the script resumes instead of restarting from the initial conditions |
| [kinetics.py](Physics/kinetics.py) | Two-particle straight-line drift test under `ConstantForce2d` (zero gravity) and `Kinetics2d`, animated as a scatter plot |
| [magTest.py](Physics/magTest.py) | Single charged particle orbiting in a time-oscillating uniform magnetic field via `MagField3d`, animated as a projected 2D scatter plot |
| [mic_test.py](Physics/mic_test.py) | Feeds `CantinaBand.wav` through a `Speaker` into a single-cell `WaveEquation2d`, records the response with a `Microphone`, and writes it back out to `mic.wav` |
| [nbody_test.py](Physics/nbody_test.py) | Direct N-body gravity simulation of 3 randomly placed bodies via `NBodyGravity2d` and `Kinematics2d`, animated as a scatter plot |
| [noh.py](Physics/noh.py) | Planar Noh problem: two cold streams collide head-on, driving a pair of stagnation shocks; run to a fixed `tstop` and checked/plotted against the exact solution (`AnalyticSolutions.NohSolution`) |
| [nohRZ.py](Physics/nohRZ.py) | Cylindrical `(r,z)` Noh problem: cold gas converges radially onto the axis; checked/plotted against the same exact solution, including the convergent pre-shock `1 + v0 t/r` profile unique to non-planar geometry |
| [orbit_test.py](Physics/orbit_test.py) | Two-body Keplerian orbit test under `PointSourceGravity2d`, comparing the numerically integrated trajectory and specific energy against the analytic ellipse |
| [rayleigh-taylor.py](Physics/rayleigh-taylor.py) | Rayleigh–Taylor instability with a sinusoidally perturbed heavy-over-light interface in hydrostatic equilibrium under `ConstantForce2d`, evolved with `GridHydroKT2d` and dumped to Silo |
| [sedov.py](Physics/sedov.py) | Sedov–Taylor point-blast test: a Gaussian energy deposition in a uniform-density 2D Cartesian grid (cylindrical symmetry, `GridHydroHLLE2d`, reflecting boundaries), run to a fixed `tstop` and checked/plotted against the closed-form self-similar solution (`AnalyticSolutions.SedovSolution`) |
| [sedovRZ.py](Physics/sedovRZ.py) | Cylindrical `(r,z)` Sedov point blast at the axis/equatorial-wall corner, representing a true spherical explosion; run to a fixed `tstop` and checked/plotted against the same self-similar solution, alongside finite/expansion/mass-conservation invariants |
| [sod.py](Physics/sod.py) | Classic Sod shock-tube test: a density/energy step evolved with `GridHydroKT2d` under reflecting boundaries, run to a fixed `tstop` and checked/plotted against the exact Riemann solution (`AnalyticSolutions.SodSolution`) |
| [stellarEvolution_test.py](Physics/stellarEvolution_test.py) | Builds a toy Sun-sized/massed 1D hydrostatic model with `StellarEvolution` and `IdealGasEOS`, evolves it under toy pp-chain burning, checks mass convergence/monotonicity/finiteness/luminosity invariants (PASS/FAIL), and plots the radial density/pressure/temperature/luminosity profile plus surface luminosity vs. time |
| [treeGravityTest.py](Physics/treeGravityTest.py) | Barnes–Hut `TreeGravity2d` simulation of 20 randomly placed/velocity-seeded bodies, animated as a scatter plot |
| [waveBox.py](Physics/waveBox.py) | `WaveEquation2d` inside a `DirichletGridBoundary2d`-defined box-with-slit maze, driven by a `HarmonicOscillator` source and animated |

## Top level

| Item | Description |
|---|---|
| [yggdrasil.py](yggdrasil.py) | Shared bootstrap module that adds the build directories to `sys.path` and re-exports the core C++/Python bindings (`DataBase`, `State`, `Integrators`, `Controller`, `IO`, etc.) imported via `from yggdrasil import *` at the top of most other test scripts. Resolves its build directory relative to its own (real, symlink-resolved) location, not the caller's cwd. |
| [gcodeCNC.py](gcodeCNC.py) | Standalone G-code parser and kinematic path simulator/plotter (G0/G1/G90/G91/F support); does not use the Yggdrasil physics stack |
| [poisson.py](poisson.py) | Standalone Poisson-disk (blue noise) sampling implementation and scatter-plot visualization; does not use the Yggdrasil stack |
