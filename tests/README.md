# Tests

A mix of lightweight unit-style checks for core data structures (mesh, tree, integrator, database) and full physics-package acceptance tests that set up a simulation, run it, and inspect or visualize the result. Most scripts are standalone and runnable directly with `python <script>.py`; some depend on data fixtures in this directory (`CantinaBand.wav` is a short audio clip fed through the wave-equation speaker/microphone tests; `puzzle.nc` is a NetCDF data file not currently read by any script here).

## Contents

| Script | Description |
|---|---|
| [boundary_test.py](boundary_test.py) | Pass/fail regression suite for `GridHydroKT2d` boundary conditions — reflecting, outflow, periodic, and combined multi-face setups, checked via ghost-cell velocity/density values after one RK4 step |
| [bowshock.py](bowshock.py) | 2D grid-hydro simulation of a dense, fast-moving blob plunging into an exponential ambient density field (bow shock formation), animated with `GridHydroKT2d` and periodic boundaries |
| [buoyancy.py](buoyancy.py) | 2D grid-hydro Rayleigh–Taylor-style buoyancy setup: a dense circular region embedded in a hydrostatic two-density atmosphere under `ConstantForce2d` gravity |
| [centroid_test.py](centroid_test.py) | Builds a `VoronoiMesh2d` from random seed points and prints each cell's centroid, catching any exceptions |
| [cncTest.py](cncTest.py) | Exercises `CNCPathPhysics2d` by programming a rectangular toolpath, then samples and plots the path kinematically in pure Python for comparison |
| [conduction.py](conduction.py) | `ThermalConduction2d` diffusion of a hot circular region on a reflecting grid, using `IdealGasEOS` and `ConstantOpacity`, animated as a temperature field |
| [constraint_test.py](constraint_test.py) | HCP-lattice ball drop under `ConstantForce2d`/`Kinetics2d` with a `MotionConstraint2d` pinning a subset of nodes, animated as a scatter plot |
| [dataBase_test.py](dataBase_test.py) | Minimal smoke test that constructs a `WaveEquation2d` on a `NodeList`/`Grid2d` and confirms the `phi` field is reachable both by attribute and by `getFieldDouble` |
| [double-mach-reflection.py](double-mach-reflection.py) | Woodward–Colella double Mach reflection benchmark: a Mach-10 oblique shock striking a reflecting wedge, with a custom time-varying inflow boundary class and Silo dumps |
| [eosView.py](eosView.py) | Evaluates a selected EOS (ideal gas, Helmholtz, Mie–Gruneisen, Tillotson, polytropic, isothermal) over a log-spaced density/energy grid and 3D-surface-plots the resulting pressure |
| [gcodeCNC.py](gcodeCNC.py) | Standalone G-code parser and kinematic path simulator/plotter (G0/G1/G90/G91/F support); does not use the Yggdrasil physics stack |
| [grid_test.py](grid_test.py) | Prints `Grid2d` indexing, boundary flags, edge-node lists, and cell positions, and exercises `DirichletGridBoundary2d` domain/neighbor queries |
| [integrator_test.py](integrator_test.py) | Compares Euler, RK2, RK4, and Crank–Nicolson integrators against an analytic ODE solution using `ImplicitPhysics2d`, plotting solutions and residuals |
| [kelvin-helmholtz.py](kelvin-helmholtz.py) | Kelvin–Helmholtz shear-layer instability: three density/velocity bands with a sinusoidal velocity perturbation, evolved with `GridHydroKT2d` and periodic boundaries, dumped to Silo under `--vizDir`. Also drops a new `RestartWriter2d` checkpoint (named `<rootName>-restart-cycle<N>.ygr`) into `--restartDir` every `--restartCycle` steps, and on startup picks up the highest-numbered checkpoint found there automatically -- or an exact one via `--restoreCycle=N` -- so re-running the script resumes instead of restarting from the initial conditions |
| [kinetics.py](kinetics.py) | Two-particle straight-line drift test under `ConstantForce2d` (zero gravity) and `Kinetics2d`, animated as a scatter plot |
| [magTest.py](magTest.py) | Single charged particle orbiting in a time-oscillating uniform magnetic field via `MagField3d`, animated as a projected 2D scatter plot |
| [mic_test.py](mic_test.py) | Feeds `CantinaBand.wav` through a `Speaker` into a single-cell `WaveEquation2d`, records the response with a `Microphone`, and writes it back out to `mic.wav` |
| [nbody_test.py](nbody_test.py) | Direct N-body gravity simulation of 3 randomly placed bodies via `NBodyGravity2d` and `Kinematics2d`, animated as a scatter plot |
| [nodeGeneratorCompare2d.py](nodeGeneratorCompare2d.py) | Compares point-distribution uniformity (area-weight standard deviation) across 2D disk node generators — Fibonacci, constant-dtheta, Poisson, CVT |
| [nodeGeneratorCompare3d.py](nodeGeneratorCompare3d.py) | Compares point-distribution uniformity across 3D spherical-surface node generators — RPR, parameterized spiral, Fibonacci, SEAGen |
| [nodeGeneratorTest.py](nodeGeneratorTest.py) | Generates points on a spherical surface with a selectable generator and 3D-scatter-plots them colored by local area deviation from uniform |
| [orbit_test.py](orbit_test.py) | Two-body Keplerian orbit test under `PointSourceGravity2d`, comparing the numerically integrated trajectory and specific energy against the analytic ellipse |
| [outflowTest.py](outflowTest.py) | `WaveEquation2d` on a grid with `OutflowGridBoundary2d`, driven by a localized pulse, animated to check non-reflecting boundary behavior |
| [poisson.py](poisson.py) | Standalone Poisson-disk (blue noise) sampling implementation and scatter-plot visualization; does not use the Yggdrasil stack |
| [restart_test.py](restart_test.py) | Pass/fail regression suite for `RestartWriter2d`/`RestartReader2d` — checks that checkpointing a `GridHydroKT2d` run partway through and resuming from a fresh construction reproduces a continuous run's cycle/time/fields exactly, and that restoring into a mismatched-size `NodeList` raises rather than silently corrupting state |
| [rayleigh-taylor.py](rayleigh-taylor.py) | Rayleigh–Taylor instability with a sinusoidally perturbed heavy-over-light interface in hydrostatic equilibrium under `ConstantForce2d`, evolved with `GridHydroKT2d` and dumped to Silo |
| [sedov.py](sedov.py) | Sedov–Taylor point-blast test: a Gaussian energy deposition in a uniform-density grid, evolved with `GridHydroKT2d` under reflecting boundaries and plotted along the mid-line |
| [silo_test.py](silo_test.py) | Minimal `SiloMeshWriter2d` smoke test writing a `WaveEquation2d` field to a `.silo` file (contains a pre-existing `nodeList`/`myNodeList` naming bug) |
| [sod.py](sod.py) | Classic Sod shock-tube test: a density/energy step evolved with `GridHydroKT2d` under reflecting boundaries, animated and plotted along the mid-line |
| [treeGravityTest.py](treeGravityTest.py) | Barnes–Hut `TreeGravity2d` simulation of 20 randomly placed/velocity-seeded bodies, animated as a scatter plot |
| [tree_test.py](tree_test.py) | Builds a `KDTree3d` over a `Grid3d`'s node positions and looks up nearest neighbors of a central point |
| [voronoi_test.py](voronoi_test.py) | Builds a `VoronoiMesh2d` from a selectable 2D node generator, computes cell areas, and generates the dual `FEMesh` |
| [vtk_test.py](vtk_test.py) | Minimal `VTKMeshWriter2d` smoke test writing a `WaveEquation2d` field to a `.vtk` file |
| [waveBox.py](waveBox.py) | `WaveEquation2d` inside a `DirichletGridBoundary2d`-defined box-with-slit maze, driven by a `HarmonicOscillator` source and animated |
| [yggdrasil.py](yggdrasil.py) | Shared bootstrap module that adds the build directories to `sys.path` and re-exports the core C++/Python bindings (`DataBase`, `State`, `Integrators`, `Controller`, `IO`, etc.) imported via `from yggdrasil import *` at the top of most other test scripts |

### Data files

- `CantinaBand.wav` — audio fixture used by `mic_test.py` as a `Speaker` source signal.
- `puzzle.nc` — NetCDF data file; not currently referenced by any script in this directory.
