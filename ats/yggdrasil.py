import os
import sys

# libgomp sizes its thread pool from the *logical* CPU count, which oversubscribes
# hyperthreads; these kernels are memory-bound, so that costs ~3x at 200^2 (see
# PERFORMANCE.md). Must run before the compiled modules load, since libgomp reads
# OMP_NUM_THREADS when it initializes. An explicit setting is left alone.
if "OMP_NUM_THREADS" not in os.environ:
    import glob
    _cores = set()
    for _p in glob.glob("/sys/devices/system/cpu/cpu[0-9]*/topology/thread_siblings_list"):
        try:
            with open(_p) as _f:
                _cores.add(_f.read().strip())
        except OSError:
            pass
    os.environ["OMP_NUM_THREADS"] = str(len(_cores) or os.cpu_count() or 1)

buildir = "../build/src/"

for dir in ["Math",
            "DataBase",
            "Mesh",
            "State",
            "Trees",
            "Physics",
            "EOS",
            "Type",
            "Integrators",
            "Utilities",
            "Boundaries",
            "Calculators",
            "IO",
            "Generation",
            "Materials"]:
    sys.path.append(buildir+dir)
    
from Mesh import *
from DataBase import *
from LinearAlgebra import *
from Physics import *
from State import *
from EOS import *
from Opac import *
from PhysicalConstants import *
from Units import *
from Integrators import *
from Controller import *
from CommandLineArgs import *
from Boundaries import *
from Calculators import *
from Trees import *
from IO import *
from Utilities import *
from TillotsonMaterials import *
from MieGruneisenMaterials import *
from Patch import *