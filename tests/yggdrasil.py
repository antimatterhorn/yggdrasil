import os
import sys

# Resolve relative to this file's real location (not the caller's cwd) so that
# importing "yggdrasil" works identically whether a test script lives directly
# in tests/ or in one of its subdirectories (Machinery/, Physics/), including
# via the yggdrasil.py symlinks placed in those subdirectories.
_testsDir = os.path.dirname(os.path.realpath(__file__))
buildir = os.path.join(_testsDir, "..", "build", "src") + "/"

for dir in ["Math",
            "DataBase",
            "Mesh",
            "State",
            "Trees",
            "Physics",
            "EOS",
            "Materials",
            "Type",
            "Integrators",
            "Utilities",
            "Boundaries",
            "Calculators",
            "IO",
            "Generation"]:
    sys.path.append(buildir+dir)

from CodeVersion import *
#from Mesh import *
from DataBase import *
from LinearAlgebra import *
#from Physics import *
from State import *
#from EOS import *
#from Opac import *
from PhysicalConstants import *
from Units import *
from Integrators import *
from Controller import *
from CommandLineArgs import *
#from Boundaries import *
#from Calculators import *
#from Trees import *
from IO import *
from Utilities import findLatestRestart, restartFileName, restoreIfAvailable
#from TillotsonMaterials import *
#from MieGruneisenMaterials import *
from Patch import *