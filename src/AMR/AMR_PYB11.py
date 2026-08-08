# Copyright (C) 2026  Cody Raskin

import os
import sys

# Boundary/FluxObserver PYB11 descriptors live in other modules' dirs; PYB11Generator_add_module has no PYTHONPATH passthrough for this.
_here = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(_here, "..", "Boundaries"))
sys.path.append(os.path.join(_here, "..", "Physics"))

from PYB11Generator import *
PYB11includes = ['"patchNeighborBoundary.cc"',
                 '"coarseFineBoundary.cc"',
                 '"restrictionOperator.cc"',
                 '"fluxRegister.cc"']

from patchNeighborBoundary import *
from coarseFineBoundary import *
from restrictionOperator import *
from fluxRegister import *
