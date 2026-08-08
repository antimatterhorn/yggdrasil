# Copyright (C) 2026  Cody Raskin

import os
import sys

# Boundary's PYB11 descriptor lives in Boundaries/; PYB11Generator_add_module has no PYTHONPATH passthrough for this.
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "Boundaries"))

from PYB11Generator import *
PYB11includes = ['"patchNeighborBoundary.cc"',
                 '"coarseFineBoundary.cc"',
                 '"restrictionOperator.cc"']

from patchNeighborBoundary import *
from coarseFineBoundary import *
from restrictionOperator import *
