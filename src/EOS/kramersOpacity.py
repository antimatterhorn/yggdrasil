# Copyright (C) 2026  Cody Raskin

from PYB11Generator import *
from opacityModel import *
from OpacAbstractMethods import *

class KramersOpacity(OpacityModel):
    def pyinit(self,kappa0="double",kappaES="double",constants="PhysicalConstants&"):
        return

#-------------------------------------------------------------------------------
# Add the virtual interface
#-------------------------------------------------------------------------------
PYB11inject(OpacAbstractMethods, OpacityModel, pure_virtual=True)
