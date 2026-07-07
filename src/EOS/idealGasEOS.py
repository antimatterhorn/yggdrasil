# Copyright (C) 2026  Cody Raskin

from PYB11Generator import *
from equationOfState import *
from EOSAbstractMethods import *

class IdealGasEOS(EquationOfState):
    def pyinit(self,specificHeatRatio="double",constants="PhysicalConstants&"):
        return
    
    gamma = PYB11property("double", getter="getGamma", doc="Specific heat ratio.")

#-------------------------------------------------------------------------------
# Add the virtual interface
#-------------------------------------------------------------------------------
PYB11inject(EOSAbstractMethods, EquationOfState, pure_virtual=True)