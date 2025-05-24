from PYB11Generator import *
from opacityModel import *
from OpacAbstractMethods import *

class ConstantOpacity(OpacityModel):
    def pyinit(self,k0="double",constants="PhysicalConstants&"):
        return

#-------------------------------------------------------------------------------
# Add the virtual interface
#-------------------------------------------------------------------------------
PYB11inject(OpacAbstractMethods, OpacityModel, pure_virtual=True)
    
