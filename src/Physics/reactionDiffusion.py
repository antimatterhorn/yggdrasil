from PYB11Generator import *
from physics import *

@PYB11template()
@PYB11template_dict({"dim" : "2"})
class ReactionDiffusion(Physics):
    def pyinit(self,
               nodeList="NodeList*",
               constants="PhysicalConstants&",
               grid="Mesh::Grid<%(dim)s>*",
               A="double",D="double"):
        return
    def getCell(self,i="int",j="int",fieldName="std::string"):
        return "std::array<double, 3>"