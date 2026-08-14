from yggdrasil import *
from Physics import WaveEquation2d
from Mesh import Grid2d

nx = ny = 10
constants = MKS()

grid = Grid2d(nx,ny,1,1)

myNodeList = NodeList(grid.size())

waveEqn = WaveEquation2d(nodeList=myNodeList,
                            constants=constants,
                            grid=grid,C=1)

phi = myNodeList.phi
phi2 = myNodeList.getFieldDouble("phi")

print(phi)
print(phi2)
print(phi.name,phi2.name)