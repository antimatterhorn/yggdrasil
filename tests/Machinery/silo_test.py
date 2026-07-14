from yggdrasil import *
from Physics import WaveEquation2d
from Mesh import Grid2d

constants = PhysicalConstants(1,1,1.0,1.0,1.0) 
nx = 100
ny = 100
dx = 0.1
dy = 0.1

myGrid = Grid2d(nx,ny,dx,dy)
myNodeList = NodeList(nx*ny)

waveEqn = WaveEquation2d(nodeList=myNodeList,
                            constants=constants,
                            grid=myGrid,C=1.0)

meshWriter = SiloMeshWriter2d(baseName="testMesh",nodeList=myNodeList,fieldNames=["phi","xi"],grid=myGrid)
meshWriter.write("00.silo")

print(myGrid.nx,myGrid.ny,myGrid.nz,myGrid.dx,myGrid.dy,myGrid.dz)
print(myGrid.position(0))