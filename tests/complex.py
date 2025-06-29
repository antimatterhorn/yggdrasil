from yggdrasil import *
from Mesh import Grid2d

commandLine = CommandLineArguments(nx = 20,
                                   ny = 20,
                                   dx = 0.1,
                                   dy = 0.1)

numNodes = nx*ny

nodeList = NodeList(numNodes)
myGrid = Grid2d(nx,ny,dx,dy)

myGrid.assignPositions(nodeList)

nodeList.insertFieldComplex("complexPosition")

positions = nodeList.getFieldVector2d("position")
cpos      = nodeList.getFieldComplex("complexPosition")

for i in range(numNodes):
    cpos.setValue(i,complex(positions[i].x,positions[i].y))
    print(positions[i],cpos[i])


print("numNodes =",nodeList.numNodes)
print("field names =",nodeList.fieldNames)