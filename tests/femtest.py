from yggdrasil import *
from Mesh import FEMesh2d

commandLine = CommandLineArguments()

mesh = FEMesh2d()
mesh.buildFromObj("example.obj",axes="(x,z)")

print(mesh.getNodes())
#print(mesh.getElements())

#mesh.writeVTK("example.vtk")