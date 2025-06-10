from yggdrasil import *
from Mesh import FEMesh2d

commandLine = CommandLineArguments()

mesh = FEMesh2d()
mesh.buildFromObj("example.obj",axes="(x,z)")

mesh.writeVTK("example.vtk")