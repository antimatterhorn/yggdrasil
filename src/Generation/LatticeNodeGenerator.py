from DataBase import *
from LinearAlgebra import *

class Lattice2d:
    def __init__(self, nx, ny):
        self.numNodes = nx*ny
        self.positions = []
        self.bounds = [[-1, -1], [1, 1]]
        from Mesh import Grid2d

        mesh = Grid2d(nx=nx,ny=ny,sx=2/nx,sy=2/ny)
        mesh.setOrigin(Vector2d(1,1))
        
        for i in range(self.numNodes):
            self.positions.append((mesh.position(i).x,mesh.position(i).y))
        

class Lattice3d:
    def __init__(self, nx,ny,nz):
        self.numNodes = nx*ny*nz
        self.positions = []
        self.bounds = [[-1, -1], [1, 1]]
        from Mesh import Grid3d

        mesh = Grid3d(nx=nx,ny=ny,nz=nz,sx=2/nx,sy=2/ny,sz=2/nz)
        mesh.setOrigin(Vector2d(1,1,1))
        
        for i in range(self.numNodes):
            self.positions.append((mesh.position(i).x,mesh.position(i).y,mesh.position(i).z))