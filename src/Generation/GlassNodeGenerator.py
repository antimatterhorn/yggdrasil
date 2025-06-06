from DataBase import *
from LinearAlgebra import *
from Trees import *
from numpy import sqrt
import random

class GlassDisk2d:
    def __init__(self, numNodes, overSample = 1):
        self.numNodes = numNodes
        self.positions = []
        self.bounds = [[-1, -1], [1, 1]]
        self.trials = overSample*numNodes
        
        pfield = FieldofVector2d("trials")
        minx,miny = self.bounds[0]
        maxx,maxy = self.bounds[1]
        
        area = (maxx - minx) * (maxy - miny)
        mean_spacing = sqrt(area / numNodes)
        dmin = mean_spacing * 0.8
        radius = mean_spacing * 2.0

        x = random.uniform(self.bounds[0][0], self.bounds[1][0])
        y = random.uniform(self.bounds[0][1], self.bounds[1][1])
        pfield.addValue(Vector2d(x,y)) # initial seed point
        self.positions.append([x,y])
        for i in range(numNodes-1):
            tree = KDTree2d(pfield)
            nbrs = numNodes
            next = (0,0)
            decay = (1.0 - i / numNodes)
            local_radius = radius * decay**0.5
            dist = 0
            for j in range(self.trials):
                x = random.uniform(self.bounds[0][0], self.bounds[1][0])
                y = random.uniform(self.bounds[0][1], self.bounds[1][1])
                if sqrt(x*x+y*y) > 1: continue
                tn = tree.findNearestNeighbors(Vector2d(x, y), dmin)
                if tn: continue  # too close to an existing point, skip
                tn_big = tree.findNearestNeighbors(Vector2d(x, y), local_radius)
                if len(tn_big) > 0:                    
                    for nbr in tn_big:
                        d = (Vector2d(x,y)-pfield[nbr]).magnitude()
                        if d > dist:
                            next = (x,y)
                else:  
                    next = (x,y)
                    continue   # didn't find any neighbors
            pfield.addValue(Vector2d(next[0],next[1]))
            self.positions.append([next[0],next[1]])

class GlassNodeGenerator2d:
    def __init__(self, numNodes, overSample = 1):
        self.numNodes = numNodes
        self.positions = []
        self.bounds = [[-1, -1], [1, 1]]
        self.trials = overSample*numNodes
        
        pfield = FieldofVector2d("trials")
        minx,miny = self.bounds[0]
        maxx,maxy = self.bounds[1]
        
        area = (maxx - minx) * (maxy - miny)
        mean_spacing = sqrt(area / numNodes)
        dmin = mean_spacing * 0.8
        radius = mean_spacing * 2.0

        x = random.uniform(self.bounds[0][0], self.bounds[1][0])
        y = random.uniform(self.bounds[0][1], self.bounds[1][1])
        pfield.addValue(Vector2d(x,y)) # initial seed point
        self.positions.append([x,y])
        for i in range(numNodes-1):
            tree = KDTree2d(pfield)
            nbrs = numNodes
            next = (0,0)
            decay = (1.0 - i / numNodes)
            local_radius = radius * decay**0.5
            for j in range(self.trials):
                x = random.uniform(self.bounds[0][0], self.bounds[1][0])
                y = random.uniform(self.bounds[0][1], self.bounds[1][1])
                tn = tree.findNearestNeighbors(Vector2d(x, y), dmin)
                if tn: continue  # too close to an existing point, skip
                tn_big = tree.findNearestNeighbors(Vector2d(x, y), local_radius)
                if len(tn_big) < nbrs:
                    nbrs = len(tn)
                    next = (x,y)
            pfield.addValue(Vector2d(next[0],next[1]))
            self.positions.append([next[0],next[1]])

class GlassNodeGenerator3d:
    def __init__(self, numNodes):
        self.numNodes = numNodes
        self.positions = []
        self.bounds = [[-1, -1, -1], [1, 1, 1]]
        self.trials = 20 * numNodes

        pfield = FieldofVector3d("trials")
        minx, miny, minz = self.bounds[0]
        maxx, maxy, maxz = self.bounds[1]

        volume = (maxx - minx) * (maxy - miny) * (maxz - minz)
        mean_spacing = (volume / numNodes) ** (1.0 / 3.0)
        dmin = mean_spacing * 0.8
        radius = mean_spacing * 2.0

        # Initial seed point
        x = random.uniform(minx, maxx)
        y = random.uniform(miny, maxy)
        z = random.uniform(minz, maxz)
        pfield.addValue(Vector3d(x, y, z))
        self.positions.append([x, y, z])

        for i in range(numNodes - 1):
            tree = KDTree3d(pfield)
            nbrs = numNodes
            next = (0, 0, 0)
            decay = (1.0 - i / numNodes)
            local_radius = radius * decay**0.5

            for j in range(self.trials):
                x = random.uniform(minx, maxx)
                y = random.uniform(miny, maxy)
                z = random.uniform(minz, maxz)
                trial = Vector3d(x, y, z)

                # Reject points too close to existing ones
                if tree.findNearestNeighbors(trial, dmin):
                    continue

                # Prefer points with fewest neighbors in a larger radius
                tn_big = tree.findNearestNeighbors(trial, local_radius)
                if len(tn_big) < nbrs:
                    nbrs = len(tn_big)
                    next = (x, y, z)

            pfield.addValue(Vector3d(next[0], next[1], next[2]))
            self.positions.append([next[0], next[1], next[2]])