from DataBase import *
from LinearAlgebra import *
from Trees import *
from numpy import sqrt

class GlassNodeGenerator2d:
    def __init__(self, numNodes):
        self.numNodes = numNodes
        self.positions = []
        self.bounds = [[-1, -1], [1, 1]]
        self.trials = 20*numNodes
        
        pfield = FieldofVector2d("position")
        minx = self.bounds[0][0]
        maxx = self.bounds[1][0]
        miny = self.bounds[0][1]
        maxy = self.bounds[1][1]
        area = (maxx - minx) * (maxy - miny)
        mean_spacing = sqrt(area / numNodes)
        dmin = mean_spacing * 0.8
        radius = mean_spacing * 2.0
        
        import random
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

