from DataBase import *
from LinearAlgebra import *
from Trees import *
from numpy import sqrt
import random
from progressBar import ProgressBar

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

        rm = 2
        while rm > 1:
            x = random.uniform(self.bounds[0][0], self.bounds[1][0])
            y = random.uniform(self.bounds[0][1], self.bounds[1][1])
            r = Vector2d(x,y)
            rm = r.magnitude
        pfield.addValue(r) # initial seed point
        self.positions.append([x,y])
        for i in range(numNodes - 1):
            ProgressBar((i+1)/numNodes,"Generating Glass Disk Nodes")
            tree = KDTree2d(pfield)
            next = None
            decay = 1.0 - i / numNodes
            local_radius = radius * decay**0.5
            dist = 0

            for _ in range(self.trials):
                x = random.uniform(self.bounds[0][0], self.bounds[1][0])
                y = random.uniform(self.bounds[0][1], self.bounds[1][1])
                candidate = Vector2d(x, y)
                if candidate.magnitude > 1: continue
                if tree.findNearestNeighbors(candidate, dmin):
                    continue  # too close

                tn_big = tree.findNearestNeighbors(candidate, local_radius)
                if tn_big:
                    for nbr in tn_big:
                        d = (candidate - pfield[nbr]).magnitude
                        if d > dist:
                            next = [candidate.x,candidate.y]
                            dist = d
                else:
                    next = [candidate.x,candidate.y]
                    break  # accept immediately if no nearby points

            if next is not None:
                pfield.addValue(Vector2d(next[0],next[1]))
                self.positions.append(next)
            else:
                # fallback: generate random until it fits
                while True:
                    x = random.uniform(self.bounds[0][0], self.bounds[1][0])
                    y = random.uniform(self.bounds[0][1], self.bounds[1][1])
                    candidate = Vector2d(x, y)
                    if not tree.findNearestNeighbors(candidate, dmin):
                        pfield.addValue(candidate)
                        self.positions.append([x, y])
                        break

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
        for i in range(numNodes - 1):
            ProgressBar((i+1)/numNodes,"Generating Glass Nodes")
            tree = KDTree2d(pfield)
            next = None
            decay = 1.0 - i / numNodes
            local_radius = radius * decay**0.5
            dist = 0

            for _ in range(self.trials):
                x = random.uniform(self.bounds[0][0], self.bounds[1][0])
                y = random.uniform(self.bounds[0][1], self.bounds[1][1])
                candidate = Vector2d(x, y)

                if tree.findNearestNeighbors(candidate, dmin):
                    continue  # too close

                tn_big = tree.findNearestNeighbors(candidate, local_radius)
                if tn_big:
                    for nbr in tn_big:
                        d = (candidate - pfield[nbr]).magnitude
                        if d > dist:
                            next = [candidate.x,candidate.y]
                            dist = d
                else:
                    next = [candidate.x,candidate.y]
                    break  # accept immediately if no nearby points

            if next is not None:
                pfield.addValue(Vector2d(next[0],next[1]))
                self.positions.append(next)
            else:
                # fallback: generate random until it fits
                while True:
                    x = random.uniform(self.bounds[0][0], self.bounds[1][0])
                    y = random.uniform(self.bounds[0][1], self.bounds[1][1])
                    candidate = Vector2d(x, y)
                    if not tree.findNearestNeighbors(candidate, dmin):
                        pfield.addValue(candidate)
                        self.positions.append([x, y])
                        break

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