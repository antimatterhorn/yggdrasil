# Hopfield memory simulator
# stores a handful of binary images to memory and then attempts to retrieve one of them
# from a randomized field.
# this is an extremely crude but illustrative example of neural memory

import numpy as np
from random import randint
from yggdrasil import *
from Animation import AnimateGrid2d

# === Hopfield Memory Simulator ===
class HopfieldSimulator:
    def __init__(self, nx, memories):
        self.nx = nx
        self.ny = nx
        self.num_nodes = nx * nx
        self.state = np.zeros(self.num_nodes)
        self.randomize()
        self.t = 0
        self.conv = 0
        self.title = "Hopfield Simulation"

        # Build weight matrix
        self.weights = sum(self.memorize(mem) for mem in memories)
        self.weights = self.weights / np.linalg.norm(self.weights)

    def memorize(self, mat):
        wtm = np.outer(mat, mat)
        np.fill_diagonal(wtm, 0)
        return wtm

    def randomize(self):
        for i in range(len(self.state)):
            val = randint(-1, 1)
            self.state[i] = -1 if val == 0 else val

    def remember(self):
        j = randint(0, self.num_nodes - 1)
        s = sum(self.weights[j][i] * self.state[i] for i in range(self.num_nodes) if i != j)
        self.state[j] = 1 if s >= 0 else -1

    def module_stepper(self):
        prev = np.copy(self.state)
        self.remember()
        self.t += 1
        if np.array_equal(prev, self.state):
            self.conv += 1
        else:
            self.conv = 0
        if self.conv == 200:
            print(f"Converged at step {self.t - 200}")
            exit()

    def module_call(self, i, j, fieldName):
        index = j * self.nx + i
        return 0.5 * (self.state[index] + 1.0)  # map [-1,1] to [0,1]

    def module_title(self):
        return f"Step {self.t}"
    
    def __call__(self, i, j):
        return self.module_call(i, j, None)

# === Define memories ===
memories = [
    np.array([-1,-1,-1,-1,1,1,1,-1,-1,-1,
              -1,-1,-1,1,1,1,1,-1,-1,-1,
              -1,-1,1,1,1,1,1,-1,-1,-1,
              -1,1,1,1,1,1,1,-1,-1,-1,
              -1,1,1,1,1,1,1,-1,-1,-1,
              -1,-1,-1,-1,1,1,1,-1,-1,-1,
              -1,-1,-1,-1,1,1,1,-1,-1,-1,
              -1,-1,-1,-1,1,1,1,-1,-1,-1,
              -1,1,1,1,1,1,1,1,1,-1,
              -1,1,1,1,1,1,1,1,1,-1]),
    
    np.array([1,1,1,1,1,1,1,1,1,1,
              1,1,1,1,1,1,1,1,1,1,
              1,1,-1,-1,-1,-1,-1,-1,1,1,
              1,1,-1,-1,-1,-1,-1,-1,1,1,
              1,1,-1,-1,-1,-1,-1,-1,1,1,
              1,1,-1,-1,-1,-1,-1,-1,1,1,
              1,1,-1,-1,-1,-1,-1,-1,1,1,
              1,1,-1,-1,-1,-1,-1,-1,1,1,
              1,1,1,1,1,1,1,1,1,1,
              1,1,1,1,1,1,1,1,1,1]),

    np.array([1,1,-1,-1,-1,-1,-1,-1,1,1,
              1,1,1,-1,-1,-1,-1,1,1,1,
              -1,1,1,1,-1,-1,1,1,1,-1,
              -1,-1,1,1,1,1,1,1,-1,-1,
              -1,-1,-1,1,1,1,1,-1,-1,-1,
              -1,-1,-1,1,1,1,1,-1,-1,-1,
              -1,-1,1,1,1,1,1,1,-1,-1,
              -1,1,1,1,-1,-1,1,1,1,-1,
              1,1,1,-1,-1,-1,-1,1,1,1,
              1,1,-1,-1,-1,-1,-1,-1,1,1])
]

# === Run animation ===
sim = HopfieldSimulator(10, memories)
AnimateGrid2d((10, 10), sim, frames=300, interval=50, cmap='gray', extremis=(0, 1))
