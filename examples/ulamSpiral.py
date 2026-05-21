from yggdrasil import *
from Mesh import Grid2d
from progressBar import ProgressBar
from random import random

import numpy as np
import matplotlib.pyplot as plt

commandLine = CommandLineArguments(n=200,
                                   randomize=False)

size = 2 * n + 1
numNodes = size * size

print(f"Generating Ulam Spiral with {numNodes} nodes ({size}×{size})...")
print("Checking all primes up to", numNodes)

sieve = np.ones(numNodes + 1, dtype=bool)
sieve[0] = sieve[1] = False
for i in range(2, int(numNodes**0.5) + 1):
    if sieve[i]:
        sieve[i*i::i] = False

print("Finished generating primes. Now creating grid and assigning prime values...")

myGrid = Grid2d(size, size, 1, 1)
myNodes = NodeList(numNodes)
myNodes.insertFieldDouble("prime")
prime_field = myNodes.prime

for k in range(numNodes):
    prime_field[k] = 0.0

# Build a parallel integer map: integer_map[k] = the integer at grid cell k
integer_map = np.zeros(numNodes, dtype=int)

directions = [(1, 0), (0, 1), (-1, 0), (0, -1)]
x, y = n, n
num = 1

prime_field[myGrid.index(x, y, 0)] = 0.0
integer_map[myGrid.index(x, y, 0)] = 1

def prob(j):
    if j < 5:
        return 0.5
    else:
        return np.sum(sieve[:j]) / j

dir_idx = 0
while num < numNodes:
    ProgressBar(num / numNodes, "Generating Ulam Spiral")
    dx, dy = directions[dir_idx % 4]
    steps = (dir_idx // 2) + 1
    for _ in range(steps):
        x += dx
        y += dy
        num += 1
        if num > numNodes:
            break
        idx = myGrid.index(x, y, 0)
        if not randomize:
            prime_field[idx] = 1.0 if sieve[num] else 0.0
        else: 
            prime_field[idx] = 1.0 if random() < prob(idx) else 0.0
        integer_map[idx] = num
    if num >= numNodes:
        break
    dir_idx += 1

data = np.array([prime_field[k] for k in range(numNodes)])
data2D = data.reshape((size, size))
integer2D = integer_map.reshape((size, size))

fig, ax = plt.subplots(figsize=(10, 10))
ax.imshow(data2D, origin='lower', cmap='binary_r', interpolation='nearest')
ax.set_title(f"Ulam Spiral  —  {size}×{size}  ({sieve.sum()} primes)")
ax.axis('off')

def format_coord(x, y):
    col, row = int(round(x)), int(round(y))
    if 0 <= col < size and 0 <= row < size:
        val = integer2D[row, col]
        is_prime = "prime" if sieve[val] else "composite"
        return f"({col}, {row})  n={val}  [{is_prime}]"
    return ""

ax.format_coord = format_coord

plt.tight_layout()
plt.show()