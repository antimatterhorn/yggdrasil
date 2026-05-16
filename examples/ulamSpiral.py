from yggdrasil import *
from Mesh import Grid2d

import numpy as np
import matplotlib.pyplot as plt

commandLine = CommandLineArguments(n=200)

# Grid is (2n+1) x (2n+1) so the center cell is exactly at (n, n)
size = 2 * n + 1
numNodes = size * size

# Sieve of Eratosthenes
sieve = np.ones(numNodes + 1, dtype=bool)
sieve[0] = sieve[1] = False
for i in range(2, int(numNodes**0.5) + 1):
    if sieve[i]:
        sieve[i*i::i] = False

myGrid = Grid2d(size, size, 1, 1)
myNodes = NodeList(numNodes)
myNodes.insertFieldDouble("prime")
prime_field = myNodes.prime

for k in range(numNodes):
    prime_field[k] = 0.0

# Ulam spiral: direction sequence R, U, L, D with step counts 1, 1, 2, 2, 3, 3, ...
directions = [(1, 0), (0, 1), (-1, 0), (0, -1)]
x, y = n, n  # start at center
num = 1

# 1 is not prime, so center gets 0
prime_field[myGrid.index(x, y, 0)] = 0.0

dir_idx = 0
while num < numNodes:
    dx, dy = directions[dir_idx % 4]
    steps = (dir_idx // 2) + 1
    for _ in range(steps):
        x += dx
        y += dy
        num += 1
        if num > numNodes:
            break
        prime_field[myGrid.index(x, y, 0)] = 1.0 if sieve[num] else 0.0
    if num >= numNodes:
        break
    dir_idx += 1

# Reshape: grid stores index = i + j*nx, so reshape to (ny, nx) gives array[j, i]
data = np.array([prime_field[k] for k in range(numNodes)])
data2D = data.reshape((size, size))

plt.figure(figsize=(10, 10))
plt.imshow(data2D, origin='lower', cmap='binary_r', interpolation='nearest')
plt.title(f"Ulam Spiral  —  {size}×{size}  ({sieve.sum()} primes)")
plt.axis('off')
plt.tight_layout()
plt.show()
