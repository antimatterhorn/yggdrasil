import random
import math

def poisson_disk_sample(width, height, r, k=30):
    """Generate Poisson disk samples in a width × height domain with min distance r."""
    cell_size = r / math.sqrt(2)
    grid_width = int(width / cell_size) + 1
    grid_height = int(height / cell_size) + 1

    grid = [[None for _ in range(grid_height)] for _ in range(grid_width)]
    samples = []
    active = []

    # Initial point
    x = random.uniform(0, width)
    y = random.uniform(0, height)
    samples.append((x, y))
    active.append((x, y))
    gx, gy = int(x / cell_size), int(y / cell_size)
    grid[gx][gy] = (x, y)

    while active:
        idx = random.randint(0, len(active) - 1)
        ox, oy = active[idx]
        found = False

        for _ in range(k):
            theta = random.uniform(0, 2 * math.pi)
            rr = random.uniform(r, 2 * r)
            nx = ox + rr * math.cos(theta)
            ny = oy + rr * math.sin(theta)

            if not (0 <= nx < width and 0 <= ny < height):
                continue

            gx, gy = int(nx / cell_size), int(ny / cell_size)

            # Check neighboring cells
            ok = True
            for i in range(max(0, gx-2), min(grid_width, gx+3)):
                for j in range(max(0, gy-2), min(grid_height, gy+3)):
                    neighbor = grid[i][j]
                    if neighbor:
                        dx = neighbor[0] - nx
                        dy = neighbor[1] - ny
                        if dx*dx + dy*dy < r*r:
                            ok = False
                            break
                if not ok:
                    break

            if ok:
                samples.append((nx, ny))
                active.append((nx, ny))
                grid[gx][gy] = (nx, ny)
                found = True
                break

        if not found:
            active.pop(idx)

    return samples


samples = poisson_disk_sample(2, 2, 0.05)
print(f"Generated {len(samples)} samples")

# Optional visualization
import matplotlib.pyplot as plt
x, y = zip(*samples)
plt.figure(figsize=(6, 6))
plt.scatter(x, y, s=5)
plt.gca().set_aspect('equal')
plt.title("Poisson Disk Sample (Blue Noise)")
plt.show()