import random
import math

class PoissonDisk2d:
    def __init__(self, npoints, k=30):
        """
        npoints: approximate number of points to generate
        k: number of candidate attempts per active point
        """
        self.positions = []
        self.k = k
        self._generate_points(npoints)

    def _generate_points(self, npoints):
        radius = 1.0  # unit disk
        area = math.pi * radius**2
        target_spacing = math.sqrt(area / npoints)
        r = target_spacing
        cell_size = r / math.sqrt(2)
        grid_dim = int(2 * radius / cell_size) + 1
        grid = [[None for _ in range(grid_dim)] for _ in range(grid_dim)]

        def in_disk(x, y):
            return x**2 + y**2 <= radius**2

        def get_cell_coords(x, y):
            return int((x + radius) / cell_size), int((y + radius) / cell_size)

        def is_valid(x, y):
            gx, gy = get_cell_coords(x, y)
            for i in range(max(0, gx-2), min(grid_dim, gx+3)):
                for j in range(max(0, gy-2), min(grid_dim, gy+3)):
                    neighbor = grid[i][j]
                    if neighbor:
                        dx = neighbor[0] - x
                        dy = neighbor[1] - y
                        if dx*dx + dy*dy < r*r:
                            return False
            return True

        # Initial seed point
        while True:
            x = random.uniform(-radius, radius)
            y = random.uniform(-radius, radius)
            if in_disk(x, y):
                break

        samples = [(x, y)]
        active = [(x, y)]
        gx, gy = get_cell_coords(x, y)
        grid[gx][gy] = (x, y)

        while active and len(samples) < npoints:
            idx = random.randint(0, len(active) - 1)
            ox, oy = active[idx]
            found = False

            for _ in range(self.k):
                theta = random.uniform(0, 2 * math.pi)
                rr = random.uniform(r, 2 * r)
                nx = ox + rr * math.cos(theta)
                ny = oy + rr * math.sin(theta)

                if not in_disk(nx, ny):
                    continue

                if is_valid(nx, ny):
                    samples.append((nx, ny))
                    active.append((nx, ny))
                    gx, gy = get_cell_coords(nx, ny)
                    grid[gx][gy] = (nx, ny)
                    found = True
                    break

            if not found:
                active.pop(idx)

        self.positions = samples

class PoissonNodeGenerator2d:
    def __init__(self, npoints, k=30):
        """
        npoints: approximate number of points to generate
        k: number of candidate attempts per active point
        """
        self.positions = []
        self.k = k
        self._generate_points(npoints)

    def _generate_points(self, npoints):
        radius = 1.0  # unit disk
        area = math.pi * radius**2
        target_spacing = math.sqrt(area / npoints)
        r = target_spacing
        cell_size = r / math.sqrt(2)
        grid_dim = int(2 * radius / cell_size) + 1
        grid = [[None for _ in range(grid_dim)] for _ in range(grid_dim)]

        def in_disk(x, y):
            return (-1 < x < 1) and (-1 < y < 1)

        def get_cell_coords(x, y):
            return int((x + radius) / cell_size), int((y + radius) / cell_size)

        def is_valid(x, y):
            gx, gy = get_cell_coords(x, y)
            for i in range(max(0, gx-2), min(grid_dim, gx+3)):
                for j in range(max(0, gy-2), min(grid_dim, gy+3)):
                    neighbor = grid[i][j]
                    if neighbor:
                        dx = neighbor[0] - x
                        dy = neighbor[1] - y
                        if dx*dx + dy*dy < r*r:
                            return False
            return True

        # Initial seed point
        while True:
            x = random.uniform(-radius, radius)
            y = random.uniform(-radius, radius)
            if in_disk(x, y):
                break

        samples = [(x, y)]
        active = [(x, y)]
        gx, gy = get_cell_coords(x, y)
        grid[gx][gy] = (x, y)

        while active and len(samples) < npoints:
            idx = random.randint(0, len(active) - 1)
            ox, oy = active[idx]
            found = False

            for _ in range(self.k):
                theta = random.uniform(0, 2 * math.pi)
                rr = random.uniform(r, 2 * r)
                nx = ox + rr * math.cos(theta)
                ny = oy + rr * math.sin(theta)

                if not in_disk(nx, ny):
                    continue

                if is_valid(nx, ny):
                    samples.append((nx, ny))
                    active.append((nx, ny))
                    gx, gy = get_cell_coords(nx, ny)
                    grid[gx][gy] = (nx, ny)
                    found = True
                    break

            if not found:
                active.pop(idx)

        self.positions = samples
