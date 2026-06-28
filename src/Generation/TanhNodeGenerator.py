import math

def _tanh_1d(n, beta, mode):
    # Returns n values in [0, 1] stretched by the chosen tanh map.
    if n == 1:
        return [0.5]
    tb = math.tanh(beta)
    tb2 = math.tanh(beta / 2.0)
    pts = []
    for i in range(n):
        t = i / (n - 1)
        if mode == 'uniform':
            xi = t
        elif mode == 'right':
            xi = math.tanh(beta * t) / tb
        elif mode == 'left':
            xi = 1.0 - math.tanh(beta * (1.0 - t)) / tb
        elif mode == 'both':
            xi = 0.5 + math.tanh(beta * (t - 0.5)) / (2.0 * tb2)
        elif mode == 'center':
            if t <= 0.5:
                xi = math.tanh(2.0 * beta * t) / (2.0 * tb)
            else:
                xi = 1.0 - math.tanh(2.0 * beta * (1.0 - t)) / (2.0 * tb)
        else:
            raise ValueError(f"mode must be 'left', 'right', 'both', 'center', or 'uniform', got '{mode}'")
        pts.append(xi)
    return pts


class TanhNodeGenerator1d:
    """Tanh-stretched 1D node distribution.

    mode:
      'left'    — cluster near x0 (boundary layer at left wall)
      'right'   — cluster near x1 (boundary layer at right wall)
      'both'    — cluster near both ends (channel flow, two-wall BL)
      'center'  — cluster near center (shock tube interface)
      'uniform' — no stretching (same as Lattice1d)

    beta controls stretching intensity; larger beta = stronger concentration.
    beta=0 is degenerate (use 'uniform' instead); beta~3 is moderate, beta~8 is aggressive.
    """
    def __init__(self, n, x0=-1.0, x1=1.0, beta=3.0, mode='both'):
        xs = _tanh_1d(n, beta, mode)
        L = x1 - x0
        self.numNodes = n
        self.positions = [x0 + L * xi for xi in xs]
        self.bounds = [x0, x1]


class TanhNodeGenerator2d:
    """Tanh-stretched 2D node distribution (tensor product of two 1D stretches).

    Each axis has its own beta and mode, allowing e.g. boundary-layer clustering
    in y with uniform spacing in x.

    See TanhNodeGenerator1d for mode descriptions.
    """
    def __init__(self, nx, ny,
                 x0=-1.0, x1=1.0, beta_x=3.0, mode_x='both',
                 y0=-1.0, y1=1.0, beta_y=3.0, mode_y='both'):
        xs = _tanh_1d(nx, beta_x, mode_x)
        ys = _tanh_1d(ny, beta_y, mode_y)
        Lx, Ly = x1 - x0, y1 - y0
        self.numNodes = nx * ny
        self.positions = [
            [x0 + Lx * xi, y0 + Ly * yj]
            for yj in ys
            for xi in xs
        ]
        self.bounds = [[x0, y0], [x1, y1]]


class TanhNodeGenerator3d:
    """Tanh-stretched 3D node distribution (tensor product of three 1D stretches).

    See TanhNodeGenerator1d for mode descriptions.
    """
    def __init__(self, nx, ny, nz,
                 x0=-1.0, x1=1.0, beta_x=3.0, mode_x='both',
                 y0=-1.0, y1=1.0, beta_y=3.0, mode_y='both',
                 z0=-1.0, z1=1.0, beta_z=3.0, mode_z='both'):
        xs = _tanh_1d(nx, beta_x, mode_x)
        ys = _tanh_1d(ny, beta_y, mode_y)
        zs = _tanh_1d(nz, beta_z, mode_z)
        Lx, Ly, Lz = x1 - x0, y1 - y0, z1 - z0
        self.numNodes = nx * ny * nz
        self.positions = [
            [x0 + Lx * xi, y0 + Ly * yj, z0 + Lz * zk]
            for zk in zs
            for yj in ys
            for xi in xs
        ]
        self.bounds = [[x0, y0, z0], [x1, y1, z1]]
