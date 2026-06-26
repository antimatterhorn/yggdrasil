import numpy as np
from scipy.spatial import Voronoi


class CVTNodeGenerator2d:
    def __init__(self, numPoints, n_iterations=10):
        from ConstantDThetaGenerator import ConstantDThetaDisk2d
        points = np.array(ConstantDThetaDisk2d(numPoints).positions, dtype=float)

        print("performing CVT")

        for _ in range(n_iterations):
            # Mirror across all 4 walls of [-1,1]^2 so every boundary cell is
            # bounded and has a well-defined centroid.
            mirrors = np.vstack([
                points,
                np.column_stack([-2.0 - points[:, 0], points[:, 1]]),
                np.column_stack([ 2.0 - points[:, 0], points[:, 1]]),
                np.column_stack([points[:, 0], -2.0 - points[:, 1]]),
                np.column_stack([points[:, 0],  2.0 - points[:, 1]]),
            ])
            vor = Voronoi(mirrors)
            n = len(points)
            for i in range(n):
                region = vor.regions[vor.point_region[i]]
                if -1 not in region and len(region) >= 3:
                    centroid = vor.vertices[region].mean(axis=0)
                    points[i] = np.clip(centroid, -1.0, 1.0)

        self.positions = points.tolist()


class CVTDiskGenerator2d:
    def __init__(self, numPoints, n_iterations=10):
        from PoissonNodeGenerator import PoissonDisk2d
        points = np.array(PoissonDisk2d(numPoints).positions, dtype=float)

        print("performing CVT in disk")

        for _ in range(n_iterations):
            vor = Voronoi(points)
            for i, region_idx in enumerate(vor.point_region):
                region = vor.regions[region_idx]
                if -1 not in region and len(region) >= 3:
                    centroid = vor.vertices[region].mean(axis=0)
                    r2 = centroid[0]**2 + centroid[1]**2
                    if r2 > 1.0:
                        centroid /= np.sqrt(r2)
                    points[i] = centroid

        self.positions = points.tolist()
