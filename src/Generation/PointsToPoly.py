import numpy as np
from scipy.spatial import Delaunay

class PointsToPoly2d:
    def __init__(self, positions, angleTol=25.0):
        self.positions = positions
        self.cells = []
        self.angleTol = angleTol

        self.points = np.array(self.positions)

        tri = Delaunay(self.points)
        triangles = tri.simplices  # (ntri, 3)

        # ---------------------------------------------------------------
        # 1. Edge -> incident triangle adjacency
        # ---------------------------------------------------------------
        edge_map = {}
        for t_idx, t in enumerate(triangles):
            for a, b in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])]:
                key = frozenset((a, b))
                edge_map.setdefault(key, []).append(t_idx)

        interior_edges = [(edge, tris) for edge, tris in edge_map.items() if len(tris) == 2]


        # ---------------------------------------------------------------
        # 2. Build the ordered quad polygon for a candidate triangle pair
        # ---------------------------------------------------------------
        def opposite_vertex(triangle, shared_edge):
            """Return the vertex of `triangle` not in `shared_edge`."""
            for v in triangle:
                if v not in shared_edge:
                    return v
            raise ValueError("shared edge not found in triangle")


        def ordered_quad(edge, t1, t2, triangles):
            a, b = tuple(edge)
            c = opposite_vertex(triangles[t1], edge)
            d = opposite_vertex(triangles[t2], edge)
            # walking the quad boundary: a -> c -> b -> d -> a
            # (skip the shared diagonal a-b entirely)
            return [a, c, b, d]


        def polygon_angles(pts):
            """Interior angles (degrees) at each vertex of a polygon, given in order."""
            n = len(pts)
            angles = []
            for i in range(n):
                p_prev = pts[i - 1]
                p_cur = pts[i]
                p_next = pts[(i + 1) % n]
                v1 = p_prev - p_cur
                v2 = p_next - p_cur
                cos_a = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
                angles.append(np.degrees(np.arccos(np.clip(cos_a, -1.0, 1.0))))
            return np.array(angles)

        def is_convex(pts):
            """Check all cross products of consecutive edges have the same sign."""
            n = len(pts)
            signs = []
            for i in range(n):
                p0 = pts[i - 1]
                p1 = pts[i]
                p2 = pts[(i + 1) % n]
                e1 = p1 - p0
                e2 = p2 - p1
                cross = e1[0] * e2[1] - e1[1] * e2[0]
                signs.append(np.sign(cross))
            return all(s == signs[0] for s in signs) and signs[0] != 0


        # ---------------------------------------------------------------
        # 3. Score every candidate quad merge, keep only those passing the
        #    angle-quality gate (90 +/- ANGLE_TOL), and rank by quality
        # ---------------------------------------------------------------
        candidates = []
        for edge, (t1, t2) in interior_edges:
            quad_verts = ordered_quad(edge, t1, t2, triangles)
            pts = self.points[quad_verts]

            if not is_convex(pts):
                continue

            angles = polygon_angles(pts)
            if np.all(np.abs(angles - 90.0) <= self.angleTol):
                score = np.sum((angles - 90.0) ** 2)  # lower = closer to a perfect square/rect
                candidates.append((score, t1, t2, quad_verts))

        candidates.sort(key=lambda c: c[0])  # best quality first

        # ---------------------------------------------------------------
        # 4. Greedy merge: claim triangles for the best-scoring quads first
        # ---------------------------------------------------------------
        used = set()
        quads = []
        for score, t1, t2, quad_verts in candidates:
            if t1 in used or t2 in used:
                continue
            quads.append(quad_verts)
            used.add(t1)
            used.add(t2)

        leftover_triangles = [t for i, t in enumerate(triangles) if i not in used]

        print(f"Original triangles: {len(triangles)}")
        print(f"Merged into quads:  {len(quads)}  (consuming {2*len(quads)} triangles)")
        print(f"Leftover triangles: {len(leftover_triangles)}")

        # quad_verts entries and Delaunay's triangle rows are numpy int64,
        # not plain int -- convert both to the same plain-int-list shape
        # every other cell entry in this module uses (and what addCell's
        # std::vector<size_t> binding expects).
        self.cells = [[int(v) for v in q] for q in quads] + \
                     [[int(v) for v in t] for t in leftover_triangles]
        self.ncells = len(self.cells)