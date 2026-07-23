from yggdrasil import *
from Mesh import ALEMesh2d, Geometry
import math

# Isolated test of the CylindricalRZ ALEMesh metric (no physics), mirroring
# tests/Machinery/metricRZ.py for the structured Grid case. Verifies, on a
# hand-built quad-strip ALEMesh (unstructured machinery, not a Grid):
#   1. A Cartesian ALEMesh is unchanged (cellVolume == plain polygon area).
#   2. An RZ ALEMesh weights cell volume by centroid radius (V = r*A, exact
#      for any polygon by definition of centroid -- not an approximation).
#   3. Face area is weighted by the average of its own two endpoint node
#      radii (exact for a straight edge -- r varies linearly along it).
#   4. Summed RZ cell volume matches the analytic half-disk*length: 0.5*R^2*L
#      (per radian), same closed form as the Grid metric test since this is
#      the same rectangular (r,z) region, just built from explicit quads.

def approx(a, b, tol=1e-9):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def buildStrip(geometry, nr, nz, dr, dz):
    mesh = ALEMesh2d(geometry)

    def nodeId(i, j):
        return i + j * (nr + 1)

    for j in range(nz + 1):
        for i in range(nr + 1):
            mesh.addNode(Vector2d(float(i) * dr, float(j) * dz))

    for j in range(nz):
        for i in range(nr):
            mesh.addCell([nodeId(i, j), nodeId(i + 1, j), nodeId(i + 1, j + 1), nodeId(i, j + 1)])

    mesh.computeFaces()
    return mesh


def cellId(i, j, nr):
    return i + j * nr


if __name__ == "__main__":
    nr, nz = 40, 10
    dr, dz = 0.05, 0.1          # R = nr*dr = 2.0, L = nz*dz = 1.0
    R = nr * dr
    L = nz * dz

    # 1. Cartesian: cellVolume == plain polygon area (dr*dz for a unit quad).
    cart = buildStrip(Geometry.Cartesian, nr, nz, dr, dz)
    assert cart.geometry() == Geometry.Cartesian
    for i, j in [(0, 0), (nr // 2, nz // 2), (nr - 1, nz - 1)]:
        idx = cellId(i, j, nr)
        assert approx(cart.cellVolume(idx), dr * dz), \
            f"cartesian cellVolume[{idx}]={cart.cellVolume(idx)} != {dr*dz}"
    print(f"[ok] Cartesian ALEMesh: uniform cellVolume = {dr*dz}")

    # 2. RZ: volume weighted by centroid radius, exact for a quad too.
    rz = buildStrip(Geometry.CylindricalRZ, nr, nz, dr, dz)
    assert rz.geometry() == Geometry.CylindricalRZ
    for i in range(nr):
        idx = cellId(i, 0, nr)
        r_expect = (i + 0.5) * dr
        r_centroid = rz.cellCentroid(idx).x
        assert approx(r_centroid, r_expect), \
            f"cellCentroid.x[{i}]={r_centroid} != {r_expect}"
        assert approx(rz.cellVolume(idx), r_expect * dr * dz), \
            f"RZ cellVolume[{i}]={rz.cellVolume(idx)} != {r_expect*dr*dz}"
    print("[ok] RZ ALEMesh: cellVolume = r_centroid * plainArea")

    # 3. Face area == 0.5*(r0+r1)*edgeLength, exact -- check both an
    # r-normal face (constant z, varying r -- length dz, endpoints share the
    # same r) and a z-normal face (constant r, varying z -- length dr,
    # endpoints at r and r+dr).
    faces = rz.getFaces()
    nodes = rz.getNodes()
    checked_r_face = checked_z_face = False
    for f in faces:
        n0, n1 = f.nodeIndices[0], f.nodeIndices[1]
        r0, r1 = nodes[n0].x, nodes[n1].x
        p0, p1 = nodes[n0], nodes[n1]
        edgeLength = math.hypot(p1.x - p0.x, p1.y - p0.y)
        expectArea = 0.5 * (r0 + r1) * edgeLength
        assert approx(f.area, expectArea), \
            f"face area {f.area} != 0.5*(r0+r1)*L = {expectArea} (r0={r0}, r1={r1}, L={edgeLength})"
        if abs(r0 - r1) < 1e-12 and not checked_r_face:
            checked_r_face = True
        if abs(r0 - r1) > 1e-12 and not checked_z_face:
            checked_z_face = True
    assert checked_r_face and checked_z_face, "expected to see both face orientations in the strip"
    print("[ok] RZ ALEMesh: every face area = 0.5*(r0+r1)*edgeLength")

    # 4. Summed volume = analytic 0.5*R^2*L (per radian).
    total = sum(rz.cellVolume(cellId(i, j, nr)) for j in range(nz) for i in range(nr))
    analytic = 0.5 * R * R * L
    rel = abs(total - analytic) / analytic
    assert rel < 1e-9, f"summed RZ volume {total} vs analytic {analytic} (rel {rel})"
    print(f"[ok] RZ ALEMesh total volume {total:.6f} == 0.5*R^2*L {analytic:.6f} (per radian)")

    print("\nAll ALE RZ metric checks passed.")
