from yggdrasil import *
from Mesh import Grid2d, Geometry
import math

# Isolated test of the CylindricalRZ grid metric (no physics). Verifies:
#   1. A Cartesian grid is unchanged (uniform cell volume dx*dy).
#   2. An RZ grid weights cell volume by radius (V = r*dr*dz, per radian).
#   3. Summed RZ cell volume matches the analytic half-disk*length: 0.5*R^2*L.
#   4. cellRadius is strictly positive everywhere (no p/r singularity).

def approx(a, b, tol=1e-9):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))

if __name__ == "__main__":
    nr, nz = 40, 10
    dr, dz = 0.05, 0.1          # R = nr*dr = 2.0, L = nz*dz = 1.0
    R = nr * dr
    L = nz * dz

    # 1. Cartesian: uniform cell volume
    cart = Grid2d(nr, nz, dr, dz)                       # default Geometry.Cartesian
    assert cart.geometry() == Geometry.Cartesian
    for (i, j) in ((0, 0), (nr // 2, nz // 2), (nr - 1, nz - 1)):
        idx = cart.index(i, j, 0)
        assert approx(cart.cellVolume(idx), dr * dz), \
            f"cartesian cellVolume[{idx}]={cart.cellVolume(idx)} != {dr*dz}"
    print(f"[ok] Cartesian: uniform cellVolume = {dr*dz}")

    # 2. RZ: volume weighted by radius; cellRadius = (i+0.5)*dr on axis 0
    rz = Grid2d(nr, nz, dr, dz, Geometry.CylindricalRZ)
    assert rz.geometry() == Geometry.CylindricalRZ
    for i in range(nr):
        idx = rz.index(i, 0, 0)
        r_expect = (i + 0.5) * dr
        assert approx(rz.cellRadius(idx), r_expect), \
            f"cellRadius[{i}]={rz.cellRadius(idx)} != {r_expect}"
        assert approx(rz.cellVolume(idx), r_expect * dr * dz), \
            f"RZ cellVolume[{i}]={rz.cellVolume(idx)} != {r_expect*dr*dz}"
        assert rz.cellRadius(idx) > 0.0, "cellRadius must be > 0 (no p/r singularity)"
    print("[ok] RZ: cellVolume = r*dr*dz and cellRadius > 0 everywhere")

    # 3. Summed volume = analytic 0.5*R^2*L (per radian, i.e. 2*pi omitted).
    # Logical cells only -- rz.size() includes the ghost halo, whose r<0 side
    # would corrupt the sum (cellVolume is r*dr*dz, and r extrapolates negative
    # there).
    total = sum(rz.cellVolume(rz.index(i, j, 0)) for j in range(nz) for i in range(nr))
    analytic = 0.5 * R * R * L
    rel = abs(total - analytic) / analytic
    assert rel < 1e-12, f"summed RZ volume {total} vs analytic {analytic} (rel {rel})"
    print(f"[ok] RZ total volume {total:.6f} == 0.5*R^2*L {analytic:.6f} (per radian)")

    print("\nAll RZ metric checks passed.")
