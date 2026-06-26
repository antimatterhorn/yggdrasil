from yggdrasil import *
from Utilities import *
import matplotlib.pyplot as plt
from Physics import CNCPathPhysics2d


def build_cnc_physics(nodeList, constants):
    """
    Create a CNCPathPhysics2d instance and add a simple rectangular toolpath.
    """
    physics = CNCPathPhysics2d(nodeList, constants)

    # Define a simple rectangular path in XY (units whatever your system uses, e.g. mm).
    # Start is implicitly whatever the NodeList position is (we'll assume (0,0)).
    # We'll also keep a mirror of these moves in Python so we can sample them.
    feed = 5.0  # units/sec

    # Corner points (absolute positions)
    corners = [
        Vector2d(10.0, 0.0),
        Vector2d(10.0, 10.0),
        Vector2d(0.0, 10.0),
        Vector2d(0.0, 0.0),
    ]

    # Program these into the C++ physics object
    for p in corners:
        physics.addLinearMove(p, feed, rapid=False)

    return physics, corners, feed


def sample_path(points, feed, dt):
    """
    Pure-Python kinematic sampler:
    given a list of *corner points* [p1, p2, ..., pN],
    a feed rate, and timestep dt, generate (x,y) samples
    along the polyline.
    """
    # Start at origin -> first point
    current = Vector2d(0.0, 0.0)
    xs = [current.x]
    ys = [current.y]

    for p in points:
        segment = p - current
        seg_len2 = segment.x**2 + segment.y**2
        if seg_len2 == 0.0:
            continue
        seg_len = seg_len2**0.5

        # Unit direction along the segment
        ux = segment.x / seg_len
        uy = segment.y / seg_len

        # Time and number of steps for this segment
        seg_time = seg_len / feed
        n_steps = max(1, int(seg_time / dt))

        # Step along the segment
        for _ in range(n_steps):
            current = Vector2d(current.x + ux * feed * dt,
                               current.y + uy * feed * dt)
            xs.append(current.x)
            ys.append(current.y)

        # Snap exactly to the endpoint
        current = Vector2d(p.x, p.y)
        xs.append(current.x)
        ys.append(current.y)

    return xs, ys


if __name__ == "__main__":
    # Basic command-line style config (tweak to taste)
    commandLine = CommandLineArguments(
        animate=True,
        siloDump=False,
        cycles=3000,
        nx=100,
        ny=20,
        dx=1.0,
        dy=1.0,
        dtmin=0.001,
    )

    constants = CGS()

    # One-node NodeList representing the toolhead
    nodeList = NodeList(1)

    # Optional: if you have a way to set the initial "position" field in Python,
    # you can do it here. Otherwise we'll assume it defaults to (0,0).
    # For example, if NodeList exposes a position field:
    # pos_field = nodeList.getFieldVector2d("position")
    # pos_field[0] = Vector2d(0.0, 0.0)

    # Build CNC physics and mirrored Python path description
    physics, corners, feed = build_cnc_physics(nodeList, constants)

    # For now, *just* sample the path kinematically in Python using dtmin
    dt = dtmin
    xs, ys = sample_path(corners, feed, dt)

    # --- Visualization ---
    fig, ax = plt.subplots()
    ax.plot(xs, ys, "-", linewidth=2)
    ax.scatter([0.0], [0.0], marker="o", label="Start")
    ax.scatter([corners[-1].x], [corners[-1].y], marker="x", label="End")

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_title("CNC Test Path (Rectangle)")
    ax.legend()
    ax.grid(True)

    plt.show()
