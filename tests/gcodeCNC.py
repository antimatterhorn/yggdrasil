#!/usr/bin/env python3
"""
gcodeCNC.py

Read a G-code file, simulate the CNC toolpath kinematically, and plot it.
Supports:
  - G0 / G00 : rapid move
  - G1 / G01 : linear feed move
  - G90      : absolute positioning
  - G91      : relative positioning
  - F        : feed rate (assumed in units/min, converted to units/sec)

We only care about X/Y here; Z, arcs, etc. are ignored for now.
"""

import argparse
import math
import re

import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

class Move:
    def __init__(self, start, end, feed, rapid=False):
        self.start = start  # (x, y)
        self.end = end      # (x, y)
        self.feed = feed    # units/sec
        self.rapid = rapid  # bool


# ---------------------------------------------------------------------------
# G-code parsing helpers
# ---------------------------------------------------------------------------

def strip_comments(line: str) -> str:
    """Remove ';' and '(...)' style comments from a line."""
    # Remove ; comments
    line = line.split(';', 1)[0]
    # Remove (...) comments
    line = re.sub(r'\(.*?\)', '', line)
    return line.strip()


def parse_gcode(path: str):
    """
    Parse a G-code file into a list of Move objects.

    Returns:
        moves: list[Move]
    """
    moves = []

    # State
    x = 0.0
    y = 0.0
    z = 0.0
    absolute_mode = True  # G90 default
    feed_rate = 1000.0    # units/min (default)
    feed_rate_sec = feed_rate / 60.0

    with open(path, "r") as f:
        for raw_line in f:
            line = strip_comments(raw_line)
            if not line:
                continue

            # Make parsing case-insensitive
            uline = line.upper()

            # --- Extract all G/M codes like G0, G00, G1, G90, G91, M03, etc. ---
            code_matches = re.findall(r'([GM])(\d+)', uline)
            g_codes = [int(num) for (gm, num) in code_matches if gm == 'G']
            m_codes = [int(num) for (gm, num) in code_matches if gm == 'M']

            # --- Extract all coordinate / feed words X, Y, Z, F ---
            # Matches things like X17.733, Y-1.0, Z3, F1000.0
            word_matches = re.findall(r'([XYZF])\s*([-+]?\d*\.?\d+)', uline)

            x_val = None
            y_val = None
            z_val = None
            f_val = None

            for var, val in word_matches:
                v = float(val)
                if var == 'X':
                    x_val = v
                elif var == 'Y':
                    y_val = v
                elif var == 'Z':
                    z_val = v
                elif var == 'F':
                    f_val = v

            # --- Handle modal changes (G90/G91, units, etc.) ---
            for g in g_codes:
                if g == 90:      # G90 absolute
                    absolute_mode = True
                elif g == 91:    # G91 relative
                    absolute_mode = False
                elif g == 20:    # G20 inches (we'll just treat as "units", no scaling)
                    pass
                elif g == 21:    # G21 mm (same: just "units")
                    pass

            # Update feed if given
            if f_val is not None:
                feed_rate = f_val          # units/min
                feed_rate_sec = feed_rate / 60.0

            # --- Determine if this line commands a motion (G0/G1) ---
            motion_type = None
            for g in g_codes:
                if g == 0:        # G0 / G00 rapid
                    motion_type = 'G0'
                    break
                elif g == 1:      # G1 / G01 feed
                    motion_type = 'G1'
                    break

            if motion_type is None:
                # No motion command on this line we care about
                continue

            # --- Compute target position ---
            old_x, old_y, old_z = x, y, z

            if absolute_mode:
                if x_val is not None:
                    x = x_val
                if y_val is not None:
                    y = y_val
                if z_val is not None:
                    z = z_val
            else:
                if x_val is not None:
                    x += x_val
                if y_val is not None:
                    y += y_val
                if z_val is not None:
                    z += z_val

            new_x, new_y, new_z = x, y, z

            # If the XY position didn't change, nothing to do for our 2D path
            if new_x == old_x and new_y == old_y:
                continue

            rapid = (motion_type == 'G0')
            # Some G-code may specify G0 with no explicit feed; we still need a nonzero speed
            move_feed = feed_rate_sec if feed_rate_sec > 0.0 else 1.0

            moves.append(Move(start=(old_x, old_y),
                              end=(new_x, new_y),
                              feed=move_feed,
                              rapid=rapid))

    return moves


# ---------------------------------------------------------------------------
# Kinematic sampling
# ---------------------------------------------------------------------------

def simulate_moves(moves, dt):
    """
    Given a list of Move objects and timestep dt,
    generate a sampled trajectory (xs, ys).

    dt is in seconds; feed is in units/sec.
    """
    if not moves:
        return [0.0], [0.0]

    xs = []
    ys = []

    # Start from the start of the first move
    current_x, current_y = moves[0].start
    xs.append(current_x)
    ys.append(current_y)

    for move in moves:
        sx, sy = move.start
        ex, ey = move.end

        dx = ex - sx
        dy = ey - sy
        seg_len2 = dx*dx + dy*dy
        if seg_len2 == 0.0:
            continue

        seg_len = math.sqrt(seg_len2)
        ux = dx / seg_len
        uy = dy / seg_len

        speed = move.feed  # units/sec

        seg_time = seg_len / speed
        n_steps = max(1, int(seg_time / dt))

        # Step along the segment
        for _ in range(n_steps):
            current_x += ux * speed * dt
            current_y += uy * speed * dt
            xs.append(current_x)
            ys.append(current_y)

        # Snap exactly to the endpoint
        current_x, current_y = ex, ey
        xs.append(current_x)
        ys.append(current_y)

    return xs, ys


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_path(moves, xs, ys, title=None):
    fig, ax = plt.subplots()

    # Plot sampled path
    ax.plot(xs, ys, "-", linewidth=2, label="Simulated path")

    if moves:
        sx, sy = moves[0].start
        ex, ey = moves[-1].end
        ax.scatter([sx], [sy], marker="o", label="Start")
        ax.scatter([ex], [ey], marker="x", label="End")

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.grid(True)
    if title:
        ax.set_title(title)
    ax.legend()

    plt.show()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Simulate and plot a G-code toolpath."
    )
    parser.add_argument("gcode",
                        help="Path to the G-code file.")
    parser.add_argument("--dt", type=float, default=0.01,
                        help="Timestep for kinematic sampling (sec). Default: 0.01")
    args = parser.parse_args()

    moves = parse_gcode(args.gcode)
    if not moves:
        print("No motion commands (G0/G1) found in G-code.")
        return

    xs, ys = simulate_moves(moves, args.dt)
    title = f"G-code Toolpath: {args.gcode}"
    plot_path(moves, xs, ys, title=title)


if __name__ == "__main__":
    main()
