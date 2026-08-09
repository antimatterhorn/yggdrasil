# Copyright (C) 2026  Cody Raskin

import sys
import math

from Colors import BRIGHT_GREEN, BRIGHT_WHITE, DIM, GREEN, colorize

def ProgressBar(pct, text):
    barWidth = 20
    filledLength = int(round(pct * barWidth))

    bar = colorize("[", DIM)
    bar += colorize("=" * filledLength, GREEN)
    bar += colorize(">", BRIGHT_GREEN)
    bar += colorize("·" * (barWidth - filledLength - 1), DIM)
    bar += colorize("]", DIM)

    bar += colorize(f" {pct * 100:.1f}%", BRIGHT_WHITE)
    bar += f" {text}"

    # Print the progress bar and overwrite the previous line
    sys.stdout.write("\r" + bar)
    sys.stdout.flush()

    if pct >= 1.0:
        print()  # Move to the next line when done
