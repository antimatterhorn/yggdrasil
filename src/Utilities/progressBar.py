import sys
import math

def ProgressBar(pct, text):
    barWidth = 20
    filledLength = int(round(pct * barWidth))

    bar = "["
    bar += "=" * filledLength
    bar += ">"
    bar += "·" * (barWidth - filledLength - 1)
    bar += "]"

    bar += f" {pct * 100:.1f}%"
    bar += f" {text}"

    # Print the progress bar and overwrite the previous line
    sys.stdout.write("\r" + bar)
    sys.stdout.flush()

    if pct >= 1.0:
        print()  # Move to the next line when done