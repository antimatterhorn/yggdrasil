# Copyright (C) 2026  Cody Raskin

from Colors import DIM, colorize, rgb

_logo = ["               _             _ _ ",
         " _ _ ___ ___ _| |___ ___ ___|_| |",
         "| | | . | . | . |  _| .'|_ -| | |",
         "|_  |_  |_  |___|_| |__,|___|_|_|",
         "|___|___|___|"]

_version = "v0.9.5"

# Vertical gradient down the ASCII logo: canopy green to root brown.
_top = (120, 200, 110)
_bottom = (97, 41, 18)

for i, line in enumerate(_logo):
    t = i / (len(_logo) - 1)
    tint = rgb(*[int(round(a + (b - a) * t)) for a, b in zip(_top, _bottom)])
    tail = colorize("       %s       " % _version, DIM) if i == len(_logo) - 1 else ""
    print(colorize(line, tint) + tail)
print("\n")
