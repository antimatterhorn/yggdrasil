# Copyright (C) 2026  Cody Raskin

from Colors import DIM, colorize, rgb

_logo = ["               _             _ _ ",
         " _ _ ___ ___ _| |___ ___ ___|_| |",
         "| | | . | . | . |  _| .'|_ -| | |",
         "|_  |_  |_  |___|_| |__,|___|_|_|",
         "|___|___|___|"]

_version = "v0.9.0"

# Vertical gradient down the ASCII logo: canopy green to root cyan.
_top = (120, 200, 110)
_bottom = (60, 170, 200)

for i, line in enumerate(_logo):
    t = i / (len(_logo) - 1)
    tint = rgb(*[int(round(a + (b - a) * t)) for a, b in zip(_top, _bottom)])
    tail = colorize("       %s       " % _version, DIM) if i == len(_logo) - 1 else ""
    print(colorize(line, tint) + tail)
print("\n")
