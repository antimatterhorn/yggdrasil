# Copyright (C) 2026  Cody Raskin

import os
import sys

RESET     = "\033[0m"
BOLD      = "\033[1m"
DIM       = "\033[2m"
ITALIC    = "\033[3m"
UNDERLINE = "\033[4m"

BLACK   = "\033[30m"
RED     = "\033[31m"
GREEN   = "\033[32m"
YELLOW  = "\033[33m"
BLUE    = "\033[34m"
MAGENTA = "\033[35m"
CYAN    = "\033[36m"
WHITE   = "\033[37m"

BRIGHT_BLACK   = "\033[90m"
BRIGHT_RED     = "\033[91m"
BRIGHT_GREEN   = "\033[92m"
BRIGHT_YELLOW  = "\033[93m"
BRIGHT_BLUE    = "\033[94m"
BRIGHT_MAGENTA = "\033[95m"
BRIGHT_CYAN    = "\033[96m"
BRIGHT_WHITE   = "\033[97m"


def _detect():
    """Colors go to interactive terminals only, so redirecting a run into a log
    file yields plain text rather than embedded escape sequences. FORCE_COLOR
    overrides the check; NO_COLOR (https://no-color.org) disables it outright."""
    if os.environ.get("NO_COLOR") is not None:
        return False
    if os.environ.get("FORCE_COLOR") is not None:
        return True
    if os.environ.get("TERM") in ("dumb", None):
        return False
    return sys.stdout.isatty()


_enabled = _detect()


def enabled():
    return _enabled


def setEnabled(flag):
    global _enabled
    _enabled = bool(flag)


def rgb(r, g, b):
    """24-bit truecolor foreground code."""
    return "\033[38;2;%d;%d;%dm" % (r, g, b)


def rgbBackground(r, g, b):
    return "\033[48;2;%d;%d;%dm" % (r, g, b)


def color256(n):
    return "\033[38;5;%dm" % n


def colorize(text, *codes):
    """Wrap text in the given SGR codes, or return it untouched when color is off."""
    if not _enabled or not codes:
        return text
    return "".join(codes) + text + RESET


def strip(text):
    """Remove SGR sequences, for measuring the visible width of a colored string."""
    out = []
    i = 0
    while i < len(text):
        if text[i] == "\033":
            j = text.find("m", i)
            if j < 0:
                break
            i = j + 1
        else:
            out.append(text[i])
            i += 1
    return "".join(out)


def visibleLength(text):
    return len(strip(text))


def gradient(text, start, end):
    """Interpolate a truecolor foreground across the characters of text."""
    if not _enabled:
        return text
    n = len(text)
    if n <= 1:
        return colorize(text, rgb(*start))
    out = []
    for i, ch in enumerate(text):
        t = i / (n - 1)
        out.append(rgb(*[int(round(s + (e - s) * t)) for s, e in zip(start, end)]))
        out.append(ch)
    out.append(RESET)
    return "".join(out)
