// Copyright (C) 2026  Cody Raskin
#pragma once

#include <cstdio>
#include <cstdlib>
#include <string>

#ifndef _WIN32
#include <unistd.h>
#endif

namespace Color {

// Colors go to interactive terminals only, so redirecting a run into a log file
// yields plain text rather than embedded escape sequences. FORCE_COLOR overrides
// the check; NO_COLOR (https://no-color.org) disables it outright.
inline bool detect() {
    if (std::getenv("NO_COLOR") != nullptr) return false;
    if (std::getenv("FORCE_COLOR") != nullptr) return true;
    const char* term = std::getenv("TERM");
    if (term == nullptr || std::string(term) == "dumb") return false;
#ifdef _WIN32
    return false;
#else
    return isatty(fileno(stdout)) != 0;
#endif
}

inline bool& enabledFlag() {
    static bool flag = detect();
    return flag;
}

inline bool enabled() { return enabledFlag(); }
inline void setEnabled(bool flag) { enabledFlag() = flag; }

inline const char* code(const char* seq) { return enabled() ? seq : ""; }

inline const char* reset()     { return code("\033[0m"); }
inline const char* bold()      { return code("\033[1m"); }
inline const char* dim()       { return code("\033[2m"); }
inline const char* underline() { return code("\033[4m"); }

inline const char* red()     { return code("\033[31m"); }
inline const char* green()   { return code("\033[32m"); }
inline const char* yellow()  { return code("\033[33m"); }
inline const char* blue()    { return code("\033[34m"); }
inline const char* magenta() { return code("\033[35m"); }
inline const char* cyan()    { return code("\033[36m"); }
inline const char* white()   { return code("\033[37m"); }

inline const char* brightRed()     { return code("\033[91m"); }
inline const char* brightGreen()   { return code("\033[92m"); }
inline const char* brightYellow()  { return code("\033[93m"); }
inline const char* brightBlue()    { return code("\033[94m"); }
inline const char* brightMagenta() { return code("\033[95m"); }
inline const char* brightCyan()    { return code("\033[96m"); }
inline const char* brightWhite()   { return code("\033[97m"); }

inline std::string rgb(int r, int g, int b) {
    if (!enabled()) return std::string();
    return "\033[38;2;" + std::to_string(r) + ";" + std::to_string(g) + ";" + std::to_string(b) + "m";
}

inline std::string colorize(const std::string& text, const std::string& seq) {
    if (!enabled()) return text;
    return seq + text + "\033[0m";
}

}
