# ats.py

import importlib
import scripts
import difflib

# ANSI color codes
RED    = "\033[91m"
GREEN  = "\033[92m"
YELLOW = "\033[93m"
RESET  = "\033[0m"
CYAN   = "\033[96m"

# Two kinds of test, distinguished by the dict each run() returns:
#
#   {"mode": "analytic", "checks": [(label, ok_bool, detail_str), ...]}
#       Self-verifying: the test compares itself against an exact/analytical
#       solution (with its own tolerances) and reports one pass/fail per check.
#       No stored reference is needed -- these prove correctness, not just
#       "unchanged".
#
#   {"mode": "snapshot", "values": [float, ...], "tol": rel_tol}
#       Regression: the returned values are diffed against ats_reference.txt.
#       For problems with no closed form (chaotic N-body, pattern formation),
#       this only proves the output hasn't drifted. Compared with a relative
#       tolerance (default below) rather than exact equality, since OpenMP
#       reductions are not bit-reproducible run to run.
#
# A bare list return is treated as a snapshot for backward compatibility.

SNAPSHOT_RTOL = 1e-6
SNAPSHOT_ATOL = 1e-9


def load_reference_file(filename="ats_reference.txt"):
    references = {}
    current_test = None
    current_output = []
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                if current_test is not None:
                    references[current_test] = current_output
                current_test = line[2:]
                current_output = []
            elif line:
                current_output.append(eval(line))  # assumes trusted source
        if current_test is not None:
            references[current_test] = current_output
    return references


def snapshot_matches(expected, actual, rtol=SNAPSHOT_RTOL, atol=SNAPSHOT_ATOL):
    if expected is None or len(expected) != len(actual):
        return False
    return all(abs(a - e) <= atol + rtol * abs(e) for e, a in zip(expected, actual))


def print_diff(expected, actual):
    diff = difflib.unified_diff([repr(x) for x in expected], [repr(x) for x in actual],
                                fromfile="expected", tofile="actual", lineterm="")
    print(f"{CYAN}--- Diff ---{RESET}")
    for line in diff:
        if line.startswith('+'):
            print(GREEN + line + RESET)
        elif line.startswith('-'):
            print(RED + line + RESET)
        elif line.startswith('@@'):
            print(YELLOW + line + RESET)
        else:
            print(line)


def main():
    references = load_reference_file()
    failures = []

    for name in scripts.tests:
        print(f"Testing {name}...")
        module = importlib.import_module(name)
        if not (hasattr(module, 'run') and callable(module.run)):
            print(f"{YELLOW}⚠️ {name} has no callable run(){RESET}")
            failures.append((name, "snapshot", None, None))
            continue

        result = module.run()
        mode = result.get("mode", "snapshot") if isinstance(result, dict) else "snapshot"

        if mode == "analytic":
            checks = result["checks"]
            for label, ok, detail in checks:
                mark = f"{GREEN}[PASS]{RESET}" if ok else f"{RED}[FAIL]{RESET}"
                print(f"   {mark} {label}: {detail}")
            if all(ok for _, ok, _ in checks):
                print(f"{GREEN}✅ {name} passed.{RESET}")
            else:
                print(f"{RED}❌ {name} failed.{RESET}")
                failures.append((name, "analytic", None, None))
        else:  # snapshot
            actual = result["values"] if isinstance(result, dict) else result
            expected = references.get(name)
            if snapshot_matches(expected, actual):
                print(f"{GREEN}✅ {name} passed (snapshot).{RESET}")
            else:
                print(f"{RED}❌ {name} failed (snapshot).{RESET}")
                failures.append((name, "snapshot", expected, actual))

    print("\nSummary:")
    if not failures:
        print(f"{GREEN}🎉 All tests passed.{RESET}")
    else:
        for name, mode, expected, actual in failures:
            print(f"\n{RED}- {name} failed.{RESET}")
            if mode == "snapshot" and expected is not None and actual is not None:
                print_diff(expected, actual)
        print(f"\n{RED}{len(failures)} test(s) failed.{RESET}")


if __name__ == "__main__":
    main()
