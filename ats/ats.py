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

def print_diff(expected, actual):
    expected_lines = [repr(x) for x in expected]
    actual_lines = [repr(x) for x in actual]
    diff = difflib.unified_diff(
        expected_lines,
        actual_lines,
        fromfile="expected",
        tofile="actual",
        lineterm=""
    )
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

        if hasattr(module, 'run') and callable(module.run):
            actual = module.run()
            expected = references.get(name)

            if actual == expected:
                print(f"{GREEN}✅ {name} passed.{RESET}")
            else:
                print(f"{RED}❌ {name} failed.{RESET}")
                failures.append((name, expected, actual))
        else:
            print(f"{YELLOW}⚠️ {name} has no callable run(){RESET}")
            failures.append((name, None, None))

    print("\nSummary:")
    if not failures:
        print(f"{GREEN}🎉 All tests passed.{RESET}")
    else:
        for name, expected, actual in failures:
            print(f"\n{RED}- {name} failed.{RESET}")
            if expected is not None and actual is not None:
                print_diff(expected, actual)
        print(f"\n{RED}{len(failures)} test(s) failed.{RESET}")

if __name__ == "__main__":
    main()
