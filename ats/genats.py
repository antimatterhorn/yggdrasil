import importlib
import scripts


def generate():
    """Regenerate ats_reference.txt. Only snapshot-mode tests are recorded --
    analytic tests verify themselves against an exact solution and need no
    stored reference (see ats.py)."""
    results = {}
    for name in scripts.tests:
        print(f"Running {name}...")
        module = importlib.import_module(name)
        if not (hasattr(module, 'run') and callable(module.run)):
            print(f"Warning: {name} has no callable run() method")
            continue
        result = module.run()
        mode = result.get("mode", "snapshot") if isinstance(result, dict) else "snapshot"
        if mode == "analytic":
            print(f"  (analytic, self-verifying -- no reference stored)")
            continue
        results[name] = result["values"] if isinstance(result, dict) else result

    with open("ats_reference.txt", "w") as f:
        for name, output in results.items():
            f.write(f"# {name}\n")
            for item in output:
                f.write(f"{repr(item)}\n")
            f.write("\n")

    print("Done. Snapshot reference outputs saved to ats_reference.txt")


if __name__ == "__main__":
    generate()
