import importlib
import scripts 

def generate():
    results = {}
    for name in scripts.tests:
        print(f"Running {name}...")
        module = importlib.import_module(name)
        if hasattr(module, 'run') and callable(module.run):
            output = module.run()
            results[name] = output
        else:
            print(f"Warning: {name} has no callable run() method")

    # Save results to a file
    with open("ats_reference.txt", "w") as f:
        for name, output in results.items():
            f.write(f"# {name}\n")
            for item in output:
                f.write(f"{repr(item)}\n")
            f.write("\n")  # newline to separate test cases

    print("Done. Reference outputs saved to ats_reference.txt")

if __name__ == "__main__":
    generate()
