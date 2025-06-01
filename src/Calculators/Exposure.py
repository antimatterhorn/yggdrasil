import math
from fractions import Fraction

def brightness(mag, ang_diam):
    """Compute surface brightness from magnitude and angular diameter."""
    return math.pow(10, 0.4 * (9.5 - mag)) / (ang_diam * ang_diam)

def f_number_prime(focal_length, diameter):
    """Calculate the f-ratio."""
    return focal_length / diameter

def exposure_duration(f_ratio, iso, sb):
    """Compute and return exposure time in a human-readable format."""
    ex = f_ratio * f_ratio / (iso * sb)
    
    if ex > 3599:
        return f"{ex / 3600:.2f} hours"
    elif ex > 59:
        return f"{ex / 60:.2f} minutes"
    elif ex < 0.25:
        frac = Fraction(ex).limit_denominator(1000)
        return f"{frac.numerator}/{frac.denominator} seconds"
    else:
        return f"{ex:.2f} seconds"

def main():
    print("\nTelescope Properties")
    focal_length = float(input("Focal Length (mm): "))
    diameter = float(input("Diameter (mm): "))
    f_ratio = f_number_prime(focal_length, diameter)
    print(f"f{f_ratio:.2f}")

    print("\nObject Properties")
    mag = float(input("Apparent Mag: "))
    ang_diam = float(input("Angular Diameter (arc-sec): "))
    sb = brightness(mag, ang_diam)

    print("\nExposure -------")
    print("ISO\tTime")
    for i in range(5):
        iso = 100 * (2 ** i)
        duration_str = exposure_duration(f_ratio, iso, sb)
        print(f"{iso}\t{duration_str}")

if __name__ == "__main__":
    main()
