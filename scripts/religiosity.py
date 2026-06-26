import random
import numpy as np
from dataclasses import dataclass


@dataclass
class Person:
    id: int
    birth_year: int
    lifespan: int       # drawn from normal(80, 10)
    religiosity: int    # 0 or 1
    has_reproduced: bool = False


def child_religiosity(r1: int, r2: int) -> int:
    """Inherit religiosity from two parents according to transmission rules."""
    if r1 == 0 and r2 == 0:
        return 0 if random.random() < 0.95 else 1
    elif r1 == 1 and r2 == 1:
        return 1 if random.random() < 0.65 else 0
    else:  # mixed
        return random.randint(0, 1)


def random_lifespan() -> int:
    """Draw a lifespan from N(80, 10), minimum 21 so reproduction can occur."""
    return max(21, round(random.gauss(80, 10)))


def simulate(
    initial_population: int = 2000,
    years: int = 80,
    initial_religiosity: float = 0.6,
    seed: int = None,
) -> list[dict]:
    """
    Simulate population evolution with religiosity transmission.

    Reproduction rules:
    - People pair up at age 20 and produce 2 or 3 children (equal probability → mean 2.5).
    - Each person reproduces only once.
    - Lifespan is drawn from N(80, 10).

    Religiosity transmission:
    - Both parents 0  → child is 0 with p=0.95, else 1
    - Both parents 1  → child is 1 with p=0.65, else 0
    - Mixed parents   → child is 0 or 1 with equal probability
    """
    if seed is not None:
        random.seed(seed)

    next_id = 0
    population: list[Person] = []

    # Seed population with a spread of ages so reproduction begins immediately.
    for _ in range(initial_population):
        age = random.randint(0, 79)
        rel = 1 if random.random() < initial_religiosity else 0
        p = Person(
            id=next_id,
            birth_year=-age,
            lifespan=random_lifespan(),
            religiosity=rel,
            has_reproduced=(age >= 20),
        )
        next_id += 1
        population.append(p)

    history = []

    for year in range(years):
        # --- Deaths ---
        population = [p for p in population if (year - p.birth_year) < p.lifespan]

        # --- Reproduction at age 20 ---
        reproducers = [
            p for p in population
            if (year - p.birth_year) >= 20 and not p.has_reproduced
        ]
        random.shuffle(reproducers)

        new_children: list[Person] = []
        for i in range(0, len(reproducers) - 1, 2):
            p1 = reproducers[i]
            p2 = reproducers[i + 1]
            p1.has_reproduced = True
            p2.has_reproduced = True

            num_children = random.choice([2, 3])   # mean = 2.5
            for _ in range(num_children):
                child = Person(
                    id=next_id,
                    birth_year=year,
                    lifespan=random_lifespan(),
                    religiosity=child_religiosity(p1.religiosity, p2.religiosity),
                )
                next_id += 1
                new_children.append(child)

        population.extend(new_children)

        # --- Record stats ---
        n = len(population)
        if n > 0:
            n_religious = sum(p.religiosity for p in population)
            history.append({
                "year": year,
                "population": n,
                "n_religious": n_religious,
                "n_secular": n - n_religious,
                "religiosity_fraction": n_religious / n,
            })

    return history


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    print("Running simulation...")
    history = simulate(initial_population=2000, years=200, initial_religiosity=0.6, seed=42)

    years_      = [h["year"]                for h in history]
    pop_        = [h["population"]          for h in history]
    rel_frac    = [h["religiosity_fraction"] for h in history]

    # Theoretical equilibrium for religiosity fraction.
    # Solving f = 0.05*(1-f)^2 + 0.65*f^2 + 0.5*2f(1-f) gives f ~ 0.274
    theoretical_eq = (-0.1 + (0.07 ** 0.5)) / 0.6
    print(f"Theoretical religiosity equilibrium: {theoretical_eq:.3f}")
    print(f"Year   0 — population: {pop_[0]:>6},  religiosity: {rel_frac[0]:.3f}")
    print(f"Year {years_[-1]:>3} — population: {pop_[-1]:>6},  religiosity: {rel_frac[-1]:.3f}")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    ax1.plot(years_, pop_, color="steelblue")
    ax1.set_ylabel("Population size")
    ax1.set_title("Population over time")
    ax1.grid(True, alpha=0.3)

    ax2.plot(years_, rel_frac, color="firebrick", label="Simulated")
    ax2.axhline(theoretical_eq, color="gray", linestyle="--", label=f"Equilibrium ≈ {theoretical_eq:.3f}")
    ax2.set_xlabel("Year")
    ax2.set_ylabel("Fraction religiosity = 1")
    ax2.set_title("Religiosity fraction over time")
    ax2.set_ylim(0, 1)
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    out = "religiosity_simulation.png"
    #plt.savefig(out, dpi=150)
    print(f"Plot saved to {out}")
    plt.show()
