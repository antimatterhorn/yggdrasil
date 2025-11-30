import itertools
import random

# Define possible attributes
genders = ["Boy", "Girl"]
birthdays = ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"]
birthmonths = ["January", "February", "March", "April", "May", "June",
               "July", "August", "September", "October", "November", "December"]

# Function to simulate many two-child families
def simulate(num_trials=100000):
    count_valid = 0
    count_other_is_girl = 0

    for _ in range(num_trials):
        # Randomly generate two children
        child1 = (random.choice(genders), random.choice(birthdays), random.choice(birthmonths))
        child2 = (random.choice(genders), random.choice(birthdays), random.choice(birthmonths))

        children = [child1, child2]

        # Check if at least one is a boy born on Tuesday
        if any(g == "Boy" and d == "Tuesday" and m == "April" for g, d, m in children):
            count_valid += 1
            # Check if the other child is a girl
            if child1[0] == "Girl" or child2[0] == "Girl":
                count_other_is_girl += 1

    probability = count_other_is_girl / count_valid if count_valid > 0 else 0
    print(f"Out of {num_trials} valid families, the probability the other child is a girl is approximately {probability:.4f}")

simulate()
