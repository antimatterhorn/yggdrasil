from math import sqrt
a = b = 0
for i in range(1, 101):
    a = i**3
    b = (i+1)**3
    print(f"Sum of cubes from {a} + {b} is {a+b} = {sqrt(a+b)}^2")