import numpy as np

n_inputs = 5
n_outputs = 3

# Ground truth (hidden from solver)
At = np.random.randn(n_outputs, n_inputs)  # True transformation matrix
bt = np.random.randn(n_outputs)            # True bias vector

# Black box function using fixed At and bt
def black_box(x):
    return At @ x + bt

# Step 1: Estimate bias vector b
x0 = np.zeros(n_inputs)
y0 = black_box(x0)
b = y0

# Step 2: Estimate matrix A
A = np.zeros((n_outputs, n_inputs))
for i in range(n_inputs):
    xi = np.zeros(n_inputs)
    xi[i] = 1.0
    yi = black_box(xi)
    A[:, i] = yi - b

# Show results
print("Estimated A:\n", A)
print("Estimated b:\n", b)

# Compare to true At and bt
print("\nTrue A (At):\n", At)
print("True b (bt):\n", bt)
