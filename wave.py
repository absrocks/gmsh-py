import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Given constants
x0, x1, xr = 30, 104, 0
h0, h1 = 1.5, 9.0
n_tau = 1.5
n_epsilon = 0.75 * (n_tau + 1)
beta = (n_epsilon - 1) / 2
exponent = 1 + beta
alpha = (h1 ** exponent - h0 ** exponent) / (x1 - x0)

# Generate data using the original model
x_vals = np.linspace(x0, x1, 300)
h_vals = ((x_vals * alpha) + h0 ** exponent) ** (1 / exponent)

# Define fitting function: h = A * (x - b)^n
def power_law(x, A, b, n):
    return A * np.power(x - b, n)

# Initial guess: A, b, n
initial_guess = [1.0, 0.0, 0.3]
bounds = ([0, -50, 0.25], [np.inf, 50, 0.41])  # constrain n in the given range

# Perform the curve fit
popt, _ = curve_fit(power_law, x_vals, h_vals, p0=initial_guess, bounds=bounds)

# Extract fitted parameters
A_fit, b_fit, n_fit = popt
print(f"Best fit parameters:\nA = {A_fit:.4f}, b = {b_fit:.4f}, n = {n_fit:.4f}")

# Plot
plt.plot(x_vals, h_vals, label='Original Model', linewidth=2)
plt.plot(x_vals, power_law(x_vals, *popt), '--', label=f'Fitted: A(x - b)^n\nn={n_fit:.3f}', linewidth=2)
plt.xlabel('x')
plt.ylabel('h(x)')
plt.title('Fitting A(x - b)^n to Original Profile')

plt.gca().invert_yaxis()
plt.grid(True)
plt.legend()
plt.show()
