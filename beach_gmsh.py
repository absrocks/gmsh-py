import numpy as np
import gmsh
import sys
from scipy.optimize import curve_fit
gmsh.initialize(sys.argv)
gmsh.model.add("beach_profile")

# Parameters
x1, x2, x0 = 76, 104, 1
h0, h1 = 1.5, 9.0
n_tau = 2
n_epsilon = 0.75 * (n_tau + 1)
beta = (n_epsilon - 1) / 2
exponent = 1 + beta
alpha = (h1 ** exponent - h0 ** exponent) / (x1 - x0)

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
h_vals = power_law(x_vals, *popt)

x_min, h_min = min(x_vals), min(h_vals)
x_max, h_max = max(x_vals), max(h_vals)
point_id = 1  # start unique point ID counter
lc = 1.0
profile_pts = []

# Beach profile points
for x, z in zip(x_vals, h_vals):
    if x == x1 and z == h_max or x == x0 and z == h_min:
        gmsh.model.geo.addPoint(x, 0, z, lc / 4, point_id)

    else:
        gmsh.model.geo.addPoint(x, 0, z, lc, point_id)
    profile_pts.append(point_id)
    point_id += 1

# Flat bottom section (shoreline to left)
for x in np.arange(x1, x2, lc/4):
    gmsh.model.geo.addPoint(x, 0, h_max, lc / 4, point_id)
    profile_pts.append(point_id)
    point_id += 1

# Base rectangle points (not added to profile_pts!)
gmsh.model.geo.addPoint(x0, 0, 0, lc, point_id)
base_left = point_id
point_id += 1

gmsh.model.geo.addPoint(x2, 0, 0, lc, point_id)
base_right = point_id
point_id += 1

# Profile lines
line_tags = []

for i in range(len(profile_pts) - 1):
    gmsh.model.geo.addLine(profile_pts[i], profile_pts[i+1], i+1)
    line_tags.append(i+1)

top_line = len(line_tags) + 1
gmsh.model.geo.addLine(base_right, base_left, top_line)

right_line = top_line + 1
gmsh.model.geo.addLine(base_right, profile_pts[-1], right_line)

left_line = right_line + 1
#print(profile_pts[-1], base_left, left_line)
gmsh.model.geo.addLine(profile_pts[0], base_left, left_line)
curve_loop = line_tags + [left_line, top_line, right_line]
#gmsh.model.geo.addCurveLoop(curve_loop, 1)
#gmsh.model.geo.addPlaneSurface([1], 1)

gmsh.option.setNumber("Geometry.PointNumbers", 1)
#gmsh.option.setNumber("Geometry.LineNumbers", 1)
gmsh.model.geo.synchronize()
gmsh.model.addPhysicalGroup(2, [1], 1)
#gmsh.model.mesh.generate(2)
gmsh.write("beach_profile.msh")


if "-nopopup" not in sys.argv:
    gmsh.fltk.run()
gmsh.finalize()
