import numpy as np
import gmsh
import sys

gmsh.initialize(sys.argv)
gmsh.model.add("beach_profile")

# Parameters
x0, x1 = 0, 70
h0, h1 = 1.5, 9.0
n_tau = 1.5
n_epsilon = 0.75 * (n_tau + 1)
beta = (n_epsilon - 1) / 2
exponent = 1 + beta
alpha = (h1 ** exponent - h0 ** exponent) / (x1 - x0)

x_vals = np.linspace(x0, x1, 300)
h_vals = ((x_vals * alpha) + h0 ** exponent) ** (1 / exponent)

lc = 1.0
profile_pts = []
for i, (x, z) in enumerate(zip(x_vals, h_vals)):
    if x == x1 and z == h1:
        gmsh.model.geo.addPoint(x1, 0, h1, lc / 4, i+1)
    elif x == x0 and z == h0:
        gmsh.model.geo.addPoint(x0, 0, h0, lc / 4, i+1)
    else:
        gmsh.model.geo.addPoint(x, 0, z, lc, i+1)
    profile_pts.append(i+1)

# Add bottom rectangle corners
base_left = len(profile_pts) + 1
base_right = base_left + 1
gmsh.model.geo.addPoint(x0, 0, h1+4, lc, base_left)
gmsh.model.geo.addPoint(x1, 0, h1+4, lc, base_right)

# Profile lines
line_tags = []
for i in range(len(profile_pts) - 1):
    gmsh.model.geo.addLine(profile_pts[i], profile_pts[i+1], i+1)
    line_tags.append(i+1)

bottom_line = len(line_tags) + 1
gmsh.model.geo.addLine(base_right, base_left, bottom_line)

right_line = bottom_line + 1
gmsh.model.geo.addLine(profile_pts[-1], base_right, right_line)

left_line = right_line + 1
gmsh.model.geo.addLine(base_left, profile_pts[0], left_line)
curve_loop = line_tags + [right_line, bottom_line, left_line]
gmsh.model.geo.addCurveLoop(curve_loop, 1)
gmsh.model.geo.addPlaneSurface([1], 1)

gmsh.model.geo.synchronize()
gmsh.model.addPhysicalGroup(2, [1], 1)
gmsh.model.mesh.generate(2)
gmsh.write("beach_profile.msh")

if "-nopopup" not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
