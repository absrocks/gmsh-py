import numpy as np
import gmsh
import sys
from scipy.optimize import curve_fit
gmsh.initialize(sys.argv)
gmsh.model.add("beach_profile")

def dim_tag(entities, n):
    dim = []
    tag = []
    for e in entities:
        # Dimension and tag of the entity:
        dim = e[0]
        if dim == n:
            tag.append(e[1])
            print("tag is", tag, "for corresponding dimension", dim)

    return tag
# Parameters
x1, x2, x0 = 76, 104, 1
h0, h1 = 1.5, 9.0
y = 0.8
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

lc = 1.0
profile_pts = []
curve_pts_front = []
curve_pts_back = []
point_id = 1  # start unique point ID counter

# Beach profile points
for x, z in zip(x_vals, h_vals):

    curve_pts_back.append(gmsh.model.geo.addPoint(x, 0, z, lc / 5))
    curve_pts_front.append(gmsh.model.geo.addPoint(x, y, z, lc / 5))

# Define Points
p_bbrc = gmsh.model.geo.addPoint(x2, 0, h_max, lc/5) # Bottom back right corner
#p_bfrc = gmsh.model.geo.addPoint(x2, y, h_max, lc/4) # Bottom front right corner
#p_tbrc = gmsh.model.geo.addPoint(x2, 0, h_min, lc) # Top back right corner
#p_tfrc = gmsh.model.geo.addPoint(x2, y, h_min, lc) # Top front right corner
#p_tblc = gmsh.model.geo.addPoint(x0, 0, h_min, lc) # Top back left corner
#p_tflc = gmsh.model.geo.addPoint(x0, y, h_min, lc) # Top front left corner

# Draw spline based on the points
spl_back = gmsh.model.geo.addSpline(curve_pts_back)

# Draw the other lines on bottom surface
l_bs = gmsh.model.geo.addLine(curve_pts_back[-1], p_bbrc) # Connecting line between spline end point and bottom back right corner

# Extrude to Y
spl_ext = gmsh.model.geo.extrude([(1, spl_back), (1, l_bs)], 0, y, 0, recombine=True)

# Make surface transfinite and recombine for quadrangles
gmsh.model.geo.synchronize()
surface = dim_tag(gmsh.model.getEntities(2), 2)
for s in surface:
    gmsh.model.geo.mesh.setTransfiniteSurface(s)
    gmsh.model.geo.mesh.setRecombine(2, s)

# Add boundary layer
N = 12 # number of layers
r = 1.1 # ratio
s = [-0.02] # thickness of first layer
d = [0.02]
for i in range(1, N): d.append(d[-1] + d[0] * r**i)
for i in range(1, N): s.append(s[-1] + s[0] * r**i)
print("height of the last layer of bl: ", max(d))
print("sum", sum(d))
extbl = gmsh.model.geo.extrudeBoundaryLayer(gmsh.model.getEntities(2), [1] * N, s, True)

# get "top" surfaces created by extrusion
print("extbl", extbl)

# Extrude the top surface of boundary layer
gmsh.model.geo.synchronize()
top = []
for i in range(1, len(extbl)):
    if extbl[i][0] == 3:
        top.append(extbl[i-1])
print("top", top)
ext_3d = gmsh.model.geo.extrude([(2, 54)], 0, 0, 5, recombine=True)

'''
left_top = gmsh.model.geo.addPoint(x0, 0, -2, lc, point_id)
base_left = point_id
point_id += 1

right_top = gmsh.model.geo.addPoint(x2, 0,-2, lc, point_id)
base_right = point_id
point_id += 1


bottom_beach_right = gmsh.model.geo.addPoint(x2, 0, h_max, lc/4, point_id)
base_bottom = point_id
point_id += 1

N = 6 # number of layers
r = 1.2 # ratio
d = [-100] # thickness of first layer
for i in range(1, N): d.append(d[-1] - (-d[0]) * r**i)

extbl = gmsh.model.geo.extrudeBoundaryLayer([(1,spl), (1,bottom_line)],
                                            [1] * N, d, True)

gmsh.model.addPhysicalGroup(2, [1], 1)
gmsh.model.setPhysicalName(2, 1, "back")

'''






#gmsh.model.addPhysicalGroup(2, [dim_tag(spl_ext , 2) ,b1], name="bottom")



'''
h = 0.8
#ov = gmsh.model.geo.extrude([(2, 1)], 0, h, 0, [2, 2], [0.5, 1])

ov2 = gmsh.model.geo.extrude([(2, 1)], 0, h, 0)
#print([ov[5][1]],len(ov), ov[0], ov)
#print(ov2)
#ov = gmsh.model.geo.extrude([(2, 1)], 0, 0, h)
gmsh.model.geo.synchronize()
# Add physical groups for all 6 surfaces
gmsh.model.addPhysicalGroup(2, [ov[-3][1]], tag=1202)  # Top surface
gmsh.model.setPhysicalName(2, 1202, "right")

gmsh.model.addPhysicalGroup(2, [ov[0][1]], tag=1203)  # Side 1
gmsh.model.setPhysicalName(2, 1203, "front")

gmsh.model.addPhysicalGroup(2, [ov[-2][1]], tag=1205)  # Side 3
gmsh.model.setPhysicalName(2, 1205, "top")

gmsh.model.addPhysicalGroup(2, [ov[-1][1]], tag=1206)  # Side 4
gmsh.model.setPhysicalName(2, 1206, "left")
# Collect all relevant surface tags into a list
bottom_surfaces = [ov[i][1] for i in range(2, len(ov)-3)]

# Add them all to a single physical group
gmsh.model.addPhysicalGroup(2, bottom_surfaces, tag=1204)
gmsh.model.setPhysicalName(2, 1204, "bottom")

#gmsh.model.addPhysicalGroup(3, [1, 2, ov[1][1]], 1001)
'''
mesh = '2D'
if mesh:
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)
    gmsh.write("beach_profile.msh")
else:
    gmsh.model.mesh.generate(3)
    gmsh.write("beach_3d.msh")


if "-nopopup" not in sys.argv:
    gmsh.fltk.run()
gmsh.finalize()
