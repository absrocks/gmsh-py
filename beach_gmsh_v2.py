import numpy as np
import gmsh
import sys
from scipy.optimize import curve_fit

gmsh.initialize(sys.argv)
gmsh.model.add("beach_profile")


def dim_tag(entities, n):
    tag = []
    for e in entities:
        # Dimension and tag of the entity:
        dim = e[0]
        if dim == n:
            tag.append(e[1])
            print("tag is", tag, "for corresponding dimension", dim)

    return tag


def get_coords(line):
    point_ids = gmsh.model.getAdjacencies(1, line)
    if len(point_ids):
        print(" - For line", line)
        coords0 = gmsh.model.getValue(0, point_ids[1][0], [])
        print(f"Coordinates of point1 {point_ids[1][0]}: x={coords0[0]}, y={coords0[1]}, z={coords0[2]}")
        point_id0 = point_ids[1][0]
        coords1 = gmsh.model.getValue(0, point_ids[1][1], [])
        print(f"Coordinates of point2 {point_ids[1][1]}: x={coords1[0]}, y={coords1[1]}, z={coords1[2]}")
        point_id1 = point_ids[1][1]

    return coords0, coords1, point_id0, point_id1


# Parameters
x1, x2, x0 = 76, 104, 1
h1, h2, h0 = 1.5, 9.0, 0
y = 0.8
n_tau = 2
n_epsilon = 0.75 * (n_tau + 1)
beta = (n_epsilon - 1) / 2
exponent = 1 + beta
alpha = (h2 ** exponent - h1 ** exponent) / (x1 - x0)

x_vals = np.linspace(x0, x1, 300)
h_vals = ((x_vals * alpha) + h1 ** exponent) ** (1 / exponent)


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
p_bbrc = gmsh.model.geo.addPoint(x2, 0, h_max, lc / 5)  # Bottom back right corner
# p_bfrc = gmsh.model.geo.addPoint(x2, y, h_max, lc/4) # Bottom front right corner
p_tbrc = gmsh.model.geo.addPoint(x2, 0, h0, lc)  # Top back right corner
p_tfrc = gmsh.model.geo.addPoint(x2, y, h0, lc)  # Top front right corner
p_tblc = gmsh.model.geo.addPoint(x0 + 0.866, 0, h0, lc / 4)  # Top back left corner
p_tflc = gmsh.model.geo.addPoint(x0 + 0.866, y, h0, lc / 4)  # Top front left corner

# Draw spline based on the points
spl_back = gmsh.model.geo.addSpline(curve_pts_back)

# Draw the other lines on bottom surface
l_bs = gmsh.model.geo.addLine(curve_pts_back[-1],
                              p_bbrc)  # Connecting line between spline end point and bottom back right corner

# Extrude to Y
spl_ext = gmsh.model.geo.extrude([(1, spl_back), (1, l_bs)], 0, y, 0, recombine=True)

# Make surface transfinite and recombine for quadrangles
gmsh.model.geo.synchronize()
bot_surface = dim_tag(gmsh.model.getEntities(2), 2)
print("bot", bot_surface)
# for s in surface:
#    gmsh.model.geo.mesh.setTransfiniteSurface(s)
#    gmsh.model.geo.mesh.setRecombine(2, s)

# Add boundary layer
N = 12  # number of layers
r = 1.1  # ratio
s = [-0.02]  # thickness of the first layer
d = [0.02]
for i in range(1, N): d.append(d[-1] + d[0] * r ** i)
for i in range(1, N): s.append(s[-1] + s[0] * r ** i)
print("height of the last layer of bl: ", max(d))
print("sum", sum(d))
extbl = gmsh.model.geo.extrudeBoundaryLayer(gmsh.model.getEntities(2), [1] * N, s, True)

# Extrude the top surface of the boundary layer
gmsh.model.geo.synchronize()
bl_top = []
for i in range(1, len(extbl)):
    if extbl[i][0] == 3:
        bl_top.append(extbl[i - 1])

# Draw the other lines on top surface
l_tf = gmsh.model.geo.addLine(p_tflc, p_tfrc)  # Connecting line between top front left and right corner
l_tr = gmsh.model.geo.addLine(p_tfrc, p_tbrc)  # Connecting line between top right front and back corner
l_tb = gmsh.model.geo.addLine(p_tbrc, p_tblc)  # Connecting line between top back right and left corner
l_tl = gmsh.model.geo.addLine(p_tblc, p_tflc)  # Connecting line between top left front and back corner

# Make loop and surface
ll_top = gmsh.model.geo.addCurveLoop([l_tf, l_tr, l_tb, l_tl])
s_top = gmsh.model.geo.addPlaneSurface([ll_top])

# Extract lines and coordinates from top surface of boundary layer

bl_top_tag = dim_tag(bl_top, 2)
print("bl_top", bl_top_tag)
gmsh.model.geo.synchronize()
for i in range(len(bl_top)):
    up, lines = gmsh.model.getAdjacencies(2, bl_top_tag[i])
    print("Adjacent lines", lines, "for surface", bl_top_tag[i])
    if len(lines):
        for j in range(len(lines)):
            coords0, coords1, point_id0, point_id1 = get_coords(lines[j])
            if coords0[0] == x2 and coords1[0] == x2:
                if coords0[1] == y:
                    l_trf = gmsh.model.geo.addLine(point_id0,
                                                   p_tfrc)  # Connecting line between front right top and bottom point
                    l_trb = gmsh.model.geo.addLine(p_tbrc,
                                                   point_id1)  # Connecting line between back right top and bottom point
                elif coords1[1] == y:
                    l_trf = gmsh.model.geo.addLine(point_id1, p_tfrc)
                    l_trb = gmsh.model.geo.addLine(p_tbrc, point_id0)
                l_fbr = lines[j]  # Connecting line between bottom right front and back corner
            elif coords0[1] == 0 and coords1[1] == 0:
                if coords0[0] == x0 or coords1[0] == x0:
                    l_bls = lines[j]  # Connecting line between back bottom left and right corner
                elif coords0[0] == x2 or coords1[0] == x2:
                    l_bsr = lines[j]  # Connecting line between back bottom left and right corner
            elif coords0[1] == y and coords1[1] == y:
                if coords0[0] == x0 or coords1[0] == x0:
                    l_fls = lines[j]  # Connecting line between back bottom left and right corner
                elif coords0[0] == x2 or coords1[0] == x2:
                    l_fsr = lines[j]  # Connecting line between back bottom left and right corner
            elif coords0[0] == x0 and coords1[0] == x0:
                l_fbl = lines[j]  # Connecting line between bottom left front and back corner
                p_fbl0 = point_id0
                p_fbl1 = point_id1
                if coords0[1] == y:
                    l_tlf = gmsh.model.geo.addLine(p_tflc, point_id0
                                                   )  # Connecting line between front right top and bottom point
                    l_tlb = gmsh.model.geo.addLine(point_id1, p_tblc
                                                   )  # Connecting line between back right top and bottom point
                elif coords1[1] == y:
                    l_tlf = gmsh.model.geo.addLine(point_id1, p_tfrc)
                    l_tlb = gmsh.model.geo.addLine(p_tblc, point_id0)

# Create curve and surface loops
ll_right = gmsh.model.geo.addCurveLoop([l_tr, l_trf, l_fbr, l_trb])
s_right = gmsh.model.geo.addPlaneSurface([ll_right])
ll_left = gmsh.model.geo.addCurveLoop([l_tl, l_tlf, l_fbl, l_tlb])
s_left = gmsh.model.geo.addPlaneSurface([ll_left])
ll_back = gmsh.model.geo.addCurveLoop([l_tb, -l_tlb, l_bls, l_bsr, -l_trb])
s_back = gmsh.model.geo.addPlaneSurface([ll_back])
ll_front = gmsh.model.geo.addCurveLoop([l_tf, -l_tlf, l_fls, l_fsr, -l_trf])
s_front = gmsh.model.geo.addPlaneSurface([ll_front])
ss_loop = gmsh.model.geo.addSurfaceLoop([s_front, bl_top_tag[0], bl_top_tag[1], s_right, s_top, s_left, s_back])

# Add physical name for B.C
gmsh.model.addPhysicalGroup(2, bot_surface, 500, "bottom")
gmsh.model.addPhysicalGroup(2, [s_front, s_back, 27, 19, 41, 49], 501, "frontandback")
gmsh.model.addPhysicalGroup(2, [31, s_left], 502, "inlet")
gmsh.model.addPhysicalGroup(2, [s_top], 503, "top")
gmsh.model.addPhysicalGroup(2, [45, s_right], 504, "outlet")

# Add volume
gmsh.model.geo.synchronize()
vol = gmsh.model.geo.addVolume([ss_loop])

mesh = '3D'
if mesh == '2D':
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)
    gmsh.write("beach_profile.msh")
else:
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)
    gmsh.write("beach_3d.msh")

if "-nopopup" not in sys.argv:
    gmsh.fltk.run()
gmsh.finalize()
