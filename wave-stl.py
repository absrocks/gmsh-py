import numpy as np
from stl import mesh

# Given parameters
x0, x1 = 0, 70
h0, h1 = 1.5, 9.0
n_tau = 2
n_epsilon = 0.75 * (n_tau + 1)
beta = (n_epsilon - 1) / 2
exponent = 1 + beta
alpha = (h1 ** exponent - h0 ** exponent) / (x1 - x0)

# Profile resolution
x_vals = np.linspace(x0, x1, 300)
h_vals = ((x_vals * alpha) + h0 ** exponent) ** (1 / exponent)

# Thickness in y-direction
y_thickness = 1.0
vertices = []
faces = []

def add_face(v0, v1, v2, v3):
    base_idx = len(vertices)
    vertices.extend([v0, v1, v2, v3])
    faces.append([base_idx, base_idx+1, base_idx+2])
    faces.append([base_idx, base_idx+2, base_idx+3])

# Front and back surfaces (sides of the extrusion)
for i in range(len(x_vals) - 1):
    x1_, x2_ = x_vals[i], x_vals[i+1]
    z1_, z2_ = h_vals[i], h_vals[i+1]

    # Front face (y = 0)
    v0 = [x1_, 0, z1_]
    v1 = [x2_, 0, z2_]
    v2 = [x2_, y_thickness, z2_]
    v3 = [x1_, y_thickness, z1_]
    add_face(v0, v1, v2, v3)

# Top face
v0 = [x_vals[-1], 0, h_vals[-1]]
v1 = [x_vals[-1], y_thickness, h_vals[-1]]
v2 = [x_vals[0], y_thickness, h_vals[0]]
v3 = [x_vals[0], 0, h_vals[0]]
add_face(v0, v1, v2, v3)

# Bottom face
v0 = [x_vals[0], 0, 0]
v1 = [x_vals[0], y_thickness, 0]
v2 = [x_vals[-1], y_thickness, 0]
v3 = [x_vals[-1], 0, 0]
add_face(v0, v1, v2, v3)

# Left cap (start x = 0)
v0 = [x_vals[0], 0, 0]
v1 = [x_vals[0], 0, h_vals[0]]
v2 = [x_vals[0], y_thickness, h_vals[0]]
v3 = [x_vals[0], y_thickness, 0]
add_face(v0, v1, v2, v3)

# Right cap (end x = 70)
v0 = [x_vals[-1], 0, 0]
v1 = [x_vals[-1], y_thickness, 0]
v2 = [x_vals[-1], y_thickness, h_vals[-1]]
v3 = [x_vals[-1], 0, h_vals[-1]]
add_face(v0, v1, v2, v3)

vertices = np.array(vertices)
faces = np.array(faces)

# Create and save the closed mesh
profile_mesh = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
for i, f in enumerate(faces):
    for j in range(3):
        profile_mesh.vectors[i][j] = vertices[f[j]]

profile_mesh.save("profile_wall_closed.stl")
