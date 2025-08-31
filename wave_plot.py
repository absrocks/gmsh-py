import numpy as np
import matplotlib.pyplot as plt
import os
from collections import defaultdict
from scipy.spatial import KDTree
import matplotlib.tri as tri
import copy
import pyvista as pv

#full_path = '/home/am2455/Documents/research/dam_break_kim'
full_path = '/run/media/am2455/Backup/mac_work/dam_break/Kim'
cases = ['Case_I', 'Case_II', 'Case_III', 'Case_IV']  #
case_path = []


def read_scalar_field(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # Skip header
    start_idx = next(i for i, line in enumerate(lines) if line.strip().isdigit())
    n_values = int(lines[start_idx].strip())
    values = []
    for line in lines[start_idx+2 : start_idx+2+n_values]:
        values.append(float(line.strip()))
    return np.array(values)

def span_average(Cx, Cy, U):

    # Round coordinates to collapse floating-point noise

    U_avg = np.mean(U, axis=0)
    x_rounded = np.round(Cx, 6)
    y_rounded = np.round(Cy, 6)
    xy_keys = np.stack((x_rounded, y_rounded), axis=1)

    # Step 1: Group indices by (x, y)
    xy_to_indices = defaultdict(list)
    for i, key in enumerate(map(tuple, xy_keys)):
        xy_to_indices[key].append(i)

    # Step 2: Initialize output array
    U_spanwise_avg = np.zeros_like(U_avg)  # shape (N, 3)

    # Step 3: Fill each (x, y) group with the averaged vector
    for indices in xy_to_indices.values():
        avg_vec = np.mean(U_avg[indices], axis=0)
        U_spanwise_avg[indices] = avg_vec  # assign to all matching cells

    return U_spanwise_avg

def alpha_mask(alpha_range, alpha_water, Cx, Cy, Cz):

    alpha_mask = (alpha_water >= alpha_range).flatten()

    # Apply mask to U and coordinates

    Cx_filtered = Cx[alpha_mask]
    Cy_filtered = Cy[alpha_mask]
    Cz_filtered = Cz[alpha_mask]
    return alpha_mask, Cx_filtered, Cy_filtered, Cz_filtered

def read_vector_field(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()

    start_idx = next(i for i, line in enumerate(lines) if line.strip().isdigit())
    n_values = int(lines[start_idx].strip())
    vectors = []
    for line in lines[start_idx+2 : start_idx+2+n_values]:
        vec = line.strip().strip('()').split()
        vectors.append([float(v) for v in vec])
    return np.array(vectors)

def filter_mask(Cz, target_z, tolerance, Cx, x_min, x_max):
    #slice_mask = np.abs(Cz - target_z) < tolerance  # or Cy for y-slice

    x_range_mask = (Cx >= x_min) & (Cx <= x_max)

    # Combine both masks
    #mask = slice_mask & x_range_mask
    mask = x_range_mask
    x = Cx[mask]
    y = Cy[mask]

    return x, y, mask

def write_vtk(Cx, Cy, Cz, TKE, U_spanwise_avg):
    # Block definitions (order and sizes must match how Cx, Cy, Cz are laid out!)
    blocks = [
        {"dims": (960, 70, 90)},  # Block 1
        {"dims": (960, 92, 90)},  # Block 2
        {"dims": (1140, 92, 90)}  # Block 3
    ]
    offset_pts = 0
    offset_cells = 0
    block_grids = []

    for block in blocks:
        nx, ny, nz = block["dims"]
        npts = (nx + 1) * (ny + 1) * (nz + 1)
        ncells = nx * ny * nz

        # Extract point coordinates
        cx_block = Cx[offset_pts:offset_pts + npts].reshape((nx + 1, ny + 1, nz + 1))
        cy_block = Cy[offset_pts:offset_pts + npts].reshape((nx + 1, ny + 1, nz + 1))
        cz_block = Cz[offset_pts:offset_pts + npts].reshape((nx + 1, ny + 1, nz + 1))

        # Create StructuredGrid
        grid = pv.StructuredGrid(cx_block, cy_block, cz_block)

        # Assign TKE as cell data
        tke_block = TKE[offset_cells:offset_cells + ncells].reshape((nx, ny, nz))
        grid.cell_data["TKE"] = tke_block.ravel()

        # Track offsets and collect grid
        block_grids.append(grid)
        offset_pts += npts
        offset_cells += ncells

    # Combine all blocks into a single dataset
    combined = pv.MultiBlock(block_grids).combine()
    combined.save("combined_blocks.vtk")

    return
def build_grid_from_flat(Cx, Cy, Cz, arrs, decimals=12):
    """
    Map flat arrays (len N) to rectilinear 3D grids using unique coordinates.
    - Scalars (N,) -> (Nx,Ny,Nz)
    - Vectors (N,3) -> (Nx,Ny,Nz,3)
    """
    Cx = np.asarray(Cx); Cy = np.asarray(Cy); Cz = np.asarray(Cz)
    # stabilize float noise before unique()
    Cx_r = np.round(Cx, decimals=decimals)
    Cy_r = np.round(Cy, decimals=decimals)
    Cz_r = np.round(Cz, decimals=decimals)

    X = np.unique(Cx_r); Y = np.unique(Cy_r); Z = np.unique(Cz_r)
    Nx, Ny, Nz = len(X), len(Y), len(Z)

    # index of each sample along each axis
    ix = np.searchsorted(X, Cx_r)
    iy = np.searchsorted(Y, Cy_r)
    iz = np.searchsorted(Z, Cz_r)

    grids = []
    for A in arrs:
        A = np.asarray(A)
        if A.ndim == 1:
            G = np.empty((Nx, Ny, Nz), dtype=A.dtype)
            G[ix, iy, iz] = A
        elif A.ndim == 2 and A.shape[1] == 3:
            G = np.empty((Nx, Ny, Nz, 3), dtype=A.dtype)
            G[ix, iy, iz, :] = A
        else:
            raise ValueError(f"Array must be (N,) or (N,3); got {A.shape}")
        grids.append(G)
    return X, Y, Z, grids

def epsilon(UPrime, Cx, Cy, Cz, nut, nu, cases):
    UPrime = np.array(UPrime)
    nut = np.array(nut)
    s_contract = []
    s_contract_nut = []
    print("uprime.shape", UPrime.shape)
    print("nut.shape", nut.shape)
    for icase in range(len(cases)):
        X, Y, Z, (U3, nut3) = build_grid_from_flat(Cx, Cy, Cz, [UPrime[icase, :, :], nut[icase, :, :]])
        print("U3.shape", U3.shape)
        print("nut3.shape", nut3.shape)
        Ux, Uy, Uz = U3[..., 0], U3[..., 1], U3[..., 2]


        # Velocity gradients (structured grid)
        dUx_dx, dUx_dy, dUx_dz = np.gradient(Ux, Cx, Cy, Cz, edge_order=2)
        dUy_dx, dUy_dy, dUy_dz = np.gradient(Uy, Cx, Cy, Cz, edge_order=2)
        dUz_dx, dUz_dy, dUz_dz = np.gradient(Uz, Cx, Cy, Cz, edge_order=2)

        # Symmetric strain-rate components s'_{ij}
        sxx = dUx_dx
        syy = dUy_dy
        szz = dUz_dz
        sxy = 0.5 * (dUx_dy + dUy_dx)
        sxz = 0.5 * (dUx_dz + dUz_dx)
        syz = 0.5 * (dUy_dz + dUz_dy)

        # Contraction s'_{ij}s'_{ij} = diag^2 + 2*offdiag^2
        s_contract.append((sxx ** 2 + syy ** 2 + szz ** 2) + 2.0 * (sxy ** 2 + sxz ** 2 + syz ** 2))
        s_contract_nut.append(((sxx ** 2 + syy ** 2 + szz ** 2) + 2.0 * (sxy ** 2 + sxz ** 2 + syz ** 2))* nut[icase])
    # Span average
    s_contract_avg = span_average(s_contract)
    eps_bar = 2.0 * nu * s_contract_avg
    # Component-wise multiplication by ν_tk BEFORE averaging, then scale by 2
    # Broadcast s_contract to (..., 1) to multiply each component
    eps_vec = 2.0 * span_average(s_contract_nut)
    eps = eps_bar + eps_vec

    return eps

# === CONFIGURATION ===
case_dir = '3.5'
target_z = 0.225     # z-plane value to extract
tolerance = 0.0025   # tolerance to match z
nu = 1e-6
U = []
nut = []

# Add X range filter
x_min = 6.9
x_max = 9.9
y_min= 0.2
y_max = 0.5
alpha = []
U_alpha = []

print("Reading the data from OpenFOAM\n")
for icase in range(len(cases)):
    case_path.append(full_path + '/' + cases[icase] + '/' + case_dir)
    # === READ FIELDS ===
    U.append(read_vector_field(os.path.join(case_path[icase], 'U')))
    alpha.append(read_vector_field(os.path.join(case_path[icase], 'alpha1')))
    U_alpha.append(U[icase]*alpha[icase])
    nut.append(read_vector_field(os.path.join(case_path[icase], 'nuSgs')))
#U_alpha = U * alpha
print("U_alpha.shape:", U_alpha[0].shape)
Cx = read_scalar_field(os.path.join(full_path, 'Cx'))
Cy = read_scalar_field(os.path.join(full_path, 'Cy'))
Cz = read_scalar_field(os.path.join(full_path, 'Cz'))
#C = read_vector_field(os.path.join(case_path[0], 'C'))
change_indices = np.where(np.diff(Cz) != 0)[0] + 1
#print("First change occurs at index:", change_indices[0])
#print("Coordinate at Cz0:", C[0])
#print("U.shape:", U[0].shape)

#alpha = read_vector_field(os.path.join(case_path[0], 'alpha1'))
print("Data Reading Completed\n")

U_spanwise_avg = span_average(Cx, Cy, U_alpha)
alpha_spanwise_avg = span_average(Cx, Cy, alpha)
#print("U_spanwise_avg.shape:", U_spanwise_avg.shape)

UPrime2square = []
UPrime = []
#print("UPrime2square[0].shape:", UPrime2square[0].shape)

for icase in range(len(cases)):
    UPrime2square.append((U_alpha[icase] - U_spanwise_avg)**2)
    UPrime.append(U[icase] - U_spanwise_avg)
#print("UPrime.shape:", UPrime.shape)
#print("UPrime[0].shape:", UPrime[0].shape)
#print("UPrime2square[0].shape:", UPrime2square[0].shape)

#raise ValueError("This is an intentional ValueError.")

UPrime2square_spanwise_avg = span_average(Cx, Cy, UPrime2square)
print("UPrime2square_spanwise_avg.shape:", UPrime2square_spanwise_avg.shape)
#SUPrimeMag2 = SU(UPrime, Cx, Cy, Cz, cases)
#SUPrimeMag2_nut = copy.deepcopy(SUPrimeMag2)
#for icase in range(len(cases)):
#    SUPrimeMag2_nut[icase] = SUPrimeMag2_all[icase] * nut[icase]

eps = epsilon(UPrime, Cx, Cy, Cz, nut, nu, cases)

TKE = 0.5 * np.sum(UPrime2square_spanwise_avg, axis=1)
print("TKE.shape:", TKE.shape)
# Write VTK
#VTK = write_vtk(Cx, Cy, Cz, TKE, U_spanwise_avg)

# Filter based on alpha.water
alpha_mask, Cx, Cy, Cz= alpha_mask(0.5,alpha_spanwise_avg, Cx, Cy, Cz)
TKE = TKE[alpha_mask]
eps = eps[alpha_mask]

# === FILTER CELLS AT SPECIFIED Y ===
x, y, mask = filter_mask(Cz, target_z, tolerance, Cx, x_min, x_max)
TKE = TKE[mask]
eps =eps[mask]

print("x.shape:", x.shape)
print("TKE.shape:", TKE.shape)
print("eps.shape:", eps.shape)

# === PLOT ===
# Create the figure and axes FIRST
# Triangulation once
triang = tri.Triangulation(x, y)

# Create 1 row, 2 columns of subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 7), constrained_layout=True)

# --- Left plot: TKE ---
vmin, vmax = 0, 0.05
levels = np.linspace(vmin, vmax, 100)

sc1 = ax1.tricontourf(triang, TKE, levels=levels, cmap='jet', vmin=vmin, vmax=vmax)
ax1.set_xlim(x_min, x_max)
ax1.set_ylim(y_min, y_max)
ax1.set_xlabel('x [m]', fontsize=14, fontweight='bold')
ax1.set_ylabel('z [m]', fontsize=14, fontweight='bold')
ax1.set_title('Turbulent Kinetic Energy', fontsize=15, fontweight='bold')

cbar1 = fig.colorbar(sc1, ax=ax1, orientation='vertical')
cbar1.set_label(r'$K\ [\mathrm{m}^{2}/\mathrm{s}^{2}]$', fontsize=13)

# --- Right plot: ε ---
vmin_eps, vmax_eps = 0, 0.2
levels_eps = np.linspace(vmin_eps, vmax_eps, 100)

sc2 = ax2.tricontourf(triang, eps, levels=levels_eps, cmap='jet', vmin=vmin_eps, vmax=vmax_eps)
ax2.set_xlim(x_min, x_max)
ax2.set_ylim(y_min, y_max)
ax2.set_xlabel('x [m]', fontsize=14, fontweight='bold')
ax2.set_ylabel('z [m]', fontsize=14, fontweight='bold')
ax2.set_title(r'Dissipation Rate $\varepsilon$', fontsize=15, fontweight='bold')

cbar2 = fig.colorbar(sc2, ax=ax2, orientation='vertical')
cbar2.set_label(r'$\varepsilon\ [\mathrm{m}^{2}/\mathrm{s}^{3}]$', fontsize=13)

# Beautify axes (thicker borders)
for ax in (ax1, ax2):
    ax.tick_params(axis='both', labelsize=12)
    for spine in ax.spines.values():
        spine.set_linewidth(2.0)

plt.show()
