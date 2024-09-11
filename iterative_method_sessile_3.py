import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Define the contact radius R and height h0 of the droplet
R = 1.0  # Contact radius (can be adjusted as needed)
h0 = 0.3  # Maximum height of the droplet (can be adjusted as needed)
theta = 2 * np.arctan(h0 / R)  # Contact angle (in radians)
H = 0.4  # Humidity value

# Define the domain
R_max = 4 * R
Z_max = 5 * h0
r_max = 1
z_max = 1
C_max = 2.32e-8  # Saturation concentration of water vapor in air in gm/mm^3

# Define the number of grid points
N = 300  # Number of points in r direction
M = 500  # Number of points in z direction
del_r = 1 / N
del_z = 1 / M
K = ((R_max / Z_max) ** 2) * (del_r / del_z)

# Create the grid
r = np.linspace(0, r_max, N + 1)
z = np.linspace(0, z_max, M + 1)
R_grid, Z_grid = np.meshgrid(r, z)

# Define the droplet height function
def h(r, R, theta):
    r = np.clip(r, 0, R / R_max)  # Limit r to R/R_max
    return (R_max / Z_max) * (np.sqrt((R / np.sin(theta)) ** 2 - r ** 2) - (R / np.tan(theta)))

# Use np.where to define the droplet region
interior = np.where(Z_grid <= h(R_grid, R / R_max, theta), 1, 0)

# Initialize C with humidity value and set interior points to -1
C = np.where(interior == 1, -1, H * np.ones((M + 1, N + 1)))

# Assign boundary conditions for the droplet cap
for i in range(N):
    if r[i] <= R / R_max:  # Ensure h(r) is used only for r <= R
        boundary_indices = np.where(z <= h(r[i], R / R_max, theta))[0]
        if boundary_indices.size > 0:
            boundary_index = boundary_indices[-1]  # Choose the last boundary point for each r
            C[boundary_index, i] = 1

# Set boundary conditions at r = r_max and z = z_max
C[:, -1] = H  # r = r_max
C[-1, :] = H  # z = z_max

# Initiate the iterative method
iter_max = 18000
tolerance = 0.000001
C_old = C.copy()
rms_error_list = []

for iter_num in range(iter_max + 1):
    C_new = C.copy()
    for i in range(1, N):
        for j in range(1, M):
            if interior[j, i] == 0:  # Only calculating for the exterior region
                C_new[j, i] = ((K * C_old[j - 1, i]) + C_old[j, i - 1] + (C_old[j, i + 1] * (1 + del_r / r[i])) + K * C_old[j + 1, i]) / (2 + (2 * K) + (del_r / r[i]))
        if i > int((R / R_max) / del_r):
            C_new[0, i:] = C_new[1, i:]  # No penetration condition at z=0 for r>R
    for j in range(1, M):
        if j > int((h0 / Z_max) / del_z):
            C_new[j:, 0] = C_new[j:, 1]  # Symmetry condition at r=0 and z>h0
    rms_error = np.abs(np.sqrt(np.mean((C_new - C_old)**2)))
    rms_error_list.append(rms_error)
    C_old = C_new.copy()
    if round(rms_error, 6) <= tolerance:
        print('Convergence Reached')
        break
    print(f"Iteration={iter_num} : RMS_error={rms_error:.6f}")

if round(rms_error, 6) > tolerance:
    print('Convergence not reached after maximum iterations')

# Plot of RMS error vs iteration number
plt.figure()
plt.plot(rms_error_list)
plt.xlabel('Iteration Number')
plt.ylabel('RMS Error')
plt.title('RMS Error vs Iteration Number')
plt.savefig('rms_error_plot.jpg')  # Save the plot as JPG

# Mask the interior region for the contour plot
C_masked = np.ma.masked_where(interior == 1, C_new)

# Plot the final concentration contour for the exterior region
plt.figure(figsize=(10, 6))
plt.contourf(R_grid, Z_grid, C_masked, levels=500, cmap='rainbow', alpha=0.6)
plt.colorbar(label='Relative Concentration')
plt.plot(r, h(r, R / R_max, theta), 'r', label='Droplet Surface')
plt.xlabel('r in normalized form')
plt.ylabel('z in normalized form')
plt.title('Sessile Droplet Profile and Interior Region')
plt.ylim(0, 1)
plt.xlim(0, 1)
plt.legend()
plt.savefig('final_concentration_contour.jpg')  # Save the plot as JPG

# Define the range for r from 0 to R/R_max
r_limit = R / R_max
interior_indices = np.where(r <= r_limit)[0]
indices = np.linspace(0, len(interior_indices) - 1, 20, dtype=int)
r_values = r[interior_indices][indices]
h_values = h(r_values, R / R_max, theta)
h_indices = np.searchsorted(z, h_values)
pairs = np.column_stack((interior_indices[indices], h_indices))

# Calculate concentration gradient at boundary points
C_grad_r = []
C_grad_z = []

for r_index, h_index in pairs:
    if r_index + 1 < len(r) and h_index + 1 < len(z):
        grad_r = -(C_max / R_max) * ((C_new[h_index, r_index + 1] - C_new[h_index, r_index]) / del_r)
        grad_z = -(C_max / Z_max) * ((C_new[h_index + 1, r_index] - C_new[h_index, r_index]) / del_z)
        C_grad_r.append(grad_r)
        C_grad_z.append(grad_z)

# Convert to NumPy arrays
C_grad_r = np.array(C_grad_r)
C_grad_z = np.array(C_grad_z)

# Calculate magnitudes and normalize
magnitudes = np.sqrt(C_grad_r**2 + C_grad_z**2)
if magnitudes.size > 0:
    max_magnitude = np.max(magnitudes)
    C_grad_r /= max_magnitude
    C_grad_z /= max_magnitude

# Plotting the vector gradient field
plt.figure(figsize=(10, 6))
plt.plot(r, h(r, R / R_max, theta), 'r', label='Droplet Surface')
plt.quiver(r[interior_indices][indices], z[h_indices], C_grad_r, C_grad_z, color='k', scale=10**-5)
plt.xlabel('r in normalized form')
plt.ylabel('z in normalized form')
plt.title('Sessile Droplet Profile and Gradient Vector Field at the Cap Boundary')
plt.ylim(0, 1)
plt.xlim(0, 1)
plt.legend()
plt.savefig('concentration_gradient_vector_field.jpg')  # Save the plot as JPG

# Prepare data for CSV
gradient_data = []
for i in range(len(pairs)):
    r_index, h_index = pairs[i]
    gradient_data.append([z[h_index], r[r_index], C_grad_r[i], C_grad_z[i], magnitudes[i], np.arctan2(C_grad_z[i], C_grad_r[i])])

gradient_df = pd.DataFrame(gradient_data, columns=['z_position', 'r_position', 'grad_r', 'grad_z', 'magnitude', 'direction'])

# Save to CSV
gradient_df.to_csv('gradient_vectors.csv', index=False)
