import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Define the contact radius R and height h0 of the droplet
R = 1.0  # Contact radius (can be adjusted as needed)
h0 = 0.3  # Maximum height of the droplet (can be adjusted as needed)
theta = 2 * np.arctan(h0 / R)  # Contact angle (in radians)
H = 0.4  # Humidity value

# Define the domain
R_max = 5 * R
Z_max = 8 * h0
### Converting into the normalized form
r_max = 1
z_max = 1
C_max = 2.32e-8  ### in gm/mm3 the saturation concentration of water vapor in air

# Define the number of grid points
N = 500  # Number of points in r direction
M = 800  # Number of points in z direction
del_r = 1 / N
del_z = 1 / M
K = ((R_max / Z_max) ** 2) * (del_r / del_z)

# Create the grid
r = np.linspace(0, r_max, N + 1)
z = np.linspace(0, z_max, M + 1)
R_grid, Z_grid = np.meshgrid(r, z)

# Define the droplet height function
def h(r, R, theta):
    return (R_max / Z_max) * (np.sqrt((R / np.sin(theta)) ** 2 - r ** 2) - (R / np.tan(theta)))

# Use np.where to define the droplet region
interior = np.where(Z_grid <= h(R_grid, R / R_max, theta), 1, 0)
### Inside, let's set the concentration as -1
C = H * np.ones((M + 1, N + 1))  # Initialize C with humidity value
C = np.where(interior == 1, -0.1, C)  # Set interior points to -1

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

#### Initiating iterative method ###
iter_max = 20000
tolerance = 0.000001  ### 0.000001 ###
C_old = C.copy()
rms_error_list = []

for iter_num in range(iter_max + 1):
    C_new = C.copy()  ### so that each time the loop runs the C_new will take the form of C and then as per the following lines the elements will be replaced 
    for i in range(1, N):
        for j in range(1, M):
            if interior[j, i] == 0:  ### only calculating for the exterior region only
                C_new[j, i] = ((K * C_old[j - 1, i]) + C_old[j, i - 1] + (C_old[j, i + 1] * (1 + del_r / r[i])) + K * C_old[j + 1, i]) / (2 + (2 * K) + (del_r / r[i]))
        if i > int((R / R_max) / del_r):
            C_new[0, i:] = C_new[1, i:]  ### No penetration condition so no flux at z=0 for r>R so partial derivative of C wrt z is zero
    for j in range(1, M):
        if j > int((h0 / Z_max) / del_z):
            C_new[j:, 0] = C_new[j:, 1]  ### For symmetry condition at r=0 and z>ho, no flux partial derivative of C wrt r is zero
    rms_error = np.abs(np.sqrt(np.mean((C_new - C_old)**2)))
    rms_error_list.append(rms_error)
    C_old = C_new.copy()
    if round(rms_error, 6) <= tolerance:
        print('Convergence Reached')
        break
    print(f"Iteration={iter_num} : RMS_error={rms_error:.6f}")

if round(rms_error, 6) > tolerance:
    print('Convergence not reached after maximum iterations')

# Plot the RMS error vs iteration number
plt.figure()
plt.plot(rms_error_list)
plt.xlabel('Iteration Number')
plt.ylabel('RMS Error')
plt.title('RMS Error vs Iteration Number')
plt.savefig('RMS_Error_vs_Iteration_Number_14:00.jpg')
plt.show()

# Mask the interior region for the contour plot
C_masked = np.ma.masked_where(interior == 1, C_new)
# Plot the final concentration contour for the exterior region
plt.figure(figsize=(10, 6))
contour = plt.contourf(R_grid, Z_grid, C_masked, levels=np.linspace(C_masked.min(), C_masked.max(),200), cmap='coolwarm',alpha=0.5)
plt.colorbar(contour, label='Relative Concentration')
plt.plot(r, h(r, R / R_max, theta), 'r', label='Droplet Surface')
plt.xlabel('r in normalized form')
plt.ylabel('z in normalized form')
plt.title('Sessile droplet profile and Concentration Contour')
plt.ylim(0, 1)
plt.legend()
plt.savefig('Sessile_droplet_profile_and_Concentration_Contour_14:00.jpg')
plt.show()

# Calculate concentration gradient at boundary points
C_grad_r = np.zeros_like(C_new)
C_grad_z = np.zeros_like(C_new)

# Forward differences for left, bottom, and interior points
C_grad_r[:, :-1] = -(C_max / R_max) * (C_new[:, 1:] - C_new[:, :-1]) / del_r
C_grad_z[:-1, :] = -(C_max / Z_max) * (C_new[1:, :] - C_new[:-1, :]) / del_z

# Backward differences for right and top wall
C_grad_r[:, -1] = -(C_max / R_max) * (C_new[:, -1] - C_new[:, -2]) / del_r
C_grad_z[-1, :] = -(C_max / Z_max) * (C_new[-1, :] - C_new[-2, :]) / del_z

# Find boundary positions within the droplet cap surface
boundary_positions = []
for i in range(N):
    if r[i] <= R / R_max:
        boundary_indices = np.where(z <= h(r[i], R / R_max, theta))[0]
        if boundary_indices.size > 0:
            boundary_index = boundary_indices[-1]
            boundary_positions.append((boundary_index, i))

boundary_positions = np.array(boundary_positions)

# Extract gradient components at boundary positions
C_grad_R_boundary = C_grad_r[boundary_positions[:, 0], boundary_positions[:, 1]]
C_grad_Z_boundary = C_grad_z[boundary_positions[:, 0], boundary_positions[:, 1]]
magnitudes = np.sqrt(C_grad_R_boundary**2 + C_grad_Z_boundary**2)
direction_boundary = np.arctan2(C_grad_Z_boundary, C_grad_R_boundary)  # Angle in radians
max_magnitude=np.max(magnitudes)
C_grad_r=C_grad_R_boundary/max_magnitude
C_grad_z=C_grad_Z_boundary/max_magnitude
magnitude_boundary=np.sqrt(C_grad_r**2 + C_grad_z**2)
#...

step = 4  # adjust this value to control the density of the plot

# Plotting the gradient vector field at boundary positions
fig, ax = plt.subplots(figsize=(10, 8))
Q = ax.quiver(r[boundary_positions[:, 1]][::step], z[boundary_positions[:, 0]][::step], 
              C_grad_r[::step], C_grad_z[::step], 
              magnitude_boundary[::step], scale=1, cmap='rainbow')
fig.colorbar(Q, ax=ax)
ax.set_xlabel('r in normalized form')
ax.set_ylabel('z in normalized form')
ax.set_title('Concentration Gradient Vector Field at Droplet Surface')
ax.set_xlim(0, r_max)
ax.set_ylim(0, z_max)
plt.savefig('Concentration_gradient_vector_field_at_boundary_scale_1_22:00.jpg')
plt.show()

step = 6  # adjust this value to control the density of the plot

# Plotting the gradient vector field at boundary positions
fig, ax = plt.subplots(figsize=(10, 8))
Q = ax.quiver(r[boundary_positions[:, 1]][::step], z[boundary_positions[:, 0]][::step], 
              C_grad_r[::step], C_grad_z[::step], 
              magnitude_boundary[::step], scale=0.1, cmap='rainbow')
fig.colorbar(Q, ax=ax)
ax.set_xlabel('r in normalized form')
ax.set_ylabel('z in normalized form')
ax.set_title('Concentration Gradient Vector Field at Droplet Surface')
ax.set_xlim(0, r_max)
ax.set_ylim(0, z_max)
plt.savefig('Concentration_gradient_vector_field_at_boundary_scale_0.1_22:00.jpg')
plt.show()

step = 8  # adjust this value to control the density of the plot

# Plotting the gradient vector field at boundary positions
fig, ax = plt.subplots(figsize=(10, 8))
Q = ax.quiver(r[boundary_positions[:, 1]][::step], z[boundary_positions[:, 0]][::step], 
              C_grad_r[::step], C_grad_z[::step], 
              magnitude_boundary[::step], scale=0.01, cmap='rainbow')
fig.colorbar(Q, ax=ax)
ax.set_xlabel('r in normalized form')
ax.set_ylabel('z in normalized form')
ax.set_title('Concentration Gradient Vector Field at Droplet Surface')
ax.set_xlim(0, r_max)
ax.set_ylim(0, z_max)
plt.savefig('Concentration_gradient_vector_field_at_boundary_scale_0.01_22:00.jpg')
plt.show()

step = 8  # adjust this value to control the density of the plot

# Plotting the gradient vector field at boundary positions
fig, ax = plt.subplots(figsize=(10, 8))
Q = ax.quiver(r[boundary_positions[:, 1]][::step], z[boundary_positions[:, 0]][::step], 
              C_grad_r[::step], C_grad_z[::step], 
              magnitude_boundary[::step], scale=2, cmap='rainbow')
fig.colorbar(Q, ax=ax)
ax.set_xlabel('r in normalized form')
ax.set_ylabel('z in normalized form')
ax.set_title('Concentration Gradient Vector Field at Droplet Surface')
ax.set_xlim(0, r_max)
ax.set_ylim(0, z_max)
plt.savefig('Concentration_gradient_vector_field_at_boundary_scale_2_22:00.jpg')
plt.show()
# Prepare data for CSV
gradient_data = []
for i in range(boundary_positions.shape[0]):
    z_index, r_index = boundary_positions[i]
    gradient_data.append([z[z_index], r[r_index], C_grad_R_boundary[i], C_grad_Z_boundary[i], magnitude_boundary[i], direction_boundary[i]])

gradient_df = pd.DataFrame(gradient_data, columns=['z_position', 'r_position', 'grad_r', 'grad_z', 'magnitude', 'direction'])

# Save to CSV
gradient_df.to_csv('gradient_vectors_22:00.csv', index=False)