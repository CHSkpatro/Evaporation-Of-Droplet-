import numpy as np
import matplotlib.pyplot as plt

# Input values
H = float(input('Enter the relative humidity value: '))
Ro = float(input('Enter the contact line radius for the sessile droplet in mm: '))

# Parameters
theta = np.pi / 4.5
ho = (Ro / np.sin(theta)) - (Ro / np.tan(theta))

# Domain of computation
R1 = 5 * Ro
R = Ro / R1
Z1 = 5 * ho
r_max = 1
z_max = 1
N = 1000
M = 1000
del_r = 1 / N
K = ((R1 / Z1) ** 2) * ((M / N) ** 2)

# Grid setup
r = np.linspace(0, r_max, N)
z = np.linspace(0, z_max, M)
r_grid, z_grid = np.meshgrid(r, z)  # Mesh grid formation

# Function for the spherical cap boundary with error handling
def h(r, R, theta):
    argument = (R / np.sin(theta)) ** 2 - r ** 2
    argument = np.maximum(argument, 0)  # Ensure argument is non-negative
    return np.sqrt(argument) - (R / np.tan(theta))

# Differentiate the interior from the exterior region
interior = np.where(z_grid <= h(r_grid, R, theta), 1, 0)

# Initializing concentration field
C = H * np.ones_like(r_grid)

# Set the initial values inside the boundary to different random values
random_concn_values = -1 - 9 * np.random.random(C.shape)
C = np.where(interior == 1, random_concn_values, C)

# Set the cap boundary concentration value
for i in range(N):
    boundary_index = np.where(z_grid[i] <= h(r_grid[i], R, theta))[0]
    if boundary_index.size > 0:
        boundary_index = boundary_index[-1]
        C[i, boundary_index] = 1  # Concentration value at boundary of the cap

# Iterative method initialization
tolerance = 2e-6
iter_max = 2000
C_old = C.copy()
rms_error_list = []

# Plot setup
fig, ax = plt.subplots()
ax.set_xlabel('Iteration Number')
ax.set_ylabel('RMS error')

# Turn on interactive mode
plt.ion()

# Iterative solver
for iter_num in range(iter_max):
    C_new = C.copy()
    for i in range(1, N-1):  # Radial direction
        for j in range(1, M-1):  # Vertical direction
            if interior[i, j] == 0:  # Compute at external grid points of the sessile droplet
                C_new[i, j] = (
                    K * C_old[i, j-1] +        # Below neighbor in the vertical direction
                    C_old[i-1, j] +            # Left neighbor in the radial direction
                    (1 + del_r / r[i]) * C_old[i+1, j] +  # Right neighbor in the radial direction
                    K * C_old[i, j+1]          # Above neighbor in the vertical direction
                ) / (2 + 2 * K + del_r / r[i])

                # Ensure concentration values remain non-negative
                C_new[i, j] = np.maximum(C_new[i, j], 0)

    # Handling edge conditions for the boundaries
    C_new[0, :] = C_new[1, :]  # Bottom boundary (z = 0)
    C_new[:, 0] = C_new[:, 1]  # Left boundary (r = 0)

    rms_error = np.abs(np.sqrt(np.mean(C_new ** 2)) - np.sqrt(np.mean(C_old ** 2)))
    rms_error_list.append(rms_error)
    C_old = C_new.copy()

    if round(rms_error, 6) <= tolerance:
        break

    print(f"Iteration={iter_num} : RMS_error={rms_error:.6f}")

    # Update the iteration error plot
    ax.plot(range(len(rms_error_list)), rms_error_list, color='blue')
    plt.pause(0.01)  # Pause for a short time to update the plot

plt.ioff()  # Turn off interactive mode
plt.show()

# Plotting the concentration contours
fig, ax = plt.subplots(figsize=(10, 6))
contourf_plot = ax.contourf(r_grid, z_grid, C_new, levels=np.linspace(C_new.min(), C_new.max(), 200), cmap='viridis')
ax.set_title('Concentration contour for a sessile droplet')
ax.set_xlabel('Radial distance in normalized form: r/R_max')
ax.set_ylabel('Height in normalized form: z/Z_max')
fig.colorbar(contourf_plot, ax=ax, label='Relative Concentration Values: C/Cv')
ax.set_aspect('equal', 'box')
plt.show()

