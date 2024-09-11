import numpy as np
import matplotlib.pyplot as plt

H = float(input('The relative humidity value: '))
Ro = float(input('Contact line radius for the sessile droplet in mm: '))
theta = np.pi / 4.5
ho = (Ro / np.sin(theta)) - (Ro / np.tan(theta))

#### domain of computation ###
R1 = 5 * Ro
R = Ro / R1
Z1 = 5 * ho
r_max = 1
z_max = 1
N = 1000
del_r = 1 / N
M = 1000
K = ((R1 / Z1) ** 2) * ((M / N) ** 2)
r = np.linspace(0, r_max, N)
z = np.linspace(0, z_max, M)
r_grid, z_grid = np.meshgrid(r, z)  ### mesh grid is formed ###

#### differentiate the interior from exterior region using spherical cap shaped boundary ###
def h(r, R, theta):
    return np.sqrt((R / np.sin(theta)) ** 2 - r **2) - (R / np.tan(theta))

interior = np.where(z_grid <= h(r_grid, R, theta), 1, 0)

##### initializing concentration field ###
C = H * np.ones_like(r_grid)

### set initial values inside the boundary to different random values ###
random_concn_values = -1 - 9 * np.random.random(C.shape)
C = np.where(interior == 1, random_concn_values, C)

##### set cap boundary ###
for j in range(M):
    boundary_index = np.where(r_grid[:, j] <= h(r_grid[:, j], R, theta))[0]
    if boundary_index.size > 0:
        boundary_index = boundary_index[-1]
        C[boundary_index, j] = 1  #### concentration value at boundary of cap ####

##### initiating iterative method ###
tolerance = 2e-6
iter_max = 2000
C_old = C.copy()
rms_error_list = []

#### to plot iteration vs RMS error values in real time ###
fig, ax = plt.subplots()
ax.set_xlabel('Iteration Number')
ax.set_ylabel('RMS error')
### Turn on interactive Mode ###
plt.ion()

for iter_num in range(iter_max):
    C_new = C.copy()
    for i in range(1, N-1):
        for j in range(1, M-1):
            if interior[i, j] == 0:
                C_new[i, j] = (
                    K * C_old[i, j-1] + 
                    C_old[i-1, j] + 
                    (1 + del_r / r[i]) * C_old[i+1, j] + 
                    K * C_old[i, j+1]
                ) / (2 + 2 * K + del_r / r[i])
                
    # Handling edge conditions for the boundaries
    C_new[:, 0] = C_new[:, 1]  # Handle left boundary (r=0)
    C_new[0, :] = C_new[1, :]  # Handle bottom boundary (z=0)
    
    rms_error = np.abs(np.sqrt(np.mean(C_new ** 2)) - np.sqrt(np.mean(C_old ** 2)))
    rms_error_list.append(rms_error)
    C_old = C_new.copy()
    
    if round(rms_error, 6) <= tolerance:
        break
    
    print(f"Iteration = {iter_num}: RMS_error = {rms_error:.6f}")
    ax.plot(range(len(rms_error_list)), rms_error_list)
    plt.pause(0.01)

plt.ioff()
plt.show()

### plotting the concentration contours ###
fig, ax = plt.subplots(figsize=(10, 6))
contourf_plot = ax.contourf(r_grid, z_grid, C_new, levels=np.linspace(C_new.min(), C_new.max(), 200), cmap='viridis')
ax.set_title('Concentration contour for a sessile droplet')
ax.set_xlabel('Radial distance in normalized form: r/R_max')
ax.set_ylabel('Height in normalized form: z/Z_max')
fig.colorbar(contourf_plot, ax=ax, label='Relative Concentration Values: C/Cv')
ax.set_aspect('equal', 'box')
plt.show()
