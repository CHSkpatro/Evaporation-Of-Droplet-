import numpy as np
import matplotlib.pyplot as plt
R = 1.0 
h0 = 0.3  
theta = 2*np.arctan(h0 / R)  
H = 0.4  
R_max = 4*R
Z_max = 5*h0
r_max = 1
z_max = 1
C_max = 2.32e-8  ### in gm/mm3 the saturation concentration of water vapor in air
# Define the number of grid points
N = 3000  # Number of points in r direction
M = 3000  # Number of points in z direction
del_r = 1 / N
del_z = 1 / M
K = ((R_max / Z_max) ** 2) * (del_r / del_z)
# Create the grid
r = np.linspace(0, r_max, N + 1)
z = np.linspace(0, z_max, M + 1)
R_grid, Z_grid = np.meshgrid(r, z)

def h(r, R, theta):
    return (R_max / Z_max) * (np.sqrt((R / np.sin(theta)) ** 2 - r ** 2) - (R / np.tan(theta)))

####differentiating interior and exterior of droplet 
interior = np.where(Z_grid <= h(R_grid, R / R_max, theta), 1, 0)
### Inside, let's set the concentration as -1
C = np.zeros((M + 1, N + 1))  #### Initializing C with humidity value
C = np.where(interior == 1, -1, C)  ##### Setting interior points to -1
#### on cap surfacce setting the concentration value as 1 
for i in range(N):
    if r[i] <= R / R_max:  
        boundary_indices = np.where(z <= h(r[i], R / R_max, theta))[0]
        if boundary_indices.size > 0:
            boundary_index = boundary_indices[-1]  # Choose the last boundary point for each r
            C[boundary_index, i] = 1

C[:, -1] = H  #### BC FOR  r = r_max
C[-1, :] = H  #### bc for z = z_max

#### Initiating iterative method ###
iter_max = 1500
tolerance = 0.000001 
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

#### Plot of RMS error vs iteration number####
plt.figure()
plt.plot(rms_error_list)
plt.xlabel('Iteration Number')
plt.ylabel('RMS Error')
plt.title('RMS Error vs Iteration Number')
plt.savefig('rms_error_plot_9million_11.jpg')  ##### Save the plot as JPG
plt.show()

#### Masking the interior region for the contour plot####
C_masked = np.ma.masked_where(interior == 1, C_new)

##### Plot the final concentration contour for the exterior region
plt.figure(figsize=(10, 6))
plt.contourf(R_grid, Z_grid, C_masked, levels=200, cmap='rainbow', alpha=0.6)
plt.colorbar(label='Relative Concentration')
plt.plot(r, h(r, R / R_max, theta), 'r', label='Droplet Surface')
plt.xlabel('r in normalized form')
plt.ylabel('z in normalized form')
plt.title('Sessile Droplet Profile and Interior Region ')
plt.ylim(0, 1)
plt.legend()
plt.savefig('final_concentration_contour_9million_11.jpg')  # Save the plot as JPG
plt.show()
#### Now I want to plot the vector gradient firld exactly at the cap boundary ###
###So 1st lets access the boundary points for 20 r  values in 0 to R/Rmax  and find the gradient vector field at the cap boundary ###
# Define the range for r from 0 to R/R_max
# Define the range for r from 0 to R/R_max
r_limit = R / R_max
interior_indices = np.where(r <= r_limit)[0]
indices = np.linspace(0, len(interior_indices) - 1, 20, dtype=int)
r_values = r[interior_indices][indices]
h_values = h(r_values, R / R_max, theta)
h_indices = np.searchsorted(z, h_values)
pairs = np.column_stack((indices, h_indices))

unique_z_indices, unique_indices = np.unique(h_indices, return_index=True)
unique_pairs = pairs[unique_indices]
r_unique_indices = unique_pairs[:, 0]
h_unique_indices = unique_pairs[:, 1]

C_grad_r = []
C_grad_z = []

for r_index, h_index in unique_pairs:
    if r_index + 1 < len(r) and h_index + 1 < len(z):
        grad_r = -(C_max/R_max) *((C_new[h_index, r_index + 1] - 1) / del_r)
        grad_z = -(C_max/Z_max)*((C_new[h_index + 1, r_index] - 1) / del_z)
        C_grad_r.append(grad_r)
        C_grad_z.append(grad_z)
C_grad_r = np.array(C_grad_r)
C_grad_z = np.array(C_grad_z)
##magnitudes = np.sqrt(C_grad_r**2 + C_grad_z**2)

##C_grad_r_normalized = np.where(magnitudes != 0, C_grad_r / magnitudes, 0)
##C_grad_z_normalized = np.where(magnitudes != 0, C_grad_z / magnitudes, 0)

#### Plotting the vector gradient field
plt.figure(figsize=(10, 6))
plt.plot(r, h(r, R / R_max, theta), 'r', label='Droplet Surface')
plt.quiver(r[unique_pairs[:, 0]], z[unique_pairs[:, 1]],C_grad_r, C_grad_z, color='k', scale=21)
plt.xlabel('r in normalized form')
plt.ylabel('z in normalized form')
plt.title('Sessile Droplet Profile and Gradient Vector Field at the Cap Boundary')
plt.ylim(0, 1)
plt.xlim(0,1)
plt.legend()
plt.savefig('concentration_gradient_vector_field_9million_11.jpg')
plt.show()