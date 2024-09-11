import numpy as np
import matplotlib.pyplot as plt

# Define the grid size
Nx, Ny = 200, 200

# Create a 2D grid
x = np.linspace(0, 1, Nx)
y = np.linspace(0, 1, Ny)
X, Y = np.meshgrid(x, y)

# Define the boundary conditions
u_left = 1.0
u_right = 0.0
u_top = 0.0
u_bottom = 0.0

# Create a 2D array to store the solution
u = np.zeros((Nx, Ny))

# Apply the boundary conditions
u[:, 0] = u_left
u[:, -1] = u_right
u[0, :] = u_top
u[-1, :] = u_bottom

# Define the finite difference coefficients
dx = x[1] - x[0]
dy = y[1] - y[0]
alpha = 1.0 / (dx**2)
beta = 1.0 / (dy**2)

# Create a matrix to store the coefficients
A = np.zeros((Nx*Ny, Nx*Ny))

# Assemble the matrix
for i in range(Nx):
    for j in range(Ny):
        idx = i*Ny + j
        if i == 0 or i == Nx-1 or j == 0 or j == Ny-1:
            A[idx, idx] = 1.0
        else:
            A[idx, idx] = -2*(alpha + beta)
            A[idx, idx-1] = beta
            A[idx, idx+1] = beta
            A[idx, idx-Ny] = alpha
            A[idx, idx+Ny] = alpha

# Create a vector to store the right-hand side
b = np.zeros(Nx*Ny)

# Assemble the right-hand side
for i in range(Nx):
    for j in range(Ny):
        idx = i*Ny + j
        if i == 0:
            b[idx] = u_left
        elif i == Nx-1:
            b[idx] = u_right
        elif j == 0:
            b[idx] = u_top
        elif j == Ny-1:
            b[idx] = u_bottom

# Solve the system of linear equations
u_flat = np.linalg.solve(A, b)

# Reshape the solution to a 2D array
u = u_flat.reshape((Nx, Ny))

# Plot the 2D pseudocolor plot
plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
plt.pcolor(X, Y, u)
plt.colorbar()
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Laplace Equation Solution (2D)')

# Plot the 3D surface plot
from mpl_toolkits.mplot3d import Axes3D
ax = plt.subplot(1, 2, 2, projection='3d')
ax.plot_surface(X, Y, u, cmap='viridis')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('u')
ax.set_title('Laplace Equation Solution (3D)')

plt.show()