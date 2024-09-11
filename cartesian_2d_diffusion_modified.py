import numpy as np
import matplotlib.pyplot as mpt
import pandas as pd
from scipy.sparse import diags
from scipy.sparse.linalg import inv

# Inputs
L = int(input('Length = '))
W = int(input('Width = '))
N = int(input('No. of divisions of L = '))
C1 = float(input('Concentration value for x=0 : '))
C2 = float(input('Concentration value for y=0 : '))
C3 = float(input('Concentration value for x=L : '))
C4 = float(input('Concentration value for y=W : '))

# Derived value
M = round(W * N / L)
print(M)

# Constructing the sparse matrix A
main_diag = -4 * np.ones((N-1)*(M-1))
off_diag1 = np.ones((N-1)*(M-1)-1)
off_diag1[np.arange(1, len(off_diag1)+1) % (N-1) == 0] = 0
off_diag2 = np.ones((N-1)*(M-1)-(N-1))

diagonals = [main_diag, off_diag1, off_diag1, off_diag2, off_diag2]
offsets = [0, 1, -1, (N-1), -(N-1)]
A = diags(diagonals, offsets).tocsc()

# Inverting the sparse matrix A
A_inv = inv(A).toarray()

# For B
B = np.zeros(((N-1)*(M-1), 1))
B[0, 0] = -C1 - C2
for i in range(1, N-2):
    B[i, 0] = -C2 
B[N-2, 0] = -C2 - C3
B[-1, 0] = -C4 - C3
for i in range(-(N-2), -1):
    B[i, 0] = -C4
B[-(N-1), 0] = -C1 - C4
for k in range(1, M-2):
    B[k*(N-1), 0] = -C1
    B[((k+1)*(N-1))-1, 0] = -C3

# Ensure matrix multiplication dimensions match
if A_inv.shape[0] != B.shape[0]:
    raise ValueError(f"Matrix dimension mismatch: A_inv.shape[0] = {A_inv.shape[0]}, B.shape[0] = {B.shape[0]}")

C = np.dot(A_inv, B)  # C values obtained

# Prepare C_actual for plotting
C_actual = np.zeros((M+1, N+1))
C_actual[:, 0] = C1
C_actual[0, :] = C2
C_actual[:, -1] = C3
C_actual[-1, :] = C4
C_actual[1:M, 1:N] = np.round(C.reshape((M-1, N-1)), 3)

# Plotting
x = np.linspace(0, L, N+1)
y = np.linspace(0, W, M+1)
X, Y = np.meshgrid(x, y)

fig, axs = mpt.subplots(1, 2, figsize=(12, 8))

# 2D contour plot
pc = axs[0].pcolor(X, Y, C_actual, cmap='rainbow')
fig.colorbar(pc, ax=axs[0])
axs[0].set_xlabel('X: Length')
axs[0].set_ylabel('Y: Width values')
axs[0].set_title('2D contour plot of concentration values')

# 3D surface plot
axs[1] = fig.add_subplot(122, projection='3d')
surf = axs[1].plot_surface(X, Y, C_actual, cmap='rainbow')
fig.colorbar(surf, ax=axs[1])
axs[1].set_title('3D surface plot of concentration')
axs[1].set_xlabel('X: Length')
axs[1].set_ylabel('Y')
axs[1].set_zlabel('C')

mpt.tight_layout()  # Adjust the layout to avoid overlapping
mpt.show()
