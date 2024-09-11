import numpy as np
import pandas as pd
import matplotlib.pyplot as mpt
H = float(input('The value of relative humidity is = '))
tolerance = float(input('The convergence tolerance(e.g.:1e-6) = '))
###Initialization of variables###
previous_C = None
N = 2 
while True:
    A = np.zeros((N - 1, N - 1))
    np.fill_diagonal(A, -2)  
    np.fill_diagonal(A[:, 1:], 1)  ### Fill the super diagonal with 1###
    np.fill_diagonal(A[1:, :], 1)  ### Fill the sub diagonal  with 1###
    A_inv = np.linalg.inv(A)
### A*c = B ###
    B = np.zeros((N - 1, 1))
    B[0, 0] = -1  
    B[-1, 0] = -H  ###We are calculating relative concentration values###
    C = np.dot(A_inv, B)
    if previous_C is not None and previous_C.shape == C.shape:
        difference = C-previous_C
        diff = np.linalg.norm(difference)###Difference in magnitude of the vector column###
        if diff < tolerance:
            print("Result converged for N =",N)
            break

    previous_C = np.copy(C)
    N += 50

df_C = pd.DataFrame(C, columns=['Relative Concentration'])
output_file = 'C_values.xlsx'
df_C.to_excel(output_file, index=False)
print(f"Saved C values to {output_file}")
x=np.linspace(0,1,N+1)###The domain is converted into non dimensonal form ###
###Matrix C does not contain the boundary values###
C0=1
C1=H
C_entire=np.concatenate(([C0],C,[C1]))
mpt.subplot(1,2,1)
mpt.scatter(x,C_entire.flatten(), marker='*', color='r')
mpt.show()
mpt.subplot(1,2,2)
mpt.plot(x,B.flatten(),'-',color='black',label='Computational analysis Solution')
mpt.legend()
mpt.show()
mpt.subplot(1,2,2)
mpt.plot(x,1+x*(H-1),marker='--',color='y',label='Analytical Solution')
mpt.legend()
mpt.show()