import numpy as np
import matplotlib.pyplot as mpt 
N=int(input('The number of divisions = '))
H=float(input('The value of relative humidity is = '))
A=np.zeros((N-1,N-1))
np.fill_diagonal(A,-2)
np.fill_diagonal(A[:,1:],1)##for filling the super diagonal by 1##
np.fill_diagonal(A[1:,:],1)## for filling the sub diagonals by 1##
print(A.shape)
A_inv=np.linalg.inv(A)
### A*C = B ### where C is the concentration values to be found 
B=np.zeros((N-1,1))
B[0,0]=-1
B[N-2,0]=-H
C=np.dot(A_inv,B)
C0=1
C1=H
C_actual=np.concatenate(([C0], C.flatten(), [C1]))
x=np.linspace(0,1,N+1)
mpt.subplot(1,2,1)
mpt.scatter(x,C_actual,marker='*', color='r')
mpt.title(" Relative Concentration Profile")
mpt.xlabel("Position in non dimensional form")
mpt.ylabel("Relative Concentration")
mpt.xlim(0,1)
mpt.grid(True)
mpt.subplot(1, 2, 2)
mpt.plot(x, C_actual, '-', color='black', label='Computational Solution')
mpt.plot(x, 1 + x * (H - 1), '--', color='y', label='Analytical Solution')
mpt.title("Comparison with Analytical Solution")
mpt.xlabel("Position in non dimensional form")
mpt.ylabel("Relative Concentration")
mpt.grid(True)
mpt.legend()
mpt.show() 



