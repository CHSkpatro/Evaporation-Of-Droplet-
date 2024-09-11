import numpy as np
import matplotlib.pyplot as mpt
import pandas as pd
L=np.int16(input('Length ='))
W=np.int16(input('Width ='))
N=int(input('No. of divisions of L = '))
C1=float(input('Concentration value for x=0 : '))
C2=float(input('Concentration value for y=0 : '))
C3=float(input('Concentration value for x=L : '))
C4=float(input('Concentration value for y=W : '))
M=round(W*N/L)
print(M)
###print(L/N)
###print(W/M)
###A*C=B###
A=np.zeros(((N-1)*(M-1),(N-1)*(M-1)))
###filling up diagonal elements###
for i in range ((N-1)*(M-1)):
    A[i,i]=-4 
    if i< (((N-1)*(M-1))-1):  
        if i==0 or (i+1)%(N-1)!=0:
          A[i,i+1]=1###super diagonal elements###but there are still slective rows in which the super diagonal element is 0
        if i==0 or (i+1)%(N-1)!=0:
          A[i+1,i]=1###sub diagonal elements ####but also there are selective rows in which the sub diagonal element is 0 
for i in range((N-1),((N-1)*(M-1))):
    A[i,i-(N-1)]=1  ###filling the second sub diagonal elements###
    A[i-(N-1),i] =1 ###filling up the 2nd super diagonal elements###
A_inv=np.linalg.inv(A)
##df = pd.DataFrame(A)
#df.to_csv('matrix_A.csv', index=False)
##print("Matrix A has been saved to 'matrix_A.csv'.")
####For B #####
B=np.zeros(((N-1)*(M-1),1))
B[0,0]=-C1-C2
for i in range(1,N-2):
    B[i,0]=-C2 
B[N-2,0]=-C2-C3
B[-1,0]=-C4-C3
for i in range(-(N-2),-1):
      B[i,0]=-C4
B[-(N-1),0] =-C1-C4
for k in range(1,M-2):
    B[k*(N-1),0]=-C1
    B[((k+1)*(N-1))-1,0]=-C3
###df = pd.DataFrame(B)
###df.to_csv('matrix_B.csv', index=False)
###print("Matrix B has been saved to 'matrix_B.csv'.")   
C=np.dot(A_inv,B) ###C values obtained ###(C11,C21,C31,C41,.............) THE VALUES ARE IN THIS ORDER
C_actual=np.zeros((M+1,N+1))
C_actual[:,0]=C1
C_actual[0,:]=C2
C_actual[:,-1]=C3
C_actual[-1,:]=C4
C_actual[1:M,1:N]=np.round(C.reshape((M-1,N-1)),3)
x=np.linspace(0,L,N+1)
y=np.linspace(0,W,M+1)
X,Y=np.meshgrid(x,y)
###2d contour plot ###
fig, axs = mpt.subplots(1, 2, figsize=(12, 8))
pc = axs[0].pcolor(X, Y, C_actual, cmap='rainbow')
fig.colorbar(pc, ax=axs[0]) 
axs[0].set_xlabel('X :Length ')
axs[0].set_ylabel('Y:Width values ')
axs[0].set_title('2D contour plot of concentration values ')

### 3d surface contour plot###
axs[1] = fig.add_subplot(122, projection='3d')  
surf = axs[1].plot_surface(X, Y, C_actual, cmap='rainbow')
fig.colorbar(surf, ax=axs[1])  
axs[1].set_title('3D surface plot of concentration')
axs[1].set_xlabel('X:Length')
axs[1].set_ylabel('Y')
axs[1].set_zlabel('C')

mpt.tight_layout()  ###adjust the layout to avoid overlapping
mpt.show()
