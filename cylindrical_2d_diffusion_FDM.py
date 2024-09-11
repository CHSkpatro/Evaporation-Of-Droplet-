import numpy as np
import matplotlib.pyplot as mpt
import pandas as pd
r_min=float(input('The minimun radius value in mm: '))
r_max=float(input('The maximum radius value in mm: '))
N=int(input('No. of divisions of r axis : '))
del_r=(r_max-r_min)/N
Z_min=float(input('The minimun z value in mm: '))
Z_max=float(input('The maximum z value in mm: '))
M=round((Z_max-Z_min)/(del_r))
print(M)
C1=float(input('Concentration value for r=r_min : '))
C2=float(input('Concentration value for Z=Z_min : '))
C3=float(input('Concentration value for r=r_max : '))
C4=float(input('Concentration value for Z=Z_max : '))
##print(del_r)
##print((Z_max-Z_min)/M)
r=np.linspace(r_min,r_max,N+1)
##print(r)
z=np.linspace(Z_min,Z_max,M+1)
###  A*C=  B ####
A=np.zeros(((N-1)*(M-1),(N-1)*(M-1)))
###filling up diagonal elements###
for i in range ((N-1)*(M-1)):
    A[i,i]=round(-(4 +(del_r/r[(i%(N-1)) + 1])),3)
    if i< (((N-1)*(M-1))-1):  
        if i==0 or (i+1)%(N-1)!=0:
          A[i,i+1]=round(1+(del_r/r[(i%(N-1)) + 1]),3)###super diagonal elements###but there are still slective rows in which the super diagonal element is 0
        if i==0 or (i+1)%(N-1)!=0:
          A[i+1,i]=1###sub diagonal elements ####but also there are selective rows in which the sub diagonal element is 0 
for i in range((N-1),((N-1)*(M-1))):
    A[i,i-(N-1)]=1  ###filling the second sub diagonal elements###
    A[i-(N-1),i] =1 ###filling up the 2nd super diagonal elements###
#df = pd.DataFrame(A)
#df.to_csv('matrix_A_cyl.csv', index=False)
#print("Matrix A cyl has been saved to 'matrix_A.csv'.")    
A_inv=np.linalg.inv(A)   
###For B ###
B=np.zeros(((N-1)*(M-1),1))
##1st N-1 rows 
B[0,0]=-C1-C2
for i in range(1,N-2):
    B[i,0]=-C2 
B[N-2,0]=-C2-round(C3*(1+(del_r/r[N-1])),3)
##last N-1 rows
B[-1,0]=-C4-round(C3*(1+(del_r/r[N-1])),3)
for i in range(-(N-2),-1):
      B[i,0]=-C4
B[-(N-1),0] =-C1-C4
##the middle rows
for k in range(1,M-2):
    B[k*(N-1),0]=-C1
    B[((k+1)*(N-1))-1,0]=-round(C3*(1+(del_r/r[N-1])),3)
##df = pd.DataFrame(B)
##df.to_csv('matrix_B_cyl.csv', index=False)
##print("Matrix B has been saved to 'matrix_B_cyl.csv'.") 
C=np.dot(A_inv,B) ###C values obtained ###(C11,C21,C31,C41,.............) THE VALUES ARE IN THIS ORDER
C_actual=np.zeros((M+1,N+1))
C_actual[:,0]=C1
C_actual[0,:]=C2
C_actual[:,-1]=C3
C_actual[-1,:]=C4
C_actual[1:M,1:N]=np.round(C.reshape((M-1,N-1)),3) 
###conversion from cylindrical to cartesian ###
theta=np.linspace(0,2*np.pi,N+1)
R,Theta=np.meshgrid(r,theta) ###instead of 1d array 2d arrays are created###
X=R*np.cos(Theta)
Y=R*np.sin(Theta) 
###2d plot in xz plane ###
mpt.figure()
mpt.contourf(r, z, C_actual, cmap='rainbow')  ###taking X here is causing shape mismatch between X and C_actual grid ?
mpt.colorbar()
mpt.xlabel('X: radial distance in mm')
mpt.ylabel('Z : height in mm')
mpt.title('2D contour plot of concentration values in XZ plane')
mpt.show()
###2d plot in YZ plane###
mpt.figure()
mpt.contourf(r, z, C_actual, cmap='rainbow')###taking Y here is causing shape mismatch between Y and C_actual grid ?
mpt.colorbar()
mpt.xlabel('Y : radial distance in mm')
mpt.ylabel('Z : height in mm')
mpt.title('2D contour plot of concentration values in YZ plane')
mpt.show()
###2d plot in XY plane###at middle of Z##
middle_z_index = len(z) // 2  ####Middle index for the Z values
mpt.figure()
C_middle_z = np.tile(C_actual[middle_z_index, :], (len(theta), 1))###np.tile recreates the c_actual middle to the same shape as required by repeating the the middle_z_index th row for len(theta) times aqnd omly one column
mpt.contourf(X, Y, C_middle_z, cmap='rainbow')
mpt.colorbar()
mpt.xlabel('X : radial distance in mm')
mpt.ylabel('Y : radial distance in mm')
mpt.gca().set_aspect('equal',adjustable='box')
mpt.title('2D contour plot of concentration values in XY plane at middle Z')
mpt.show()