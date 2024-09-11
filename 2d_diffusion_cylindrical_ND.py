import numpy as np
import matplotlib.pyplot as mpt
import pandas as pd
R_min=float(input('The minimun radius value in mm: '))
R_max=float(input('The maximum radius value in mm: '))
N=int(input('No. of divisions of r axis : '))
r_min=(R_min/R_max)
del_r_nd=(1 - r_min)/N
Z_min=float(input('The minimun z value in mm: '))
Z_max=float(input('The maximum z value in mm: '))
M=int(input('No. of divisions of Z axis :'))
z_min=(Z_min/Z_max)
del_z_nd=(1- z_min)/M
print('Boundary_condition_choices :')
print('1) All 4 walls are assigned constant value.')
print('2) All walls except the bottom wall are assigned constant value and below wall has zero flux condition.')
print('3) Left wall has zero flux condition.')
print('4) Both bottom and top walls have zero flux condition.')
print('5) Both left and right walls have zero flux condition.')
choice=input('Choose number from the above for required boundary condition : ')
###### Dirichlet's Boundary condition at all 4 walls ###
if choice =='1':
 C1=float(input('Concentration value for r=r_min : '))
 C2=float(input('Concentration value for Z=Z_min : '))
 C3=float(input('Concentration value for r=r_max : '))
 C4=float(input('Concentration value for Z=Z_max : '))
 C_max=max(C1,C2,C3,C4)
 C1_nd=C1/C_max
 C2_nd=C2/C_max
 C3_nd=C3/C_max
 C4_nd=C4/C_max
 r=np.linspace(r_min,1,N+1)
 ##print(r)
 z=np.linspace(z_min,1,M+1)
 K=round(((R_max*del_r_nd)/(Z_max*del_z_nd))**2, 3)
 ###  A*C=  B ####
 A=np.zeros(((N-1)*(M-1),(N-1)*(M-1)))
 ###filling up diagonal elements###
 for i in range ((N-1)*(M-1)):
    A[i,i]=-round((2+(2*K)+(del_r_nd/r[(i%(N-1)) + 1])),3)
    if i< (((N-1)*(M-1))-1):  
        if i==0 or (i+1)%(N-1)!=0:
          A[i,i+1]=round(1+(del_r_nd/r[(i%(N-1)) + 1]),3)###super diagonal elements###but there are still slective rows in which the super diagonal element is 0
        if i==0 or (i+1)%(N-1)!=0:
          A[i+1,i]=1###sub diagonal elements ####but also there are selective rows in which the sub diagonal element is 0 
 for i in range((N-1),((N-1)*(M-1))):
    A[i,i-(N-1)]=K  ###filling the second sub diagonal elements###
    A[i-(N-1),i]=K ###filling up the 2nd super diagonal elements###
 #df = pd.DataFrame(A)
 #df.to_csv('matrix_A_cyl.csv', index=False)
 #print("Matrix A cyl has been saved to 'matrix_A.csv'.")    
 A_inv=np.linalg.inv(A)   
 ###For B ###
 B=np.zeros(((N-1)*(M-1),1))
 ##1st N-1 rows 
 B[0,0]=-C1_nd-(K*C2_nd)
 for i in range(1,N-2):
    B[i,0]=-(K*C2_nd)
 B[N-2,0]=-(K*C2_nd)-round(C3_nd*(1+(del_r_nd/r[N-1])),3)
 ##last N-1 rows
 B[-1,0]=-(K*C4_nd)-round(C3_nd*(1+(del_r_nd/r[N-1])),3)
 for i in range(-(N-2),-1):
      B[i,0]=-(K*C4_nd)
 B[-(N-1),0] =-C1_nd-(K*C4_nd)
  ##the middle rows
 for k in range(1,M-2):
    B[k*(N-1),0]=-C1_nd
    B[((k+1)*(N-1))-1,0]=-round(C3_nd*(1+(del_r_nd/r[N-1])),3)
 ##df = pd.DataFrame(B)
 ##df.to_csv('matrix_B_cyl.csv', index=False)
 ##print("Matrix B has been saved to 'matrix_B_cyl.csv'.") 
 C=np.dot(A_inv,B) ###C values obtained ###(C11,C21,C31,C41,.............) THE VALUES ARE IN THIS ORDER
 C_actual=np.zeros((M+1,N+1))
 C_actual[:,0]=C1_nd
 C_actual[0,:]=C2_nd
 C_actual[:,-1]=C3_nd
 C_actual[-1,:]=C4_nd
 C_actual[1:M,1:N]=np.round(C.reshape((M-1,N-1)),3) 
 ###conversion from cylindrical to cartesian ###
 theta=np.linspace(0,2*np.pi,N+1)
 R,Theta=np.meshgrid(r,theta) ###instead of 1d array 2d arrays are created###
 X=R*np.cos(Theta)
 Y=R*np.sin(Theta) 
 ###2d plot in xz plane ###
 mpt.figure()
 mpt.contourf(r, z, C_actual,levels=200, cmap='rainbow')  ###taking X here is causing shape mismatch between X and C_actual grid ?
 mpt.colorbar()
 mpt.xlabel('X: radial distance in non dimensional form,(r/R_max)*cos(theta)')
 mpt.ylabel('Z: height in non dimensional form z/Z_max')
 mpt.xlim(0,1)
 mpt.ylim(0,1)
 mpt.title('2D contour plot of relative concentration values C/C_max in XZ plane')
 mpt.show()
 ###2d plot in YZ plane###
 mpt.figure()
 mpt.contourf(r, z, C_actual,levels=200, cmap='rainbow')###taking Y here is causing shape mismatch between Y and C_actual grid ?
 mpt.colorbar()
 mpt.xlabel('Y : radial distance in non dimensional form,(r/R_max)*sin(theta)')
 mpt.ylabel('Z : height in non dimensional form z/Z_max')
 mpt.xlim(0,1)
 mpt.ylim(0,1)
 mpt.title('2D contour plot of relative concentration values C/C_max in YZ plane')
 mpt.show()
 ###2d plot in XY plane###at middle of Z##
 middle_z_index = len(z) // 2  ####Middle index for the Z values
 mpt.figure()
 C_middle_z = np.tile(C_actual[middle_z_index, :], (len(theta), 1))###np.tile recreates the c_actual middle to the same shape as required by repeating the the middle_z_index th row for len(theta) times aqnd omly one column
 mpt.contourf(X, Y, C_middle_z,levels=200, cmap='rainbow')
 mpt.colorbar()
 mpt.xlabel('X : radial distance in non dimensional form,(r/R_max)*cos(theta)')
 mpt.ylabel('Y : radial distance in non dimensional form,(r/R_max)*sin(theta)')
 mpt.gca().set_aspect('equal',adjustable='box')
 mpt.title('2D contour plot of relative concentration values in XY plane at middle Z')
 mpt.show()
 ###Now to draw gradient vector field###
 C_grad_r=np.zeros_like(C_actual)
 C_grad_z=np.zeros_like(C_actual)
 ###now defining them for left ,bottom and all interior points by forward difference method###
 C_grad_r[:,:-1]=-(C_max/R_max)*(C_actual[:,1:]-C_actual[:,:-1])/(del_r_nd)
 C_grad_z[:-1,:]=-(C_max/Z_max)*(C_actual[1:,:]-C_actual[:-1,:])/(del_z_nd)
 ###Now defining for right and top wall by backward difference method###
 C_grad_r[:,-1]=-(C_max/R_max)*(C_actual[:,-1]-C_actual[:,-2])/(del_r_nd)
 C_grad_z[-1,:]=-(C_max/Z_max)*(C_actual[-1,:]-C_actual[-2,:])/(del_z_nd)
 magnitude = np.sqrt(C_grad_r**2 + C_grad_z**2)
 zero_magnitude = (magnitude == 0)
 C_grad_R = np.where(zero_magnitude, 0, C_grad_r / np.where(magnitude == 0, 1, magnitude))
 C_grad_Z = np.where(zero_magnitude, 0, C_grad_z / np.where(magnitude == 0, 1, magnitude))
 ### quiver to plot the gradient vector field ###
 fig, ax = mpt.subplots(figsize=(10, 8))
 step = 3  ###adjust this value to change the density of the arrows
 r_quiver = r[::step]
 z_quiver = z[::step]
 R_quiver,Z_quiver=np.meshgrid(r_quiver,z_quiver)
 C_grad_R_quiver = C_grad_R[::step, ::step]
 C_grad_Z_quiver = C_grad_Z[::step, ::step]
 magnitude_quiver = magnitude[::step, ::step]
 print(R_quiver.shape,Z_quiver.shape,C_grad_R_quiver.shape,C_grad_Z_quiver.shape,magnitude_quiver.shape)
 Q = ax.quiver(R_quiver, Z_quiver, C_grad_R_quiver, C_grad_Z_quiver, magnitude_quiver, scale=40, cmap='rainbow')
 fig.colorbar(Q, ax=ax)
 ax.set_xlabel('X : radial distance in normalized form,(r/R_max)*cos(theta)')
 ax.set_ylabel('Z : height in normalized form z/Z_max')
 ax.set_title('Concentration gradient Vector field in non dimensional form ')
 ax.set_xlim(0,1)
 ax.set_ylim(0,1)
 mpt.show()
 ###streamline plot###
 fig, ax = mpt.subplots(figsize=(10, 8))
 strm = ax.streamplot(r, z, C_grad_R, C_grad_Z, color=magnitude,cmap='rainbow', density=2)
 ax.set_xlabel('X : radial distance in normalized form,(r/R_max)*cos(theta)')
 ax.set_ylabel('Z : height in normalized form z/Z_max')
 ax.set_title('Concentration gradient Vector field in non dimensional form Streamline plot ')
 ax.set_xlim(0,1)
 mpt.show()
#############################################################################################################
####### Only Bottom wall with no penetration condition #####
if choice =='2':
 C1=float(input('Concentration value for r=r_min : '))
 C3=float(input('Concentration value for r=r_max : '))
 C4=float(input('Concentration value for Z=Z_max : '))
 C_max=max(C1,C3,C4)
 C1_nd=C1/C_max
 C3_nd=C3/C_max
 C4_nd=C4/C_max
 r=np.linspace(r_min,1,N+1)
 ##print(r)
 z=np.linspace(z_min,1,M+1)
 K=((R_max*del_r_nd)/(Z_max*del_z_nd))**2
 ###  A*C=  B ####
 A=np.zeros(((N-1)*(M-1),(N-1)*(M-1)))
 ###filling up diagonal elements###
 for i in range ((N-1)*(M-1)):
    A[i,i]=round(-(2+(2*K)+(del_r_nd/r[(i%(N-1)) + 1])),3)
    if i< (((N-1)*(M-1))-1):  
        if i==0 or (i+1)%(N-1)!=0:
          A[i,i+1]=round(1+(del_r_nd/r[(i%(N-1)) + 1]),3)###super diagonal elements###but there are still slective rows in which the super diagonal element is 0
        if i==0 or (i+1)%(N-1)!=0:
          A[i+1,i]=1###sub diagonal elements ####but also there are selective rows in which the sub diagonal element is 0 
 for i in range(N-1): ### diagonal elements of the 1st N-1 
    A[i,i]=round(-(2+(K)+(del_r_nd/r[(i%(N-1)) + 1])),3)
 for i in range((N-1),((N-1)*(M-1))):
    A[i,i-(N-1)]=K  ###filling the second sub diagonal elements###
    A[i-(N-1),i]=K ###filling up the 2nd super diagonal elements###
 #df = pd.DataFrame(A)
 #df.to_csv('matrix_A_cyl.csv', index=False)
 #print("Matrix A cyl has been saved to 'matrix_A.csv'.")    
 A_inv=np.linalg.inv(A)   
 ###For B ###
 B=np.zeros(((N-1)*(M-1),1))
 ##1st N-1 rows 
 B[0,0]=-C1_nd
 for i in range(1,N-2):
    B[i,0]=0
 B[N-2,0]=-round(C3_nd*(1+(del_r_nd/r[N-1])),3)
 ##last N-1 rows
 B[-1,0]=-(K*C4_nd)-round(C3_nd*(1+(del_r_nd/r[N-1])),3)
 for i in range(-(N-2),-1):
      B[i,0]=-(K*C4_nd)
 B[-(N-1),0] =-C1_nd-(K*C4_nd)
  ##the middle rows
 for k in range(1,M-2):
    B[k*(N-1),0]=-C1_nd
    B[((k+1)*(N-1))-1,0]=-round(C3_nd*(1+(del_r_nd/r[N-1])),3)
 ##df = pd.DataFrame(B)
 ##df.to_csv('matrix_B_cyl.csv', index=False)
 ##print("Matrix B has been saved to 'matrix_B_cyl.csv'.") 
 C=np.dot(A_inv,B) ###C values obtained ###(C11,C21,C31,C41,.............) THE VALUES ARE IN THIS ORDER
 C_actual=np.zeros((M+1,N+1))
 C_actual[:,0]=C1_nd
 
 C_actual[:,-1]=C3_nd
 C_actual[-1,:]=C4_nd
 C_actual[1:M,1:N]=np.round(C.reshape((M-1,N-1)),3) 
 C_actual[0,:]=C_actual[1,:]
 ###conversion from cylindrical to cartesian ###
 theta=np.linspace(0,2*np.pi,N+1)
 R,Theta=np.meshgrid(r,theta) ###instead of 1d array 2d arrays are created###
 X=R*np.cos(Theta)
 Y=R*np.sin(Theta) 
 ###2d plot in xz plane ###
 mpt.figure()
 mpt.contourf(r, z, C_actual,levels=50, cmap='rainbow')  ###taking X here is causing shape mismatch between X and C_actual grid ?
 mpt.colorbar()
 mpt.xlabel('X: radial distance in non dimensional form,(r/R_max)*cos(theta)')
 mpt.ylabel('Z: height in non dimensional form z/Z_max')
 mpt.xlim(0,1)
 mpt.ylim(0,1)
 mpt.title('2D contour plot of relative concentration values C/C_max in XZ plane')
 mpt.show()
 ###2d plot in YZ plane###
 mpt.figure()
 mpt.contourf(r, z, C_actual,levels=50, cmap='rainbow')###taking Y here is causing shape mismatch between Y and C_actual grid ?
 mpt.colorbar()
 mpt.xlabel('Y : radial distance in non dimensional form,(r/R_max)*sin(theta)')
 mpt.ylabel('Z : height in non dimensional form z/Z_max')
 mpt.xlim(0,1)
 mpt.ylim(0,1)
 mpt.title('2D contour plot of relative concentration values C/C_max in YZ plane')
 mpt.show()
 ###2d plot in XY plane###at middle of Z##
 middle_z_index = len(z) // 2  ####Middle index for the Z values
 mpt.figure()
 C_middle_z = np.tile(C_actual[middle_z_index, :], (len(theta), 1))###np.tile recreates the c_actual middle to the same shape as required by repeating the the middle_z_index th row for len(theta) times aqnd omly one column
 mpt.contourf(X, Y, C_middle_z,levels=50, cmap='rainbow')
 mpt.colorbar()
 mpt.xlabel('X : radial distance in non dimensional form,(r/R_max)*cos(theta)')
 mpt.ylabel('Y : radial distance in non dimensional form,(r/R_max)*sin(theta)')
 mpt.gca().set_aspect('equal',adjustable='box')
 mpt.title('2D contour plot of relative concentration values in XY plane at middle Z')
 mpt.show()
 ###Now to draw gradient vector field###
 C_grad_r=np.zeros_like(C_actual)
 C_grad_z=np.zeros_like(C_actual)
 ###now defining them for left ,bottom and all interior points by forward difference method###
 C_grad_r[:,:-1]=-(C_max/R_max)*(C_actual[:,1:]-C_actual[:,:-1])/(del_r_nd)
 C_grad_z[:-1,:]=-(C_max/Z_max)*(C_actual[1:,:]-C_actual[:-1,:])/(del_z_nd)
 ###Now defining for right and top wall by backward difference method###
 C_grad_r[:,-1]=-(C_max/R_max)*(C_actual[:,-1]-C_actual[:,-2])/(del_r_nd)
 C_grad_z[-1,:]=-(C_max/Z_max)*(C_actual[-1,:]-C_actual[-2,:])/(del_z_nd)
 magnitude = np.sqrt(C_grad_r**2 + C_grad_z**2)
 zero_magnitude = (magnitude == 0)
 C_grad_R = np.where(zero_magnitude, 0, C_grad_r / np.where(magnitude == 0, 1, magnitude))
 C_grad_Z = np.where(zero_magnitude, 0, C_grad_z / np.where(magnitude == 0, 1, magnitude))
 ### quiver to plot the gradient vector field ###
 fig, ax = mpt.subplots(figsize=(10, 8))
 step = 2  ###adjust this value to change the density of the arrows
 r_quiver = r[::step]
 z_quiver = z[::step]
 R_quiver,Z_quiver=np.meshgrid(r_quiver,z_quiver)
 C_grad_R_quiver = C_grad_R[::step, ::step]
 C_grad_Z_quiver = C_grad_Z[::step, ::step]
 magnitude_quiver = magnitude[::step, ::step]
 print(R_quiver.shape,Z_quiver.shape,C_grad_R_quiver.shape,C_grad_Z_quiver.shape,magnitude_quiver.shape)
 Q = ax.quiver(R_quiver, Z_quiver, C_grad_R_quiver, C_grad_Z_quiver, magnitude_quiver, scale=40, cmap='rainbow')
 fig.colorbar(Q, ax=ax)
 ax.set_xlabel('X : radial distance in normalized form,(r/R_max)*cos(theta)')
 ax.set_ylabel('Z : height in normalized form z/Z_max')
 ax.set_title('Concentration gradient Vector field in non dimensional form ')
 ax.set_xlim(0,1)
 ax.set_ylim(0,1)
 mpt.show()
 ###streamline plot###
 fig, ax = mpt.subplots(figsize=(10, 8))
 strm = ax.streamplot(r, z, C_grad_R, C_grad_Z, color=magnitude,cmap='rainbow', density=2)
 ax.set_xlabel('X : radial distance in normalized form,(r/R_max)*cos(theta)')
 ax.set_ylabel('Z : height in normalized form z/Z_max')
 ax.set_title('Concentration gradient Vector field in non dimensional form Streamline plot ')
 ax.set_xlim(0,1)
 mpt.show()
#######################################################################################################################################
### Left wall has zero flux condition ###
if choice=='3':
 C3=float(input('Concentration value for r=r_max : '))
 C2=float(input('Concentration value for z=z_min : '))
 C4=float(input('Concentration value for Z=Z_max : '))
 C_max=max(C3,C2,C4)
 C3_nd=C3/C_max
 C2_nd=C2/C_max
 C4_nd=C4/C_max
 r=np.linspace(r_min,1,N+1)
 ##print(r)
 z=np.linspace(z_min,1,M+1)
 K=round(((R_max*del_r_nd)/(Z_max*del_z_nd))**2,3)
 ###  A*C=  B ####
 A=np.zeros(((N-1)*(M-1),(N-1)*(M-1)))
 ###filling up diagonal elements###
 for i in range ((N-1)*(M-1)):
    A[i,i]=round(-(2+(2*K)+(del_r_nd/r[(i%(N-1)) + 1])),3)
    if (i % (N-1)) ==0: ###the 1st row daigonal for every N-1 rows is different
       A[i,i]=round(-(1+(2*K)+(del_r_nd/r[(i%(N-1)) + 1])),3)
    if i< (((N-1)*(M-1))-1):  
        if i==0 or (i+1)%(N-1)!=0:
          A[i,i+1]=round(1+(del_r_nd/r[(i%(N-1)) + 1]),3)###super diagonal elements###but there are still slective rows in which the super diagonal element is 0
        if i==0 or (i+1)%(N-1)!=0:
          A[i+1,i]=1###sub diagonal elements ####but also there are selective rows in which the sub diagonal element is 0 
 for i in range((N-1),((N-1)*(M-1))):
    A[i,i-(N-1)]=K  ###filling the second sub diagonal elements###
    A[i-(N-1),i]=K ###filling up the 2nd super diagonal elements###
 #df = pd.DataFrame(A)
 #df.to_csv('matrix_A_cyl.csv', index=False)
 #print("Matrix A cyl has been saved to 'matrix_A.csv'.")    
 A_inv=np.linalg.inv(A)   
 ###For B ###
 B=np.zeros(((N-1)*(M-1),1))
 ###1st N-1 rows 
 B[0,0]=-(K*C2_nd)    ###INDEPENDENT OF C1
 for i in range(1,N-2):
    B[i,0]=-(K*C2_nd)
 B[N-2,0]=-(K*C2_nd)-round(C3_nd*(1+(del_r_nd/r[N-1])),3)
 ###last N-1 rows
 B[-1,0]=-(K*C4_nd)-round(C3_nd*(1+(del_r_nd/r[N-1])),3)
 for i in range(-(N-2),-1):
      B[i,0]=-(K*C4_nd)
 B[-(N-1),0] =-(K*C4_nd) ###Independent of C1
  ###the middle rows
 for k in range(1,M-2):
    B[k*(N-1),0]=0 ###independent of C1
    B[((k+1)*(N-1))-1,0]=-round(C3_nd*(1+(del_r_nd/r[N-1])),3)
 ##df = pd.DataFrame(B)
 ##df.to_csv('matrix_B_cyl.csv', index=False)
 ##print("Matrix B has been saved to 'matrix_B_cyl.csv'.") 
 C=np.dot(A_inv,B) ###C values obtained ###(C11,C21,C31,C41,.............) THE VALUES ARE IN THIS ORDER
 C_actual=np.zeros((M+1,N+1))
 
 C_actual[0,:]=C2_nd
 C_actual[:,-1]=C3_nd
 C_actual[-1,:]=C4_nd
 C_actual[1:M,1:N]=np.round(C.reshape((M-1,N-1)),3) 
 C_actual[:,0]=C_actual[:,1]   ###1st two columns are equal###
 ###conversion from cylindrical to cartesian ###
 theta=np.linspace(0,2*np.pi,N+1)
 R,Theta=np.meshgrid(r,theta) ###instead of 1d array 2d arrays are created###
 X=R*np.cos(Theta)
 Y=R*np.sin(Theta) 
 ###2d plot in xz plane ###
 mpt.figure()
 mpt.contourf(r, z, C_actual,levels=50, cmap='rainbow')  ###taking X here is causing shape mismatch between X and C_actual grid ?
 mpt.colorbar()
 mpt.xlabel('X: radial distance in non dimensional form,(r/R_max)*cos(theta)')
 mpt.ylabel('Z: height in non dimensional form z/Z_max')
 mpt.xlim(0,1)
 mpt.ylim(0,1)
 mpt.title('2D contour plot of relative concentration values C/C_max in XZ plane')
 mpt.show()
 ###2d plot in YZ plane###
 mpt.figure()
 mpt.contourf(r, z, C_actual,levels=50, cmap='rainbow')###taking Y here is causing shape mismatch between Y and C_actual grid ?
 mpt.colorbar()
 mpt.xlabel('Y : radial distance in non dimensional form,(r/R_max)*sin(theta)')
 mpt.ylabel('Z : height in non dimensional form z/Z_max')
 mpt.xlim(0,1)
 mpt.ylim(0,1)
 mpt.title('2D contour plot of relative concentration values C/C_max in YZ plane')
 mpt.show()
 ###2d plot in XY plane###at middle of Z##
 middle_z_index = len(z) // 2  ####Middle index for the Z values
 mpt.figure()
 C_middle_z = np.tile(C_actual[middle_z_index, :], (len(theta), 1))###np.tile recreates the c_actual middle to the same shape as required by repeating the the middle_z_index th row for len(theta) times aqnd omly one column
 mpt.contourf(X, Y, C_middle_z,levels=50, cmap='rainbow')
 mpt.colorbar()
 mpt.xlabel('X : radial distance in non dimensional form,(r/R_max)*cos(theta)')
 mpt.ylabel('Y : radial distance in non dimensional form,(r/R_max)*sin(theta)')
 mpt.gca().set_aspect('equal',adjustable='box')
 mpt.title('2D contour plot of relative concentration values in XY plane at middle Z')
 mpt.show()
 ###Now to draw gradient vector field###
 C_grad_r=np.zeros_like(C_actual)
 C_grad_z=np.zeros_like(C_actual)
 ###now defining them for left ,bottom and all interior points by forward difference method###
 C_grad_r[:,:-1]=-(C_max/R_max)*(C_actual[:,1:]-C_actual[:,:-1])/(del_r_nd)
 C_grad_z[:-1,:]=-(C_max/Z_max)*(C_actual[1:,:]-C_actual[:-1,:])/(del_z_nd)
 ###Now defining for right and top wall by backward difference method###
 C_grad_r[:,-1]=-(C_max/R_max)*(C_actual[:,-1]-C_actual[:,-2])/(del_r_nd)
 C_grad_z[-1,:]=-(C_max/Z_max)*(C_actual[-1,:]-C_actual[-2,:])/(del_z_nd)
 magnitude = np.sqrt(C_grad_r**2 + C_grad_z**2)
 zero_magnitude = (magnitude == 0)
 C_grad_R = np.where(zero_magnitude, 0, C_grad_r / np.where(magnitude == 0, 1, magnitude))
 C_grad_Z = np.where(zero_magnitude, 0, C_grad_z / np.where(magnitude == 0, 1, magnitude))
 ### quiver to plot the gradient vector field ###
 fig, ax = mpt.subplots(figsize=(10, 8))
 step = 3  ###adjust this value to change the density of the arrows
 r_quiver = r[::step]
 z_quiver = z[::step]
 R_quiver,Z_quiver=np.meshgrid(r_quiver,z_quiver)
 C_grad_R_quiver = C_grad_R[::step, ::step]
 C_grad_Z_quiver = C_grad_Z[::step, ::step]
 magnitude_quiver = magnitude[::step, ::step]
 print(R_quiver.shape,Z_quiver.shape,C_grad_R_quiver.shape,C_grad_Z_quiver.shape,magnitude_quiver.shape)
 Q = ax.quiver(R_quiver, Z_quiver, C_grad_R_quiver, C_grad_Z_quiver, magnitude_quiver, scale=40, cmap='rainbow')
 fig.colorbar(Q, ax=ax)
 ax.set_xlabel('X : radial distance in normalized form,(r/R_max)*cos(theta)')
 ax.set_ylabel('Z : height in normalized form z/Z_max')
 ax.set_title('Concentration gradient Vector field in non dimensional form ')
 ax.set_xlim(0,1)
 ax.set_ylim(0,1)
 mpt.show()
 ###streamline plot###
 fig, ax = mpt.subplots(figsize=(10, 8))
 strm = ax.streamplot(r, z, C_grad_R, C_grad_Z, color=magnitude,cmap='rainbow', density=2)
 ax.set_xlabel('X : radial distance in normalized form,(r/R_max)*cos(theta)')
 ax.set_ylabel('Z : height in normalized form z/Z_max')
 ax.set_title('Concentration gradient Vector field in non dimensional form Streamline plot ')
 ax.set_xlim(0,1)
 mpt.show()
#######################################################################################################################################
##########  Both bottom and top wall with no flux condition #########
if choice =='4':
 C1=float(input('Concentration value for r=r_min : '))
 C3=float(input('Concentration value for r=r_max : '))
 C_max=max(C1,C3)
 C1_nd=C1/C_max
 C3_nd=C3/C_max
 
 r=np.linspace(r_min,1,N+1)
 ##print(r)
 z=np.linspace(z_min,1,M+1)
 K=round(((R_max*del_r_nd)/(Z_max*del_z_nd))**2 ,3)
 ###  A*C=  B ####
 A=np.zeros(((N-1)*(M-1),(N-1)*(M-1)))
 ###filling up diagonal elements###
 for i in range ((N-1)*(M-1)):
    A[i,i]=round(-(2+(2*K)+(del_r_nd/r[(i%(N-1)) + 1])),3)
    if i< (((N-1)*(M-1))-1):  
        if i==0 or (i+1)%(N-1)!=0:
          A[i,i+1]=round(1+(del_r_nd/r[(i%(N-1)) + 1]),3)###super diagonal elements###but there are still slective rows in which the super diagonal element is 0
        if i==0 or (i+1)%(N-1)!=0:
          A[i+1,i]=1###sub diagonal elements ####but also there are selective rows in which the sub diagonal element is 0 
 for i in range(N-1): ### diagonal elements of the 1st N-1 rows
    A[i,i]=round(-(2+(K)+(del_r_nd/r[(i%(N-1)) + 1])),3)
 for i in range (-(N-1),0):###diagonal elements of last N-1 rows
    A[i,i]=round(-(2+(K)+(del_r_nd/r[i+N])),3)
 for i in range((N-1),((N-1)*(M-1))):
    A[i,i-(N-1)]=K  ###filling the second sub diagonal elements###
    A[i-(N-1),i]=K ###filling up the 2nd super diagonal elements###
 #df = pd.DataFrame(A)
 #df.to_csv('matrix_A_cyl.csv', index=False)
 #print("Matrix A cyl has been saved to 'matrix_A.csv'.")   
 #print(K) 
 A_inv=np.linalg.inv(A)   
 ###For B ###
 B=np.zeros(((N-1)*(M-1),1))
 ##1st N-1 rows 
 B[0,0]=-C1_nd
 for i in range(1,N-2):
    B[i,0]=0
 B[N-2,0]=-round(C3_nd*(1+(del_r_nd/r[N-1])),3)
 ##last N-1 rows
 B[-1,0]=-round(C3_nd*(1+(del_r_nd/r[N-1])),3)
 for i in range(-(N-2),-1):
      B[i,0]=0
 B[-(N-1),0] =-C1_nd
  ##the middle rows
 for k in range(1,M-2):
    B[k*(N-1),0]=-C1_nd
    B[((k+1)*(N-1))-1,0]=-round(C3_nd*(1+(del_r_nd/r[N-1])),3)
 ##df = pd.DataFrame(B)
 ##df.to_csv('matrix_B_cyl.csv', index=False)
 ##print("Matrix B has been saved to 'matrix_B_cyl.csv'.") 
 C=np.dot(A_inv,B) ###C values obtained ###(C11,C21,C31,C41,.............) THE VALUES ARE IN THIS ORDER
 C_actual=np.zeros((M+1,N+1))
 C_actual[:,0]=C1_nd
 C_actual[:,-1]=C3_nd
 C_actual[1:M,1:N]=np.round(C.reshape((M-1,N-1)),3) 
 C_actual[0,:]=C_actual[1,:]
 C_actual[-1,:]=C_actual[-2,:]
 ###conversion from cylindrical to cartesian ###
 theta=np.linspace(0,2*np.pi,N+1)
 R,Theta=np.meshgrid(r,theta) ###instead of 1d array 2d arrays are created###
 X=R*np.cos(Theta)
 Y=R*np.sin(Theta) 
 ###2d plot in xz plane ###
 mpt.figure()
 mpt.contourf(r, z, C_actual,levels=50, cmap='rainbow')  ###taking X here is causing shape mismatch between X and C_actual grid ?
 mpt.colorbar()
 mpt.xlabel('X: radial distance in non dimensional form,(r/R_max)*cos(theta)')
 mpt.ylabel('Z: height in non dimensional form z/Z_max')
 mpt.xlim(0,1)
 mpt.ylim(0,1)
 mpt.title('2D contour plot of relative concentration values C/C_max in XZ plane')
 mpt.show()
 ###2d plot in YZ plane###
 mpt.figure()
 mpt.contourf(r, z, C_actual,levels=50, cmap='rainbow')###taking Y here is causing shape mismatch between Y and C_actual grid ?
 mpt.colorbar()
 mpt.xlabel('Y : radial distance in non dimensional form,(r/R_max)*sin(theta)')
 mpt.ylabel('Z : height in non dimensional form z/Z_max')
 mpt.xlim(0,1)
 mpt.ylim(0,1)
 mpt.title('2D contour plot of relative concentration values C/C_max in YZ plane')
 mpt.show()
 ###2d plot in XY plane###at middle of Z##
 middle_z_index = len(z) // 2  ####Middle index for the Z values
 mpt.figure()
 C_middle_z = np.tile(C_actual[middle_z_index, :], (len(theta), 1))###np.tile recreates the c_actual middle to the same shape as required by repeating the the middle_z_index th row for len(theta) times aqnd omly one column
 mpt.contourf(X, Y, C_middle_z,levels=50, cmap='rainbow')
 mpt.colorbar()
 mpt.xlabel('X : radial distance in non dimensional form,(r/R_max)*cos(theta)')
 mpt.ylabel('Y : radial distance in non dimensional form,(r/R_max)*sin(theta)')
 mpt.gca().set_aspect('equal',adjustable='box')
 mpt.title('2D contour plot of relative concentration values in XY plane at middle Z')
 mpt.show()
 ###Now to draw gradient vector field###
 C_grad_r=np.zeros_like(C_actual)
 C_grad_z=np.zeros_like(C_actual)
 ###now defining them for left ,bottom and all interior points by forward difference method###
 C_grad_r[:,:-1]=-(C_max/R_max)*(C_actual[:,1:]-C_actual[:,:-1])/(1/N)
 C_grad_z[:-1,:]=-(C_max/Z_max)*(C_actual[1:,:]-C_actual[:-1,:])/(1/M)
 ###Now defining for right and top wall by backward difference method###
 C_grad_r[:,-1]=-(C_max/R_max)*(C_actual[:,-1]-C_actual[:,-2])/(1/N)
 C_grad_z[-1,:]=-(C_max/Z_max)*(C_actual[-1,:]-C_actual[-2,:])/(1/M)
 magnitude = np.sqrt(C_grad_r**2 + C_grad_z**2)
 zero_magnitude = (magnitude == 0)
 C_grad_R = np.where(zero_magnitude, 0, C_grad_r / np.where(magnitude == 0, 1, magnitude))
 C_grad_Z = np.where(zero_magnitude, 0, C_grad_z / np.where(magnitude == 0, 1, magnitude))
 ### quiver to plot the gradient vector field ###
 fig, ax = mpt.subplots(figsize=(10, 8))
 step = 3  ###adjust this value to change the density of the arrows
 r_quiver = r[::step]
 z_quiver = z[::step]
 R_quiver,Z_quiver=np.meshgrid(r_quiver,z_quiver)
 C_grad_R_quiver = C_grad_R[::step, ::step]
 C_grad_Z_quiver = C_grad_Z[::step, ::step]
 magnitude_quiver = magnitude[::step, ::step]
 print(R_quiver.shape,Z_quiver.shape,C_grad_R_quiver.shape,C_grad_Z_quiver.shape,magnitude_quiver.shape)
 Q = ax.quiver(R_quiver, Z_quiver, C_grad_R_quiver, C_grad_Z_quiver, magnitude_quiver, scale=40, cmap='rainbow')
 fig.colorbar(Q, ax=ax)
 ax.set_xlabel('X : radial distance in normalized form,(r/R_max)*cos(theta)')
 ax.set_ylabel('Z : height in normalized form z/Z_max')
 ax.set_title('Concentration gradient Vector field in non dimensional form ')
 ax.set_xlim(0,1)
 ax.set_ylim(0,1)
 mpt.show()
 ###streamline plot###
 fig, ax = mpt.subplots(figsize=(10, 8))
 strm = ax.streamplot(r, z, C_grad_R, C_grad_Z, color=magnitude,cmap='rainbow', density=2)
 ax.set_xlabel('X : radial distance in normalized form,(r/R_max)*cos(theta)')
 ax.set_ylabel('Z : height in normalized form z/Z_max')
 ax.set_title('Concentration gradient Vector field in non dimensional form Streamline plot ')
 ax.set_xlim(0,1)
 mpt.show()
#######################################################################################################
############ When both Left and right wall have zero flux condition ###
if choice=='5':
 C2=float(input('Concentration value for z=z_min : '))
 C4=float(input('Concentration value for Z=Z_max : '))
 C_max=max(C2,C4)
 C2_nd=C2/C_max
 C4_nd=C4/C_max
 r=np.linspace(r_min,1,N+1)
 ##print(r)
 z=np.linspace(z_min,1,M+1)
 K=round(((R_max*del_r_nd)/(Z_max*del_z_nd))**2,3)
 ###  A*C=  B ####
 A=np.zeros(((N-1)*(M-1),(N-1)*(M-1)))
 ###filling up diagonal elements###
 for i in range ((N-1)*(M-1)):
    A[i,i]=round(-(2+(2*K)+(del_r_nd/r[(i%(N-1)) + 1])),3)
    if (i % (N-1)) ==0: ###the 1st row daigonal for every N-1 rows is different
       A[i,i]=round(-(1+(2*K)+(del_r_nd/r[1])),3)
       
    if i< (((N-1)*(M-1))-1):  
        if i==0 or (i+1)%(N-1)!=0:
          A[i,i+1]=round(1+(del_r_nd/r[(i%(N-1)) + 1]),3)###super diagonal elements###but there are still slective rows in which the super diagonal element is 0
        if i==0 or (i+1)%(N-1)!=0:
          A[i+1,i]=1###sub diagonal elements ####but also there are selective rows in which the sub diagonal element is 0 
 for X in range(1,M):
    A[(X*(N-1))-1,(X*(N-1))-1]=-1-(2*K) #####filling last row diagonal of every N-1 rows with -1-2K
 for i in range((N-1),((N-1)*(M-1))):
    A[i,i-(N-1)]=K  ###filling the second sub diagonal elements###
    A[i-(N-1),i]=K ###filling up the 2nd super diagonal elements###
 #df = pd.DataFrame(A)
 #df.to_csv('matrix_A_cyl_nd.csv', index=False)
 #print("Matrix A cyl has been saved to 'matrix_A.csv'.")    
 #print(K)
 #print(-1-(2*K))
 #print(round(-(1+(2*K)+(del_r_nd/r[1])),3))
 A_inv=np.linalg.inv(A)   
 ###For B ###
 B=np.zeros(((N-1)*(M-1),1))
 ###1st N-1 rows 
 B[0,0]=-(K*C2_nd)    ###INDEPENDENT OF C1
 for i in range(1,N-2):
    B[i,0]=-(K*C2_nd)
 B[N-2,0]=-(K*C2_nd)
 ###last N-1 rows
 B[-1,0]=-(K*C4_nd)
 for i in range(-(N-2),-1):
      B[i,0]=-(K*C4_nd)
 B[-(N-1),0] =-(K*C4_nd) ###Independent of C1
    
 #df = pd.DataFrame(B)
 #df.to_csv('matrix_B_cyl.csv', index=False)
 #print("Matrix B has been saved to 'matrix_B_cyl.csv'.") 
 C=np.dot(A_inv,B) ###C values obtained ###(C11,C21,C31,C41,.............) THE VALUES ARE IN THIS ORDER
 C_actual=np.zeros((M+1,N+1))
 
 C_actual[0,:]=C2_nd
 
 C_actual[-1,:]=C4_nd
 C_actual[1:M,1:N]=np.round(C.reshape((M-1,N-1)),3) 
 C_actual[:,0]=C_actual[:,1]   ###1st two columns are equal###
 C_actual[:,-1]=C_actual[:,-2]  #### last two columns are equal###
 ###conversion from cylindrical to cartesian ###
 theta=np.linspace(0,2*np.pi,N+1)
 R,Theta=np.meshgrid(r,theta) ###instead of 1d array 2d arrays are created###
 X=R*np.cos(Theta)
 Y=R*np.sin(Theta) 
 ###2d plot in xz plane ###
 mpt.figure()
 mpt.contourf(r, z, C_actual,levels=50, cmap='rainbow')  ###taking X here is causing shape mismatch between X and C_actual grid ?
 mpt.colorbar()
 mpt.xlabel('X: radial distance in non dimensional form,(r/R_max)*cos(theta)')
 mpt.ylabel('Z: height in non dimensional form z/Z_max')
 mpt.xlim(0,1)
 mpt.ylim(0,1)
 mpt.title('2D contour plot of relative concentration values C/C_max in XZ plane')
 mpt.show()
 ###2d plot in YZ plane###
 mpt.figure()
 mpt.contourf(r, z, C_actual,levels=50, cmap='rainbow')###taking Y here is causing shape mismatch between Y and C_actual grid ?
 mpt.colorbar()
 mpt.xlabel('Y : radial distance in non dimensional form,(r/R_max)*sin(theta)')
 mpt.ylabel('Z : height in non dimensional form z/Z_max')
 mpt.xlim(0,1)
 mpt.ylim(0,1)
 mpt.title('2D contour plot of relative concentration values C/C_max in YZ plane')
 mpt.show()
 ###2d plot in XY plane###at middle of Z##
 middle_z_index = len(z) // 2  ####Middle index for the Z values
 mpt.figure()
 C_middle_z = np.tile(C_actual[middle_z_index, :], (len(theta), 1))###np.tile recreates the c_actual middle to the same shape as required by repeating the the middle_z_index th row for len(theta) times aqnd omly one column
 mpt.contourf(X, Y, C_middle_z,levels=50, cmap='rainbow')
 mpt.colorbar()
 mpt.xlabel('X : radial distance in non dimensional form,(r/R_max)*cos(theta)')
 mpt.ylabel('Y : radial distance in non dimensional form,(r/R_max)*sin(theta)')
 mpt.gca().set_aspect('equal',adjustable='box')
 mpt.title('2D contour plot of relative concentration values in XY plane at middle Z')
 mpt.show()
 ###Now to draw gradient vector field###
 C_grad_r=np.zeros_like(C_actual)
 C_grad_z=np.zeros_like(C_actual)
 ###now defining them for left ,bottom and all interior points by forward difference method###
 C_grad_r[:,:-1]=-(C_max/R_max)*(C_actual[:,1:]-C_actual[:,:-1])/(del_r_nd)
 C_grad_z[:-1,:]=-(C_max/Z_max)*(C_actual[1:,:]-C_actual[:-1,:])/(del_z_nd)
 ###Now defining for right and top wall by backward difference method###
 C_grad_r[:,-1]=-(C_max/R_max)*(C_actual[:,-1]-C_actual[:,-2])/(del_r_nd)
 C_grad_z[-1,:]=-(C_max/Z_max)*(C_actual[-1,:]-C_actual[-2,:])/(del_z_nd)
 magnitude = np.sqrt(C_grad_r**2 + C_grad_z**2)
 zero_magnitude = (magnitude == 0)
 C_grad_R = np.where(zero_magnitude, 0, C_grad_r / np.where(magnitude == 0, 1, magnitude))
 C_grad_Z = np.where(zero_magnitude, 0, C_grad_z / np.where(magnitude == 0, 1, magnitude))
 ### quiver to plot the gradient vector field ###
 fig, ax = mpt.subplots(figsize=(10, 8))
 step = 2  ###adjust this value to change the density of the arrows
 r_quiver = r[::step]
 z_quiver = z[::step]
 R_quiver,Z_quiver=np.meshgrid(r_quiver,z_quiver)
 C_grad_R_quiver = C_grad_R[::step, ::step]
 C_grad_Z_quiver = C_grad_Z[::step, ::step]
 magnitude_quiver = magnitude[::step, ::step]
 print(R_quiver.shape,Z_quiver.shape,C_grad_R_quiver.shape,C_grad_Z_quiver.shape,magnitude_quiver.shape)
 Q = ax.quiver(R_quiver, Z_quiver, C_grad_R_quiver, C_grad_Z_quiver, magnitude_quiver, scale=40, cmap='rainbow')
 fig.colorbar(Q, ax=ax)
 ax.set_xlabel('X : radial distance in normalized form,(r/R_max)*cos(theta)')
 ax.set_ylabel('Z : height in normalized form z/Z_max')
 ax.set_title('Concentration gradient Vector field in non dimensional form ')
 ax.set_xlim(0,1)
 ax.set_ylim(0,1)
 mpt.show()
 ###streamline plot###
 fig, ax = mpt.subplots(figsize=(10, 8))
 strm = ax.streamplot(r, z, C_grad_R, C_grad_Z, color=magnitude,cmap='rainbow', density=2)
 ax.set_xlabel('X : radial distance in normalized form,(r/R_max)*cos(theta)')
 ax.set_ylabel('Z : height in normalized form z/Z_max')
 ax.set_title('Concentration gradient Vector field in non dimensional form Streamline plot ')
 ax.set_xlim(0,1)
 mpt.show()