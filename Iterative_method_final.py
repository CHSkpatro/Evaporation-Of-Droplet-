import matplotlib.pyplot as mpt
import numpy as np
H=float(input('The relative humidity value : '))
Ro=float(input('Contact line radius for the sessile droplet in mm = '))
theta=0.583 ###radian 
Ho= Ro *np.tan(theta/2) ### droplet height at center 
### Domain of computation ### 20 Ro to 20 ho ###
R_max=20*Ro
Z_max=20*Ho
##### Normalized dimensions####
R=Ro/R_max
ho=Ho/Z_max
C_max=2.32e-8  ### in gm/mm3 the saturation concentration of water vapour in air 
N=1200
M=1000
del_r=1/N
del_z=1/M
K=((R_max/Z_max)**2 )*(del_r/del_z)
r=np.linspace(0,1,N+1)### N divisions along r
z=np.linspace(0,1,M+1)### M divisions along z  ##### in contourf or whenever we are allocating 2d arrays to a meshgrid then the x and y can be 1d arrays and the function which is being allocated must be of shape (len(z),len(r))
C=H*np.ones((M+1,N+1)) #### Initialization of the entire domain as H relative concentration
#### Now we want to separate the exterior or interior of the sessile droplet using the cap shaped boundary###
R_grid,Z_grid=np.meshgrid(r,z)
def h(r,R,theta):
  term=(R / np.sin(theta))**2 - r**2
  h_value = np.sqrt(np.maximum(term, 0)) - (R / np.tan(theta))
  return h_value
h_grid=(R_max/Z_max)*h(R_grid,R,theta)
interior=np.where(Z_grid<=h_grid,1,0)
### inside let set the concentration as -1 ##### 
C=np.where(interior==1,-1,C) ### in interior points it is -1 other wise C will remain C 
for i in range(N):
    boundary_indices = np.where(z <= h(r[i], R, theta))[0]
    if boundary_indices.size > 0:
        boundary_index = boundary_indices[-1]  #### Choose the last boundary point for each r
        C[boundary_index, i] =1  
####Initiating iterative method###
iter_max=500
tolerance=1e-6 ###0.000001 ###
C_old=C.copy()
rms_error_list=[]
###To plot RMS error vs iteration number in real time ###
fig,ax=mpt.subplots()
ax.set_xlabel('Iteration Number')
ax.set_ylabel('RMS error')
mpt.ion() ###Turning the interactive mode ON ###

for iter_num in range(iter_max+1):
   C_new=C.copy() ### so that each time the loop runs the C_new will take the form of C and then as per the following lines the the elements will be replaced 
   for i in range(1,N):
      for j in range(1,M):
         if interior[j,i]==0:### only calculating for the exterior region only*****
            C_new[j,i]=((K*C_old[j-1,i]) + C_old[j,i-1] + (C_old[j,i+1]*(1+del_r/r[i])) + K*C_old[j+1,i] ) /(2 +(2*K)+ (del_r/r[i]))
         if j> int(ho/del_z):
          C_new[j:,0]=C_new[j:,1] ### for symmetry condition at r=0 and z>ho, no flux partial derivative of C wrt r is zero ###
      if i> int(R/del_r):
       C_new[0,i:]=C_new[1,i:] ### No penetration condition so no flux at z=0 for r>R so partial drerivative of C wrt z is zero ###
   rms_error=np.abs(np.sqrt(np.mean((C_new-C_old)**2)))
   rms_error_list.append(rms_error)
   C_old=C_new.copy()
   if round(rms_error,6) <= tolerance:
      print('Convergence Reached')
      break
   print(f"Iteartion={iter_num} : RMS_error={rms_error:.6f}")
   ax.plot(rms_error_list)
   mpt.pause(0.01)

mpt.ioff()
mpt.show()

mpt.figure()
mpt.contourf(r,z,C_new,levels=50,cmap='viridis')
mpt.colorbar()
mpt.contour(R_grid,Z_grid,h_grid,colors='red')
mpt.xlabel('Normalized Radius (r)')
mpt.ylabel('Normalized Height (z)')
mpt.title(' Relative Concentration Contours')
mpt.show()
