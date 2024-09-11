import numpy as np
import matplotlib.pyplot as mpt
H=float(input('The relative humidity value : '))
Ro=float(input('Contact line radius for the sessile droplet in mm = '))
theta=np.pi/(4.5)
ho= (Ro/np.sin(theta))-(Ro/np.tan(theta))
####domain of computation###
R1=5*Ro
R=Ro/R1
Z1=5*ho
r_max=1
z_max=1
N=1000
del_r=1/N
M=1000
K=((R1/Z1)**2) *((M/N)**2)
r=np.linspace(0,r_max,N)
z=np.linspace(0,z_max,M)
r_grid,z_grid=np.meshgrid(r,z)    ### mesh grid is formed ###
#### i want to differentiate the interior from exterior region , by the help of spherical cap shaped boundary ###
def h(r,R,theta):
    return np.sqrt((R/np.sin(theta))**2 -r**2) - (R/np.tan(theta))
interior=np.where(z_grid<= h(r_grid,R,theta),1,0)
##### initializing concentration field###
C=H*np.ones_like(r_grid)
### lets set the initial values inside the boundary different random values###
random_concn_values=1+ 9*np.random.random(C.shape)
C=np.where(interior==1,random_concn_values,C)
##### now for cap boundary ###
for i in range(N):
    boundary_index=np.where(z_grid[i]<=h(r_grid[i],R,theta))[0]
    if boundary_index.size > 0:
        boundary_index=boundary_index[-1]
        C[i,boundary_index]=1 ####### concentration value at boundary of cap####
####initiating iterative method###
tolerance=2e-6
iter_max=2000
C_old=C.copy()
rms_errror_list=[]
#### to plot iteration vs RMS error values in real time ###
fig,ax=mpt.subplots()
ax.set_xlabel('Iteration Number')
ax.set_ylabel('RMS error')
###Turn on interactive Mode###
mpt.ion()
for iter_num in range (iter_max):
    C_new=C.copy()
    for i in range(1,N-1):#####only calculating at the entire domain interior points ###
        for j in range(1,M-1):
            if interior[i,j]==0:  #### we want to compute at external grid points of the sessile droplets###
                C_new[i,j]=(K*C_old[i,j-1] + C_old[i-1,j] +((1+ (del_r/r[i]))*C_old[i+1,j]) + K*C_old[i,j+1])/ (2+ 2*K +(del_r/r[i]))###explicit form of the iterstive eqn###
    if i==0  and  j > int(ho/(Z1/M)):
                    C_new[i,j]=C_new[i+1,j]
    if j==0:
                    C_new[i,j]=C_new[i,j+1]
    rms_error=np.abs(np.sqrt(np.mean(C_new **2))-np.sqrt(np.mean(C_old **2)))
    rms_errror_list.append(rms_error)
    C_old=C_new.copy()
    if round(rms_error,6) <= tolerance :
            break
    print(f"Iteartion={iter_num} : RMS_error={rms_error:.6f}")
        ### update the iteration error plot ###
    ax.plot(range(len(rms_errror_list)),rms_errror_list)
    mpt.pause(0.01)### 0.01 sec pause before updating###

mpt.ioff()
mpt.show()
### plotting the concentration contours ##
fig,ax=mpt.subplots()
###contour_plot=ax.contour(r_grid,z_grid,interior,levels=[0.5],colors='red')
contourf_plot=ax.contourf(r_grid,z_grid,C_new,levels=50,cmap='rainbow')
ax.set_title('Concentration contour for a sessile droplet')
ax.set_xlabel('radial distance in normalized form : r/R_max')
ax.set_ylabel('height in normalized form : z/Z_max')
fig.colorbar(contourf_plot,ax=ax,label='Relative Concentration Values :C/Cv')
mpt.show()