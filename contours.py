import numpy as np
import matplotlib.pyplot as mpt
R=float(input('The radius of droplet in mm='))
H=float(input('Relative Humidity = '))
grid_size=500
r_min=R
r_max=10*R
theta=np.linspace(0,2*np.pi,grid_size)
r=np.linspace(r_min/R, r_max/R, grid_size)
r_grid,Theta_grid=np.meshgrid(r,theta)
X=r_grid*np.cos(Theta_grid)
Y=r_grid*np.sin(Theta_grid)
r_star=r_grid/R
c_star=((1-H)/r_star) + H
mpt.figure(figsize=(8,6))
cp=mpt.contourf(X,Y,c_star,levels=1000,cmap='rainbow')
mpt.colorbar(cp)
mpt.title('Concentration contours around an evaporating spherical water droplet')
mpt.xlabel('(r/R)')
mpt.ylabel('r/R')
mpt.gca().set_aspect('equal',adjustable='box')
mpt.xlim([-10,10])
mpt.ylim([-10,10])
#mpt.imshow(cp)
mpt.show()
