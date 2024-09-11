import matplotlib.pyplot as plt
import numpy as np
from  matplotlib.animation import FuncAnimation
Ro=float(input('Initial value of radius in mm= '))
H=float(input('Relative Humidity= '))
k = (1 - H) * 6.06 * 10**-4
tmax = (Ro**2) / (2 * k)
fig,ax=plt.subplots()
ax.set_aspect('equal')
ax.set_xlabel('X values in nondimensional form ')
ax.set_ylabel('Y values of non dimensional form')
ax.set_title('Relative concentration contours')
r=np.linspace(Ro,10*Ro,500)
theta=np.linspace(0,np.pi*2,500)
r_grid,theta_grid=np.meshgrid(r,theta)
r_star=r_grid/Ro
c_star=((1-H)/r_star)+ H
x=r_star*np.sin(theta_grid)
y=r_star*np.cos(theta_grid)
cp=plt.contourf(x,y,c_star,levels=100,cmap='rainbow')
fig.colorbar(cp)

def init():
    return cp,

def animate(t): 
    ax.clear()
    ax.set_aspect('equal')
    ax.set_xlabel('X values in nondimensional form ')
    ax.set_ylabel('Y values of non dimensional form')
    ax.set_title('Relative concentration contours')
    R=((Ro**2) - (2 * k*t))**0.5 
    if R>0:
      r=np.linspace(R,10*Ro,500)
      theta=np.linspace(0,2*np.pi,500)
      r_grid,theta_grid=np.meshgrid(r,theta)
      r_star=r_grid/(R)
      c_star=((1-H)/r_star)+ H
      x=r_star*np.sin(theta_grid)
      y=r_star*np.cos(theta_grid)
      cp = ax.contourf(x,y,c_star,levels=100,cmap='rainbow') 
      return cp,
    if R==0:
       r_star=np.linspace(0,10,500)
       theta=np.linspace(0,2*np.pi,500)
       r_grid,theta_grid=np.meshgrid(r_star,theta)
       c_star=H*np.ones_like(r_grid)##makes sure that c_star always remains 2D aRRAY ##
       x=r_grid*np.sin(theta_grid)
       y=r_grid*np.cos(theta_grid)
       cp = ax.contourf(x,y,c_star,levels=100,cmap='rainbow')
       return cp,

ani=FuncAnimation(fig,animate,init_func=init,frames=np.linspace(0,tmax,200),interval=20,blit=False)###blit has to be false here
ani.save('Evaporating_Droplet_TEST4.gif',writer='pillow')###saving as gif##
plt.show()