import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

Ro = float(input('Initial value of radius in mm= '))
H = float(input('Relative Humidity= '))
k = (1 - H) * 6.06 * 10**-4
tmax = (Ro**2) / (2 * k)

fig, ax = plt.subplots()
ax.set_aspect('equal')
ax.set_xlabel('Distance from center (mm)')
ax.set_ylabel('Distance from center (mm)')
ax.set_title('Concentration contours')

def init():
    ax.set_xlim([-10*Ro, 10*Ro])
    ax.set_ylim([-10*Ro, 10*Ro])
    return []

def animate(t): 
    ax.clear()
    ax.set_aspect('equal')
    ax.set_xlabel('Distance from center (mm)')
    ax.set_ylabel('Distance from center (mm)')
    ax.set_title('Concentration contours')
    ax.set_xlim([-10*Ro, 10*Ro])
    ax.set_ylim([-10*Ro, 10*Ro])
    
    R = ((Ro**2) - (2 * k*t))**0.5 
    if R>0:
      r=np.linspace(R,10*Ro,500)
      theta=np.linspace(0,2*np.pi,500)
      r_grid,theta_grid=np.meshgrid(r,theta)
      r_star=r_grid/R
      c_star=((1-H)/r_star)+ H
      x=r_grid*np.sin(theta_grid)
      y=r_grid*np.cos(theta_grid)
      cp = ax.contourf(x,y,c_star,levels=100,cmap='rainbow') 
    if R==0:
       r_grid,theta_grid=np.meshgrid(np.linspace(0,10,500),np.linspace(0,2*np.pi,500))
       c_star=H*np.ones_like(r_grid)##makes sure that c_star always remains 2D aRRAY ##
       x=r_grid*np.sin(theta_grid)
       y=r_grid*np.cos(theta_grid)
       cp = ax.contourf(x,y,c_star,levels=100,cmap='rainbow')
    if animate.first_call:
       fig.colorbar(cp)
       animate.first_call=False
    return cp,
animate.first_call=True
ani=FuncAnimation(fig,animate,init_func=init,frames=np.linspace(0,tmax,200),interval=20,blit=False)
ani.save('Droplet_Modified.gif',writer='pillow')
plt.show()