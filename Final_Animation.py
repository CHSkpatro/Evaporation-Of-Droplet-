import matplotlib.pyplot as plt
import numpy as np
from  matplotlib.animation import FuncAnimation
Ro=float(input('Initial value of radius in mm= '))
H=float(input('Relative Humidity= '))
k = (1 - H) * 6.06 * 10**-4
tmax = (Ro**2) / (2 * k)
fig,ax=plt.subplots()
ax.set_aspect('equal')
ax.set_xlabel('X')
ax.set_ylabel('Y')
r=np.linspace(Ro,10*Ro,500)
theta=np.linspace(0,np.pi*2,500)
r_grid,theta_grid=np.meshgrid(r,theta)
r_star=r_grid/Ro
c_star=((1-H)/r_star)+ H
x=r_star*np.sin(theta_grid)
y=r_star*np.cos(theta_grid)
cp=plt.contourf(x,y,c_star,levels=100,cmap='rainbow')

def init():
    return cp,

def animate(t):
    global cp
    ax.clear()  
    R=((Ro**2) - (2 * k*t))**0.5 
    r1=np.linspace(R,10*R,500)
    theta1=np.linspace(0,2*np.pi,500)
    r1_grid,theta1_grid=np.meshgrid(r1,theta1)
    r1_star=r1_grid/(R)
    c1_star=((1-H)/r1_star)+ H
    x1=r1_star*np.sin(theta1_grid)
    y1=r_star*np.cos(theta1_grid)
    cp = plt.contourf(x1,y1,c1_star,levels=100,cmap='rainbow') 
    return cp,

ani=FuncAnimation(fig,animate,init_func=init,frames=np.linspace(0,tmax,200),interval=20,blit=True)
plt.show()

