import matplotlib.pyplot as mpt
import numpy as np
from matplotlib.animation import FuncAnimation
Ro=float(input('Initial radius in mm = '))
H=float(input('Relative Humidity is = '))
k=(1-H)*6.06*10**-4
fig,ax=mpt.subplots()
ax.set_xlim(-10*Ro,10*Ro)
ax.set_ylim(-10*Ro,10*Ro)
ax.set_aspect('equal')
circle=mpt.Circle((0,0),radius=Ro,fill=False)
ax.add_patch(circle)
###### Animation initialization Function #####
def init():
    circle.set_radius(Ro)
    return circle,

########Animation Function######
def animate(t):
    R=(((Ro**2) - 2*t*k)**0.5)
    circle.set_radius(R)
    return circle,
tmax=(Ro**2)/(2*k)
ani=FuncAnimation(fig=fig,func=animate,init_func=init,frames=np.linspace(0,tmax,1000),interval=20,blit=True)
mpt.show()

