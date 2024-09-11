import matplotlib.pyplot as mpt
import numpy as np
from matplotlib.animation import FuncAnimation
fig, ax =mpt.subplots()
#ax.set_xlim(0,10)
ax.set_ylim(-1,1)
line,=ax.plot([],[],lw=2)
#####initialization of animation####
def init():
    line.set_data([],[])
    return (line,)

def animate(t):
    x=np.linspace(0,10,100)+ 3*t+ 4
    y=np.sin( x )
    line.set_data(x,y)
    ax.set_xlim(x[0],x[-1])
    return (line,)

ani=FuncAnimation(fig=fig,func=animate,init_func=init,frames=np.linspace(0,10,200),interval=50,blit=True)
mpt.show()

