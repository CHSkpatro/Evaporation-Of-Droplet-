import numpy as np
import matplotlib.pyplot as mpt
from mpl_toolkits.mplot3d import Axes3D
fig=mpt.figure()###fig is accessed from mpt.figure()###
x=np.linspace(-1,1,100)
y=np.linspace(-1,1,100)
###next has to form the 2d mesh grid from 1D arrays of x and y###
X,Y=np.meshgrid(x,y)
Z_p=(1 - (X**2 +Y**2))**0.5
Z_l=-(1 - (X**2 +Y**2))**0.5

ax=fig.add_subplot(111,projection='3d')
ax.plot_surface(X,Y,Z_p,cmap='viridis')
ax.plot_surface(X,Y,Z_l,cmap='viridis')
mpt.show()
