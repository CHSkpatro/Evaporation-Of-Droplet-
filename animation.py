import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

fig, ax = plt.subplots()
ax.set_xlim(0, 10)
ax.set_ylim(-1, 1)
line, = ax.plot([], [], lw=2)

# Initialization of animation
def init():
    line.set_data([], [])
    return (line,)

# Animation function
def animate(t):
    x = np.linspace(0, 10, 100)
    y = np.sin(3 * t + 4)
    line.set_data(x, y)
    return (line,)

# Keep a reference to the animation object
ani = FuncAnimation(fig, animate, init_func=init, frames=np.linspace(0, 10, 200), interval=50, blit=True)

plt.show()
