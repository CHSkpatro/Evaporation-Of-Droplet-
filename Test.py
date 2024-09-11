import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Define the contact radius R and height h0 of the droplet
R = 1.0  # Contact radius (can be adjusted as needed)
h0 = 1.0  # Maximum height of the droplet (can be adjusted as needed)
theta = 2 * np.arctan(h0 / R)  # Contact angle (in radians)
print(theta)