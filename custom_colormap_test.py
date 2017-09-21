# This is a test script that aim to bring a better understanding of the colormap in MatPlotLib to the world. #
# Nobody reads the headers of python files anyway. Cunt.
# Boris
# the 21th day of September of 2017

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import random


# Parameters
discrete_values = range(10)
colors = []

# Assigning the random color values
for i in range(len(discrete_values)):
	colorTemp = "#"
	for j in range(6):
		colorTemp = colorTemp + (random.choice('0123456789ABCDEF'))
	colors.append(colorTemp)

print colors

#Creating the colormaps
cm = LinearSegmentedColormap.from_list("LithoColorMap", colors, N=len(discrete_values)-1)

#Test set
imtest = [[0,1,2,3,4,5,6,7,8,9],[9,8,7,6,5,4,3,2,1,0],[5,4,6,7,9,8,3,1,2,0],[0,0,0,0,0,0,0,0,0,0]]

plt.imshow(imtest, interpolation = "none", cmap = cm)
plt.colorbar()

plt.show()