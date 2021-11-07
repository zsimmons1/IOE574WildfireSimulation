# This can be used as the main file where we simulate using cellular automata

# DESCRIPTION:

import os 
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import numpy as np
import rasterio
from rasterio.plot import show
import math
from astropy.visualization import make_lupton_rgb

shape = [11.4, 13.6, 13.0]

def ignite(fire, i, j):
    for a in range (-1, 2):
        # ensures index doesn't go out of bounds 
        if i+a > np.size(fire,0) - 1:
            a = a-1
        for b in range (-1, 2):
            # ensures index doesn't go out of bounds 
            if j+b > np.size(fire,1) - 1:
                b = b-1
            if fire[i+a][j+b] == 1:
                 return True
    return False

def calcNewBurn(fire, veg, i, j, del_t):
    # if it is already burning, update burn rate and % burn, or check to ignite new fire
    if fire[i][j] > 0 or ignite(fire, i, j) == True:
        burnRate = math.pow(10, np.random.weibull(shape[veg[i][j]-1]))
        burn = fire[i][j] + burnRate*del_t / 30
        if burn > 1:
            burn = 1
        return burn
    return 0    


img = rasterio.open('/Users/sprin/OneDrive/Desktop/IOE574/TermProject/Data/Other/us_210evc.tif')
map = img.read()

veg = np.floor_divide(map[0], np.ones([np.size(map, 1), np.size(map, 2)], dtype=int)*100)
den = np.mod(map[0], np.ones((np.size(map, 1), np.size(map, 2)), dtype=int)*100)
fire = np.zeros((np.size(map, 1), np.size(map, 2)), dtype=float)
burnRate = fire

starti = 28
startj = 24
del_t = 0.5 # in hours, what is the time step

fire[starti][startj] = 1 # Start fire at location 28, 24 with cell completely on fire
fireBorder = [starti, starti, startj, startj] # [x_lower_bound, x_upper_bound, y_lower_bound, y_upper_bound]

t = 0
min_i = starti
max_i = starti
min_j = startj
max_j = startj

while t < 25:
    t += del_t
    tempFire = fire
    for i in range(fireBorder[0] - 1, fireBorder[1] + 2):
        # ensure cell is in bounds
        if i < np.size(fire, 0):
            for j in range(fireBorder[2] - 1, fireBorder[3] + 2):
                # ensure cell is in bounds:
                if j < np.size(fire, 1):
                    # if the cell is unburnable, skip it
                    if veg[i][j]!=0:
                        tempFire[i][j] = calcNewBurn(fire, veg, i, j, del_t)
                        # update outer fire border
                        if tempFire[i][j] >= 0:
                            if i < fireBorder[0]:
                                min_i = i
                            if i > fireBorder[1]:
                                max_i = i
                            if j < fireBorder[2]:
                                min_j = j
                            if j > fireBorder[3]:
                                max_j = j
                    
    fire = tempFire
    fireBorder = [min_i, max_i, min_j, max_j]

# SHOW RESULTS

# Build pre-burn landscape image
rgbIMG = np.zeros([np.size(map, 1), np.size(map, 2), 3], dtype=int)
r = np.add(np.add(np.where(veg == 1, 56, 0), np.where(veg == 2, 147, 0)), np.where(veg == 3, 219, 0))
g = np.add(np.add(np.where(veg == 1, 118, 0), np.where(veg == 2, 196, 0)), np.where(veg == 3, 235, 0))
b = np.add(np.add(np.where(veg == 1, 29, 0), np.where(veg == 2, 125, 0)), np.where(veg == 3, 118, 0))
rgbIMG = np.dstack((r, g, b))

# Overlay burn
new_r = np.where(fire > 0, 255, r)
new_g = np.where(fire > 0, 0, g)
new_b = np.where(fire > 0, 0, b)
newRGB = np.dstack((new_r, new_g, new_b))

# Create custom legends
ax2legendElements = [lines.Line2D([0], [0], marker='o', color='w', label='Fire', markerfacecolor='#ff0000', markersize=10),
                     lines.Line2D([0], [0], marker='o', color='w', label='Tree', markerfacecolor='#38761d', markersize=10), 
                     lines.Line2D([0], [0], marker='o', color='w', label='Shrub', markerfacecolor='#93c47d', markersize=10), 
                     lines.Line2D([0], [0], marker='o', color='w', label='Herb', markerfacecolor='#dbeb76', markersize=10)]

# Plot figures
fig, ((ax1, ax2)) = plt.subplots(1, 2)
ax1.imshow(rgbIMG)
ax2.imshow(newRGB)
ax2.legend(handles = ax2legendElements, bbox_to_anchor=(1.5, 1.0), loc='upper right')
plt.show()

show(fire, cmap='Reds')

