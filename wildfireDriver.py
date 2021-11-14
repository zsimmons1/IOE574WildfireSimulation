# Include libraries
import os 
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import numpy as np
import rasterio
from rasterio.plot import show
import math
from astropy.visualization import make_lupton_rgb

# Read data and initialize variables
img = rasterio.open('/Users/Zack/Desktop/IOE574/TermProject/IOE574WildfireSimulation/us_210evc.tif')
# 'map' holds original vegetation raster data from TIFF file
map = img.read()
# 'veg' is a 2D matrix containing the cell vegetation type as an integer (0 = unburnable, 1 = trees, 2 = shrub, 3 = herb)
veg = np.floor_divide(map[0], np.ones([np.size(map, 1), np.size(map, 2)], dtype=int)*100)
# 'den' holds the densitity of the vegetation type as a number between 0 and 100
# TODO: den is currently unused
den = np.mod(map[0], np.ones((np.size(map, 1), np.size(map, 2)), dtype=int)*100)
# 'fire' is a 2D matrix containing the percentage on fire of each cell on the map
fire = np.zeros((np.size(map, 1), np.size(map, 2)), dtype=float)
# 'distance' is a 3D matrix containing the distance of fire spread for each cell from each neighboring cell
distance = np.zeros((np.size(map, 1), np.size(map, 2), 8), dtype=float)
starti = 28 # 'starti' is the i index (corresponding to latitude) where the fire initiates
startj = 24 # 'startj' is the j index (corresponding to longitude) where the fire initiates
t = 0 # time elapsed, in hours
del_t = 0.5 # in hours, the time step between updates of the fire status

# Ignite the fire
fire[starti][startj] = 1 # Start fire at location 28, 24 with cell 100% on fire
fireBorder = [starti, starti, startj, startj] # [x_lower_bound, x_upper_bound, y_lower_bound, y_upper_bound]

# Spread fire and build fire lines until fire is 95% contained OR fire spreads beyond map borders(?)
# TODO: Adjust to incorporate fire lines and 95% containment
while t > 25:
    t += del_t # increment time
    tempFire = fire # 'tempFire' is a temporary fire matrix to store new % of fire info. Use of this temp
        # matrix ensures that all fire spread depends on the state of spread in the previous time step, 
        # not the current time step (thus preventing spreading too quickly).
    tempFireBorder = fireBorder # 'tempFireBorder' stores the edges of the fire for the current time step.
        # Use of this temp vector ensures that the nested for loops do not reach cells outside of the border
        # for the previous time step's fire (thus preventing spreading too quickly).    
    # traverse the rectangular border around the fire edge plus one cell on each side   
    for i in range(fireBorder[0] - 1, fireBorder[1] + 2): # i is latitutde index
        if i < np.size(fire, 0): # ensures cell is in latitute bounds of map
            for j in range(fireBorder[2] - 1, fireBorder[3] + 2): # j is longitude index
                if j < np.size(fire, 1): # ensure cell is in longitude bounds of map
                    if veg[i][j]!=0: # skip all unburnable cells
                        # TODO: determine wind_speed and wind_direction
                        wind_speed = 15
                        wind_direction = 180
                        # determine which neighbors to spread fire from
                            # (see 'ignite' helper function)
                        cell_transition, spread_prob  = igniteCell(fire, i, j, wind_speed, wind_direction)    
                        # determine the amount of spread from each neighbor which has ignitied cell i,j
                            # (see 'advanceBurn' helper function)
                        distance[i][j] = advanceBurn(fire, veg, spread_prob, distance[i][j], i, j, del_t)
                        # calculate the total percentage on fire for cell i,j based on contributions from
                            # all neighbors (see 'spreadFire' helper function)
                        tempFire[i][j] = spreadFire(distance[i][j], cell_transition[1])
                        # update 'tempFireBorder' to be the fire border of the current time step
                        if tempFire[i][j] >= 0:
                            if i < fireBorder[0]: tempFireBorder[0] = i
                            elif i > fireBorder[1]: tempFireBorder[1] = i
                            if j < fireBorder[2]: tempFireBorder[2] = i
                            elif j > fireBorder[3]: tempFireBorder[3] = i
                    
    fire = tempFire # update fire matrix 
    fireBorder = [min_i, max_i, min_j, max_j] # update fire border

# Show results
showResults(fire, veg)