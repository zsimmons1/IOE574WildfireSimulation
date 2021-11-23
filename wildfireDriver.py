# Include libraries
import os 
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import numpy as np
import rasterio
from rasterio.plot import show
from array import array
import math
from wildfireHelpers import *
import copy

# Initialize simulation tracking variables
# N is the number of replications to run
N = 1
# totBurnArea is the total area burned before fire containment
totBurnArea = []
# burnTime is the total burn time (in hours) of the fire before containment
burnTime = []
# linesJumped is the total number of fire lines jumped
linesJumped = []

# Read data and initialize constant inputs
# img holds the vegetation raster data from the LANDFIRE database
img = rasterio.open('/Users/sprin/OneDrive/Desktop/IOE574/TermProject/IOE574WildfireSimulation/us_210evc.tif')
# 'map' holds original vegetation raster data from TIFF file
map = img.read()
# 'veg' is a 2D matrix containing the cell vegetation type as an integer (0 = unburnable, 1 = trees, 2 = shrub, 3 = herb, 4 = fire boarder)
veg = np.floor_divide(map[0], np.ones([np.size(map, 1), np.size(map, 2)], dtype=int)*100)
# 'den' holds the densitity of the vegetation type as a number between 0 and 100
# TODO: den is currently unused
den = np.mod(map[0], np.ones((np.size(map, 1), np.size(map, 2)), dtype=int)*100)

# Run one replication
for n in range(N):
    # 'fire' is a 2D matrix containing the percentage on fire of each cell on the map
    fire = np.zeros((np.size(map, 1), np.size(map, 2)), dtype=float)
    # 'distance' is a 3D matrix containing the distance of fire spread for each cell from each neighboring cell
    distance = np.zeros((np.size(map, 1), np.size(map, 2), 8), dtype=float)
    # 'fire_timeline' is a 3D matrix containing the burn state of each cell at each time step in the simulation
    fire_timeline = np.zeros((np.size(map, 1), np.size(map, 2), 20), dtype=float)
    starti = 28 # 'starti' is the i index (corresponding to latitude) where the fire initiates
    startj = 24 # 'startj' is the j index (corresponding to longitude) where the fire initiates
    t = 0 # time elapsed, in hours
    del_t = 0.5 # in hours, the time step between updates of the fire status

    # Ignite the fire
    fire[starti][startj] = 1 # Start fire at location 28, 24 with cell 100% on fire
    fireBorder = [starti, starti, startj, startj] # [x_lower_bound, x_upper_bound, y_lower_bound, y_upper_bound]

    # Initialize wind speed and wind direction
    wind_speed = initializeWind()
    wind_direction = 0 # TODO: Sample wind_direction from data too?

    # Spread fire and build fire lines until fire is 95% contained OR fire spreads beyond map borders(?)
    # TODO: Adjust to incorporate fire lines and 95% containment
    count = -1
    while t < 15:
        count += 1
        t += del_t # increment time
        tempFire = copy.deepcopy(fire) # 'tempFire' is a temporary fire matrix to store new % of fire info. Use of this temp
            # matrix ensures that all fire spread depends on the state of spread in the previous time step, 
            # not the current time step (thus preventing spreading too quickly).
        tempFireBorder = copy.deepcopy(fireBorder) # 'tempFireBorder' stores the edges of the fire for the current time step.
            # Use of this temp vector ensures that the nested for loops do not reach cells outside of the border
            # for the previous time step's fire (thus preventing spreading too quickly).    
        
        # Drawing boarder for fire lines at time = 2 hours, this is done as soon as t = 2
        if t == 2:
            # Setting the x,y lower and upper bounds for the fire line boarder
            x_lower = tempFireBorder[0] - 5
            x_upper = tempFireBorder[1] + 5
            y_lower = tempFireBorder[2] - 5
            y_upper = tempFireBorder[3] + 5
            # Setting the vegitation type to 4 for the fire lines
            for i in range(x_lower, x_upper+1):
                canJump = jumpFireLine()
                if canJump == False:
                    veg[i][y_lower] = 4 # Representative of fire boarder
                    veg[i][y_upper] = 4
            for j in range(y_lower, y_upper+1):
                canJump = jumpFireLine()
                if canJump == False:
                    veg[x_lower][j] = 4
                    veg[x_upper][j] = 4
            
        # traverse the rectangular border around the fire edge plus one cell on each side
        for i in range(fireBorder[0] - 1, fireBorder[1] + 2): # i is latitutde index
            if i < np.size(fire, 0)-1: # ensures cell is in latitute bounds of map
                for j in range(fireBorder[2] - 1, fireBorder[3] + 2): # j is longitude index
                    if j < np.size(fire, 1)-1: # ensure cell is in longitude bounds of map
                        if (veg[i][j]!=0) and (fire[i][j]!=1 and veg[i][j]!=4): # skip all unburnable cells or all completely on fire cells                        
                            # determine which neighbors to spread fire from (see 'ignite' helper function)  
                            cell_transition, spread_prob  = igniteCell(fire, i, j, distance[i][j], wind_speed, wind_direction)
                            # determine the amount of spread from each neighbor which has ignitied cell i,j
                                # (see 'advanceBurn' helper function)
                            distance[i][j] = advanceBurn(veg[i][j], cell_transition, distance[i][j], del_t)
                            # calculate the total percentage on fire for cell i,j based on contributions from
                                # all neighbors (see 'spreadFire' helper function)
                            # spread_prob = [1, 1, 1, 1, 1, 1, 1, 1]    
                            # print("before spreadFire tempFire[",i,"][",j,"]: ", tempFire[i][j])
                            tempFire[i][j] = spreadFire(distance[i][j], spread_prob)
                            # print("after spreadFire tempFire[",i,"][",j,"]: ", tempFire[i][j])
                            # update 'tempFireBorder' to be the fire border of the current time step
                            if tempFire[i][j] > 0:
                                if i < fireBorder[0]: 
                                    tempFireBorder[0] = i
                                elif i > fireBorder[1]: 
                                    tempFireBorder[1] = i
                                if j < fireBorder[2]: 
                                    tempFireBorder[2] = j
                                elif j > fireBorder[3]: 
                                    tempFireBorder[3] = j 
                            if count > 0:
                                if ((fire_timeline[i][j][count-1] > 0) & (tempFire[i][j]==0)):   
                                    print("ERROR! Reseting to zero")                
        # update the wind speed and direction across all cells based for next time step based on current time step
            # wind speed and direction
        wind_speed, wind_direction = updateWind(wind_speed, wind_direction)
        # showResults(fire, veg)
        # showResults(tempFire, veg)
        fire = copy.deepcopy(tempFire) # update fire matrix 
        fireBorder = copy.deepcopy(tempFireBorder) # update fire border
        # fire_timeline[:][:][count] = tempFire # store fire status into timeline matrix
        fire_timeline = np.insert(fire_timeline, count, tempFire, axis=2)
        # showResults(fire, veg)

    # show results in map
    showResults(fire, veg)

    # TODO: store simulation results
    # totBurnArea
    # burnTime
    # linesJumped

# Finish
print("done")    