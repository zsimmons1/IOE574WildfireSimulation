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
# linesBuilt is the number of cells of fire lines (of any type) built during the simulation run
linesBuilt = []
# linesEngaged is the number of cells of fire lines which engaged with the fire at any point during the simulation run
linesEngaged = []

# Read data and initialize constant inputs
# img holds the vegetation raster data from the LANDFIRE database
img = rasterio.open('/Users/Zack/Desktop/IOE574/TermProject/IOE574WildfireSimulation/us_210evc.tif')
# 'map' holds original vegetation raster data from TIFF file
map = img.read()
# 'veg' is a 2D matrix containing the cell vegetation type as an integer (0 = unburnable, 1 = trees, 2 = shrub, 3 = herb, 4 = fire border)
veg = np.floor_divide(map[0], np.ones([np.size(map, 1), np.size(map, 2)], dtype=int)*100)
# 'den' holds the densitity of the vegetation type as a number between 0 and 100
# TODO: 'den' is currently unused
den = np.mod(map[0], np.ones((np.size(map, 1), np.size(map, 2)), dtype=int)*100)

# Run one replication
for n in range(N):
    # Initialize replication-specific variables
    # 'fire' is a 2D matrix containing the percentage on fire of each cell on the map
    fire = np.zeros((np.size(map, 1), np.size(map, 2)), dtype=float)
    # contained is a 2D matrix containing the status of the cell's containment (1 if is is within a fire-line, 0 otherwise)
    contained = np.zeros((np.size(map, 1), np.size(map, 2)), dtype=float)
    # igniteTime is a 2D matrix containing the time (in hours) during which each cell ignites
    igniteTime = np.zeros((np.size(map, 1), np.size(map, 2)), dtype=float)
    # fullBurnTime is a 2D matrix containing the time (in hours) during which the cell reaches 100% on fire
    fullburnTime = np.zeros((np.size(map, 1), np.size(map, 2)), dtype=float)
    # 'distance' is a 3D matrix containing the distance of fire spread for each cell from each neighboring cell
    distance = np.zeros((np.size(map, 1), np.size(map, 2), 8), dtype=float)
    # 'fire_timeline' is a 3D matrix containing the burn state of each cell at each time step in the simulation
    fire_timeline = np.zeros((np.size(map, 1), np.size(map, 2), 20), dtype=float)
    starti = 28 # 'starti' is the i index (corresponding to latitude) where the fire initiates
    startj = 24 # 'startj' is the j index (corresponding to longitude) where the fire initiates
    t = 0 # time elapsed, in hours
    del_t = 0.5 # in hours, the time step between updates of the fire status
    breachProb = 0.05 # the probability that the fire jumps any given fire line
    primaryBuffer = 4 # the number of cells away from the active fire border where the primary lines are built
    responseRadius = 3 # the number of cells away from the breach where the fire line will be built
    responseTime = 2 # the number of hours before initial contingency lines are built
    concentricContingency = False # a boolean variable indicating if we will build a proactive concentric contingency line
    fireLineShape = "rectangle" # a string variable indicating whether the fire lines will be rectangular or cirucular
    contingencyBuffer = primaryBuffer + 3

    # Ignite the fire
    fire[starti][startj] = 1 # Start fire at location 28, 24 with cell 100% on fire
    fireBorder = [starti, starti, startj, startj] # [x_lower_bound, x_upper_bound, y_lower_bound, y_upper_bound]

    # Initialize wind speed and wind direction
    wind_speed = initializeWind()
    wind_direction = 0 # TODO: Sample wind_direction from data too?

    # Spread fire and build fire lines until fire is 95% contained OR fire spreads beyond map borders(?)
    count = -1 # 'count' tracks the number of time steps though the while loop
    zeroSpread = False # Initially set 'zeroSpread' to False so we ititiate the while loop
    borderReached = False # Initially set 'borderReached' to False so we initiate the while loop
    while (zeroSpread == False and borderReached == False):
        
        count += 1 # increment the count
        t += del_t # increment time
        tempFire = copy.deepcopy(fire) # 'tempFire' is a temporary fire matrix to store new % of fire info. Use of this temp
            # matrix ensures that all fire spread depends on the state of spread in the previous time step, 
            # not the current time step (thus preventing spreading too quickly).
        tempFireBorder = copy.deepcopy(fireBorder) # 'tempFireBorder' stores the edges of the fire for the current time step.
            # Use of this temp vector ensures that the nested for loops do not reach cells outside of the border
            # for the previous time step's fire (thus preventing spreading too quickly).    
        
        # Draw all primary fire lines as soon as t == responseTime
        if t == responseTime: 
            buildProactiveLines(i, j, contained, veg, primaryBuffer, breachProb, tempFireBorder, fireLineShape)
            # This will be used as either the radius or the new upper/lower bounds for the contingency lines
            if concentricContingency:
                buildProactiveLines(i, j, contained, veg, contingencyBuffer, breachProb, tempFireBorder, fireLineShape)
            # else:
            #     # To build the spoke contingency lines we can reuse the primary lines but pass in a larger buffer value
            #     buildPrimaryLines(i, j, contained, veg, contingencyBuffer, breachProb, tempFireBorder)
            show(contained, cmap='Blues')
        # traverse the rectangular border around the fire edge plus one cell on each side
        for i in range(fireBorder[0] - 1, fireBorder[1] + 2): # i is latitutde index
            if i < np.size(fire, 0)-1: # ensures cell is in latitute bounds of map
                for j in range(fireBorder[2] - 1, fireBorder[3] + 2): # j is longitude index
                    if j < np.size(fire, 1)-1: # ensure cell is in longitude bounds of map
                        if (veg[i][j]!=0) and (fire[i][j]!=1 and veg[i][j]!=4): # skip all unburnable cells or all completely on fire cells                        
                            # determine which neighbors to spread fire from (see 'ignite' helper function)  
                            cell_transition, spread_prob  = igniteCell(fire, i, j, distance[i][j], wind_speed, wind_direction)
                            # determine if the cell is newly ignited and a fire line breach
                            if (fire[i][j] == 0) and (np.sum(cell_transition) > 0 and contained[i][j]== -1):
                                # if so, build a response line!
                                buildResponseLine(i,j, contained, veg, responseRadius, breachProb, fireLineShape)    
                            # determine the amount of spread from each neighbor which has ignitied cell i,j
                                # (see 'advanceBurn' helper function)
                            distance[i][j] = advanceBurn(veg[i][j], cell_transition, distance[i][j], del_t)
                            # calculate the total percentage on fire for cell i,j based on contributions from
                                # all neighbors (see 'spreadFire' helper function)
                            tempFire[i][j] = spreadFire(distance[i][j], spread_prob)
                            # update 'tempFireBorder' to be the fire border of the current time step
                            if tempFire[i][j] > 0:
                                if i < fireBorder[0]: tempFireBorder[0] = i
                                elif i > fireBorder[1]: tempFireBorder[1] = i
                                if j < fireBorder[2]: tempFireBorder[2] = j
                                elif j > fireBorder[3]: tempFireBorder[3] = j 
                    else: borderReached = True
            else: borderReached = True

        zeroSpread = np.sum(fire) == np.sum(tempFire) # there is zero spread of fire when the total percent
            # on fire of all cells is equal between the previous time step and current time step
        fire = copy.deepcopy(tempFire) # update fire matrix 
        fire_timeline = np.insert(fire_timeline, count, tempFire, axis=2) # store fire status into timeline matrix
        fireBorder = copy.deepcopy(tempFireBorder) # update fire border
        wind_speed, wind_direction = updateWind(wind_speed, wind_direction) # update the wind speed and direction across all cells based for next time step based on current time step

    # show results in map
    showResults(fire, veg)

    # Store simulation results
    # Total area burned (m^2)
    totBurnArea.append(np.sum(fire) * 900)
    print(totBurnArea)

    # Total burn time 
    burnTime.append(t)
    print(burnTime)
    

# Finish
print("Simulation Complete")    