# Include libraries
import os 
import matplotlib.pyplot as plt
import matplotlib.lines as lines
from matplotlib.colors import *
import numpy as np
import rasterio
from rasterio.plot import show
from array import array
import math
import copy
from wildfireHelpers import *

def runOneRep(n, responseTime, fireLineShape, responseRadius, primaryBuffer, concentricContingency, contingencyBuffer, spokes, totBurnArea, burnTime, totLinesBuilt, map, cumulativeFire, starti, startj, del_t, breachProb, windSpeeds, windDirs, breachProbs, cellTransitionProbs, allRates):
    # Initialize replication-specific variables
    t = 0 # time elapsed, in hours
    linesBuilt = 0 # the number of fire lines built (of all types)
    int(linesBuilt)
    c = 0 # the number of cellTransitionProbs drawn
    int(c)
    tSpread = 0 # the number of spread rates sampled for trees
    int(tSpread)
    sSpread = 0 # the number of spread rates sampled for shurbs
    int(sSpread)
    hSpread = 0 # the number of spread rates sampled for herbs
    int(hSpread)
    rateCounts = [tSpread, sSpread, hSpread]
    # 'veg' is a 2D matrix containing the cell vegetation type as an integer (0 = unburnable, 1 = trees, 2 = shrub, 3 = herb, 4 = fire border)
    veg = np.floor_divide(map[0], np.ones([np.size(map, 1), np.size(map, 2)], dtype=int)*100)
    # 'fire' is a 2D matrix containing the percentage on fire of each cell on the map
    fire = np.zeros((np.size(map, 1), np.size(map, 2)), dtype=float)
    # contained is a 2D matrix containing the status of the cell's containment (1 if is is within a fire-line, 0 otherwise)
    contained = np.zeros((np.size(map, 1), np.size(map, 2)), dtype=float)
    # 'distance' is a 3D matrix containing the distance of fire spread for each cell from each neighboring cell
    distance = np.zeros((np.size(map, 1), np.size(map, 2), 8), dtype=float)
    
    # Ignite the fire
    fire[starti][startj] = 1 # Start fire at location 28, 24 with cell 100% on fire
    fireBorder = [starti, starti, startj, startj] # [x_lower_bound, x_upper_bound, y_lower_bound, y_upper_bound]

    # Get wind speed and wind direction for this replication
    wind_speed = windSpeeds[n]
    wind_direction = windDirs[n]

    # Spread fire and build fire lines until there is zero fire spread OR fire spreads beyond map borders
    zeroSpread = False # Initially set 'zeroSpread' to False so we ititiate the while loop
    borderReached = False # Initially set 'borderReached' to False so we initiate the while loop
    while (zeroSpread == False and borderReached == False):
        t += del_t # increment time
        tempFire = copy.deepcopy(fire) # 'tempFire' is a temporary fire matrix to store new % of fire info. Use of this temp
            # matrix ensures that all fire spread depends on the state of spread in the previous time step, 
            # not the current time step (thus preventing spreading too quickly).
        tempFireBorder = copy.deepcopy(fireBorder) # 'tempFireBorder' stores the edges of the fire for the current time step.
            # Use of this temp vector ensures that the nested for loops do not reach cells outside of the border
            # for the previous time step's fire (thus preventing spreading too quickly).    
        
        # Draw all proactive fire lines as soon as t == responseTime
        if t == responseTime: 
            linesBuilt = buildProactiveLines(i, j, contained, veg, primaryBuffer, breachProbs, linesBuilt, breachProb, tempFireBorder, fireLineShape, spokes)
            # If the policy calls for contingency lines, draw contigency lines
            if concentricContingency:
                linesBuilt = buildProactiveLines(i, j, contained, veg, contingencyBuffer, linesBuilt, breachProbs, breachProb, tempFireBorder, fireLineShape, spokes)

        # Traverse the rectangular border around the fire edge plus one cell on each side
        for i in range(fireBorder[0] - 1, fireBorder[1] + 2): # i is latitutde index
            if i < np.size(fire, 0)-1: # ensures cell is in latitute bounds of map
                for j in range(fireBorder[2] - 1, fireBorder[3] + 2): # j is longitude index
                    if j < np.size(fire, 1)-1: # ensure cell is in longitude bounds of map
                        if (veg[i][j]!=0) and (fire[i][j]!=1 and veg[i][j]!=4): # skip all unburnable cells or all completely on fire cells                        
                            # determine which neighbors to spread fire from (see 'ignite' helper function)  
                            cell_transition, spread_prob, c  = igniteCell(fire, i, j, distance[i][j], wind_speed, wind_direction, cellTransitionProbs, c)
                            # determine if the cell is newly ignited and a fire line breach
                            if (fire[i][j] == 0) and (np.sum(cell_transition) > 0 and contained[i][j]== -1):
                                # if so, build a response line!
                                linesBuilt = buildResponseLine(i,j, contained, veg, responseRadius, breachProbs, linesBuilt, breachProb, fireLineShape)    
                            # determine the amount of spread from each neighbor which has ignitied cell i,j
                                # (see 'advanceBurn' helper function)
                            distance[i][j], rateCounts = advanceBurn(veg[i][j], cell_transition, distance[i][j], del_t, allRates, rateCounts)
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
        fireBorder = copy.deepcopy(tempFireBorder) # update fire border
        wind_speed, wind_direction = updateWind(wind_speed, wind_direction) # update the wind speed and direction across all cells based for next time step based on current time step
        
    # Store simulation results
    totBurnArea.append((np.sum(fire) * 900)/1000) # Total area burned (m^2 * 10^-3)
    burnTime.append(t) # Total burn time 
    totLinesBuilt.append((linesBuilt * 30)) # Total fire lines built (m)
    cumulativeFire = np.add(cumulativeFire, fire) # Add to cumulative burn
    print("Replication " + str(n+1) + ": " + str(t) + " hours to burn " + str(np.sum(fire) * 900) +" square meters")
    return totBurnArea, burnTime, linesBuilt, cumulativeFire