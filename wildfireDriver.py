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

# Initialize simulation tracking variables
# N is the number of replications to run
N = 5
# totBurnArea is the total area burned before fire containment
totBurnArea = []
# burnTime is the total burn time (in hours) of the fire before containment
burnTime = []
# linesBuilt is the number of cells of fire lines (of any type) built during the simulation run
totLinesBuilt = []
# linesEngaged is the number of cells of fire lines which engaged with the fire at any point during the simulation run
linesEngaged = []

# Read data and initialize constant inputs
# img holds the vegetation raster data from the LANDFIRE database
img = rasterio.open('/Users/sprin/OneDrive/Desktop/IOE574/TermProject/IOE574WildfireSimulation/us_210evc.tif')
# 'map' holds original vegetation raster data from TIFF file
map = img.read()
# 'cumulativeFire' is a 2D matrix containing the number of simulations in which each cell was on fire
cumulativeFire = np.zeros((np.size(map, 1), np.size(map, 2)), dtype=float)

# Initiative other input constants
starti = 28 # 'starti' is the i index (corresponding to latitude) where the fire initiates
startj = 24 # 'startj' is the j index (corresponding to longitude) where the fire initiates
del_t = 0.5 # in hours, the time step between updates of the fire status
breachProb = 0.05 # the probability that the fire jumps any given fire line

# Initialize wind speed and wind direction for all replications
windSpeeds = []
windDirs = []
for n in range(N):
    windSpeeds.append(initializeWind())
    windDirs.append(0) # TODO: Sample wind_direction from data too?
print("windSpeeds and windDirs initialized")

# Initialize enough breach probabilities for each replications
breachProbs = []
for i in range(1000):
    breachProbs.append(np.random.uniform())
print("breachProbs initialized")

# Initialize enough cellTransitionProbs for each replication  
cellTransitionProbs = []
for i in range(999999):
    cellTransitionProbs.append(np.random.uniform())
print("cellTransitionProbs initialized")

# Initialize enough spreadRates for each vegetation type for each replication
shape = [11.4, 13.6, 13.0, 0] 
# trees
treeRates = []
for i in range(10000):
    treeRates.append(np.random.weibull(11.4))
# shrubs
shrubRates = []
for i in range(10000):
    shrubRates.append(np.random.weibull(13.6)) 
# herbs
herbRates = []
for i in range(10000):
    herbRates.append(np.random.weibull(13.0))        
allRates = [treeRates, shrubRates, herbRates]
print("allRates initialized")

# Run one replication
for n in range(N):
    # Establish policies
    responseTime = 2 # the number of hours before proctive lines are planned/built
    fireLineShape = "rectangle" # a string variable indicating whether the fire lines will be rectangular or cirucular
    responseRadius = 4 # the number of cells away from the breach where the response lines are built
    primaryBuffer = 4 # the number of cells away from the active fire border where the primary lines are built
    concentricContingency = False # a boolean variable indicating if we will build a proactive concentric contingency line
    contingencyBuffer = primaryBuffer + 3
    spokes = False # a boolean variable indicating if we will build spokes for the contingency lines

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
    showOneRep(fire, veg)
    print("Replication " + str(n+1) + ": " + str(t) + " hours to burn " + str(np.sum(fire) * 900) +" square meters")
    

# Finish
# Generate empirical CDFs
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
y = np.linspace(1/N, 1, num=N)
x1 = np.sort(totBurnArea)
ax1.scatter(x1, y)
ax1.set_title('Total Area Burned (squared meters * 10^-3)')
x2 = np.sort(burnTime)
ax2.scatter(x2, y)
ax2.set_title('Total Burn Time (hours)')
x3 = np.sort(totLinesBuilt)
ax3.scatter(x3, y)
ax3.set_title('Fire Lines Built (meters)')
veg = np.floor_divide(map[0], np.ones([np.size(map, 1), np.size(map, 2)], dtype=int)*100)
rgbIMG = np.zeros([np.size(veg, 0), np.size(veg, 1), 3], dtype=int)
r = np.add(np.add(np.where(veg == 1, 56, 0), np.where(veg == 2, 147, 0)), np.where(veg == 3, 219, 0))
g = np.add(np.add(np.where(veg == 1, 118, 0), np.where(veg == 2, 196, 0)), np.where(veg == 3, 235, 0))
b = np.add(np.add(np.where(veg == 1, 29, 0), np.where(veg == 2, 125, 0)), np.where(veg == 3, 118, 0))
vegRGB = np.dstack((r, g, b))
ax4.imshow(vegRGB)
alphaMat = np.ceil(cumulativeFire / N)
cumulativeFire[cumulativeFire==0]=['nan']
plot4 = ax4.imshow(cumulativeFire, cmap='autumn_r')
cbar = plt.colorbar(plot4, ax = ax4, orientation='vertical')
cbar.set_label('# of replications burned')
ax4.set_title('Burn Map')
plt.show()

print("Burn Area -- Average: " + str(np.average(totBurnArea)) + " thousand square meters")
print("Burn Time -- Average: " + str(np.average(burnTime)) + " hours")
print("Fire Lines Built -- Average: " + str(np.average(totLinesBuilt)) + " meters")  


