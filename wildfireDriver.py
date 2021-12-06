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
from runOneRep import *

# Initialize simulation tracking variables
# N is the number of replications to run
N = 5
# totBurnArea is the total area burned before fire containment
totBurnArea = []
# burnTime is the total burn time (in hours) of the fire before containment
burnTime = []
# linesBuilt is the number of cells of fire lines (of any type) built during the simulation run
totLinesBuilt = []

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

# Run N replications of each policy type
for n in range(N):
    # Establish policies
    responseTime = 2 # the number of hours before proctive lines are planned/built
    fireLineShape = "rectangle" # a string variable indicating whether the fire lines will be rectangular or cirucular
    responseRadius = 4 # the number of cells away from the breach where the response lines are built
    primaryBuffer = 4 # the number of cells away from the active fire border where the primary lines are built
    concentricContingency = False # a boolean variable indicating if we will build a proactive concentric contingency line
    contingencyBuffer = primaryBuffer + 3
    spokes = False # a boolean variable indicating if we will build spokes for the contingency lines
    # run replication
    totBurnArea, burnTime, linesBuilt, cumulativeFire = runOneRep(n, responseTime, fireLineShape, responseRadius, primaryBuffer, concentricContingency, contingencyBuffer, spokes, totBurnArea, burnTime, totLinesBuilt, map, cumulativeFire, starti, startj, del_t, breachProb, windSpeeds, windDirs, breachProbs, cellTransitionProbs, allRates)
    
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


