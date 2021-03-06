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
# compiledResults is a vector of vectors of results for each policy type
compiledReults = []

# Read data and initialize constant inputs
# img holds the vegetation raster data from the LANDFIRE database
img = rasterio.open('/Users/sprin/OneDrive/Desktop/IOE574/TermProject/IOE574WildfireSimulation/northCali210EVC.tif')
# 'map' holds original vegetation raster data from TIFF file
map = img.read()
# crop map
map = map[:, 1000:1150, 1000:1200]
# 'cumulativeFire' is a 2D matrix containing the number of simulations in which each cell was on fire
cumulativeFire = np.zeros((np.size(map, 1), np.size(map, 2)), dtype=float)

veg = np.floor_divide(map[0], np.ones([np.size(map, 1), np.size(map, 2)], dtype=int)*100)
showOneRep(cumulativeFire, veg)

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
for i in range(9999999):
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
# Visualize and output results for each policy
showResults(totBurnArea, burnTime, totLinesBuilt, cumulativeFire, map, N)
