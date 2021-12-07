# Include libraries
import pandas
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
numPolicies = 2 # numPolicies is the total number of policies being compared
N = 5 # N is the number of replications to run
# totBurnArea is the total area burned before fire containment
totBurnArea = []*numPolicies
# burnTime is the total burn time (in hours) of the fire before containment
burnTime = []*numPolicies
# linesBuilt is the number of cells of fire lines (of any type) built during the simulation run
totLinesBuilt = []*numPolicies

# Read data and initialize constant inputs
# img holds the vegetation raster data from the LANDFIRE database
img = rasterio.open('/Users/Zack/Desktop/IOE574/TermProject/IOE574WildfireSimulation/finalVegetationData.tif')
# 'map' holds original vegetation raster data from TIFF file
map = img.read()
map = map[:, 200:400, 150:350]
# 'cumulativeFire' is a 2D matrix containing the number of simulations in which each cell was on fire
cumulativeFire = np.zeros((np.size(map, 1), np.size(map, 2)), dtype=float)

# Initiative other input constants
starti = 100 # 'starti' is the i index (corresponding to latitude) where the fire initiates
startj = 100 # 'startj' is the j index (corresponding to longitude) where the fire initiates
del_t = 0.5 # in hours, the time step between updates of the fire status
breachProb = 0.05 # the probability that the fire jumps any given fire line

# Initialize wind speed and wind direction for all replications
windSpeeds = []
windDirs = []
for n in range(N):
    windSpeeds.append(initializeWindSpeed())
    windDirs.append(initializeWindDir()) 
print("windSpeeds and windDirs initialized")

# Initialize enough breach probabilities for each replications
breachProbs = []
for i in range(1000):
    breachProbs.append(np.random.uniform())
print("breachProbs initialized")

# Initialize enough cellTransitionProbs for each replication  
cellTransitionProbs = []
for i in range(99999999):
    cellTransitionProbs.append(np.random.uniform())
print("cellTransitionProbs initialized")

# Initialize enough spreadRates for each vegetation type for each replication
shape = [11.4, 13.6, 13.0, 0] 
# trees
treeRates = []
for i in range(100000):
    treeRates.append(np.random.weibull(11.4))
# shrubs
shrubRates = []
for i in range(100000):
    shrubRates.append(np.random.weibull(13.6)) 
# herbs
herbRates = []
for i in range(100000):
    herbRates.append(np.random.weibull(13.0))        
allRates = [treeRates, shrubRates, herbRates]
print("allRates initialized")

# Initialize enough speedNoise and dirNoise for each replication
speedNoise = []
dirNoise = []
for i in range(1000):
    speedNoise.append(np.random.uniform(0.8, 1.2))
    dirNoise.append(np.random.uniform(-11.25, 11.25))

# Run N replications of each policy type
for n in range(N):
    # Establish policies
    responseTime = 9 # the number of hours before proctive lines are planned/built
    fireLineShape = "circle" # a string variable indicating whether the fire lines will be rectangular or cirucular
    responseRadius = 6 # the number of cells away from the breach where the response lines are built
    primaryBuffer = 12 # the number of cells away from the active fire border where the primary lines are built
    concentricContingency = True # a boolean variable indicating if we will build a proactive concentric contingency line
    contingencyBuffer = primaryBuffer + 3
    spokes = True # a boolean variable indicating if we will build spokes for the contingency lines
    # run replication
    totBurnArea, burnTime, linesBuilt, cumulativeFire = runOneRep(n, responseTime, fireLineShape, responseRadius, primaryBuffer, concentricContingency, contingencyBuffer, spokes, totBurnArea, burnTime, totLinesBuilt, map, cumulativeFire, starti, startj, del_t, breachProb, windSpeeds, windDirs, breachProbs, cellTransitionProbs, allRates, speedNoise, dirNoise)
    
# Finish
# Visualize and output results for each policy
showResults(totBurnArea, burnTime, totLinesBuilt, cumulativeFire, map, N)



