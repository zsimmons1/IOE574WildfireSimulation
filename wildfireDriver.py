# Include libraries
import pandas
import os 
import math
import copy
import csv
import rasterio
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import numpy as np
from rasterio.plot import show
from array import array
from matplotlib.colors import *
from wildfireHelpers import *
from runOneRep import *

# Create Policy class
class Policy:
    def __init__(self, fireLineShape, concentric, spokes, name):
        self.fireLineShape = fireLineShape # a string variable indicating whether the fire lines will be rectangular or cirucular
        self.concentric = concentric # a boolean variable indicating if we will build a proactive concentric contingency line
        self.spokes = spokes # a boolean variable indicating if we will build spokes for the contingency lines
        self.name = name # a string with the letter name of the policy

# Establish policies
policyA = Policy("rectangle", False, False, "policyA_0.15")
policyB = Policy("rectangle", True, False, "policyB_0.15")
policyC = Policy("rectangle", True, True, "policyC_0.15")
policyD = Policy("circle", False, False, "policyD_0.15")
policyE = Policy("circle", True, False, "policyE_0.15")
policyF = Policy("circle", True, True, "policyF_0.15")
policies = [policyA, policyB, policyC, policyD, policyE, policyF]

# Create a header line in each csv file
for p in range(6):
    with open(policies[p].name + ".csv", "a") as f_out:
        writer = csv.writer(f_out, lineterminator = "\n")
        writer.writerow(["AreaBurned", "BurnTime", "LinesBuilt", "AvgWind", "numBreaches", "numResponses"])

# Initialize simulation tracking variables
N = 100 # N is the number of replications to run
# totBurnArea is the total area burned before fire containment
totBurnArea = []
# burnTime is the total burn time (in hours) of the fire before containment
burnTime = []
# linesBuilt is the number of cells of fire lines (of any type) built during the simulation run
totLinesBuilt = []

# Read data and initialize constant inputs
# img holds the vegetation raster data from the LANDFIRE database
img = rasterio.open('/Users/sprin/OneDrive/Desktop/IOE574/TermProject/IOE574WildfireSimulation/finalVegetationData.tif')
# 'map' holds original vegetation raster data from TIFF file
map = img.read()
# crop map to area of interest
map = map[:, 200:400, 150:350]
# 'cumulativeFire' is a 3D matrix containing the number of simulations in which each cell was on fire
cumulativeFire = np.zeros((6, np.size(map, 1), np.size(map, 2)), dtype=float)

# Initiative other input constants
starti = 100 # 'starti' is the i index (corresponding to latitude) where the fire initiates (100 for final)
startj = 100 # 'startj' is the j index (corresponding to longitude) where the fire initiates (100 for final)
del_t = 0.5 # in hours, the time step between updates of the fire status
breachProb = 0.15 # the probability that the fire jumps any given fire line

# Establish sensitivity analysis variables
responseTime = 5 # the number of hours before proctive lines are planned/built
responseRadius = 5 # the number of cells away from the breach where the response lines are built
primaryBuffer = 5 # the number of cells away from the active fire border where the primary lines are built
contingencyBuffer = primaryBuffer + 3

# Run N replications of each policy type
for n in range(N):
    # Generate all random numbers needed for replicaiton n
    initWindSpeed, initWindDir, breachProbs, cellTransitionProbs, allRates, speedNoise, dirNoise = getRVs()

    # Run one replication for each policy
    for p in range(6):
        totAreaBurned, burnTime, totLinesBuilt, avgWind, numBreached, responseLineCounter, cumulativeFire[p] = runOneRep(n, responseTime, policies[p].fireLineShape, responseRadius, primaryBuffer, policies[p].concentric, contingencyBuffer, policies[p].spokes, totBurnArea, burnTime, totLinesBuilt, map, cumulativeFire[p], starti, startj, del_t, breachProb, initWindSpeed, initWindDir, breachProbs, cellTransitionProbs, allRates, speedNoise, dirNoise, policies[p].name)
        with open(policies[p].name + ".csv", "a") as f_out:
            writer = csv.writer(f_out, lineterminator = "\n")
            writer.writerow([totAreaBurned, burnTime, totLinesBuilt, avgWind, numBreached, responseLineCounter])
# Finish
# Visualize and output results for each policy
showMaps(map, cumulativeFire)
print("Simulation complete")
