# Include libraries
import math
import pandas
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as lines
from rasterio.plot import show
from scipy.stats import weibull_min
from array import array

# getRVs:
def getRVs():
    # Initialize wind speed and wind direction for all replications
    initWindSpeed = initializeWindSpeed()
    initWindDir = initializeWindDir()
    print("windSpeeds and windDirs initialized")

    # Initialize enough breach probabilities for each replications
    breachProbs = []
    for i in range(10000):
        breachProbs.append(np.random.uniform())
    print("breachProbs initialized")

    # Initialize enough cellTransitionProbs for each replication  
    cellTransitionProbs = []
    for i in range(9999999):
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

    return initWindSpeed, initWindDir, breachProbs, cellTransitionProbs, allRates, speedNoise, dirNoise

# initializeWind: Based on input weather data, determine the appropriate Weibull distribution for wind speed
# and draw from this distribution to determine the initial wind speed across all cells
def initializeWindSpeed():
    data = pandas.read_csv("weather.csv")
    messy_wind_speed = data["HourlyWindSpeed"]
    # remove nan values
    wind_speed = [x for x in messy_wind_speed if np.isnan(x) == False]
    # determine the Weibull distrbution curve that fits the weather data
    shape, loc, scale = weibull_min.fit(wind_speed, floc = 1)
    # generate a random wind speeds from a weibull distribution using parameters generated from actual data
    initial_wind = weibull_min.rvs(shape, loc, scale, size=1)
    # return the generated initial wind as a single value (not an array) in km/hr
    return np.random.choice(initial_wind)*3.6

def initializeWindDir():
    data = pandas.read_csv("weather.csv")
    messy_wind_dir = data["HourlyWindDirection"]
    wh = messy_wind_dir.index[messy_wind_dir == 'VRB']
    wind_dir_1 = messy_wind_dir.drop(wh)
    wind_dir_1 = wind_dir_1[~wind_dir_1.isnull()]
    wind_direction = wind_dir_1.astype(int)
    return np.random.choice(wind_direction)
    
# updateWind: Determines the wind direction and velocity for cell i,j
    # 'wind_speed' is the wind speed (in km/hour) across all cells from the previous time period
    # 'wind_direction' is the angle (in degrees, where East is 0 degrees) of the wind across all cells from the 
        # previous time period
def updateWind(wind_speed, wind_direction, sNoise, dNoise):
    # Apply multaplicative uniform noise to the wind speed U[0.8, 1.2] per Trucchia et al
    wind_speed = wind_speed*sNoise
    # Apply multaplicative uniform noise to the wind direction U[???11.25???, 11.25???] per Trucchia et al
    wind_direction = wind_direction + dNoise
    # Ensure the angle is between 0 and 360
    if wind_direction < 0:
        wind_direction += 360
    elif wind_direction > 360:
        wind_direction -= 360    
    # Return the updated wind speed and wind direction
    return wind_speed, wind_direction

# igniteCell: Determines which of the neighbors of cell i,j will contribute to the fire spread of cell i,j 
# during the current timestep
    # 'fire'
    # 'i' is the latitude index of the current cell of interest
    # 'j' is the longitutde index of the current cell of interest
    # 'distance' is the 8-element vector of spread distances (from the last time period) from each neighbor 
        # from the 3-D distance matrix at cell i,j 
    # 'wind_speed' is the wind speed (in km/hour) across all cells for the current time period
    # 'wind_direction' is the angle (in degrees, where East is 0 degrees) of the wind across all cells for the 
        # current time step
def igniteCell(fire, i, j, distance, wind_speed, wind_direction, cellTransitionProbs, c):
    cell_transition = [0, 0, 0, 0, 0, 0, 0, 0]
    spread_angles = [315, 270, 225, 0, 180, 45, 90, 135]
    # Set baseline spread probability for each neighbor as the percent on fire of that neighbor
    spread_prob = [fire[i-1][j-1], fire[i-1][j], fire[i-1][j+1],
                   fire[i][j-1], fire[i][j+1],
                   fire[i+1][j-1], fire[i+1][j], fire[i+1][j+1]]        
    # Determine which neighbor cells contribute to the fire in cell i,j
    for k in range(len(distance)):
        # If there is a spread distance from neighbor k, then neighbor k ignited cell i,j in a previous time step
        if distance[k] > 0:
            cell_transition[k] = 1
        # Else, determine if neighbor k will ignite cell i,j in this time step    
        else: 
            # Adjust 'spread_prob' based on wind speed and direction
            angle = spread_angles[k] - wind_direction
            if angle > 180:
                angle -= 360
            elif angle < -180:
                angle += 360
            # TODO: Create function(s) to determine alpha_w
            # set function based on wind speed
            if wind_speed >= 0 and wind_speed < 10 :
                w = 1
                stdev = 100000000
            elif wind_speed >= 10 and wind_speed < 20 :
                w = 1.2
                stdev = 250
            elif wind_speed >= 20 and wind_speed < 30 :
                w = 1.8
                stdev = 160
            elif wind_speed >= 30 and wind_speed < 50 :
                w = 2.5
                stdev = 70
            elif wind_speed >= 50 and wind_speed < 60 :
                w = 3
                stdev = 80
            elif wind_speed >= 60 and wind_speed < 70 :
                w = 3.2
                stdev = 60
            elif wind_speed >= 70 and wind_speed < 80 :
                w = 3.4
                stdev = 50
            elif wind_speed >= 80 and wind_speed < 90 :
                w = 3.6
                stdev = 38
            elif wind_speed >= 90 and wind_speed < 100 :
                w = 3.75
                stdev = 25
            else :
                w = 4
                stdev = 10
            # input angle into funciton to get alpha_w facor
            alpha_w = w*math.exp(-0.5*(angle**2)/stdev**2)
            spread_prob[k] = spread_prob[k]*alpha_w
            # Ensure all probabilities are no more than 1
            if spread_prob[k] > 1:
                spread_prob[k] = 1
                cell_transition[k] = 1
            else:    
            # Draw from uniform distribution to determine if neighbor k will spread to cell i,j  
            # print("c: ", c)
            # print("size of cellTransitionProbs: ", len(cellTransitionProbs))
                u = cellTransitionProbs[c] 
                c += 1 
                if u < spread_prob[k]:
                    cell_transition[k] = 1          
            # Else, neighbor k does not contribute to the fire in cell i,j, so 'cell_transition[k]' remains 0
    # Return the cell_transition and spread_prob
    return cell_transition, spread_prob, c

# advanceBurn: Determine the distance of spread of the fire from each contributing neighbor to cell i,j. These
# distances are determined from the previous distance and velocity of fire spread sampled from a Weibull distribution
    # 'veg' is a 2D matrix containing the cell vegetation type as an integer 
        # (0 = unburnable, 1 = trees, 2 = shrub, 3 = herb)
    # 'cell_transition' is an 8-element vector for the neighbors of cell i,j where 0 means neighbor k does not
        # contribute to the fire spread in cell i,j during the current time period and 1 means neighbor k does
        # contribute to the fire spread in cell i,j during the current time period
    # 'distance' is the 8-element vector of spread distances (from the last time period) from each neighbor 
        # from the 3-D distance matrix at cell i,j  
    # 'del_t' is the size of each time step
def advanceBurn(veg, cell_transition, distance, del_t, allRates, rateCounts):
    # shape stores the Weibull distribution shape for each vegetation type
    shape = [11.4, 13.6, 13.0, 0]        
    orthogonal = [1, 3, 4, 6]
    diagonal = [0, 2, 5, 7]
    for x in orthogonal:
        # If this neighbor is contributing, then we will calculate the distance
        # diagonal neighbors can have spread distance up to 30 m
        if shape[veg-1] == 0:
            distance[x] = 0
        if cell_transition[x] == 1 and shape[veg-1] != 0:
            del_x = math.pow(10, allRates[veg-1][rateCounts[veg-1]])*del_t
            rateCounts[veg-1] += 1
            if distance[x] + del_x > 30:
                distance[x] = 30
            else: 
                distance[x] += del_x
    for x in diagonal:
        # If this neighbor is contributing, then we will calculate the distance
        if shape[veg-1] == 0:
            distance[x] = 0
        if cell_transition[x] == 1 and shape[veg-1] != 0:
            del_x = math.pow(10, allRates[veg-1][rateCounts[veg-1]])*del_t
            # diagonal neighbors can have spread distance up to 30*square root of 2 m
            if distance[x] + del_x > 30*pow(2, 0.5):
                distance[x] = 30*pow(2, 0.5)
            else: 
                distance[x] += del_x
    return distance, rateCounts

# spreadFire: Calculates the area of burn in cell i,j resulting from the spread distance in each direction and the 
# intensity of spread in each direction. Area of spread is then used to calculate total proportion of 
# the cell which is on fire.
    # 'distance' is an 8-element vector of the distance of spread in each direction, from neighbors in the
        # following order: NW, N, NE, W, E, SW, S, SE
    # 'spread_prob' is an 8-element vector of the spread probability from each neighbor to cell i,j with 
        # neighbor order same as distance vector
def spreadFire(distance, spread_prob):
    orthogonal = [1, 3, 4, 6] # set the indices for orthogonal neighbors (N, W, E, S)
    diagonal = [0, 2, 5, 7] # set the indices for diagonal neighbors (NW, NE, SW, SE)
    burnArea = 0 # initialize the burn area to be 0
    # Add burn contribution from orthogonal neighbors
    for a in orthogonal:
        burnArea += spread_prob[a]*30*distance[a]    
    # Add burn contribution from diagonal neighbors
    for a in diagonal:
        if distance[a] < 15*pow(2, 0.5):
            # Add burn contribution from diagonal neighbor when diagonal distance ??? half
            burnArea += spread_prob[a]*pow(distance[a], 2)
        else:
            # Add burn contribution from diagonal neighbor when diagonal distance ??? half 
            burnArea += spread_prob[a]*(900-pow((30*pow(2, 0.5)-distance[a]), 2))
    # Return the burn percentage 
    if (burnArea > 900): 
        burnArea = 900       
    return burnArea / 900

# breachFireLine: Determines whether any given fire line will be breached based on a uniform random probability
    # 'breachProb' is the probability that a fire line will be breached (independent of any characteristics of
        # the fire)
def breachFireLine(breachProbs, linesBuilt, breachProb, responseLineCounter):
    # print("response line counter: ", responseLineCounter)
    if responseLineCounter > 20:
        return False
    elif responseLineCounter > 10:
        if breachProbs[linesBuilt] > (1-(breachProb/2)):
            return True
        else:
            return False        
    if breachProbs[linesBuilt] > (1-breachProb):
        return True
    else:
        return False

# rectangleLine: Constructs a rectangular fire line in accordance with upper and lower bounds for rows and columns
    # 'veg' is a 2D matrix containing the cell vegetation type as an integer (4 indicates an unbreached fire border)
    # 'fireLineBounds' is a vector of [rL, rU, cL, cU] of the rectangle fire line to be built
    # 'contained' is a 2D matrix containing the status of the cell's containment 
    #   (1 if is is within a fire-line, 0 otherwise)
    # 'breachProb' is the probability that the fire jumps any given fire line
def rectangleLine(veg, fireLineBounds, contained, breachProbs, linesBuilt, breachProb, responseLineCounter, proactive, numBreached):
    # Unpack fire line bounds
    rL = fireLineBounds[0] # row lower bound (northmost)
    rU = fireLineBounds[1] # row upper bound (southmost)
    cL = fireLineBounds[2] # column lower bound (westmost)
    cU = fireLineBounds[3] # column upper bound (eastmost)

    # Set the vegitation type to 4 for the fire lines
    for i in range(rL, rU+1): # Build vertical fire lines
        if contained[i][cL] == 0: # if the potential West fire line is not already contained
            if breachFireLine(breachProbs, linesBuilt, breachProb, responseLineCounter) == False: veg[i][cL] = 4
            else: 
                contained[i][cL] = -1
                numBreached += 1
                if proactive == False: responseLineCounter += 1
            linesBuilt += 1 
        if contained[i][cU] == 0: # if the potential East fire line is not already contained 
            if breachFireLine(breachProbs, linesBuilt, breachProb, responseLineCounter) == False: veg[i][cU] = 4
            else: 
                contained[i][cU] = -1
                numBreached += 1
                if proactive == False:
                    responseLineCounter += 1
            linesBuilt += 1    
    for j in range(cL, cU+1): # Build horizontal fire lines
        if contained[rL][j] == 0: # if the potential North fire line is not already contained
            if breachFireLine(breachProbs, linesBuilt, breachProb, responseLineCounter) == False: veg[rL][j] = 4 
            else: 
                contained[rL][j] = -1
                numBreached += 1
                if proactive == False: responseLineCounter += 1
            linesBuilt += 1 
        if contained[rU][j] == 0: # if the potential South fire line is not already contained   
            if breachFireLine(breachProbs, linesBuilt, breachProb, responseLineCounter) == False: veg[rU][j] = 4
            else: 
                contained[rU][j] = -1
                numBreached += 1
                if proactive == False: responseLineCounter += 1
            linesBuilt += 1       
    # Create the contained zone from the new fire line
    for i in range(rL+1, rU):
        for j in range(cL+1, cU):
            contained[i][j] = 1            
    return linesBuilt, responseLineCounter, numBreached


def distCalculator(x1, y1, x2, y2):
    dist = ((x2 - x1)**(2) + (y2-y1)**(2))**(0.5)
    return dist

def circularLines(veg, fireLineBounds, contained, breachProbs, linesBuilt, breachProb, responseLineCounter, proactive, numBreached):
    # newLinesBuilt = 0
    # Unpack fire line bounds which will help calculate the radius of the circle
    rL = fireLineBounds[0]
    rU = fireLineBounds[1]
    cL = fireLineBounds[2]
    cU = fireLineBounds[3]

    # First we have to find the center point
    cC = (cL + cU) / 2 # Find the column center location
    rC = (rL + rU) / 2 # Find the row center location

    # Now the radius will be the distance from (rU, cU) to (rC, cC)
    rad =  distCalculator(rC, cC, rU, cU)

    for i in range(rL-7, rU+7): # Check a range well outside of the max radius
        for j in range(cL-7, cU+7):
            # Setting the x,y lower and upper bounds for the fire line boarder to ensure they do not exceed the 
            # bounds of the map
            if i < 0: i = 0 # left fire line boundary
            if i > np.size(veg, 0)-1: i = np.size(veg, 0)-1 # right fire line boundary
            if j < 0: j = 0 # lower fire line boundary      
            if j > np.size(veg, 1)-1: j = np.size(veg, 1)-1 # upper fire line boundary

            if contained[i][j] == 0: # if the potential South fire line is not already contained  
                currDist = distCalculator(rC, cC, i, j)
                if currDist - 1.5 < rad:
                    if rad <= currDist:
                        if breachFireLine(breachProbs, linesBuilt, breachProb, responseLineCounter) == False: veg[i][j] = 4
                        else:
                            contained[i][j] = -1  
                            numBreached += 1 
                            if proactive == False: responseLineCounter += 1
                        linesBuilt += 1     
                    else:
                        contained[i][j] = 1        
    return linesBuilt, responseLineCounter, numBreached

def addSpokes(veg, fireLineBounds, contained, breachProbs, linesBuilt, breachProb, fireLineShape, responseLineCounter, numBreached):
    # Unpack fire line bounds which will help calculate the radius of the circle
    rL = fireLineBounds[0]
    rU = fireLineBounds[1]
    cL = fireLineBounds[2]
    cU = fireLineBounds[3]

    # First we have to find the center point
    cC = (cL + cU) / 2 # Find the column center location
    rC = (rL + rU) / 2 # Find the row center location

    rC = int(rC)
    cC = int(cC)

    # Find slope for each
    slope1 = (rU - rC) / (cU - cC)
    slope2 = (rU - rC) / (cL - cC)

    # Find dif between rL and cL to define diagonal line #1
    b1 = rU - (slope1*cU)
    # Find dif between rL and cU to define other diagonal line #2
    b2 = rU - (slope2*cL)

    # Find the max distance the spokes can be when circular
    maxDist = distCalculator(cL, rC, cC, rC)

    for i in range(rL, rU + 1):
        for j in range(cL, cU + 1):
            # Diagonal spokes for rectangle border
            if fireLineShape == "rectangle":
                if i == np.floor((j*slope1)+b1):
                    if contained[i][j] == 0:
                        if j == cL or j == cU:
                            if breachFireLine(breachProbs, linesBuilt, breachProb, responseLineCounter) == False: veg[i][j] = 4
                            else: numBreached += 1
                            linesBuilt += 1   
                        else:
                            if breachFireLine(breachProbs, linesBuilt, breachProb, responseLineCounter) == False: veg[i][j] = 4
                            else: numBreached += 1
                            linesBuilt += 1   
                            if breachFireLine(breachProbs, linesBuilt, breachProb, responseLineCounter) == False: veg[i+1][j] = 4
                            else: numBreached += 1
                            linesBuilt += 1
                if i == np.floor((j*slope2)+b2):
                    if contained[i][j] == 0:
                        if j == cL or j == cU:
                            if breachFireLine(breachProbs, linesBuilt, breachProb, responseLineCounter) == False: veg[i][j] = 4
                            else: numBreached += 1
                            linesBuilt += 1   
                        else:
                            if breachFireLine(breachProbs, linesBuilt, breachProb, responseLineCounter) == False: veg[i][j] = 4
                            else: numBreached += 1
                            linesBuilt += 1   
                            if breachFireLine(breachProbs, linesBuilt, breachProb, responseLineCounter) == False: veg[i + 1][j] = 4
                            else: numBreached += 1
                            linesBuilt += 1    
            # Left to right spokes for both circular and rectangular
            if i == rC:
                if breachFireLine(breachProbs, linesBuilt, breachProb, responseLineCounter) == False:
                    if contained[i][j] == 0: veg[rC][j] = 4
                else: numBreached += 1
                linesBuilt += 1
            # Top to bottom spokes for both circular and rectangular
            if j == cC:
                if breachFireLine(breachProbs, linesBuilt, breachProb, responseLineCounter) == False:
                    if contained[i][j] == 0: veg[i][cC] = 4
                else: numBreached += 1
                linesBuilt += 1
    
    # Diagonal spokes for circular
    if fireLineShape == "circle":

        currRowPos1 = rC
        currColPos1 = cC
        currRowPos2 = rC
        currColPos2 = cC

        currRowNeg1 = rC
        currColNeg1 = cC
        currRowNeg2 = rC
        currColNeg2 = cC

        distFromCenter = 0

        dist1 = 0
        dist2 = 0
        dist3 = 0
        dist4 = 0
        slope1 = 1
        slope2 = -1

        while distFromCenter <= maxDist:
            currRowPos1 += slope1
            currColPos1 += slope1

            if contained[currRowPos1][currColPos1] == 0:
                dist1 = distCalculator(currColPos1, currRowPos1, cC, rC)
                if dist1 <= maxDist:
                    if breachFireLine(breachProbs, linesBuilt, breachProb, responseLineCounter) == False: veg[currRowPos1][currColPos1] = 4
                    else: numBreached += 1
                    linesBuilt += 1
                    if breachFireLine(breachProbs, linesBuilt, breachProb, responseLineCounter) == False: veg[currRowPos1][currColPos1-1] = 4
                    else: numBreached += 1
                    linesBuilt += 1
                    if breachFireLine(breachProbs, linesBuilt, breachProb, responseLineCounter) == False: veg[currRowPos1+1][currColPos1] = 4
                    else: numBreached += 1
                    linesBuilt += 1

            currRowPos2 -= slope1
            currColPos2 -= slope1

            if contained[currRowPos2][currColPos2] == 0:
                dist2 = distCalculator(currColPos2, currRowPos2, cC, rC)
                if dist2 <= maxDist:
                    if breachFireLine(breachProbs, linesBuilt, breachProb, responseLineCounter) == False: veg[currRowPos2][currColPos2] = 4
                    else: numBreached += 1
                    linesBuilt += 1
                    if breachFireLine(breachProbs, linesBuilt, breachProb, responseLineCounter) == False: veg[currRowPos2+1][currColPos2] = 4
                    else: numBreached += 1
                    linesBuilt += 1
                    if breachFireLine(breachProbs, linesBuilt, breachProb, responseLineCounter) == False: veg[currRowPos2][currColPos2-1] = 4
                    else: numBreached += 1
                    linesBuilt += 1

            currRowNeg1 -= slope2
            currColNeg1 += slope2

            if contained[currRowNeg1][currColNeg1] == 0:
                dist3 = distCalculator(currColNeg1, currRowNeg1, cC, rC)
                if dist3 <= maxDist:
                    if breachFireLine(breachProbs, linesBuilt, breachProb, responseLineCounter) == False: veg[currRowNeg1][currColNeg1] = 4
                    else: numBreached += 1
                    linesBuilt += 1
                    if breachFireLine(breachProbs, linesBuilt, breachProb, responseLineCounter) == False: veg[currRowNeg1+1][currColNeg1] = 4
                    else: numBreached += 1
                    linesBuilt += 1
                    if breachFireLine(breachProbs, linesBuilt, breachProb, responseLineCounter) == False: veg[currRowNeg1][currColNeg1+1] = 4
                    else: numBreached += 1
                    linesBuilt += 1

            currRowNeg2 += slope2
            currColNeg2 -= slope2

            if contained[currRowNeg2][currColNeg2] == 0:
                dist4 = distCalculator(currColNeg2, currRowNeg2, cC, rC)
                if dist4 <= maxDist:
                    if breachFireLine(breachProbs, linesBuilt, breachProb, responseLineCounter) == False: veg[currRowNeg2][currColNeg2] = 4
                    else: numBreached += 1
                    linesBuilt += 1
                    if breachFireLine(breachProbs, linesBuilt, breachProb, responseLineCounter) == False: veg[currRowNeg2-1][currColNeg2] = 4
                    else: numBreached += 1
                    linesBuilt += 1
                    if breachFireLine(breachProbs, linesBuilt, breachProb, responseLineCounter) == False: veg[currRowNeg2][currColNeg2-1] = 4
                    else: numBreached += 1
                    linesBuilt += 1

            distFromCenter = max(dist1, dist2, dist3, dist4)

    return linesBuilt, numBreached

# buildPrimaryLines:
    # 'i' is the latitude index of the cell where the fire line was breached
    # 'j' is the longitude index of the cell where the fire line was breached
    # 'contained' is a 2D matrix containing the status of the cell's containment 
    # 'veg' is a 2D matrix containing the cell vegetation type as an integer (4 indicates an unbreached fire border)
    # 'fireBuffer' is the number of cells between the active fire border and location of primary fire line
    # 'breachProbs' is the probability that any fire will jump a fire line
    # 'tempFireBorder' is the current position of the active fire border in each direction
def buildProactiveLines(i, j, contained, veg, concentricContingency, primaryBuffer, contingencyBuffer, breachProbs, linesBuilt, breachProb, tempFireBorder, fireLineShape, spokes, responseLineCounter, numBreached):
    # Setting the x,y lower and upper bounds for the fire line boarder
    rL = tempFireBorder[0] - primaryBuffer # row lower bound (northmost)
    rU = tempFireBorder[1] + primaryBuffer # row upper bound (southmost)
    cL = tempFireBorder[2] - primaryBuffer # column lower bound (westmost)
    cU = tempFireBorder[3] + primaryBuffer # column upper bound (eastmost)
    # check the fireLineShape and call the appropriate helper function
    proactive = True
    if fireLineShape == "rectangle":
        linesBuilt, responseLineCounter, numBreached = rectangleLine(veg, [rL, rU, cL, cU], contained, breachProbs, linesBuilt, breachProb, responseLineCounter, proactive, numBreached) 
        # print("response line counter after Proactive: ", responseLineCounter)
    elif fireLineShape == "circle":
        linesBuilt, responseLineCounter, numBreached = circularLines(veg, [rL, rU, cL, cU], contained, breachProbs, linesBuilt, breachProb, responseLineCounter, proactive, numBreached)
    
    if concentricContingency == True:
        rL = tempFireBorder[0] - contingencyBuffer # row lower bound (northmost)
        rU = tempFireBorder[1] + contingencyBuffer # row upper bound (southmost)
        cL = tempFireBorder[2] - contingencyBuffer # column lower bound (westmost)
        cU = tempFireBorder[3] + contingencyBuffer # column upper bound (eastmost)
        # check the fireLineShape and call the appropriate helper function
        if fireLineShape == "rectangle":
            if spokes == True: linesBuilt, numBreached = addSpokes(veg, [rL, rU, cL, cU], contained, breachProbs, linesBuilt, breachProb, fireLineShape, responseLineCounter, numBreached)
            linesBuilt, responseLineCounter, numBreached = rectangleLine(veg, [rL, rU, cL, cU], contained, breachProbs, linesBuilt, breachProb, responseLineCounter, proactive, numBreached) 
            # print("response line counter after contingency: ", responseLineCounter)
        elif fireLineShape == "circle":
            if spokes == True: linesBuilt, numBreached = addSpokes(veg, [rL-3, rU+3, cL-3, cU+3], contained, breachProbs, linesBuilt, breachProb, fireLineShape, responseLineCounter, numBreached)
            linesBuilt, responseLineCounter, numBreached = circularLines(veg, [rL, rU, cL, cU], contained, breachProbs, linesBuilt, breachProb, responseLineCounter, proactive, numBreached)
    return linesBuilt, numBreached   

# buildResponseLine:
    # 'i' is the latitude index of the cell where the fire line was breached
    # 'j' is the longitude index of the cell where the fire line was breached
    # 'contained' is a 2D matrix containing the status of the cell's containment 
    # 'veg' is a 2D matrix containing the cell vegetation type as an integer (4 indicates an unbreached fire border)
    # 'responseRadius' is the radius of the squre of the reponse line
    # 'breachProb' is the probability that any fire will jump a fire line
def buildResponseLine(i, j, contained, veg, responseRadius, breachProbs, linesBuilt, breachProb, fireLineShape, responseLineCounter, numBreached):
    # Setting the x,y lower and upper bounds for the fire line boarder to ensure they do not exceed the 
        # bounds of the map
    # left fire line boundary
    if i - responseRadius < 0: rL = 0
    else: rL = i - responseRadius
    # right fire line boundary
    if i + responseRadius > np.size(veg, 0)-1: rU = 0
    else: rU = i + responseRadius 
    # lower fire line boundary
    if j - responseRadius < 0: cL = 0
    else: cL = j - responseRadius         
    # upper fire line boundary
    if j + responseRadius > np.size(veg, 1)-1: cU = 0
    else: cU = j + responseRadius

    proactive = False

    # check the fireLineShape and call the appropriate helper function
    if fireLineShape == "rectangle":
        linesBuilt, responseLineCounter, numBreached = rectangleLine(veg, [rL, rU, cL, cU], contained, breachProbs, linesBuilt, breachProb, responseLineCounter, proactive, numBreached) 
        # print("response line counter after each response line: ", responseLineCounter)
    elif fireLineShape == "circle":
        linesBuilt, responseLineCounter, numBreached = circularLines(veg, [rL, rU, cL, cU], contained, breachProbs, linesBuilt, breachProb, responseLineCounter, proactive, numBreached) 
    return linesBuilt, responseLineCounter, numBreached

# showMaps: an alternative to showResults that just shows the cumulativeFire maps
def showMaps(map, cumulativeFire):
    # get a clean vegetation map
    veg = np.floor_divide(map[0], np.ones([np.size(map, 1), np.size(map, 2)], dtype=int)*100)
    rgbIMG = np.zeros([np.size(veg, 0), np.size(veg, 1), 3], dtype=int)
    r = np.add(np.add(np.where(veg == 1, 56, 0), np.where(veg == 2, 147, 0)), np.where(veg == 3, 219, 0))
    g = np.add(np.add(np.where(veg == 1, 118, 0), np.where(veg == 2, 196, 0)), np.where(veg == 3, 235, 0))
    b = np.add(np.add(np.where(veg == 1, 29, 0), np.where(veg == 2, 125, 0)), np.where(veg == 3, 118, 0))
    vegRGB = np.dstack((r, g, b))
    # clear zero-values from cumulative fire so they are transparent on the map
    cumulativeFire[cumulativeFire==0]=['nan']
    # Make a 2x3 plot for the cumulative fire maps
    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3)
    # Plot policy A
    ax1.imshow(vegRGB)
    plot1 = ax1.imshow(cumulativeFire[0], cmap='autumn_r')
    ax1.axis("off")
    ax1.set_title('Policy A')
    # Plot policy B
    ax2.imshow(vegRGB)
    plot1 = ax2.imshow(cumulativeFire[1], cmap='autumn_r')
    ax2.axis("off")
    ax2.set_title('Policy A')
    # Plot policy C
    ax3.imshow(vegRGB)
    plot3 = ax3.imshow(cumulativeFire[2], cmap='autumn_r')
    ax3.axis("off")
    ax3.set_title('Policy C')
    # Plot policy D
    ax4.imshow(vegRGB)
    plot4 = ax4.imshow(cumulativeFire[3], cmap='autumn_r')
    ax4.axis("off")
    ax4.set_title('Policy D')
    # Plot policy E
    ax5.imshow(vegRGB)
    plot5 = ax5.imshow(cumulativeFire[4], cmap='autumn_r')
    ax5.axis("off")
    ax5.set_title('Policy E')
    # Plot policy F
    ax6.imshow(vegRGB)
    plot6 = ax6.imshow(cumulativeFire[5], cmap='autumn_r')
    ax6.axis("off")
    ax6.set_title('Policy E')

    # Create the color bar and save the plot
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.025, 0.7])
    cbar = fig.colorbar(plot1, cax=cbar_ax)
    cbar.set_label('Number of replications burned')
    plt.savefig('Cumulative Fire Maps.png')  

# showOneRep: Plots the map of the area pre-burn and post-burn (USED ONLY FOR VERIFICATION AND DEBUGGING)
    # 'fire' is a 2D matrix containing the percentage on fire of each cell on the map
    # 'veg' is a 2D matrix containing the cell vegetation type as an integer (0 = unburnable, 1 = trees, 2 = shrub, 3 = herb)
def showOneRep(fire, veg):
    # Build pre-burn landscape image
    rgbIMG = np.zeros([np.size(veg, 0), np.size(veg, 1), 3], dtype=int)
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
                        lines.Line2D([0], [0], marker='o', color='w', label='Herb', markerfacecolor='#dbeb76', markersize=10),
                        lines.Line2D([0], [0], marker='o', color='w', label='Fire Line', markerfacecolor='#333000', markersize=10)]
    # Plot figures
    fig, ((ax1, ax2)) = plt.subplots(1, 2)
    ax1.imshow(rgbIMG)
    ax2.imshow(newRGB)
    ax1.legend(handles = ax2legendElements, bbox_to_anchor=(1.5, 1.0), loc='lower center')
    plt.show()
    return 0
