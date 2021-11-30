# Include libraries
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.lines as lines
from rasterio.plot import show
from scipy.stats import weibull_min
from array import array
import pandas
import random

# initializeWind: Based on input weather data, determine the appropriate Weibull distribution for wind speed
# and draw from this distribution to determine the initial wind speed across all cells
def initializeWind():
    data = pandas.read_csv("weather.csv")
    messy_wind_speed = data["HourlyWindSpeed"]
    # remove nan values
    wind_speed = [x for x in messy_wind_speed if np.isnan(x) == False]
    # determine the Weibull distrbution curve that fits the weather data
    shape, loc, scale = weibull_min.fit(wind_speed, floc = 1)
    # generate a random wind speeds from a weibull distribution using parameters generated from actual data
    initial_wind = weibull_min.rvs(shape, loc, scale, size=1)
    # return the generated initial wind as a single value (not an array) in km/hr
    return (3600/1000)*initial_wind[0]

# updateWind: Determines the wind direction and velocity for cell i,j
    # 'wind_speed' is the wind speed (in km/hour) across all cells from the previous time period
    # 'wind_direction' is the angle (in degrees, where East is 0 degrees) of the wind across all cells from the 
        # previous time period
def updateWind(wind_speed, wind_direction):
    # Apply multaplicative uniform noise to the wind speed U[0.8, 1.2] per Trucchia et al
    wind_speed = np.random.uniform(0.8*wind_speed, 1.2*wind_speed)
    # Apply multaplicative uniform noise to the wind direction U[−11.25◦, 11.25◦] per Trucchia et al
    wind_direction = np.random.uniform(wind_direction - 11.25, wind_direction + 11.25)
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
def igniteCell(fire, i, j, distance, wind_speed, wind_direction):
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
            # input angle into funciton to get alpha_w facor
            alpha_w = 1
            spread_prob[k] = spread_prob[k]*alpha_w
            
            # Ensure all probabilities are no more than 1
            if spread_prob[k] > 1:
                spread_prob[k] = 1

            # Draw from uniform distribution to determine if neighbor k will spread to cell i,j  
            u = np.random.uniform(0,1)   
            if u < spread_prob[k]:
                cell_transition[k] = 1          
            # Else, neighbor k does not contribute to the fire in cell i,j, so 'cell_transition[k]' remains 0

    # Return the cell_transition and spread_prob
    return cell_transition, spread_prob


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
def advanceBurn(veg, cell_transition, distance, del_t):
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
            del_x = math.pow(10, np.random.weibull(shape[veg-1]))*del_t
            if distance[x] + del_x > 30:
                distance[x] = 30
            else: 
                distance[x] += del_x

    for x in diagonal:
        # If this neighbor is contributing, then we will calculate the distance
        if shape[veg-1] == 0:
            distance[x] = 0
        if cell_transition[x] == 1 and shape[veg-1] != 0:
            del_x = math.pow(10, np.random.weibull(shape[veg-1]))*del_t
            # diagonal neighbors can have spread distance up to 30*square root of 2 m
            if distance[x] + del_x > 30*pow(2, 0.5):
                distance[x] = 30*pow(2, 0.5)
            else: 
                distance[x] += del_x
    
    return distance


# spreadFire: Calculates the area of burn in cell i,j resulting from the spread distance in each direction and the 
# intensity of spread in each direction. Area of spread is then used to calculate total proportion of 
# the cell which is on fire.
    # 'distance' is an 8-element vector of the distance of spread in each direction, from neighbors in the
        # following order: NW, N, NE, W, E, SW, S, SE
    # 'spread_prob' is an 8-element vector of the spread probability from each neighbor to cell i,j with 
        # neighbor order same as distance vector
def spreadFire(distance, spread_prob):
    orthogonal = [1, 3, 4, 6]
    diagonal = [0, 2, 5, 7]
    burnArea = 0

    # Equation 1: Add burn contribution from orthogonal neighbor
    for a in orthogonal:
        burnArea += spread_prob[a]*30*distance[a]    

    # Equation 2: Add burn contribution from diagonal neighbor
    for a in diagonal:
        if distance[a] < 15*pow(2, 0.5):
            # Equation 2.1: Add burn contribution from diagonal neighbor when diagonal distance ≤ half
            burnArea += spread_prob[a]*pow(distance[a], 2)
        else:
            # Equation 2.2: Add burn contribution from diagonal neighbor when diagonal distance ≥ half 
            burnArea += spread_prob[a]*(900-pow((30*pow(2, 0.5)-distance[a]), 2))

    # Return burn percentage 
    if (burnArea > 900): 
        burnArea = 900       
    return burnArea / 900


# showResults: Plots the map of the area pre-burn and post-burn
    # 'fire' is a 2D matrix containing the percentage on fire of each cell on the map
    # 'veg' is a 2D matrix containing the cell vegetation type as an integer (0 = unburnable, 1 = trees, 2 = shrub, 3 = herb)
def showResults(fire, veg):
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
    ax2.legend(handles = ax2legendElements, bbox_to_anchor=(1.5, 1.0), loc='upper right')
    plt.show()
    show(fire, cmap='Reds')

    return 0

# showTimeline:
def showTimeline(fire_timeline, veg):
    # Build pre-burn landscape image
    rgbIMG = np.zeros([np.size(veg, 0), np.size(veg, 1), 3], dtype=int)
    r = np.add(np.add(np.where(veg == 1, 56, 0), np.where(veg == 2, 147, 0)), np.where(veg == 3, 219, 0))
    g = np.add(np.add(np.where(veg == 1, 118, 0), np.where(veg == 2, 196, 0)), np.where(veg == 3, 235, 0))
    b = np.add(np.add(np.where(veg == 1, 29, 0), np.where(veg == 2, 125, 0)), np.where(veg == 3, 118, 0))
    rgbIMG = np.dstack((r, g, b))

    fires = np.zeros((20, np.size(fire_timeline, 0), np.size(fire_timeline, 1)), dtype=int)
    
    for i in range(20):
        # Overlay burn
        new_r = np.where(fire_timeline[i] > 0, 255, r)
        new_g = np.where(fire_timeline[i] > 0, 0, g)
        new_b = np.where(fire_timeline[i] > 0, 0, b)
        fires[i] = np.dstack((new_r, new_g, new_b))

        # Create custom legends
        ax2legendElements = [lines.Line2D([0], [0], marker='o', color='w', label='Fire', markerfacecolor='#ff0000', markersize=10),
                            lines.Line2D([0], [0], marker='o', color='w', label='Tree', markerfacecolor='#38761d', markersize=10), 
                            lines.Line2D([0], [0], marker='o', color='w', label='Shrub', markerfacecolor='#93c47d', markersize=10), 
                            lines.Line2D([0], [0], marker='o', color='w', label='Herb', markerfacecolor='#dbeb76', markersize=10),
                            lines.Line2D([0], [0], marker='o', color='w', label='Fire Line', markerfacecolor='#333000', markersize=10)]

    # Plot figures
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    ax1.imshow(rgbIMG)
    ax2.imshow(fires[0])
    ax3.imshow(fires[1])
    ax4.imshow(fires[2])
    ax2.legend(handles = ax2legendElements, bbox_to_anchor=(1.5, 1.0), loc='upper right')
    plt.show()
    show(fire, cmap='Reds')

    return 0

# Simple function to tell whether or not the fire will jump the fire line, can add complexity later
def jumpFireLine(jump_prob):
    new_number = random.random()
    if new_number > (1-jump_prob):
        return True
    else:
        return False

# buildResponseLine:
    # i is the latitude index of the cell where the fire line was breached/jumped
    # j is the longitude index of the cell where the fire line was breached/jumped
    # response_radius is the radius of the squre of the reponse line
    # contingency_border is a vector of [x_lower, x_upper, y_lower, y_upper] of the fire line that was breached
    # jump_prob is the probability that any fire will jump a fire line
def buildResponseLine(i,j, veg, response_radius, contingency_border, jump_prob, numLinesJumped):
        # Setting the x,y lower and upper bounds for the fire line boarder
        if i - response_radius < 0:
            x_lower = 0
        else:
            x_lower = i - response_radius

        if i + response_radius > np.size(veg, 0)-1:
            x_upper = 0
        else:
            x_upper = i + response_radius 

        if j - response_radius < 0:
            y_lower = 0
        else:
            y_lower = j - response_radius         
        
        if j + response_radius > np.size(veg, 1)-1:
            y_upper = 0
        else:
            y_upper = j + response_radius

        # Setting the vegitation type to 4 for the fire lines
        for i in range(x_lower, x_upper+1):
                # (i < contingency_border[2] or i > contingency_border[3]) 
            if (i < np.size(veg, 0)-1) and (i > 0):
                canJump = jumpFireLine(jump_prob)
                if canJump == False:
                    veg[i][y_lower] = 4 # Representative of fire boarder
                    veg[i][y_upper] = 4
                else:
                    numLinesJumped += 1
        for j in range(y_lower, y_upper+1):
            # (j < contingency_border[0] or i > contingency_border[1])
            if (j < np.size(veg, 1)-1) and (j > 0):
                canJump = jumpFireLine(jump_prob)
                if canJump == False:
                    veg[x_lower][j] = 4 # Representative of fire boarder
                    veg[x_upper][j] = 4
                else:
                    numLinesJumped += 1      


# midPointCircleDraw: 
# Source: https://www.geeksforgeeks.org/mid-point-circle-drawing-algorithm/
    # 'x_centre' 
    # 'y_centre'
    # 'r' is the radius of the circle
def midPointCircleDraw(x_centre, y_centre, r):
    x = r
    y = 0
     
    # Printing the initial point the
    # axes after translation
    print("(", x + x_centre, ", ",
               y + y_centre, ")",
               sep = "", end = "")
     
    # When radius is zero only a single
    # point be printed
    if (r > 0) :
     
        print("(", x + x_centre, ", ",
                  -y + y_centre, ")",
                  sep = "", end = "")
        print("(", y + x_centre, ", ",
                   x + y_centre, ")",
                   sep = "", end = "")
        print("(", -y + x_centre, ", ",
                    x + y_centre, ")", sep = "")
     
    # Initialising the value of P
    P = 1 - r
 
    while x > y:
     
        y += 1
         
        # Mid-point inside or on the perimeter
        if P <= 0:
            P = P + 2 * y + 1
             
        # Mid-point outside the perimeter
        else:        
            x -= 1
            P = P + 2 * y - 2 * x + 1
         
        # All the perimeter points have
        # already been printed
        if (x < y):
            break
         
        # Printing the generated point its reflection
        # in the other octants after translation
        print("(", x + x_centre, ", ", y + y_centre,
                            ")", sep = "", end = "")
        print("(", -x + x_centre, ", ", y + y_centre,
                             ")", sep = "", end = "")
        print("(", x + x_centre, ", ", -y + y_centre,
                             ")", sep = "", end = "")
        print("(", -x + x_centre, ", ", -y + y_centre,
                                        ")", sep = "")
         
        # If the generated point on the line x = y then
        # the perimeter points have already been printed
        if x != y:
         
            print("(", y + x_centre, ", ", x + y_centre,
                                ")", sep = "", end = "")
            print("(", -y + x_centre, ", ", x + y_centre,
                                 ")", sep = "", end = "")
            print("(", y + x_centre, ", ", -x + y_centre,
                                 ")", sep = "", end = "")
            print("(", -y + x_centre, ", ", -x + y_centre,
                                            ")", sep = "")