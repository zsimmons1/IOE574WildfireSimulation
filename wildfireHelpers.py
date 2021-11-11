# Include libraries
import os 
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import numpy as np
import rasterio
from rasterio.plot import show
import math
from astropy.visualization import make_lupton_rgb

# Calculates the area of burn in cell i,j resulting from the spread distance in each direction and the 
# intensity of spread in each direction. Area of spread is then used to calculate total proportion of 
# the cell which is on fire.
    # distance is an 8-element vector of the distance of spread in each direction, from neighbors in the
    # following order: SW, W, NW, S, N, SE, E, NE
    # spread_prob is an 8-element vector of the spread probability from each neighbor to cell i,j with 
    # neighbor order same as distance vector
def spreadFire(distance, spread_prob):
    orthogonal = [1, 3, 4, 6]
    diagonal = [0, 2, 5, 7]
    burnArea = 0
    # Add burn contribution from orthogonal neighbors
    for i in orthogonal:
        burnArea +=30*spread_prob[i]*distance[i]    
    # Add burn contribution from diagonal neighbors
    for i in diagonal:
        if distance[i] < 15*pow(2, 0.5):
            # Fire area is a triangle
            burnArea += spread_prob[i]*pow(distance[i], 2)
        else:
            # Non-fire area is a triangle  
            burnArea += spread_prob[i]*(900-pow((30*pow(2, 0.5)-distance[i]), 2)) 
    # Subtract overlap from orthogonal corner overlap
    for i in [1, 6]:
        for j in [3, 4]:
            burnArea -= spread_prob[i]*spread_prob[j]*distance[i]*distance[j]
    # Subtract overlap from orthogonal parallel overlap
    if (distance[1] + distance[6]) > 30:
        burnArea -= 30*spread_prob[1]*spread_prob[6]*abs(distance[1]-distance[6]) 
    if (distance[3] + distance[4])  > 30:
        burnArea = burnArea - 30*spread_prob[3]*spread_prob[4]*abs(distance[3]-distance[4])   
    # Subtract overlap from parallel diagonal cells (0-2, 0-5, 7-2, 7-5)
    for i in [1, 7]:
        for j in [2, 5]:
            # TODO: Finish equation
            burnArea = burnArea
    # Subtract overlap from diagonal diagonal cells (0-7, 2-5)
    if (distance[0]+distance[7]) > 30*pow(2, 0.5):
        # TODO: Test equation
        burnArea = burnArea - (900 - (spread_prob[0]*pow(30*pow(2, 0.5) - distance[0], 2)
        +spread_prob[7]*pow(30*pow(2, 0.5) - distance[7], 2)))
    if (distance[2]+distance[5]) > 30*pow(2, 0.5):
        # TODO: Test equation
        burnArea = burnArea - (900 - (spread_prob[2]*pow(30*pow(2, 0.5) - distance[2], 2)
        +spread_prob[5]*pow(30*pow(2, 0.5) - distance[5], 2)))
    # Subtract overlap from diagonal and orthogonal adjacent neighbors
    # (0-1, 0-3, 2-1, 2-4, 5-3, 5-6, 7-4, 7-6)
        # TODO: Finish equation
    # Subtract overlap from diagonal and orthogonal non-adjacent neighbors   
    # (0-4, 0-6, 2-3, 2-5, 5-1, 5-4, 7-1, 7-3)     
        # TODO: Finish equation
    return burnArea / 900

# Test spreadFire:
distance = [0, 0, 0, 20, 15, 0, 0, 0]    
spread_prob = [1, 1, 1, 1, 1, 1, 1, 1]
burnPercent = spreadFire(distance, spread_prob)
print(burnPercent)