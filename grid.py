# Importing GIS data to setup a grid
# Define all cells according to vegetation type

# Storing all information related to each individual cell

import os 
import matplotlib.pyplot as plt
import geopandas as gpd
import earthpy as et
import numpy as np
import rasterio
from rasterio.plot import show

img = rasterio.open('/Users/Zack/Desktop/IOE574/TermProject/IOE574WildfireSimulation/us_210evc.tif')
show(img)

map = img.read()
fire = np.zeros(map.shape, dtype=rasterio.float32)  # Create empty matrix

fire[0][28][24] = 100
fire_locations = [28, 28, 24, 24] # [x_lower_bound, x_upper_bound, y_lower_bound, y_upper_bound]

t = 0
min_i = 28
max_i = 28
min_j = 24
max_j = 24

# 11 = Water (anything less than 100)
# 1xx = Tree, xx = % Cover, most likely to spread
# 2xx = Shrub, xx = % Cover, likely to spread
# 3xx = Herb, xx = % Cover, least likely to spread

# 100 = sparse vegetation canopy

def willTreeAreaBurn(vegetation_value, fire_counter):
    # Case 1
    if vegetation_value > 99 and vegetation_value < 200:
        if fire_counter >= 100:
            return 100
    return 0

def willShrubAreaBurn(vegetation_value, fire_counter):
    # Case 1
    if vegetation_value > 199 and vegetation_value < 300: 
        if fire_counter >= 600:
            return 100   
    return 0

def willHerbAreaBurn(vegetation_value, fire_counter):
    if vegetation_value > 299: 
        if fire_counter >= 700:
            return 100   
    return 0


# Trying a simple case with 5 time steps and starting the fire at the center of our grid
while t != 5:
    temp_fire = fire
    for i in range(fire_locations[0] - 1, fire_locations[1] + 2):
        print("i: ", i)
        for j in range(fire_locations[2] - 1, fire_locations[3] + 2):
            print("j: ", j)
            if map[0][i][j] > 99:
                fire_counter = fire[0][i+1][j] + fire[0][i-1][j] + fire[0][i][j+1] + fire[0][i][j-1] + fire[0][i+1][j+1] + fire[0][i+1][j-1] + fire[0][i-1][j+1] + fire[0][i-1][j-1]
                if map[0][i][j] > 99 and map[0][i][j] < 200:
                    temp_fire[0][i][j] = willTreeAreaBurn(map[0][i][j], fire_counter)
                elif map[0][i][j] > 199 and map[0][i][j] < 300:
                    temp_fire[0][i][j] = willShrubAreaBurn(map[0][i][j], fire_counter)
                elif map[0][i][j] > 299:
                    temp_fire[0][i][j] = willHerbAreaBurn(map[0][i][j], fire_counter)                
                
                # if fire[0][i+1][j] == 100 or fire[0][i-1][j] == 100 or fire[0][i][j+1] == 100 or fire[0][i][j-1] == 100 or fire[0][i+1][j+1] == 100 or fire[0][i+1][j-1] == 100 or fire[0][i-1][j+1] == 100 or fire[0][i-1][j-1] == 100:
                    
                #     temp_fire[0][i][j] = 100
                if i < fire_locations[0]:
                    min_i = i
                if i > fire_locations[1]:
                    print("step in x")
                    max_i = i
                if j < fire_locations[2]:
                    min_j = j
                if j > fire_locations[3]:
                    print("step in y")
                    max_j = j
    
    fire_locations = [min_i, max_i, min_j, max_j]
    fire = temp_fire
    t += 1

show(fire)