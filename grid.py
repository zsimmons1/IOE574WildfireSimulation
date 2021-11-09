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
import math

img = rasterio.open('/Users/sprin/OneDrive/Desktop/IOE574/TermProject/Data/NorthernCali/us_210evc_us_210evc.tif')
show(img)

map = img.read()
fire = np.zeros(map.shape, dtype=rasterio.float32)  # Create empty matrix


# fire[0][x][x] = [0, 1, 2, 3] [not burning, partially burned, fully burning, burnt]
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
    percentage = (vegetation_value - 100) / 100
    fire_counter = fire_counter * percentage
    if vegetation_value > 99 and vegetation_value < 200:
        if fire_counter >= 100:
            return 100
    return 0

def willShrubAreaBurn(vegetation_value, fire_counter):
    # Case 1
    percentage = vegetation_value - 200
    if vegetation_value > 199 and vegetation_value < 300: 
        if fire_counter >= 600:
            return 100   
    return 0

def willHerbAreaBurn(vegetation_value, fire_counter):
    percentage = vegetation_value - 300
    if vegetation_value > 299: 
        if fire_counter >= 700:
            return 100  
    return 0

def calcWindEffects(w_s, w_d, i_pos, j_pos): # where w_s is wind speed, w_d is wind direction, and i_pos and j_pos are -1, 0, or 1
    C1 = 0.045 # (m^-1 * s)
    C2 = 0.131 # (m^-1 * s)
    
    if i_pos == -1:
        if j_pos == -1:
            s_d = 315
        elif j_pos == 0:
            s_d = 270
        elif j_pos == 1:
            s_d = 225
        else:
            print("invalid j position")            
    elif i_pos == 0:
        if j_pos == -1:
            s_d = 0
        elif j_pos == 1: 
            s_d = 180
        else:
            print("invalid j position")           
    elif i_pos == 1:
        if j_pos == -1:
            s_d = 45
        elif j_pos == 0:
            s_d = 90
        elif j_pos == 1:
            s_d = 135
        else:
            print("invalid j position")    
    else:
        print("invalid i position")

    theta = math.pi*abs(w_d - s_d)/180 
    f_t = math.exp(w_s*C2*(math.cos(theta) - 1))
    Pw = f_t*exp(C1*w_s)
    
    return Pw

def calcBurn(map, fire, i, j, w_s, w_d):
    P_burn = 0
    if map[0][i][j] >= 100 & map[0][i][j] <= 199:
        P_veg = 0.8 # TODO: Adjust this value based on literature!
        P_den = 0.7*(veg - 100)/100 + -0.4 # Based on -0.4 to 0.3 range in Mutthulakshimi et al.
    elif map[0][i][j] >= 200 & map[0][i][j] <= 299:  
        P_veg = 0.5 # TODO: Adjust this value based on literature!
        P_den = 0.7*(veg - 200)/100 + -0.4
    elif map[0][i][j] >= 300 & map[0][i][j] <= 399:
        P_veg = 0.2 # TODO: Adjust this value based on literature!
        P_den = 0.7*(veg - 300)/100 + -0.4
    else:
        print("Invalid vegetation value")  

    C1 = 0.045 # (m^-1 * s)
    C2 = 0.131 # (m^-1 * s)

    P_w = fire[i-1][j-1]*math.exp(w_s*C2*(math.cos(math.pi*abs(w_d - 315)/180) - 1))*exp(C1*w_s)
    P_burn = 1/sqrt(2)*P_h*(1+P_den)*(1+P_veg)*P_w
    if np.random.uniform >= P_burn:
        return 1
    else:    
        P_w = fire[i-1][j]*math.exp(w_s*C2*(math.cos(math.pi*abs(w_d - 270)/180) - 1))*exp(C1*w_s)
        P_burn = P_h*(1+P_den)*(1+P_veg)*P_w
        if np.random.uniform >= P_burn:
            return 1
        else:
            P_w = fire[i-1][j+1]*math.exp(w_s*C2*(math.cos(math.pi*abs(w_d - 225)/180) - 1))*exp(C1*w_s)
            P_burn = 1/sqrt(2)*(1+P_den)*(1+P_veg)*P_w
            if np.random.uniform >= P_burn:
                return 1
            else:    
                P_w = fire[i][j+1]*math.exp(w_s*C2*(math.cos(math.pi*abs(w_d - 0)/180) - 1))*exp(C1*w_s)
                P_burn = P_h*(1+P_den)*(1+P_veg)*P_w
                if np.random.uniform >= P_burn:
                    return 1
                else:    
                    P_w = fire[i][j-1]*math.exp(w_s*C2*(math.cos(math.pi*abs(w_d - 180)/180) - 1))*exp(C1*w_s)
                    P_burn = P_h*(1+P_den)*(1+P_veg)*P_w
                    if np.random.uniform >= P_burn:
                        return 1
                    else:       
                        P_w = fire[i+1][j-1]*math.exp(w_s*C2*(math.cos(math.pi*abs(w_d - 45)/180) - 1))*exp(C1*w_s)
                        P_burn = 1/sqrt(2)*P_h*(1+P_den)*(1+P_veg)*P_w
                        if np.random.uniform >= P_burn:
                            return 1
                        else:    
                            P_w = fire[i+1][j]*math.exp(w_s*C2*(math.cos(math.pi*abs(w_d - 90)/180) - 1))*exp(C1*w_s)
                            P_burn = P_h*(1+P_den)*(1+P_veg)*P_w
                            if np.random.uniform >= P_burn:
                                return 1
                            else:    
                                P_w = fire[i+1][j+1]*math.exp(w_s*C2*(math.cos(math.pi*abs(w_d - 135)/180) - 1))*exp(C1*w_s)
                                P_burn = 1/sqrt(2)*P_h*(1+P_den)*(1+P_veg)*P_w
                                if np.random.uniform >= P_burn:
                                    return 1
                                else:
                                    return 0



# Trying a simple case with 5 time steps and starting the fire at the center of our grid
while t != 5:
    temp_fire = fire
    for i in range(fire_locations[0] - 1, fire_locations[1] + 2):
        for j in range(fire_locations[2] - 1, fire_locations[3] + 2):
            if map[0][i][j] > 99:
                w_s = 8 # TODO: Make stochastic
                w_d = 0 # TODO: Make stochastic
                P_burn = calcBurnProb(map, fire, i, j, w_s, w_d)               
                temp_fire[i][j] = P_burn
                    
                #     temp_fire[0][i][j] = 100
                if i < fire_locations[0]:
                    min_i = i
                if i > fire_locations[1]:
                    max_i = i
                if j < fire_locations[2]:
                    min_j = j
                if j > fire_locations[3]:
                    max_j = j
    
    # fire_locations = [min_i, max_i, min_j, max_j]
    fire = temp_fire
    t += 1

show(fire, cmap='Reds')