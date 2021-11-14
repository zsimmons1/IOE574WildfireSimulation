# Include libraries
import numpy as np
import math

# igniteCell: 
def igniteCell(fire, i, j, distance, wind_speed, wind_direction):
    cell_transition = [0, 0, 0, 0, 0, 0, 0, 0]
    spread_angles = [45, 0, 315, 90, 270, 135, 180, 225]

    # Set baseline spread probability for each neighbor as the percent on fire of that neighbor
    spread_prob = [fire[i-1][j+1], fire[i-1][j], fire[i-1][j-1],
                   fire[i][j+1], fire[i][j-1],
                   fire[i+1][j+1], fire[i+1][j], fire[i+1][j-1]]        

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
            alpha_w = 0.5
            spread_prob[k] = spread_prob[k]*alpha_w
            
            # Ensure all probabilities are no more than 1
            if spread_prob[k] > 1:
                spread_prob[k] = 1

            # Draw from uniform distribution to determine if neighbor k will spread to cell i,j     
            if np.random.uniform < spread_prob[k]:
                cell_transition[k] = 1          
            # Else, neighbor k does not contribute to the fire in cell i,j, so 'cell_transition[k]' remains 0

    # Return the cell_transition and spread_prob
    return cell_transition, spread_prob


# advanceBurn: 
def advanceBurn(veg, cell_transition, distance, del_t):
    # shape stores the Weibull distribution shape for each vegetation type
    # Cell transition, size 8, 0s and 1s, 0 not contributing, 1 contributing
    shape = [11.4, 13.6, 13.0]
    
    for x in range(len(distance)):
        # If this neighbor is contributing, then we will calculate the distance
        if cell_transition[x] == 1:
            del_x = math.pow(10, np.random.weibull(shape[veg-1]))*del_t
            if distance[x] + del_x > 30:
                distance[x] = 30
            else: 
                distance[x] += del_x
    
    return distance


# spreadFire: Calculates the area of burn in cell i,j resulting from the spread distance in each direction and the 
# intensity of spread in each direction. Area of spread is then used to calculate total proportion of 
# the cell which is on fire.
    # 'distance' is an 8-element vector of the distance of spread in each direction, from neighbors in the
        # following order: SW, W, NW, S, N, SE, E, NE
    # 'spread_prob' is an 8-element vector of the spread probability from each neighbor to cell i,j with 
        # neighbor order same as distance vector
def spreadFire(distance, spread_prob):
    orthogonal = [1, 3, 4, 6]
    diagonal = [0, 2, 5, 7]
    burnArea = 0

    # Equation 1: Add burn contribution from orthogonal neighbor
    for i in orthogonal:
        burnArea += 30*spread_prob[i]*distance[i]    

    # Equation 2: Add burn contribution from diagonal neighbor
    for i in diagonal:
        if distance[i] < 15*pow(2, 0.5):
            # Equation 2.1: Add burn contribution from diagonal neighbor when diagonal distance ≤ half
            burnArea += spread_prob[i]*pow(distance[i], 2)
        else:
            # Equation 2.2: Add burn contribution from diagonal neighbor when diagonal distance ≥ half 
            burnArea += spread_prob[i]*(900-pow((30*pow(2, 0.5)-distance[i]), 2)) 

    # Equation 3: Subtract overlap from opposite orthogonal neighbors (1-6, 3-4)
    for i, j in [1, 6], [3, 4]:
        # Equation 3.1: Subtract overlap from opposite orthogonal neighbors when they overlap (1-6, 3-4)
        if (distance[i] + distance[j]) > 30:
            burnArea -= 30*((spread_prob[i] + spread_prob[j])/2)*(30 - (distance[i]+distance[j])) 
        # Else, Equation 3.0: opposite orthogonal neighbors do not overlap

    # Equation 4: Subtract overlap from kitty-corner orthogonal neighbors (1-3, 1-4, 6-3, 6-4)
    for i in [1, 6]:
        for j in [3, 4]:
            burnArea -= ((spread_prob[i] + spread_prob[j])/2)*distance[i]*distance[j]

    # Equation 5: Subtract overlap from opposite diagonal neighbors (0-7, 2-5)
    for i, j in [0, 7], [2, 5]:    
        if distance[i]+distance[j] > 30*pow(2, 0.5):
            # Equation 5.1: Subtract overlap from opposite diagonal neighbors when both diagonal distances > half
            if (distance[i] > 15*pow(2, 0.5)) & (distance[j] > 15*pow(2, 0.5)):
                A = 900-pow((30*pow(2, 0.5)-distance[i]), 2)
                B = 900-pow((30*pow(2, 0.5)-distance[j]), 2)
                burnArea -= -((spread_prob[i] + spread_prob[j])/2)*(900 - A - B)
            # Equation 5.2: Subtract overlap from opposite diagonal neighbors when only one diagonal distance is > half
            else:
                if distance[i] > 15*pow(2, 0.5):
                    B = pow(distance[j], 2)
                    A_inv = pow((30*pow(2, 0.5)-distance[i]), 2) 
                else: 
                    B = pow(distance[i], 2)
                    A_inv = pow((30*pow(2, 0.5)-distance[j]), 2) 
                burnArea -= ((spread_prob[i] + spread_prob[j])/2)*(B - A_inv)  
        # Else, Equation 5.0: opposite diagonal neighbors do not overlap 
        
    # Equation 6: Subtract overlap from diagonal neighbors in same row or same column (0-2, 0-5, 7-2, 7-5)
    for i in [0, 7]:
        for j in [2, 5]:
            if (pow(2, 0.5)*distance[i] + pow(2, 0.5)*distance[j] > 30) & (distance[i]!= 0) & (distance[j]!= 0):
                # Equation 6.1: Subtract overlap from diagonal neighbors in same row or same column when both diagonal distances ≤ half
                if (distance[i] <= 15*pow(2, 0.5)) & (distance[j] <= 15*pow(2, 0.5)):
                    print(burnArea)
                    burnArea -= ((spread_prob[i] + spread_prob[j])/2)*(pow(30-(pow(2, 0.5)*distance[i] - pow(2, 0.5)*distance[j]), 2) / 4)
                    print(burnArea)
                elif (distance[i] >= 15*pow(2, 0.5)) & (distance[j] >= 15*pow(2, 0.5)):
                    # Equation 6.2: Subtract overlap from diagonal neighbors in same row or same column when both diagonal distances are > half and there is not complete overlap
                    if pow(2, 0.5)*(30*pow(2, 0.5) - distance[i]) + pow(2, 0.5)*(30*pow(2, 0.5) - distance[j]) > 30:
                        A = pow((30*pow(2, 0.5)-distance[i]), 2)
                        B = pow((30*pow(2, 0.5)-distance[j]), 2)
                        C = pow(30 - (30*pow(2, 0.5)-distance[i] + 30*pow(2, 0.5)-distance[j]), 2)
                        burnArea -= ((spread_prob[i] + spread_prob[j])/2)*(900 - A - B + C)
                    else:
                    # Equation 6.3: Subtract overlap from diagonal neighbors in same row or same column when both diagonal distances are > half and there is complete overlap
                        burnArea -= ((spread_prob[i] + spread_prob[j])/2)*(900 - pow((30*pow(2, 0.5)-distance[i]), 2) - pow((30*pow(2, 0.5)-distance[j]), 2))
                # Equation 6.4: Subtract overlap from diagonal neighbors in same row or same column when only one diagonal distance is > half and there is not complete overlap
                elif pow(2, 0.5)*distance[j] + pow(2, 0.5)*(30*pow(2, 0.5) - distance[i]) > 30:  
                    if distance[i] > 15*pow(2, 0.5):
                        B = pow(distance[j], 2)
                        B_no_C = pow((30 - (pow(2, 0.5)*distance[j] + pow(2, 0.5)*(30*pow(2, 0.5) - distance[i])))*pow(2, 0.5), 2)
                    else:
                        B = pow(distance[i], 2)
                        B_no_C = pow((30 - (pow(2, 0.5)*distance[i] + pow(2, 0.5)*(30*pow(2, 0.5) - distance[j])))*pow(2, 0.5), 2)
                    burnArea -= ((spread_prob[i] + spread_prob[j])/2)*(B - B_no_C) 
                # Equation 6.5: Subtract overlap from diagonal neighbors in same row or same column when only one diagonal distance is > half and there complete overlap
                else:
                    if distance[i] > 15*pow(2, 0.5):
                        burnArea -= ((spread_prob[i] + spread_prob[j])/2)*pow(distance[j], 2)
                    else:
                        burnArea -= ((spread_prob[i] + spread_prob[j])/2)*pow(distance[i], 2)
            # Else: Equation 6.0: Diagonal neighbors in same row or column do not overlap

    # Equation 7: Subtract overlap from adjacent diagonal and orthogonal neighbors (0-1, 0-3, 2-1, 2-4, 5-3, 5-6, 7-4, 7-6)
    for i, j in [0,1], [0,3], [2,1], [2,4], [5,3], [5,6], [7,4], [7,6]:
        if distance[i]!=0 & distance[j]!=0:
            # Equation 7.1: Subtract overlap from adjacent orthogonal and diagonal neighbors when orthogonal distance completely overlaps diagonal
            if (distance[i] <= 15*pow(2, 0.5)) & (distance[j] > pow(2, 0.5)*distance[i]):
                burnArea -= ((spread_prob[i] + spread_prob[j])/2)*pow(distance[i], 2)
            # Equation 7.2: Subtract overlap from adjacent orthogonal and diagonal neighbors when diagonal distance completely overlaps orthogonal    
            elif pow(2, 0.5)*(30*pow(2, 0.5)-distance[i]) < distance[j]:
                burnArea -= ((spread_prob[i] + spread_prob[j])/2)*30*distance[j]
            # Equation 7.3: Subtract overlap from adjacent orthogonal and diagonal neighbors when diagonal distance ≤ half and there is not complete overlap   
            elif (distance[i] <= 15*pow(2, 0.5)) & (distance[j] < pow(2, 0.5)*distance[i]):
                B = pow(distance[i], 2)
                B_no_C = pow(pow(2, 0.5)*distance[i] - distance[j], 2)/2
                burnArea -= ((spread_prob[i] + spread_prob[j])/2)*(B-B_no_C)
            # Equation 7.4: Subtract overlap from adjacent orthogonal and diagonal neighbors when diagonal distance ≥ half 
            elif pow(2, 0.5)*(30*pow(2, 0.5)-distance[i]) + distance[j] > 30:   
                A_side = 30 - (pow(2, 0.5)*(30*pow(2, 0.5)-distance[i]) + distance[j])
                burnArea -= ((spread_prob[i] + spread_prob[j])/2)*(pow(A_side, 2)/2)

    # Equation 8: Subtract overlap from non-adjacent diagonal and orthogonal  neighbors (0-4, 0-6, 2-3, 2-5, 5-1, 5-4, 7-1, 7-3)     
    for i, j in [0,4], [0,6], [2,3], [2,5], [5,1], [5,4], [7,1], [7,3]:
        if ((pow(2, 0.5)*distance[i] + distance[j] > 30) & (distance[i]!=0) & (distance[j]!=0)):
            # Equation 8.1: Subtract overlap from non-adjacent orthogonal and diagonal neighbors when the diagonal distance ≤ half 
            if((distance[i] <= 15*pow(2, 0.5)) & (pow(2, 0.5)*distance[i] + distance[j] > 30)):
                burnArea -= ((spread_prob[i] + spread_prob[j])/2)*(pow((30 - (pow(2, 0.5)*distance[i] + distance[j])), 2)/2)
            # Equation 8.2: Subtract overlap from non-adjacent orthogonal and diagonal neighbors when the diagonal distance ≥ half and does not fill the entire area
            elif((distance[i] > 15*pow(2, 0.5)) & (pow(2, 0.5)*(30*pow(2, 0.5)-distance[i])+distance[j] < 30)):
                A_inv = 900 - 30*distance[j] 
                B_inv = pow((30*pow(2, 0.5)-distance[i]), 2)
                unburned = pow(((pow(2, 0.5)*(30*pow(2, 0.5)-distance[i]))-distance[j]), 2)
                burnArea -= ((spread_prob[i] + spread_prob[j])/2)*(900 - B_inv - A_inv + unburned)
            # Equation 8.3: Subtract overlap from non-adjacent orthogonal and diagonal neighbors when the diagonal distance ≥ half and fills the entire area       
            elif((distance[i] > 15*pow(2, 0.5)) & (pow(2, 0.5)*(30*pow(2, 0.5)-distance[i])+distance[j] > 30)):
                A_inv = 900 - 30*distance[j] 
                B_inv = pow((30*pow(2, 0.5)-distance[i]), 2)
                burnArea -= ((spread_prob[i] + spread_prob[j])/2)*(900 - A_inv - B_inv)
        # Else: Equation 8.0: Non-adjacent orthogonal and diagonal neighbors do not overlap
    
    # Return burn percentage
    return burnArea / 900


# showResults: Plots the map of the area pre-burn and post-burn
    # 'fire' is a 2D matrix containing the percentage on fire of each cell on the map
    # 'veg' is a 2D matrix containing the cell vegetation type as an integer (0 = unburnable, 1 = trees, 2 = shrub, 3 = herb)
def showResults(fire, veg):
    # Build pre-burn landscape image
    rgbIMG = np.zeros([np.size(map, 1), np.size(map, 2), 3], dtype=int)
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
                        lines.Line2D([0], [0], marker='o', color='w', label='Herb', markerfacecolor='#dbeb76', markersize=10)]

    # Plot figures
    fig, ((ax1, ax2)) = plt.subplots(1, 2)
    ax1.imshow(rgbIMG)
    ax2.imshow(newRGB)
    ax2.legend(handles = ax2legendElements, bbox_to_anchor=(1.5, 1.0), loc='upper right')
    plt.show()
    show(fire, cmap='Reds')

    return 0