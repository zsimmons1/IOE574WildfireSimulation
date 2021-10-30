# Importing GIS data to setup a grid
# Define all cells according to vegetation type

# Storing all information related to each individual cell

import os 
import matplotlib.pyplot as plt
import geopandas as gpd
import earthpy as et

# Set home directory and download data
et.data.get_data("spatial-vector-lidar")
os.chdir(os.path.join(et.io.HOME, 'earth-analytics'))

# Import shapefile using geopandas
sjer_plot_locations = gpd.read_file('data/spatial-vector-lidar/california/neon-sjer-site/vector_data/SJER_plot_centroids.shp')

sjer_plot_locations.head(6)