import rasterio
from rasterio.plot import show
from matplotlib import pyplot as plt
import numpy as np

img = rasterio.open('/Users/Zack/Desktop/IOE574/TermProject/IOE574WildfireSimulation/us_210evc.tif')
show(img)

full_img = img.read()

# # To find the number of bands in an image
# num_bands = img.count
# print("Number of bands in the image = ", num_bands)

# img_band1 = img.read(1)

# fig = plt.figure(figsize=(10,10))
# ax1 = fig.add_subplot(2,2,1)
# ax1.imshow(img_band1, cmap='pink')

cellValue = full_img[0][0][0]
print(cellValue)

check = np.less_equal(full_img, 20)  # Create check raster with True/False values
water = np.where(check == 1, 1, 0) # Set values above 0 as water and otherwise leave it at 0
show(water, cmap='Blues')

fire = np.zeros(full_img.shape, dtype=rasterio.float32)  # Create empty matrix
# start fire in upper right corner:
fire[0][0][48] = 100

show(fire, cmap='Reds')

print(fire[0][56][0])
print(fire[0][3][5])
