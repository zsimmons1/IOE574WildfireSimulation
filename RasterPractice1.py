import rasterio
from rasterio.plot import show
from matplotlib import pyplot as plt

img = rasterio.open('/Users/Zack/Desktop/IOE574/TermProject/IOE574WildfireSimulation/us_210evc.tif')
show(img)

full_img = img.read()

# To find the number of bands in an image
num_bands = img.count
print("Number of bands in the image = ", num_bands)

img_band1 = img.read(1)

fig = plt.figure(figsize=(10,10))
ax1 = fig.add_subplot(2,2,1)
ax1.imshow(img_band1, cmap='pink')

# To find out the coordinate reference system
print("Coordinate reference system: ", img.crs)

# Read metadata
metadata = img.meta
print('Metadata: {metadata}\n'.format(metadata=metadata))

# Read description, if any
desc = img.descriptions
print('Raster description: {desc}\n'.format(desc=desc))

# To find out geo transform
print("Geotransform: ", img.transform)
