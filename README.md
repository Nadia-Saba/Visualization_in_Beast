# Visualization_in_Beast
This file contains the code for Multilevel Visualization in Beast .
The code aims at solving the issue faced by multilevel visualization in Beast by aggregating data when it crosses the threshold value
using the object delineation algorithm.

val maxPointsPerTile defines the threshold value. A large value can be assigned for larger datasets.
getPixelOffset() is used to identify the pixel location and rasterizeFeature() method is used for rasterizing the vector data(convert vector data to raster data).
vectorTile() builds the tiles called intermediate tiles based on the threshold. If the value is greater than the threshold set then aggregation happens by
the object delineation algorithm and if the value is lesser then the data points known as features just add one by one on the tile.

To run the code:
Add these methods mentioned above in the Intermediate Tile class in the Beast application and then run the main class.



