################################################################################
# Data generating process 
# Code adapted from Andre's paper without the treatment or spatial aggregation

# Steps
# 1. create a raster grid as values of X 
# 2. define a distribution for X 
# 3. create patterns on the grid 
# 4. create a spatial weight matrix in y 
# 5. generate y 

##############################################################################

library(raster)
library(splm) 
library(tidyverse)
library(reshape2)

# Define extent
lon_min = 0
lat_min = 0
lon_max = 25
lat_max = 25

grid_extent = extent(lon_min, lon_max, lat_min, lat_max)

# Create the raster grid and get xy coordinates 
grid_raster = raster(ext=grid_extent, resolution=1)


# 
values(grid_raster) = 1:ncell(grid_raster)

#Extract cell centre coordinates

x_centres=xFromCol(rast1)
y_centres=yFromRow(rast1)

#Select some random points
random_point_count = 0
random_point_sample_number = 100 #the number of points you want



mat1 = as.matrix(grid_raster)
dat1 = melt(mat1)
names(dat1) = c("easting","northing","value")

#Shift over the easting and northing values by half pixel (ggplot uses coordinates as midpoints, whereas we want it to plot from the mid-point minus half a pixel lenght or width)

dat1$easting = dat1$easting - 0.5
dat1$northing = dat1$northing - 0.5

p = ggplot(dat1, aes(x = easting, y = northing)) +
  geom_tile(aes(fill=value), colour="grey20") +
  scale_fill_gradientn(colours = terrain.colors(10)) +
#  geom_point(data=random_points, mapping = aes(x=EASTING, y=NORTHING), colour="black") +
  labs(x = "Easting", y = "Northing") #+ 
#theme(legend.position = "none")

p

# Y = a + bX  + r WY + e， where a = b = 1 


# Y is land use (binary: 0 if deforestation, 1 if covered )

# X variables: slope, elevation, or soil type 

# Soil type pattern, Moisture level 
# Justify why we use type 1 or type 2 X, specified in the 
# 1. X ~ U (0, 5) (spatial grid pattern) 
# 2. X = (1- piW)K  or X = W^2 K , with K ~ U (0,5) 

# Y is also a function of contiguous parcels W, where our spatial matrix come in to depict the deforestation patterns.





# OLS 
# # coefficients:
# 1. Constant 
# 2. X = 
# 3. Spatial lag 
# 4. X type 
# 5. Disaggregation SD 


# xType Options:
#   "lagged" : lagged X using 0.9
#   "wwx"    : W*W*X
#   "random" : purely random


############################################
# Creating Spatial Weight Matrices:
############################################

# Functions used 


# Spatial weights for neighbours lists
?nb2listw

# generate neighbors based on the type 
?cell2nb(nrow, ncol, type="rook", torus=FALSE)

# For the circle pattern, it would just be a spatial lag thing (i.e. spatial matrix using Queen)


# Create spatial weight matrix W 
size = 120  
p=1
eval(parse(text = paste("Wn", p, "= as(as_dgRMatrix_listw(nb2listw(cell2nb(size, size, type = \"queen\"))), \"CsparseMatrix\")", sep = "")))


# eval(parse(text = paste("Wl", p, "= as(as_dgRMatrix_listw(nb2listw(cell2nb(size, size, type = \"queen\"))))", sep = "")))



# For the fishbone, you might consider dropping fake roads on the map 
# and then having a spatial lag that was not based on a rook or queen's weights matrix
# but one that was more linear

# Y = bX + a 





# generate the X, based on X type 

nRep = 1
size =120
xType = "random"

X = runif(size^2, min = 0, max = 5)   

if (xType != "random"){
  if (xType == "lagged"){
    X = as(powerWeights(Wn, 0.9, X = as.matrix(X)), "matrix")[,1]
  }else{
    X = as(lag.listw(Wl, lag.listw(Wl, X)), "matrix")[,1]
  }
}


# generate the error term 
e = rnorm(size^2, mean = 0, sd = 1)

# generate the Y for different spatial weights in X and Y 


# Y random

trueCoeff = c(1, 1, 0.5, 0)

cat(paste("    Coefficients: const =", trueCoeff[1] ,"| X =", trueCoeff[2], "| T =", trueCoeff[3], ifelse(trueCoeff[4] > 0, paste("| Spatial Lag Y =", trueCoeff[4]), "| No Spatial Lag Y"),"\n"))

Y = powerWeights(Wn1, trueCoeff[4], X = as.matrix(trueCoeff[1] + trueCoeff[2]*X  + e))


#  Y + Lag 0.9
trueCoeff  = c(1, 1, 0.5, 0.9)
cat(paste("    Coefficients: const =", trueCoeff[1] ,"| X =", trueCoeff[2], "| T =", trueCoeff[3], ifelse(trueCoeff[4] > 0, paste("| Spatial Lag Y =", trueCoeff[4]), "| No Spatial Lag Y"),"\n"))
Y = powerWeights(Wn1, trueCoeff[4], X = as.matrix(trueCoeff[1] + trueCoeff[2]*X  + e))



image(Y)
image(Wn1)
image(Wl1)


# Simulation 

# 1000 simulations 

#The results of aggregating the observed data to different resolutions are illustrated by mapping the estimated coefficients and their standard errors against their values at the true level. 

# the vertical axis represents the magnitudes of the coefficient estimates and the horizontal axis represents the resolution chosen 

# Thus, it ranges from the pixel level, up to the highest aggregation level of 12x12 (the lowest resolution). The vertical dashed line denotes the true level, which in the case of Fig 5 case is a 6x6 resolution. 





