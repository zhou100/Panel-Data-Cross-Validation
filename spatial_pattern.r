################################################################################
# Data generating process 
# Code adapted from Andre's paper without the treatment or spatial aggregation

# Steps
# 1. create a raster grid as values of X and T
# 2. Define a distribution for X 
# 3. create patterns on the grid for T 
# 4. create a spatial weight matrix in Y   
# 5. generate y based on the values of X, T, and spatial weight matrix of Y 

##############################################################################

library(raster)
library(splm) 
library(tidyverse)
library(reshape2)

# Define extent

grid_size = 25 

lon_min = 0
lat_min = 0
lon_max = grid_size
lat_max = grid_size

grid_extent = extent(lon_min, lon_max, lat_min, lat_max)

# Create the raster grid and get xy coordinates 
grid_raster = raster(ext=grid_extent, resolution=1)

# fill with random X (such as slope, soil  with some patterns in value)
values(grid_raster) = 1:ncell(grid_raster)

slope_error = runif(grid_size^2, min = 40, max = 100)

# X with error
plot(grid_raster+slope_error)




##############################################
# 1. radial pattern: define T as circular points 
###################################################################

radial_grid = grid_raster 

 
# Consider “Treatment” (human activities) as distance-based,
# with intensity varying with distance from the center. 

xy_center = c(grid_size/2,grid_size/2)

dist_radial <- 100000/distanceFromPoints(radial_grid,xy_center ) 

radial_error = runif(grid_size^2, min = 0, max = 0.3)

plot(dist_radial+radial_error)

   


dist_radial.matrix= as.matrix(dist_radial+radial_error)

# create points with distance to center at given intervals 


##############################################
# 2. fishbone pattern: define T as roads 
###################################################################

fishbone_grid = grid_raster

# set the back ground to 0 
values(fishbone_grid)=0

# choose the index of the main road
grid_range = (grid_size*12+5) :(grid_size*12+23)

# set the main road value to be 1 
values(fishbone_grid)[grid_range] = 1


plot(fishbone_grid)


# randomly choose index where the side roads will emerge 

random_index = sample(grid_range)[1:10]

# for each index, generate north or south road

for (j in 1:length(random_index)){
  
  # randomly choose north or south, i.e. i randomly equals 1 or -1 
  i= floor(runif(1, 0, 2)) * 2 - 1
  
  values(fishbone_grid)[random_index[j]+i*grid_size]=1
  values(fishbone_grid)[random_index[j]+i*grid_size*2]=1
  values(fishbone_grid)[random_index[j]+i*grid_size*3]=1
}
plot(fishbone_grid)
 






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





