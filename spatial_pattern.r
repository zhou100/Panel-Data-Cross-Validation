Andre paper 

Data generating process 
Do we need to consider treatment? No 


http://www.econ.uiuc.edu/~lab/workshop/Spatial_in_R.html#running-spatial-regressions


Various degrees of complexity 

Y = a + bX  + r WY + e

a=b= 1 


Y is land use (binary ?) 0 if deforestation , 1 if covered 

X variables: slope, elevation, or soil type 

Soil type pattern ? 
  
  Moisture level 

Justify why we use spatially lagged X 

1. X ~ U (0, 5) (spatial grid pattern) 
2. X = (1- piW)K  or X = W^2 K , with K ~ U (0,5) 



Y is also a function of contiguous parcels W 

Clustering 


# OLS 
# coefficients:
1. Constant 
2. X = 
  3. Spatial lag 
4. X type 
5. Disaggregation SD 


# xType Options:
#   "lagged" : lagged X using 0.9
#   "wwx"    : W*W*X
#   "random" : purely random



# Creating Spatial Weight Matrices:


aggList = c(1, 2, 4, 5, 6, 10, 12)
disSteps. = agglist 
Just set a s 1 

Different level of aggregation 

size = (120 / disSteps[p])



# Spatial weights for neighbours lists
nb2listw

# generate neighbors based on the type 
cell2nb(nrow, ncol, type="rook", torus=FALSE)



# Create spatial weight matrix W 

size = 120  
p=1
eval(parse(text = paste("Wn", p, "= as(as_dgRMatrix_listw(nb2listw(cell2nb(size, size, type = \"queen\"))), \"CsparseMatrix\")", sep = "")))


Queen:   For the circle pattern, it would just be a spatial lag thing.  



For the fishbone, you might consider dropping fake roads on the map and then having a spatial lag that was not based on a rook or queen's weights matrix, but one that was more linear

Y = bX + a 

eval(parse(text = paste("Wl", p, "= nb2listw(cell2nb(size, size, type = \"queen\"))", sep = "")))




nRep = 2



# generate the X 


for (k in 1:nRep){

X = runif(size^2, min = 0, max = 5)   

if (xType != "random"){
if (xType == "lagged"){
X = as(powerWeights(Wn, 0.9, X = as.matrix(X)), "matrix")[,1]
}else{
X = as(lag.listw(Wl, lag.listw(Wl, X)), "matrix")[,1]
}
}




#. error term 
e = rnorm(size^2, mean = 0, sd = 1)

# generate the Y 

Y = powerWeights(Wn, trueCoeff[4], X = as.matrix(trueCoeff[1] + trueCoeff[2]*X  + e))





Simulation 

1000 simulations 

The results of aggregating the observed data to different resolutions are illustrated by mapping the estimated coefficients and their standard errors against their values at the true level. 

the vertical axis represents the magnitudes of the coefficient estimates and the horizontal axis represents the resolution chosen 

Thus, it ranges from the pixel level, up to the highest aggregation level of 12x12 (the lowest resolution). The vertical dashed line denotes the true level, which in the case of Fig 5 case is a 6x6 resolution. 




X = 

