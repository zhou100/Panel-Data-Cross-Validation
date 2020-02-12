####      Cross-validation strategies for data with temporal,
####       spatial, hierarchical, or phylogenetic structure
####
####  ------------------------------------------------------------
####
####  David R. Roberts, Volker Bahn, Simone Ciuti, Mark S. Boyce,
####  Jane Elith, Gurutzeta Guillera-Arroita, Severin Hauenstein,
####   Jose J. Lahoz-Monfort, Boris Schröder, Wilfried Thuiller,
####    David I. Warton, Brendan A. Wintle, Florian Hartig,
####                  and Carsten F. Dormann
####
####                       ECOGRAPHY 2017
####
####
####  BOX 1 - Spatial blocking
####  ------------------------
####
####  CONTENTS:
####  ---------
####  1. Simulate the landscapes
####  2. Measure autocorrelation distances
####  3. Modelling and cross-validation
####     3a. Resubstitution
####     3b. Ideal error estimation
####     3c. Random CV
####     3d. Spatially blocked CV
####  4. Spatially buffered leave-one-out CV
####  5. Combine results and plot RMSE for all simulations


####  PACKAGES AND FUNCTIONS

# Random Gaussian fields
#install.packages("RandomFields")
library(RandomFields)

# Random Forests
#install.packages("randomForest")
library(randomForest)

# Parallel processing
#install.packages("doParallel")
library(doParallel)

# RMSE calculation
#install.packages("hydroGOF")
library(hydroGOF)

# Autocorrelation
#install.packages("gstat")
library(gstat)
#install.packages(automap)
library(automap)

# Distance matrices
#install.packages("fields")
library(fields)

# Nearest neighbour searches
#install.packages("yaImpute")
library(yaImpute)

# Data management
#install.packages("plyr")
library(plyr)
#install.packages("reshape")
library(reshape)

# Plotting
#install.packages("ggplot2")
library(ggplot2)


## PART 1 - Simulate the landscapes
## -------------------------------
##
## These landscapes have:
##    - One long range SAC as "temp" with a RMgauss covariance model,
##    - A "precip" variable with a bit more noise and a bit shorter range, and
##    - RMexp and another intermediate RMexp as "noise"


# Set the number of simulated surfaces to run
n.sims <- 25

# Build the base grids
gridsize = c(50L, 50L)
Xvec <- seq(0, 1, len = gridsize[1])
Yvec <- seq(0, 1, len = gridsize[2])
grd <- expand.grid(Y = Yvec, X = Xvec)
lon <- grd$X
lat <- grd$Y

## Simulate landscapes
# LOOP for all simulated surfaces
for(i in 1:n.sims){
  
  ## FIRST, smulate the spatially structured environment
  
  # Some other env var is just standard (could be soil or topography or whatever)
  expCov <- RMexp(var = 0.1, scale = 0.1)
  x.1 <- as.vector(t(RFsimulate(expCov, x = Xvec, y = Yvec, spConform = FALSE)))
  x.1 <- scale(x.1)
  
  # "Precip" has a higher variance
  expCov <- RMexp(var = 0.3, scale = 0.1)
  x.2 <- as.vector(t(RFsimulate(expCov, x = Xvec, y = Yvec, spConform = FALSE)))
  x.2 <- scale(x.2)
  
  # "Temp" has a longer range
  expCov <- RMgauss(var = 0.1, scale = 0.3) #0.4 was a bit much
  x.3 <- as.vector(t(RFsimulate(expCov, x = Xvec, y = Yvec, spConform = FALSE)))
  x.3 <- scale(x.3)
  
  
  ## SECOND, simulate the biotic interaction and the disease

  # Prevalence of disease
  p.d <- 0.05
  
  # Add constant to x.2 and x.3 to keep them positive
  t.ps <- -floor(min(x.2, x.3))
  
  # Set up disease with 1 where it occurs
  x.4 <- rep(1, prod(gridsize))
  
  # Change locations to zero where the ratio of precip to temp isn't large enough to be 
  # in the percentile specified by prevalence (the top wet and warm locations have disease)
  x.4[((x.2 + t.ps)/(x.3 + t.ps)) < quantile(((x.2 + t.ps)/(x.3 + t.ps)), (1 - p.d))] <- 0
  
  # Putting in "food" as a linear combination of previous variables
  x.5 <- (x.1 + x.2 + x.3 + x.2*x.3)
  x.5 <- scale(x.5)

  #  Add more "normal" variables as additional material for overfitting
  expCov <- RMexp(var = 0.1, scale = 0.1)
  x.6 <- as.vector(t(RFsimulate(expCov, x = Xvec, y = Yvec, spConform = FALSE)))
  x.6 <- scale(x.6)
  x.7 <- as.vector(t(RFsimulate(expCov, x = Xvec, y = Yvec, spConform = FALSE)))
  x.7 <- scale(x.7)
  x.8 <- as.vector(t(RFsimulate(expCov, x = Xvec, y = Yvec, spConform = FALSE)))
  x.8 <- scale(x.8)
  expCov <- RMgauss(var = 0.1, scale = 0.3)
  x.9 <- as.vector(t(RFsimulate(expCov, x = Xvec, y = Yvec, spConform = FALSE)))
  x.9 <- scale(x.9)
  expCov <- RMgauss(var = 0.1, scale = 0.3)
  x.10 <- as.vector(t(RFsimulate(expCov, x = Xvec, y = Yvec, spConform = FALSE)))
  x.10 <- scale(x.10)
  
  
  ## THIRD, smulate the species
  
  # Species should depend linearly on food (x.5).
  # Species should be limited by an interaction between precip and temp.
  # In reality, water potential depends linearly on temp, but temp is in K, so the
  # real variance is small and relatively unimportant compared to precip. I will simulate that
  # here by turning the simulated temp into somewhat realistic temp and then doing calcs (precip/temp)
  # as if it were in K (so add 273 first and in the end standardize the variable again)
  # dependence on temp is unimodal f(x) = 1/(sqrt(2*pi)*sigma)*e^(-((x - mean)^2/(2*sigma^2)))
  # I'll have sigma = 1 here and mean = 0
  
  # "water availability"
  x.11 <- x.2/(x.3 + 273)
  x.11 <- scale(x.11)
  
  # Gaussian dependence on temperature
  x.12 <- 1/(sqrt(2*pi))*exp(-(x.3^2/4))
  x.12 <- scale(x.12)
  
  # Gaussian dependence on water
  x.13 <- 1/(sqrt(2*pi))*exp(-(x.2^2/4))
  x.13 <- scale(x.13)
  
  # x.1 = some standard var such as soil - unknown to model
  # x.5 = "food" a combination of x.1, x.2, x.3 and x.2*x.3
  # x.12 and x.13, Gaussian response of temp and precip
  y <- x.1 + x.5 + x.12 + x.13 + x.6
  y <- scale(y)
  
  # Then use the water potential as limiting factor by reducing y to water potential where 
  # water potential is lower (x.11, "water potential" derived from x.2 and x.3)
  y[y>x.11] <- x.11[y>x.11]
  y[x.4 == 1] <- min(y)
  
  # Empty list for simulated landscapes (created on first loop)
  if(i==1){sim.data <- list()}
  sim.data[[i]] <- data.frame(y, mget(paste0("x.", 1:13)))

} # LOOP for all simulated surfaces 

# Save the simulated data
save(sim.data, file="Complete simulated data.Rdata")
#load("Complete simulated data.Rdata")

# Create complete table of simulated data (for ideal validation)
for(i in 1:n.sims){
  if(i==1){complete.data <- data.frame(SIM=i, sim.data[[1]])}
  complete.data <- rbind(complete.data, data.frame(SIM=i,sim.data[[i]]))
}



## PART 2 - Measure autocorrelation
## --------------------------------

# LOOP for each simulated landscape
for(i in 1:n.sims){

  # Make empty data frame for SAC results (first loop only)
  if(i==1){sac.ranges <- data.frame(sac.dist.lm = rep(NA, n.sims), sac.dist.rf = NA)}
  
  # Define the data
  data.used <- sim.data[[i]]
  
  # Build the full data models
  mod.1 <- lm(y ~ x.2 + x.3 + x.6 + x.7 + x.8 + x.9 + x.10, data = data.used)
  mod.2 <- randomForest(y ~ x.2 + x.3 + x.6 + x.7 + x.8 + x.9 + x.10, data = data.used)
  
  # Turn residuals and coords into a "SpatialPointsDataFrame" (sp)
  # LM
  res.1 <- predict(mod.1, data.used, type="response") - data.used$y
  resid.spatial <- SpatialPointsDataFrame(grd, data.frame(resids=res.1))
  vari.foo <- autofitVariogram(resids ~ 1, resid.spatial, model=c("Sph"))
  sac.ranges[i,1] <- vari.foo$var_model$range[2]
  # RF
  res.2 <- predict(mod.2, data.used, type="response") - data.used$y
  resid.spatial <- SpatialPointsDataFrame(grd, data.frame(resids = predict(mod.2)-data.used$y))
  vari.foo <- autofitVariogram(resids ~ 1, resid.spatial, model=c("Sph"), verbose=F)
  sac.ranges[i,2] <- vari.foo$var_model$range[2]
  
} # LOOP for each simulated landscape

# Convert distances to actual cells
sac.ranges <- sac.ranges*50

# Check the autocorrelation distances
summary(sac.ranges)

# Save autocorrelation distances
write.csv(sac.ranges, "Results/Residual spatial autocorrelation ranges for LM and RM models.csv", row.names=F)

# Plot autocorrelation distances
png("Results/Residual autocorrelation plots.png")
par(mfrow=c(1,2), mar=c(3,4,0.5,0.5))
layout(matrix(c(1,2), 2), widths=c(1,1), heights=c(1.5,2.5))
layout.show(2)
boxplot(sac.ranges, horizontal=T, names=c("GLM","RF"), xlab="", las=1, ylim=range(sac.ranges,na.rm=T), col=c("pink","lightblue"))
par(mar=c(4,4,0.5,0.5))
plot(density(sac.ranges[,1], na.rm=T), main="", xlab="Semivariance", las=1, xlim=range(sac.ranges,na.rm=T), col="pink", lwd=2)
lines(density(sac.ranges[,2], na.rm=T), lty=1, col="lightblue", lwd=2)
legend("topright", col=c("lightblue","pink"), lty=1, lwd=2, legend=c("RF","LM"))
dev.off()



## PART 3 - Modelling and cross-validation
## ---------------------------------------
##
## Missing are the "food" x.5 that is directly a linear combination of the first 
## three x and thus highly related to them, x.1 that is just a low cor noise
## and x.4, which is "disease".

# Define spatial arrangements for blocked CV
block.arrange <- data.frame(x=c(5,10,15,20,25,50),y=c(5,10,15,20,25,25))

## Parallel processing
## This is a Revolution R library for parallelizable loops and the back end that can 
## handle parallel processes - also calls library(foreach)
# Set up parallel backend to use multiple processors
# Need to finish with stopCluster(cl) below
cl<-makeCluster(4)
registerDoParallel(cl)

# Produce an output vector in each loop (rbind with dopar)
cv.results <- foreach(i=1:n.sims, .combine=rbind) %dopar% {
  
  # Packages for paralleisation
  library(randomForest)
  library(hydroGOF)

  # Define simulated data set to use
  data.used <- sim.data[[i]]

  # Full data models
  mod.1 <- lm(y ~ x.2 + x.3 + x.6 + x.7 + x.8 + x.9 + x.10, data = data.used)
  mod.2 <- randomForest(y=data.used$y, x=data.used[, c(3,4,7,8,9,10,11)], ntree=250, nodesize=10)

  
  ######### PART 3a. Resubstitution all data
  
  # RMSE
  resub.1 <- rmse(predict(mod.1, data.used), data.used$y)
  resub.2 <- rmse(predict(mod.2, data.used), data.used$y)
  
  
  ######### PART 3b. Ideal = predicted to new landscapes
  
  # Make prediction data set comprised of all other simulated datasets
  ideal.pred <- complete.data[complete.data$SIM!=i,]
  # Identify no-analogs
  foo <- rbind(data.used[,c("x.2","x.3","x.6")], ideal.pred[,c("x.2","x.3","x.6")])
  test.na <- rowSums(sapply(foo, function(x) 
    x[(prod(gridsize)+1):(prod(gridsize)+nrow(ideal.pred))] < min(x[1:prod(gridsize)]) | 
      x[(prod(gridsize)+1):(prod(gridsize)+nrow(ideal.pred))] > max(x[1:prod(gridsize)])))!=0
  # RMSE
  ideal.1 <- c(rmse(predict(mod.1, ideal.pred), ideal.pred$y),                       # All data
               rmse(predict(mod.1, ideal.pred[!test.na,]), ideal.pred$y[!test.na]))  # No-analogs removed
  ideal.2 <- c(rmse(predict(mod.2, ideal.pred), ideal.pred$y),                       # All data
               rmse(predict(mod.2, ideal.pred[!test.na,]), ideal.pred$y[!test.na]))  # No-analogs removed

    
  ######### PART 3c. Random hold-out (50%)
  
  # FOLD 1 (first half)
  # Define training and testing data
  train <- rep(T, (prod(gridsize)))
  train[sample(1:(prod(gridsize)), (prod(gridsize))/2)] <- F
  test <- !train
  # Identify no-analogs
  test.na <- rowSums(sapply(data.used[,c(3,4,7)], function(x) x[test] < min(x[train]) | x[test] > max(x[train])))!=0
  # Build models
  mod.1 <- lm(y ~ x.2 + x.3 + x.6 + x.7 + x.8 + x.9 + x.10, data = data.used[train,])
  mod.2 <- randomForest(y=data.used$y[train], x=data.used[train, c(3,4,7,8,9,10,11)], ntree=250, nodesize=10)
  # RMSE
  foo.1 <- c(rmse(predict(mod.1, data.used[test,]), data.used$y[test]),                      # All data
             rmse(predict(mod.1, data.used[test,][!test.na,]), data.used$y[test][!test.na])) # No-analogs removed
  foo.2 <- c(rmse(predict(mod.2, data.used[test,]), data.used$y[test]),                      # All data
             rmse(predict(mod.2, data.used[test,][!test.na,]), data.used$y[test][!test.na])) # No-analogs removed
  
  # FOLD 2 (second half)
  # Define training and testing data
  test <- train
  train <- !train
  test.na <- rowSums(sapply(data.used[,c(3,4,7)], function(x) x[test] < min(x[train]) | x[test] > max(x[train])))!=0
  # Build models
  mod.1 <- lm(y ~ x.2 + x.3 + x.6 + x.7 + x.8 + x.9 + x.10, data = data.used[train,])
  mod.2 <- randomForest(y=data.used$y[train], x=data.used[train, c(3,4,7,8,9,10,11)], ntree=250, nodesize=10)
  # RMSE
  foo.3 <- c(rmse(predict(mod.1, data.used[test,]), data.used$y[test]),                      # All data
             rmse(predict(mod.1, data.used[test,][!test.na,]), data.used$y[test][!test.na])) # No-analogs removed
  foo.4 <- c(rmse(predict(mod.2, data.used[test,]), data.used$y[test]),                      # All data
             rmse(predict(mod.2, data.used[test,][!test.na,]), data.used$y[test][!test.na])) # No-analogs removed  
 
   # Populate results vector with averages of two halves
  random.1 <- c(mean(c(foo.1[1], foo.3[1])), mean(c(foo.1[2], foo.3[2])))
  random.2 <- c(mean(c(foo.2[1], foo.4[1])), mean(c(foo.2[2], foo.4[2])))

  
  ######### PART 3d. Spatial blocking
  
  ## LOOP for all blocking arrangements
  for(b in 1:nrow(block.arrange)){
    
    # Make empty block CV output table (first loop only)
    if(b==1){
      block.1 <- data.frame(NoAnalogs=c(T,F), matrix(NA,nrow=2,ncol=nrow(block.arrange)))
      names(block.1) <- c("NoAnalogs", apply(block.arrange, 1, function(x){paste0("RMSE.",x[1],"x",x[2])}))
      block.2 <- block.1
    }
    
    # Set the size of blocks
    cvblock.size <- c(block.arrange[b,"x"],block.arrange[b,"y"])
    dim.1 <- ceiling(gridsize[1]/cvblock.size[1])
    dim.2 <- ceiling(gridsize[2]/cvblock.size[2])
    dim.3 <- ifelse(as.logical(dim.1%%2), dim.1, dim.1 + 1)
    
    # Block matrix
    m.init <- matrix(c(T, F), dim.3, dim.2)
    if(dim.1 != dim.3)  m.init <- m.init[-(dim.1+1), ]
    
    # Expand to real size of matrix, then cut back edges
    m.half <- matrix(rep(as.vector(m.init), each=cvblock.size[1]), nrow=dim.1*cvblock.size[1], ncol=dim.2)
    m.full <- matrix(rep(as.vector(t(m.half)), each=cvblock.size[2]), ncol=dim.1*cvblock.size[1], nrow=dim.2*cvblock.size[2])
    # Cut back to size
    m <- t(m.full)[1:gridsize[1], 1:gridsize[2]]
    
    # Define training & testing data
    train <- as.vector(m)
    test <- !train
    
    # Visual check
    #image(1:gridsize[1], 1:gridsize[2], matrix(train, gridsize[1]), col=heat.colors(n.sims), las=1, asp=1); box()
    
    # FOLD 1
    # Determine which test data are no-analogs
    test.na <- rowSums(sapply(data.used[,c(3,4,7)], function(x) x[test] < min(x[train]) | x[test] > max(x[train])))!=0
    # Build models
    mod.1 <- lm(y ~ x.2 + x.3 + x.6 + x.7 + x.8 + x.9 + x.10, data = data.used[train,])
    mod.2 <- randomForest(y=data.used$y[train], x=data.used[train, c(3,4,7,8,9,10,11)], ntree=250, nodesize=10)
    # RMSE
    foo.1 <- c(rmse(predict(mod.1, data.used[test,]), data.used$y[test]),                      # All data
               rmse(predict(mod.1, data.used[test,][!test.na,]), data.used$y[test][!test.na])) # No-analogs removed
    foo.2 <- c(rmse(predict(mod.2, data.used[test,]), data.used$y[test]),                      # All data
               rmse(predict(mod.2, data.used[test,][!test.na,]), data.used$y[test][!test.na])) # No-analogs removed  
    
    # FOLD 2
    # Re-define training and testing data
    test <- train
    train <- !train
    # Determine which test data are no-analogs
    test.na <- rowSums(sapply(data.used[,c(3,4,7)], function(x) x[test] < min(x[train]) | x[test] > max(x[train])))!=0
    # Build models
    mod.1 <- lm(y ~ x.2 + x.3 + x.6 + x.7 + x.8 + x.9 + x.10, data = data.used[train,])
    mod.2 <- randomForest(y=data.used$y[train], x=data.used[train, c(3,4,7,8,9,10,11)], ntree=250, nodesize=10)
    # RMSE
    foo.3 <- c(rmse(predict(mod.1, data.used[test,]), data.used$y[test]),                      # All data
               rmse(predict(mod.1, data.used[test,][!test.na,]), data.used$y[test][!test.na])) # No-analogs removed
    foo.4 <- c(rmse(predict(mod.2, data.used[test,]), data.used$y[test]),                      # All data
               rmse(predict(mod.2, data.used[test,][!test.na,]), data.used$y[test][!test.na])) # No-analogs removed  
    
    # Populate results table
    block.1[block.1$NoAnalogs==T,paste0("RMSE.",cvblock.size[1],"x",cvblock.size[2])] <-  # LM with all data
      mean(c(foo.1[1], foo.3[1]))
    block.1[block.1$NoAnalogs==F,paste0("RMSE.",cvblock.size[1],"x",cvblock.size[2])] <-  # LM with no-analogs removed
      mean(c(foo.1[2], foo.3[2]))
    block.2[block.1$NoAnalogs==T,paste0("RMSE.",cvblock.size[1],"x",cvblock.size[2])] <-  # RF with all data
      mean(c(foo.2[1], foo.4[1]))
    block.2[block.1$NoAnalogs==F,paste0("RMSE.",cvblock.size[1],"x",cvblock.size[2])] <-  # RF with no-analogs removed
      mean(c(foo.2[2], foo.4[2]))
    
  } # LOOP for all block arrangements
  
  cv.out.1 <- data.frame(SIM=i, Model="LM",
                         RMSE.ideal=ideal.1,       # Ideal RMSE (w/ and w/o no-analogues)
                         RMSE.resub=c(resub.1,NA), # Never no-analogs in resubstitution
                         RMSE.random=random.1,     # Random CV (w/ and w/o no-analogues)
                         block.1)                  # Spatial block CV RMSE (w/ and w/o no-analogues)
  cv.out.2 <- data.frame(SIM=i, Model="RF",
                         RMSE.ideal=ideal.2,       # Ideal RMSE (w/ and w/o no-analogues)
                         RMSE.resub=c(resub.2,NA), # Never no-analogs in resubstitution
                         RMSE.random=random.2,     # Random CV RMSE (w/ and w/o no-analogues)
                         block.2)                  # Spatial block CV RMSE (w/ and w/o no-analogues)
  rbind(cv.out.1[c(1,2,6,3:5,7:12)],cv.out.2[c(1,2,6,3:5,7:12)])

} # Parallel LOOP for all simulated surfaces

# Finish the multiple clusters for parallel processing
stopCluster(cl)

# Save the results output
dir.create("Results")
write.csv(cv.results, "Results/Results - RMSE for all simulations (resub, random, block).csv", row.names=F)
#cv.results <- read.csv("Results/Results - RMSE for all simulations (resub, random, block).csv")



## PART 4 - Evaluation with spatially buffered leave-one-out
## ---------------------------------------------------------
##
## This spatially buffered LOO approach is slightly different from the one used in other boxes.
## Results, however, should be consistent between methods.

# List of LOO buffer sizes to run
buffer.list <- c(5,8,10)

# Number of centers (LOO test sites) in x and y dimension (in that sequence)
n.centers <- c(10,10)

# General way of creating a sequence of optimally distributed centers
# foo.x <- seq(1, gridsize[1], round(gridsize[1]/n.centers[1])) + ceiling((gridsize[1]/n.centers[1]-1)/2)/50
# foo.y <- seq(1, gridsize[2], gridsize[2]/n.centers[2]) + ceiling((gridsize[2]/n.centers[2]-1)/2)

# Locate the centers of the LOO
x.centers <- seq(0,1,length.out=n.centers[1]+1)[2:(n.centers[1]+1)] - 0.5/n.centers[1]
y.centers <- seq(0,1,length.out=n.centers[2]+1)[2:(n.centers[2]+1)] - 0.5/n.centers[2]

# In case XY not values between 0 and 1
foo.x <- as.numeric(quantile(Xvec, x.centers))
foo.y <- as.numeric(quantile(Yvec, y.centers))

# Make coordinates of XY centers
centers <- expand.grid(Y=foo.y, X=foo.x)

# Assign grid cells as centers (or not)
loo.select <- grd
loo.select$center <- F
loo.select[ann(as.matrix(grd[c("X","Y")]), as.matrix(centers[c("X","Y")]), k=1)$knnIndexDist[,1],"center"] <- T

# Quick look a the LOO centers
#plot(Y~X, data=loo.select, pch=15, asp=1, col=c("black","red")[as.factor(loo.select$center)])
# On a simulated surface
#image(Xvec, Yvec, matrix(y,50), asp=1, col=heat.colors(n.sims))
#points(Y~X, data=loo.select[loo.select$center==TRUE,], pch=15, cex=0.5)

#setup parallel backend to use multiple processors
# Need to finish with stopCluster(cl) below
cl <- makeCluster(4)
registerDoParallel(cl)

# LOOP for all simulations
loo.results <- foreach(i=1:n.sims, .combine=rbind) %dopar% {
  
  require(randomForest)
  library(fields)
  
  # Define the simulation data to use
  data.used <- sim.data[[i]]
  
  # LOOP for all buffer radii
  for(dist.buffer in buffer.list){
    
    # Empty result tables for each buffer size (first loop only)
    # Need one table for ALL DATA and one table for ANALOGS ONLY
    if(dist.buffer==buffer.list[1]){
      # All Data
      loo.results.buffers.1 <- data.frame(Model=c("LM","RF"), matrix(rep(NA,length(buffer.list)), nrow=1))
      names(loo.results.buffers.1) <- c("Model",paste0("RMSE.LOO.",buffer.list))
      # No-analogs removed
      loo.results.buffers.2 <- loo.results.buffers.1
    }
    
    # Scale the buffer to the grid size
    scaled.buffer <- dist.buffer*(Xvec[2]-Xvec[1])
    
    # LOOP for all LOO centers (test points)
    for(j in 1:(n.centers[1]*n.centers[2])){
      
      # Empty results table for each center (first loop only)
      if(j==1){
        loo.results.centers <- data.frame(rep(NA, 100),  NA)
        names(loo.results.centers) <- c("LM.error","RF.error")
        loo.results.centers$NoAnalog <- F
      }
      
      # Define the LOO center (test point)
      coord.center <- which(loo.select$center)[j]
      
      # Create distance matrix to center
      foo.dist <- rdist(loo.select[,1:2], loo.select[coord.center,1:2])
      
      # Define training and testing data
      train <- as.vector(foo.dist) > scaled.buffer
      test <- rep(F, gridsize[1]*gridsize[2])
      test[which(loo.select$center)[j]] <- T
      
      # Determine and record whether the center is no-analog
      if(sum(sapply(data.used[,c(3,4,7)], function(x) x[test] < min(x[train]) | x[test] > max(x[train])))!=0) loo.results.centers[j,"NoAnalog"] <- T
      
      # Build models
      mod.1 <- lm(y ~ x.2 + x.3 + x.6 + x.7 + x.8 + x.9 + x.10, data = data.used[train,])
      mod.2 <- randomForest(y=data.used$y[train], x=data.used[train,c(3,4,7,8,9,10,11)], ntree=250, nodesize=10)
      
      # Test prediction to center
      # Write residuals to temporary results table
      loo.results.centers[j,"LM.error"] <-  predict(mod.1, data.used[test,]) - data.used$y[test]
      loo.results.centers[j,"RF.error"] <-  predict(mod.2, data.used[test,]) - data.used$y[test]
      
    } # LOOP for all LOO centers (j)
    
    # Populate the LOO results table with RMSE from all centers
    # All data
    loo.results.buffers.1[loo.results.buffers.1$Model=="LM",paste0("RMSE.LOO.",dist.buffer)] <- 
      sqrt(mean(loo.results.centers$LM.error^2))
    loo.results.buffers.1[loo.results.buffers.1$Model=="RF",paste0("RMSE.LOO.",dist.buffer)] <- 
      sqrt(mean(loo.results.centers$RF.error^2))
    # No-analogs removed
    loo.results.buffers.2[loo.results.buffers.2$Model=="LM",paste0("RMSE.LOO.",dist.buffer)] <- 
      sqrt(mean(loo.results.centers[loo.results.centers$NoAnalog==F,"LM.error"]^2))
    loo.results.buffers.2[loo.results.buffers.2$Model=="RF",paste0("RMSE.LOO.",dist.buffer)] <- 
      sqrt(mean(loo.results.centers[loo.results.centers$NoAnalog==F,"RF.error"]^2))
    
  } # LOOP for all buffer radii (dist.buffer)
  
  loo.1 <- data.frame(SIM=i, NoAnalogs=T, loo.results.buffers.1) # Summary with all data
  loo.2 <- data.frame(SIM=i, NoAnalogs=F, loo.results.buffers.2) # Summary with no-analogs removed
  rbind(loo.1,loo.2)
  
} # LOOP for all simulated surfaces (i)
stopCluster(cl)

# Save the LOO results output
write.csv(loo.results, "Results/Results - RMSE for all simulations (loo).csv", row.names=F)
#loo.results <- read.csv("Results/Results - RMSE for all simulations (loo).csv")



## PART 5 - Combine results and plot RMSE for all simulations
## ----------------------------------------------------------

alldat <- merge(cv.results, loo.results, by=c("SIM","Model","NoAnalogs"), all.x=T, all.y=T, sort=F)
alldat

# Reshape data
plotdat <- melt(alldat, id.vars=c("SIM","Model","NoAnalogs"))
names(plotdat)[4:5] <- c("CV","RMSE")
plotdat$CVtype <- ifelse(plotdat$CV=="RMSE.ideal","Ideal",
                    ifelse(plotdat$CV=="RMSE.resub","Resub",
                      ifelse(plotdat$CV=="RMSE.random","Random",
                        ifelse(plotdat$CV %in% apply(block.arrange, 1, function(x){paste0("RMSE.",x[1],"x",x[2])}),"Spatial",
                          ifelse(plotdat$CV %in% paste0("RMSE.LOO.",buffer.list),"LOO",NA)))))
plotdat$Data <- ifelse(plotdat$NoAnalogs,"All data","No-analogues removed")

# Table of averages
CV.means <- ddply(plotdat, .(CV, CVtype, NoAnalogs, Data), summarise, AVG=mean(RMSE), MED=median(RMSE))

# Boxplot
png("Results/Boxplots of RMSE by CV.png", w=640)
ggplot(plotdat, aes(y=RMSE, x=CV, fill=CVtype)) + 
  geom_boxplot() +
  coord_flip() +
  facet_grid(Data~.)
dev.off()

# Denisty plot
png("Results/Density plots of RMSE by CV.png", width=640)
ggplot(plotdat, aes(x=RMSE, col=CV)) + 
  geom_vline(data=CV.means, aes(xintercept=MED, col=CV), lty=2) +
  geom_density() + 
  facet_grid(CVtype~Data)
dev.off()




