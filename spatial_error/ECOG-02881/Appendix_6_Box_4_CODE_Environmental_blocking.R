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
####  BOX 4 - Blocking for extrapolation
####  ----------------------------------
####
####  CONTENTS:
####  ---------
####  1. Import the data and set up the models
####  2. Create the full data models
####  3. Measure the spatial autocorrelation
####  4. Run the cross validations
####       a. Resubstitution CV
####       b. Random k-fold CV
####       c. Spatially blocked k-fold CV
####       d. Environmentall blocked k-fold CV with PCA cuts
####       e. Environmentall blocked k-fold CV with kmeans clustering
####       f. Spatially buffered leave-one-out CV
####   5. Compile the CV results (tables and plots)
####   6. Measure distances between model training and testing data
####       a. Geographic distances
####       b. Environmental distances
####   7. Compare distances (geog and env) to error estimates



####  PACKAGES

# Stepwise model selection
#install.packages('MASS')
library(MASS)

# Random forests
#install.packages('randomForest')
library(randomForest)

# AUC calculation
#install.packages('ROCR')
library(ROCR)

# Approximate nearest neighbour computation
#install.packages('yaImpute')
library(yaImpute)

# Spatial autocorrelation
#install.packages('ncf')
library(ncf)
#install.packages('geoR')
library(geoR)

# Data management
#install.packages('plyr')
library(plyr)
#install.packages('SDMTools')
library(SDMTools)

# Plotting
#install.packages('ggplot2')
library('ggplot2')
#install.packages('gridExtra')
library(gridExtra)

#### FUNCTIONS

# Function to calculate quick AUC
auc.f <- function(pred, obs){
  
  # pred: vector of predicted values
  # obs: vector of observed (true) values
  
  # Calculate the AUC
  library(ROCR)
  ROC_perf <- performance(prediction(pred,obs),"tpr","fpr")
  ROC_sens <- performance(prediction(pred,obs),"sens","spec")
  ROC_auc <- performance(prediction(pred,obs),"auc")
  AUC <- ROC_auc@y.values[[1]] # AUC
  
  # Mean sensitivity across all thresholds
  x.Sens <- mean(as.data.frame(ROC_sens@y.values)[,1])
  # Mean specificity across all thresholds
  x.Spec <- mean(as.data.frame(ROC_sens@x.values)[,1])
  
  # Create output table
  cbind(AUC, x.Sens, x.Spec)
}



## PART 1 - Import the data and set up the models
## ----------------------------------------------

# Import the North America Douglas-fir survey data (contains one data frame named 'tr')
load("Appendix_6_Box_4_DATA_NorthAmerica_DougFir.RData")

##  These data were compiled by:
##
##    Rehfeldt, G. E., Jaquish, B. C., Lopez-Upton, J., Saenz-Romero, C.,
##    St Clair, J. B., Leites, L. P., & Joyce, D. G. 2014. Comparative
##    genetic responses to climate for the varieties of Pinus ponderosa 
##    and Pseudotsuga menziesii: Realized climate niches. Forest Ecology
##    and Management, 324, 126-137.
##    http://dx.doi.org/10.1016/j.foreco.2014.02.035
##
##  These data are provided here with permission from G. E. Rehfeldt.
##  Many thanks to them for sharing their hard work.
##
##  Data have been trimmed to west of -95 degrees Longitude.
##  XY coordinates are in Lambert Conformal Conic projection (consistent spatial grids):
##
##  Projection    LAMBERT
##  Datum         WGS84
##  Spheroid      WGS84
##  Units         METERS
##  Zunits        NO
##  Xshift        0.0
##  Yshift        0.0
##  Parameters    
##  49  0  0.0 /* 1st standard parallel
##  77  0  0.0 /* 2nd standard parallel
##  -95  0  0.0 /* central meridian
##  0  0  0.0 /* latitude of projection's origin
##  0.0 /* false easting (meters)
##  0.0 /* false northing (meters)


# Define the model variables (PC1 - PC6)
modvars <- paste0("PC",1:6)

# Set up the linear (LIN) model
lin.modform  <- as.formula(paste("PRES~",paste(modvars,collapse="+")))

# Set up the quadratic (QUAD) model
quad.modform <- as.formula(paste("PRES~",paste(modvars,collapse="+"),"+",paste0("I(",modvars,"^2)",collapse=" + ")))

# Set up the Random Forest (RF) model (linear variables only)
rf.modform   <- as.formula(paste("as.factor(PRES)~",paste(modvars[1:6],collapse="+")))
# Set the number of RF trees
n.trees <- 1000

# Make output directories
dir.create("Models")   # For models and model summaries
dir.create("CV")       # For cross-validation results
dir.create("Plots")    # Save plots here


## PART 2 - Create the full data models
## ------------------------------------

# Create the full data LIN model using stepwise selection
lin.glm.full <- stepAIC(glm(lin.modform, data=tr, family="binomial", na.action="na.fail"), direction="both", trace=F)
# Correlation between predictors should ideally be < 0.7
summary(lin.glm.full, cor=T)
# Save the full data LIN model & model summary
save(lin.glm.full, file="Models/Linear GLM (Full Data) Model.RData")
capture.output(summary(lin.glm.full, cor=T), file="Models/Linear GLM (Full Data) Model Summary.txt")

# Create the full data QUAD model using stepwise selection
quad.glm.full <- stepAIC(glm(quad.modform, data=tr, family="binomial", na.action="na.fail"), direction="both", trace=F)
# Correlation between predictors should ideally be < 0.7
summary(quad.glm.full, cor=T)
# Save the full data QUAD model & model summary
save(quad.glm.full, file="Models/Quadratic GLM (Full Data) Model.RData")
capture.output(summary(quad.glm.full, cor=T), file="Models/Quadratic GLM (Full Data) Model Summary.txt")

# Create the full data RF model
rf.full <- randomForest(rf.modform, data=tr, ntree=n.trees, do.trace=T, importance=T)
rf.full
# Save the full data RF model & model summary
save(rf.full, file="Models/Random Forest (Full Data) Model.RData")
capture.output(rf.full, file="Models/Random Forest (Full Data) Model Summary.txt")

# Load the model back in later
# load("Models/Linear GLM (Full Data) Model.RData")
# load("Models/Quadratic GLM (Full Data) Model.RData")
# load("Models/Random Forest (Full Data) Model.RData")

# Reduced the model formulae to those variables retained in the stepwise selection
lin.modform.red  <- as.formula(paste("PRES~",paste(modvars[1:6],collapse="+")))
quad.modform.red <- as.formula(paste("PRES~",paste(modvars[1:6],collapse="+"),"+",paste0("I(",modvars[1:5],"^2)",collapse=" + "))) 
rf.modform.red  <- as.formula(paste("as.factor(PRES)~",paste(modvars[1:6],collapse="+")))



## PART 3 - Measure autocorrelation
## --------------------------------

# Create output directory
dir.create("Autocorrelation")

# Make model residuals
tr.res <- data.frame(tr[c("x","y","PRES")],
                     LinRes=round(tr$PRES-predict(lin.glm.full, tr, type="response"),5),
                     QuadRes=round(tr$PRES-predict(quad.glm.full, tr, type="response"),5),
                     RFRes=round(tr$PRES-predict(rf.full, tr, type="prob")[,2],5))

# Subsample the training data for memory limitations
t.sam <- sample(1:nrow(tr.res),5000)
dat.t <- tr.res[t.sam,]

# Convert to km
dat.t[c("x","y")] <- round(dat.t[c("x","y")]/1000)

# List of the models to process
model.list <- c("Linear GLM","Quadratic GLM","Random Forest")

## CORRELOGRAMS

# LOOP for each model
for(i in 1:length(model.list)){
  
  # Set the model name
  model <- model.list[i]
  
  # Correlogram (up to 4000km)
  ncf.t <- ncf::correlog(dat.t$x, dat.t$y, dat.t[,3+i], increment=50, resamp=1)
  ncf.t <- data.frame(Model=model, Dist=ncf.t$mean.of.class, Moran=ncf.t$correlation, p=ncf.t$p, n=ncf.t$n)
  
  # Output table
  if(i==1){cor.out <- ncf.t} else {cor.out <- rbind(cor.out,ncf.t)} 
  
}

# Save the correlogram data
write.csv(cor.out, "Autocorrelation/Correlogram data for model residuals.csv", row.names=F)

# Plot the correlograms
png(file="Plots/Autocorrelation correlograms.png")
par(mfrow=c(length(model.list),1)); layout.show(length(model.list))
for(i in model.list){
  plot(Moran~Dist, data=cor.out[cor.out$Model==i,], type="l", xlim=c(0,4000), main=i, ylab="Moran's I", xlab="Distance")
}
dev.off()

## SEMIVARIOGRAMS

# LOOP for each model
for(i in 1:length(model.list)){
  
  # Set the model name
  model <- model.list[i]
  
  # Semivariogram (up to 4000km)
  v.xy <- variog(coords=dat.t[c("x","y")], data=dat.t[,3+i], breaks=seq(0,4000,25), messages=F)
  v.xy.sum <- data.frame(Model=model, Dist=v.xy$u, SV=v.xy$v, N=v.xy$n)
  
  # Output table
  if(i==1){var.out <- v.xy.sum} else {var.out <- rbind(var.out,v.xy.sum)} 
  
}

# Save the semivariogram data
write.csv(var.out, "Autocorrelation/Semivariogram data for model residuals.csv", row.names=F)

# Plot the semivariograms
png(file="Plots/Autocorrelation semivariograms.png")
par(mfrow=c(length(model.list),1)); layout.show(length(model.list))
for(i in model.list){
  plot(SV~Dist, data=var.out[var.out$Model==i,], type="l", xlim=c(0,4000), main=i, ylab="Semivariance", xlab="Distance")
}
dev.off()



## PART 4a - Resubstitution
## ------------------------

# Set up the folds for the resubstitution (resub)
# NOTE: for the resub, this is all one single fold
f.resub <- rep(1,nrow(tr))
write.csv(data.frame(ID=tr$ID, f.resub=f.resub), "CV/Fold List - Resubstitution.csv", row.names=F)
# f.resub <- read.csv("CV/Fold List - Resubstitution.csv")$f.resub

# Quick plot to visualise the folds
plot(y~x, data=tr, asp=1, las=1, axes=F, ann=F, pch=20, cex=0.5, col=rainbow(length(unique(f.resub)))[f.resub])

# Table of Coefficients/Importance
c.resub <- data.frame(Fold=1, Variable=c("(Intercept)",modvars,paste0("I(",modvars,"^2)")))#, Linear=NA, Quadratic=NA, RF_IncMSE=NA)
c.resub[c.resub$Variable %in% names(lin.glm.full$coefficients),"Linear"]     <- as.numeric(lin.glm.full$coefficients)
c.resub[c.resub$Variable %in% names(quad.glm.full$coefficients),"Quadratic"] <- as.numeric(quad.glm.full$coefficients)
c.resub[c.resub$Variable %in% modvars,"MeanDecreaseAccuracy"]                <- importance(rf.full)[,"MeanDecreaseAccuracy"]
c.resub[c.resub$Variable %in% modvars,"MeanDecreaseGini"]                    <- importance(rf.full)[,"MeanDecreaseGini"]
write.csv(c.resub, "CV/Coefficient Table - Resubstitution.csv", row.names=F)
# c.resub <- read.csv("CV/Coefficient Table - Resubstitution.csv")

# Table of model predictions
p.resub <- data.frame(Fold=1, ID=tr$ID, PRES=tr$PRES)
p.resub["PRED_lin"]  <- round(predict(lin.glm.full, tr[modvars], type="response"),5)
p.resub["PRED_quad"] <- round(predict(quad.glm.full, tr[modvars], type="response"),5)
p.resub["PRED_rf"]   <- round(predict(rf.full, tr[modvars], type="prob")[,2],5)
write.csv(p.resub, "CV/Prediction Table - Resubstitution.csv", row.names=F)
# p.resub <- read.csv("CV/Prediction Table - Resubstitution.csv")

# Calculate the model residuals
p.resub$RES_lin  <- p.resub$PRES-p.resub$PRED_lin
p.resub$RES_quad <- p.resub$PRES-p.resub$PRED_quad
p.resub$RES_rf   <- p.resub$PRES-p.resub$PRED_rf

# Calculate a table of error estimates for resubstitution
# LOOP for each model
for(j in 1:3){
  
  # Create table (first loop only)
  if(j==1){e.resub <- data.frame(Model=c("Linear GLM","Quadratic GLM","Random Forest"), CV="resub", RMSE=NA, AUC=NA, Sens=NA, Spec=NA, LogLik=NA, LogLik_n=NA)}
  
  # Name models
  longmod <- c("Linear GLM","Quadratic GLM","Random Forest")[j]
  shortmod <- c("lin","quad","rf")[j]
  
  # Predictions and residuals
  e.resub.pred  <- p.resub[,paste0("PRED_",shortmod)]
  e.resub.res   <- p.resub[,paste0("RES_",shortmod)]
  e.resub.truth <- p.resub[,"PRES"]
  # RMSE
  e.resub[e.resub$Model==longmod,"RMSE"] <- sqrt(mean(e.resub.res^2))
  # AUC
  e.resub[e.resub$Model==longmod,c("AUC","Sens","Spec")] <- auc.f(e.resub.pred,e.resub.truth)
  
  # Adjustment for RF to calculate likelihood (removes perfect 1s and 0s)
  if(shortmod=="rf"){
    e.resub.pred[e.resub.pred==0] <- 1/n.trees
    e.resub.pred[e.resub.pred==1] <- 1-1/n.trees
  } else {
    e.resub.pred[e.resub.pred==0] <- min(e.resub.pred[e.resub.pred!=0])
    e.resub.pred[e.resub.pred==1] <- min(e.resub.pred[e.resub.pred!=1])
  }
  
  # Log likelihood (overall)
  e.resub[e.resub$Model==longmod,c("LogLik")] <- sum(dbinom(e.resub.truth, size=1, prob=e.resub.pred, log=T))
  # Log likelihood (per data point, to compare uneven folds)
  e.resub[e.resub$Model==longmod,c("LogLik_n")] <- e.resub[e.resub$Model==longmod,"LogLik"]/nrow(tr)
  # e.resub <- read.csv("CV/Validation Table - Resubstitution.csv")
  
} # LOOP for each model

# Save the error table
write.csv(e.resub, "CV/Validation Table - Resubstitution.csv", row.names=F)
# e.resub <- read.csv("CV/Validation Table - Resubstitution.csv")



## PART 4b - Random Folds
## ----------------------

# List of number of folds to run (for random CV)
ran.fold.list <- c(2,4,6,8)

# LOOP for all numbers of folds
for(n.ran in ran.fold.list){

  # Create the random folds
  f.ran <- sample(1:n.ran,nrow(tr), replace=T)
  table(f.ran)
  # f.ran <- read.csv(paste0("CV/Fold List - Random-",n.ran,".csv"))$f.ran
  
  # Quick plot to visualise the folds
  plot(y~x, data=tr, asp=1, las=1, axes=F, ann=F, pch=20, cex=0.5, col=rainbow(length(unique(f.ran)))[f.ran])
  
  # LOOP for all folds in the k-fold random CV
  for(i in 1:n.ran){
    
    # Model Building
    li <- glm(lin.modform.red, data=tr[f.ran!=i,], family="binomial")
    qu <- glm(quad.modform.red, data=tr[f.ran!=i,], family="binomial")
    rf <- randomForest(rf.modform.red, data=tr[f.ran!=i,], ntree=n.trees, do.trace=F, importance=T)
    
    # Table of Coefficients/Importance
    if(i==1){c.ran <- data.frame(Fold=rep(1:n.ran,each=length(modvars)*2+1), Variable=c("(Intercept)",modvars,paste0("I(",modvars,"^2)")), Linear=NA, Quadratic=NA, RF_IncMSE=NA)}
    c.ran[c.ran$Variable %in% names(li$coefficients) & c.ran$Fold==i,"Linear"]    <- as.numeric(li$coefficients)
    c.ran[c.ran$Variable %in% names(qu$coefficients) & c.ran$Fold==i,"Quadratic"] <- as.numeric(qu$coefficients)
    c.ran[c.ran$Variable %in% modvars & c.ran$Fold==i,"RF_IncMSE"]                <- importance(rf)[,1]
    #c.ran <- read.csv(paste0("CV/Coefficient Table - Random-",n.ran,".csv",))
    
    # Table of model predictions
    if(i==1){p.ran <- data.frame(Fold=f.ran, ID=tr$ID, PRES=tr$PRES)}
    p.ran[p.ran$Fold==i, "PRED_lin"]  <- round(predict(li, tr[p.ran$Fold==i,modvars], type="response"),5)
    p.ran[p.ran$Fold==i, "PRED_quad"] <- round(predict(qu, tr[p.ran$Fold==i,modvars], type="response"),5)
    p.ran[p.ran$Fold==i, "PRED_rf"]   <- round(predict(rf, tr[p.ran$Fold==i,modvars], type="prob")[,2],5)
    #p.ran <- read.csv(paste0("CV/Prediction Table - Random-",n.ran,".csv"))
    
  } # LOOP for each fold
  
  # Save the tables
  write.csv(data.frame(ID=tr$ID, f.ran=f.ran), paste0("CV/Fold List - Random-",n.ran,".csv"), row.names=F)
  write.csv(c.ran, paste0("CV/Coefficient Table - Random-",n.ran,".csv"), row.names=F)
  write.csv(p.ran, paste0("CV/Prediction Table - Random-",n.ran,".csv"), row.names=F)
  
  # Calculate the model residuals
  p.ran$RES_lin  <- p.ran$PRES-p.ran$PRED_lin
  p.ran$RES_quad <- p.ran$PRES-p.ran$PRED_quad
  p.ran$RES_rf   <- p.ran$PRES-p.ran$PRED_rf
  
  # Calculate a table of error estimates for random k-fold CV
  # LOOP for each model
  for(j in 1:3){
    
    # Create table (first loop only)
    if(j==1){e.ran <- data.frame(Model=c("Linear GLM","Quadratic GLM","Random Forest"), CV=paste0("ran",n.ran), RMSE=NA, AUC=NA, Sens=NA, Spec=NA, LogLik=NA, LogLik_n=NA)}
    
    # Name models
    longmod <- c("Linear GLM","Quadratic GLM","Random Forest")[j]
    shortmod <- c("lin","quad","rf")[j]
    
    # Predictions and residuals
    e.ran.pred  <- p.ran[,paste0("PRED_",shortmod)]
    e.ran.res   <- p.ran[,paste0("RES_",shortmod)]
    e.ran.truth <- p.ran[,"PRES"]
    # RMSE
    e.ran[e.ran$Model==longmod,"RMSE"] <- sqrt(mean(e.ran.res^2))
    # AUC
    e.ran[e.ran$Model==longmod,c("AUC","Sens","Spec")] <- auc.f(e.ran.pred,e.ran.truth)
    
    # Adjustment for RF to calculate likelihood (removes perfect 1s and 0s)
    if(shortmod=="rf"){
      e.ran.pred[e.ran.pred==0] <- 1/n.trees
      e.ran.pred[e.ran.pred==1] <- 1-1/n.trees
    } else {
      e.ran.pred[e.ran.pred==0] <- min(e.ran.pred[e.ran.pred!=0])
      e.ran.pred[e.ran.pred==1] <- min(e.ran.pred[e.ran.pred!=1])
    }
  
    # Log likelihood (overall)
    e.ran[e.ran$Model==longmod,c("LogLik")] <- sum(dbinom(e.ran.truth, size=1, prob=e.ran.pred, log=T))
    # Log likelihood (per data point, to compare uneven folds)
    e.ran[e.ran$Model==longmod,c("LogLik_n")] <- e.ran[e.ran$Model==longmod,"LogLik"]/nrow(tr)
  
  } # LOOP for each model
  
  # Save the error table
  write.csv(e.ran, paste0("CV/Validation Table - Random-",n.ran,".csv"), row.names=F)
  #e.ran <- read.csv(paste0("CV/Validation Table - Random-",n.ran,".csv"))
  
} # LOOP for all number of folds



## PART 4c - Spatially blocked folds
## ---------------------------------

# Function to create spatially blocked folds
grid.fun <- function(dat, dat.x, dat.y, nx, ny){
  
  # dat   = the data table to block
  # dat.x = variable name of the x coordinate
  # dat.y = variable name of the y coordinate
  # nx    = the number of spatial divisions in x
  # ny    = the number of spatial divisions in y
  
  library('yaImpute')
  
  # Calculat the block sizes based on the data extents
  # NOTE: this calculation will be sensitive to spatial outliers
  x.size <- diff(range(dat[,dat.x]))/(nx)
  y.size <- diff(range(dat[,dat.y]))/(ny)
  
  # Create centroid coordinates for the spatial blocks
  cent.x <- (min(dat[,dat.x]) + (x.size*(1:nx))) - x.size/2
  cent.y <- (min(dat[,dat.y]) + (y.size*(1:ny))) - y.size/2
  cent.grid <- expand.grid(cent.x,cent.y); names(cent.grid) <- c("x","y")
  
  # Randomly assign fold numbers to the blocks
  foldnames <- sample(1:nrow(cent.grid),replace=F)
  cent.grid <- cent.grid[order(rev(cent.grid$y),cent.grid$x),]
  cent.grid$FOLD=foldnames
  
  # Use a nearest neighbour approach to assign all points to the closest centroid
  f <- cent.grid[ann(as.matrix(cent.grid[c("x","y")]), as.matrix(data.frame(x=dat[,dat.x],y=dat[,dat.y])), k=1, verbose=F)$knnIndexDist[,1],"FOLD"]
  
  # Output the spatially blocked folds
  ft <- table(f)
  list(N_Folds=nrow(cent.grid), FoldTable=ft, GridCoords=cent.grid, Folds=f)
}

# Table of block arrangements to run
block.arrange <- data.frame(x=c(1,2,2,2,3,4,5,6,7,8,9,10,15,20),
                            y=c(2,2,3,4,6,8,10,12,14,16,18,20,30,40))

# LOOP for each spatial block arrangement
for(b in 1:nrow(block.arrange)){
  
  # Set up the blocking arrangement
  in.x <- block.arrange[b,"x"]
  in.y <- block.arrange[b,"y"]
  
  # Calculate and assign block with blocking function
  grid.b <- grid.fun(dat=tr, dat.x="x", dat.y="y", nx=in.x, ny=in.y)
  
  # Quick plot to visualise the folds
  plot(y~x, data=tr, asp=1, las=1, axes=F, ann=F, pch=20, cex=0.5, col=rainbow(grid.b$N_Folds)[grid.b$Folds])
  
  # Create output table of fold assignments
  f.b.out <- data.frame(ID=tr$ID, f.b=grid.b$Folds)
  names(f.b.out) <- c("ID",paste0("F.",in.x,"x",in.y))

  # LOOP for all folds in the spatial k-fold CV
  fold.names <- unique(grid.b$Folds)[order(unique(grid.b$Folds))]
  for(i in fold.names){
    
    # Model Building
    li <- glm(lin.modform.red, data=tr[grid.b$Folds!=i,], family="binomial")
    qu <- glm(quad.modform.red, data=tr[grid.b$Folds!=i,], family="binomial")
    rf <- randomForest(rf.modform.red, data=tr[grid.b$Folds!=i,], ntree=n.trees, do.trace=F, importance=T)
    
    # Table of Coefficients/Importance
    if(i==min(fold.names)){c.b <- data.frame(Fold=rep(1:grid.b$N_Folds,each=length(modvars)*2+1), Variable=c("(Intercept)",modvars,paste0("I(",modvars,"^2)")), Linear=NA, Quadratic=NA, RF_IncMSE=NA)}
    c.b[c.b$Variable %in% names(li$coefficients) & c.b$Fold==i,"Linear"]    <- as.numeric(li$coefficients)
    c.b[c.b$Variable %in% names(qu$coefficients) & c.b$Fold==i,"Quadratic"] <- as.numeric(qu$coefficients)
    c.b[c.b$Variable %in% modvars & c.b$Fold==i,"RF_IncMSE"]                <- importance(rf)[,1]
    
    # Table of model predictions
    if(i==1){p.b <- data.frame(Fold=grid.b$Folds, ID=tr$ID, PRES=tr$PRES)}
    p.b[p.b$Fold==i, "PRED_lin"]  <- round(predict(li, tr[p.b$Fold==i,modvars], type="response"),5)
    p.b[p.b$Fold==i, "PRED_quad"] <- round(predict(qu, tr[p.b$Fold==i,modvars], type="response"),5)
    p.b[p.b$Fold==i, "PRED_rf"]   <- round(predict(rf, tr[p.b$Fold==i,modvars], type="prob")[,2],5)
    
  } # LOOP for each fold
   
  # Save tables
  write.csv(f.b.out, paste0("CV/Fold List - Spatial Grid ",in.x,"x",in.y,".csv"), row.names=F)
  write.csv(c.b, paste0("CV/Coefficient Table - Spatial Grid ",in.x,"x",in.y,".csv"), row.names=F)
  write.csv(p.b, paste0("CV/Prediction Table - Spatial Grid ",in.x,"x",in.y,".csv"), row.names=F)
  #f.b <- read.csv(paste0("CV/Fold List - Spatial Grid ",in.x,"x",in.y,".csv"))$f.b
  #c.b <- read.csv(paste0("CV/Coefficient Table - Spatial Grid ",in.x,"x",in.y,".csv"))
  #p.b <- read.csv(paste0("CV/Prediction Table - Spatial Grid ",in.x,"x",in.y,".csv"))
  
  # Calculate the model residuals
  p.b$RES_lin  <- p.b$PRES-p.b$PRED_lin
  p.b$RES_quad <- p.b$PRES-p.b$PRED_quad
  p.b$RES_rf   <- p.b$PRES-p.b$PRED_rf
  
  # Calculate a table of error estimates for blocked k-fold CV
  # LOOP for each model
  for(j in 1:3){
    
    # Create table (first loop only)
    if(j==1){e.b <- data.frame(Model=c("Linear GLM","Quadratic GLM","Random Forest"), CV=paste0(in.x,"x",in.y), RMSE=NA, AUC=NA, Sens=NA, Spec=NA, LogLik=NA, LogLik_n=NA)}
    
    # Name models
    longmod <- c("Linear GLM","Quadratic GLM","Random Forest")[j]
    shortmod <- c("lin","quad","rf")[j]
    
    # Predictions and residuals
    e.b.pred  <- p.b[,paste0("PRED_",shortmod)]
    e.b.res   <- p.b[,paste0("RES_",shortmod)]
    e.b.truth <- p.b[,"PRES"]
    # RMSE
    e.b[e.b$Model==longmod,"RMSE"] <- sqrt(mean(e.b.res^2))
    # AUC
    e.b[e.b$Model==longmod,c("AUC","Sens","Spec")] <- auc.f(e.b.pred,e.b.truth)
    
    # Adjustment for RF to calculate likelihood (removes perfect 1s and 0s)
    if(shortmod=="rf"){
      e.b.pred[e.b.pred==0] <- 1/n.trees
      e.b.pred[e.b.pred==1] <- 1-1/n.trees
    } else {
      e.b.pred[e.b.pred==0] <- min(e.b.pred[e.b.pred!=0])
      e.b.pred[e.b.pred==1] <- min(e.b.pred[e.b.pred!=1])
    }
    
    # Log likelihood (overall)
    e.b[e.b$Model==longmod,c("LogLik")] <- sum(dbinom(e.b.truth, size=1, prob=e.b.pred, log=T), na.rm=T)
    # Log likelihood (per data point, to compare uneven folds)
    e.b[e.b$Model==longmod,c("LogLik_n")] <- e.b[e.b$Model==longmod,"LogLik"]/nrow(tr)
    
  } # LOOP for each model
  
  # Save the error table
  write.csv(e.b, paste0("CV/Validation Table - Spatial Grid ",in.x,"x",in.y,".csv"), row.names=F)
  #e.b <- read.csv(paste0("CV/Validation Table - Spatial Grid ",in.x,"x",in.y,".csv"))

} # LOOP for all spatial fold arrangements



## PART 4d - Environmentally blocked folds using PCA cuts
## ------------------------------------------------------

# PCA component variance explained for adjustment of the weightings of the variables in the splitting
# Ensures that higher components explaining less data variance do not drive the environmental splits
pca.varex <- c(0.5484876,0.2757948,0.07201832,0.05321062,0.04440524,0.00608348)
# Make adjusted PCA variables for environmental splitting 
pc.adj <- tr[modvars]*rep(pca.varex,nrow(tr))

# Set the number of PCA components to cut (n folds in CV = 2 * n cuts)
pca.cut.list <- c(1,2,3,4)

# LOOP for each level of PCA cutting
for(pc in pca.cut.list){
  
  # Cut the PCA components
  for(j in 1:pc){
    if(j==1){q <- list()}
    q[[j]] <- quantile(pc.adj[,j],c(0.25,0.75))
  }
  names(q) <- paste0("PC",1:pc)
  
  # Make a table of PCA centroids and assign fold numbers
  pc.cent <- expand.grid(q)
  pc.cent$Fold <- 1:nrow(pc.cent)
  n.cut <- nrow(pc.cent)
  
  # Use nearest neighbour calculation to assign data to folds based on centroids
  f.cut <- pc.cent$Fold[ann(as.matrix(pc.cent[modvars[1:pc]]), 
                            as.matrix(tr[modvars][1:pc]*rep(pca.varex[1:pc],nrow(tr))), 
                            k=1, verbose=F)$knnIndexDist[,1]]
  
  # Quick plot to visualise the folds
  plot(y~x, data=tr, asp=1, las=1, axes=F, ann=F, pch=20, cex=0.5, col=rainbow(n.cut)[f.cut])
  
  # LOOP for all folds in the PCA blocked CV
  for(i in 1:n.cut){
    
    # Model Building
    li <- glm(lin.modform.red, data=tr[f.cut!=i,], family="binomial")
    qu <- glm(quad.modform.red, data=tr[f.cut!=i,], family="binomial")
    rf <- randomForest(rf.modform.red, data=tr[f.cut!=i,], ntree=n.trees, do.trace=F, importance=T)
    
    # Table of Coefficients/Importance
    if(i==1){c.cut <- data.frame(Fold=rep(1:n.cut,each=length(modvars)*2+1), Variable=c("(Intercept)",modvars,paste0("I(",modvars,"^2)")), Linear=NA, Quadratic=NA, RF_IncMSE=NA)}
    c.cut[c.cut$Variable %in% names(li$coefficients) & c.cut$Fold==i,"Linear"]    <- as.numeric(li$coefficients)
    c.cut[c.cut$Variable %in% names(qu$coefficients) & c.cut$Fold==i,"Quadratic"] <- as.numeric(qu$coefficients)
    c.cut[c.cut$Variable %in% modvars & c.cut$Fold==i,"RF_IncMSE"]                <- importance(rf)[,1]
    
    # Table of model predictions
    if(i==1){p.cut <- data.frame(Fold=f.cut, ID=tr$ID, PRES=tr$PRES)}
    p.cut[p.cut$Fold==i, "PRED_lin"]  <- round(predict(li, tr[p.cut$Fold==i,modvars], type="response"),5)
    p.cut[p.cut$Fold==i, "PRED_quad"] <- round(predict(qu, tr[p.cut$Fold==i,modvars], type="response"),5)
    p.cut[p.cut$Fold==i, "PRED_rf"]   <- round(predict(rf, tr[p.cut$Fold==i,modvars], type="prob")[,2],5)
    
  } # LOOP for each fold
  
  # Save the tables
  write.csv(data.frame(ID=tr$ID, f.cut=f.cut), paste0("CV/Fold List - PCAcut-",n.cut,".csv"), row.names=F)
  write.csv(c.cut, paste0("CV/Coefficient Table - PCAcut-",n.cut,".csv"), row.names=F)
  write.csv(p.cut, paste0("CV/Prediction Table - PCAcut-",n.cut,".csv"), row.names=F)
  #f.cut <- read.csv(paste0("CV/Fold List - PCAcut-",n.cut,".csv"))$f.cut
  #c.cut <- read.csv(paste0("CV/Coefficient Table - PCAcut-",n.cut,".csv"))
  #p.cut <- read.csv(paste0("CV/Prediction Table - PCAcut-",n.cut,".csv"))
  
  # Calculate the model residuals
  p.cut$RES_lin  <- p.cut$PRES-p.cut$PRED_lin
  p.cut$RES_quad <- p.cut$PRES-p.cut$PRED_quad
  p.cut$RES_rf   <- p.cut$PRES-p.cut$PRED_rf
  
  # Calculate a table of error estimates by PCA blocked CV
  # LOOP for each model
  for(j in 1:3){
    
    # Create table (first loop only)
    if(j==1){e.cut <- data.frame(Model=c("Linear GLM","Quadratic GLM","Random Forest"), CV=paste0("cut",n.cut), RMSE=NA, AUC=NA, Sens=NA, Spec=NA, LogLik=NA, LogLik_n=NA)}
    
    # Name models
    longmod <- c("Linear GLM","Quadratic GLM","Random Forest")[j]
    shortmod <- c("lin","quad","rf")[j]
    
    # Predictions and residuals
    e.cut.pred  <- p.cut[,paste0("PRED_",shortmod)]
    e.cut.res   <- p.cut[,paste0("RES_",shortmod)]
    e.cut.truth <- p.cut[,"PRES"]
    
    # RMSE
    e.cut[e.cut$Model==longmod,"RMSE"] <- sqrt(mean(e.cut.res^2))
    # AUC
    e.cut[e.cut$Model==longmod,c("AUC","Sens","Spec")] <- auc.f(e.cut.pred,e.cut.truth)
    
    # Adjustment for RF to calculate likelihood (removes perfect 1s and 0s)
    if(shortmod=="rf"){
      e.cut.pred[e.cut.pred==0] <- 1/n.trees
      e.cut.pred[e.cut.pred==1] <- 1-1/n.trees
    } else {
      e.cut.pred[e.cut.pred==0] <- min(e.cut.pred[e.cut.pred!=0])
      e.cut.pred[e.cut.pred==1] <- min(e.cut.pred[e.cut.pred!=1])
    }
    # Log likelihood (overall)
    e.cut[e.cut$Model==longmod,c("LogLik")] <- sum(dbinom(e.cut.truth, size=1, prob=e.cut.pred, log=T))
    # Log likelihood (per data point, to compare uneven folds)
    e.cut[e.cut$Model==longmod,c("LogLik_n")] <- e.cut[e.cut$Model==longmod,"LogLik"]/nrow(tr)
    
  } # LOOP for each model
  
  # Save the error table
  write.csv(e.cut, paste0("CV/Validation Table - PCAcut-",n.cut,".csv"), row.names=F)
  #e.cut <- read.csv(paste0("CV/Validation Table - PCAcut-",n.cut,".csv"))
  
} # LOOP for each level of PCA cutting



## PART 4e - Environmentally blocked folds using k-means clustering
## ----------------------------------------------------------------

# PCA component variance explained for adjustment of the weightings of the variables in the splitting
# Ensures that higher components explaining less data variance do not drive the environmental splits
pca.varex <- c(0.5484876,0.2757948,0.07201832,0.05321062,0.04440524,0.00608348)
# Make adjusted PCA variables for environmental splitting 
pc.adj <- tr[modvars]*rep(pca.varex,nrow(tr))

# Set the number of clusters to run
cluster.list <- c(2,4,6,8)

# LOOP for all numbers of clusters (folds)
for(n.clus in cluster.list){
  
  # Create folds by environment using k-means cluster analysis
  # Start by clustering presence data (ensures each fold contains some presence data)
  clus <- kmeans(pc.adj[tr$PRES==1,], centers=n.clus, iter.max=100, nstart=n.clus*10)
  f.clus.pres <- as.numeric(clus$cluster)
  # Match absences to cluster groups using nearest neighbour
  f.clus.abs <- f.clus.pres[ann(as.matrix(pc.adj[tr$PRES==1,]), as.matrix(pc.adj[tr$PRES==0,]), k=1, verbose=F)$knnIndexDist[,1]]
  
  # Make the vecor of block/fold assignments
  f.clus <- rep(NA,nrow(tr))
  f.clus[tr$PRES==1] <- f.clus.pres
  f.clus[tr$PRES==0] <- f.clus.abs
  
  # Quick plot to visualise the folds
  plot(y~x, data=tr, asp=1, las=1, axes=F, ann=F, pch=20, cex=0.5, col=rainbow(n.clus)[f.clus])
  
  # LOOP for all folds in the cluster PCA blocked CV
  for(i in 1:n.clus){
   
    # Model Building
    if(i==1){c.clus <- data.frame(Fold=rep(1:n.clus,each=length(modvars)*2+1), Variable=c("(Intercept)",modvars,paste("I(",modvars,"^2)",sep="")), Linear=NA, Quadratic=NA, RF_IncMSE=NA)}
    li <- glm(lin.modform.red, data=tr[f.clus!=i,], family="binomial")
    qu <- glm(quad.modform.red, data=tr[f.clus!=i,], family="binomial")
    rf <- randomForest(rf.modform.red, data=tr[f.clus!=i,], ntree=n.trees, do.trace=F, importance=T)
    
    # Table of Coefficients/Importance
    if(i==1){p.clus <- data.frame(Fold=f.clus, ID=tr$ID, PRES=tr$PRES)}
    c.clus[c.clus$Variable %in% names(li$coefficients) & c.clus$Fold==i,"Linear"]    <- as.numeric(li$coefficients)
    c.clus[c.clus$Variable %in% names(qu$coefficients) & c.clus$Fold==i,"Quadratic"] <- as.numeric(qu$coefficients)
    c.clus[c.clus$Variable %in% modvars & c.clus$Fold==i,"RF_IncMSE"]                <- importance(rf)[,1]
    
    # Table of model predictions
    p.clus[p.clus$Fold==i, "PRED_lin"]  <- round(predict(li, tr[p.clus$Fold==i,modvars], type="response"),5)
    p.clus[p.clus$Fold==i, "PRED_quad"] <- round(predict(qu, tr[p.clus$Fold==i,modvars], type="response"),5)
    p.clus[p.clus$Fold==i, "PRED_rf"]   <- round(predict(rf, tr[p.clus$Fold==i,modvars], type="prob")[,2],5)
     
  } # LOOP for all folds
  
  # Save the tables
  write.csv(data.frame(ID=tr$ID, f.clus=f.clus), paste0("CV/Fold List - Cluster-",n.clus,".csv"), row.names=F)
  write.csv(c.clus, paste0("CV/Coefficient Table - Cluster-",n.clus,".csv"), row.names=F)
  write.csv(p.clus, paste0("CV/Prediction Table - Cluster-",n.clus,".csv"), row.names=F)
  #f.clus <- read.csv(paste0("CV/Fold List - Cluster-",n.ran,".csv")$f.clus
  #c.clus <- read.csv(paste0("CV/Coefficient Table - Cluster-",n.clus,".csv"))
  #p.clus <- read.csv(paste0("CV/Prediction Table - Cluster-",n.clus,".csv"))
  
  # Calculate the model residuals
  p.clus$RES_lin  <- p.clus$PRES-p.clus$PRED_lin
  p.clus$RES_quad <- p.clus$PRES-p.clus$PRED_quad
  p.clus$RES_rf   <- p.clus$PRES-p.clus$PRED_rf
  
  # Calculate a table of error estimates by PCA blocked CV
  # LOOP for each model
  for(j in 1:3){
    
    # Create table (first loop only)
    if(j==1){e.clus <- data.frame(Model=c("Linear GLM","Quadratic GLM","Random Forest"), CV=paste0("clus",n.clus), RMSE=NA, AUC=NA, Sens=NA, Spec=NA, LogLik=NA, LogLik_n=NA)}
    
    # Name models
    longmod <- c("Linear GLM","Quadratic GLM","Random Forest")[j]
    shortmod <- c("lin","quad","rf")[j]
    
    # Predictions and residuals
    e.clus.pred  <- p.clus[,paste0("PRED_",shortmod)]
    e.clus.res   <- p.clus[,paste0("RES_",shortmod)]
    e.clus.truth <- p.clus[,"PRES"]
    
    # RMSE
    e.clus[e.clus$Model==longmod,"RMSE"] <- sqrt(mean(e.clus.res^2))
    # AUC
    e.clus[e.clus$Model==longmod,c("AUC","Sens","Spec")] <- auc.f(e.clus.pred,e.clus.truth)
    
    # Adjustment for RF to calculate likelihood (removes perfect 1s and 0s)
    if(shortmod=="rf"){
      e.clus.pred[e.clus.pred==0] <- 1/n.trees
      e.clus.pred[e.clus.pred==1] <- 1-1/n.trees
    } else {
      e.clus.pred[e.clus.pred==0] <- min(e.clus.pred[e.clus.pred!=0])
      e.clus.pred[e.clus.pred==1] <- min(e.clus.pred[e.clus.pred!=1])
    }
    # Log likelihood (overall)
    e.clus[e.clus$Model==longmod,c("LogLik")] <- sum(dbinom(e.clus.truth, size=1, prob=e.clus.pred, log=T))
    # Log likelihood (per data point, to compare uneven folds)
    e.clus[e.clus$Model==longmod,c("LogLik_n")] <- e.clus[e.clus$Model==longmod,"LogLik"]/nrow(tr)
    
  } # LOOP for each model
  
  # Save the error table
  write.csv(e.clus, paste0("CV/Validation Table - Cluster-",n.clus,".csv"), row.names=F)
  #e.clus <- read.csv(paste0("CV/Validation Table - Cluster-",n.clus,".csv"))

} # LOOP for all numbers of clusters (folds)



## PART 4e - Spatially buffered leave-one-out (LOO)
## ------------------------------------------------

# Spatially buffered leave-one-out function with GLM or RF
# NOTE: LOO analysis can be comuptationally intensive for large data sets & complex models
# See below for data subsampling approaches
loo.fun <- function(dat, x_coord, y_coord, resp, rad, modform, modmethod, AUC){
  
  # dat              = complete data file
  # x_coord(y_coord) = names of x(y) coordinates in data
  # truth            = true value of response
  # ss               = sample size to process (number of LOO runs)
  # rad              = radius for the spatial buffer around test point
  # modform          = formula for the GLM
  # modmethod        = one of "GLM" or "RF"
  # AUC              = logical T/F: should AUC be calculated?
  
  library()
  
  # LOOP for all points in data (exhaustive LOO)
  for(i in 1:nrow(dat)){
    
    # Make empty vector for predictions
    if(i==1){p <- c()}
    
    # Make training data (test point & buffer removed)
    train <- dat[sqrt((dat[,x_coord]-dat[i,x_coord])^2 + (dat[,y_coord]-dat[i,y_coord])^2)>rad,]
    
    # Build the model
    # GLM
    if (modmethod=="GLM"){
      m <- glm(modform, data=train, family="binomial")
      # Predict on test point
      p <- c(p, predict(m, dat[i,], type="response"))
    }
    # Random Forest
    if (modmethod=="RF"){
      m <- randomForest(modform, data=train, ntree=n.trees, do.trace=F, importance=F)
      # Predict on test point
      p <- c(p, predict(m, dat[i,], type="prob")[,2])
    }
    
  } # for all points
  
  # Calculate complete residuals
  p.res <- dat[,resp]-p
  
  # RMSE
  p.rmse <- sqrt(mean(p.res^2))
  
  # Output list
  p.out <- list(SampleRows=as.numeric(rownames(dat)),
                Truth=dat[,resp],  # List entry for true values
                Predictions=p,     # List entry for model predictions
                Residuals=p.res,   # List entry for model residuals
                RMSE=p.rmse)       # List entry for RMSE

  # AUC (if set to TRUE in function call)
  if(AUC==T){library(ROCR); p.out$AUC <- performance(prediction(p,dat[,resp]),"auc")@y.values[[1]]}
  
  # Report table
  p.out
  
}

# List of buffer sizes to run (in km)
buffer.list <- c(100,500,1000,1500)

# Beacuse exhaustive LOO is computationally intensive, you can subsample the dataset.
# Perform a sensitivity analysis to find a sample size providing stable results.
# Use the largest buffer you intend to run in the function call (removes most data).
# Can be tested with various model structures (substitute in the LOO function call)
test.samples <- c(50,100,250,500,1000,2500,5000)
test.reps <- 5
for(i in test.samples){
  if(i==test.samples[1]){test.tab <- data.frame(n=rep(test.samples, each=5), RMSE=NA, AUC=NA)}
  for(j in 1:test.reps){
    test.tab[test.tab$n==i,c("RMSE","AUC")][j,] <- loo.fun(tr[sample(1:nrow(tr),i),], "x", "y", "PRES", max(buffer.list), lin.modform.red, "GLM", AUC=T)[c("RMSE","AUC")]
  }
}
# Save the sensitivity analysis data
write.csv(test.tab, paste0("CV/Sensitivity Table - LOO (n=",n.loo,", ",max(buffer.list)/1000,"km buffer).csv"), row.names=F)
# Plot the sensitivity analysis
test.tab2 <- ddply(test.tab, .(n), summarise,
                   x.RMSE=mean(RMSE), sd.RMSE=sd(RMSE), se.RMSE=sd.RMSE/sqrt(length(RMSE)),
                   x.AUC=mean(AUC), sd.AUC=sd(AUC), se.AUC=sd.AUC/sqrt(length(AUC)))
png(file="Plots/LOO sample size sensitivity.png")
par(mfrow=c(2,2)); layout.show(4)
boxplot(RMSE~n, data=test.tab, las=1, ylab="RMSE", xlab="Sample size", main="RMSE")
boxplot(AUC~n, data=test.tab, las=1, ylab="AUC", xlab="Sample size", main="AUC")
plot(se.RMSE~n, data=test.tab2, las=1, type="b", ylab="RMSE Std Error", xlab="Sample size", main="RMSE")
plot(se.AUC~n, data=test.tab2, las=1, type="b", ylab="AUC Std Error", xlab="Sample size", main="AUC")
dev.off()

# LOO sample size (large enough to ensure accurate error estimation for the entire data set)
n.loo <- 500

# Sample from the data
f.loo <- tr[sample(nrow(tr),n.loo),]
f.loo <- f.loo[order(f.loo$ID),]
# Save the data subsample
write.csv(f.loo[,c("ID"),drop=F], paste0("CV/Data subsample - LOO (n=",n.loo,").csv"), row.names=F)

# LOOP for all buffer sizes in the LOO (x 1000 to convert km to m)
for(buffer in buffer.list*1000){
  
  # LOO function (exhaustive for the data subsample)
  p.loo.lin  <- loo.fun(f.loo, "x", "y", "PRES", buffer, lin.modform.red, "GLM", AUC=T)
  p.loo.quad <- loo.fun(f.loo, "x", "y", "PRES", buffer, quad.modform.red, "GLM", AUC=T)
  p.loo.rf   <- loo.fun(f.loo, "x", "y", "PRES", buffer, rf.modform.red, "RF", AUC=T)
  
  # Make the prediction and residual table
  p.loo <- data.frame(ID=p.loo.lin[[1]],
                       PRED_lin=p.loo.lin$Predictions,    # LIN predictions
                       PRED_quad=p.loo.quad$Predictions,  # QUAD Predictions
                       PRED_rf=p.loo.rf$Predictions,      # RF predictions
                       RES_lin=p.loo.lin$Residuals,       # LIN residuals
                       RES_quad=p.loo.quad$Residuals,     # QUAD residuals
                       RES_rf=p.loo.rf$Residuals,         # RF residuals
                       PRES=f.loo$PRES)                   # True data values (pres/abs)
  p.loo <- p.loo[order(p.loo$ID),]
  
  # Check to ensure that the RMSE stabilised in the subsample (i.e. if sample is large enough)
  #### I tried this, which I'm not sure is 100% right. Please check and let me know what you think.
  #### https://davidrroberts.wordpress.com/2016/03/11/spatial-leave-one-out-sloo-cross-validation/
  for(i in 1:n.loo){
    p.loo[i,"CumRMSE_lin"] <- sqrt(mean(p.loo[1:i,"RES_lin"]^2))
    p.loo[i,"CumRMSE_quad"] <- sqrt(mean(p.loo[1:i,"RES_quad"]^2))
    p.loo[i,"CumRMSE_rf"] <- sqrt(mean(p.loo[1:i,"RES_rf"]^2))
  }
  par(mfrow=c(1,1))
  plot(CumRMSE_lin~c(1:n.loo), data=p.loo, type="l", las=1, ylim=c(0,max(p.loo$CumRMSE_lin)), col="blue", xlab="Data points", ylab="RMSE")
  points(CumRMSE_quad~c(1:n.loo), data=p.loo, type="l", col="red")
  points(CumRMSE_rf~c(1:n.loo), data=p.loo, type="l", col="orange")
  legend("bottomright", legend=c("Linear GLM","Quadratic GLM","Random Forest"), col=c("blue","red","orange"), lty=1, cex=0.8)
  
  # Save the table of LOO predictions
  write.csv(p.loo, paste0("CV/Prediction Table - LOO (n=",n.loo,", b=",buffer/1000,"km).csv"), row.names=F)
  #p.loo <- read.csv(paste0("CV/Prediction Table - LOO (n=",n.loo,", b=",buffer/1000,"km).csv"))
  
  # Table of error estimates for LOO
  # LOOP for each model
  for(j in 1:3){
    
    # Create table (first loop only)
    if(j==1){e.loo <- data.frame(Model=c("Linear GLM","Quadratic GLM","Random Forest"), CV=paste0("loo",buffer/1000), RMSE=NA, AUC=NA, Sens=NA, Spec=NA, LogLik=NA, LogLik_n=NA)}
    
    # Name models
    longmod <- c("Linear GLM","Quadratic GLM","Random Forest")[j]
    shortmod <- c("lin","quad","rf")[j]
    
    # Predictions and residuals
    e.loo.pred  <- p.loo[,paste("PRED_",shortmod,sep="")]
    e.loo.res   <- p.loo[,paste("RES_",shortmod,sep="")]
    e.loo.truth <- p.loo[,"PRES"]
    
    # RMSE
    e.loo[e.loo$Model==longmod,"RMSE"] <- sqrt(mean(e.loo.res^2))
    # AUC
    e.loo[e.loo$Model==longmod,c("AUC","Sens","Spec")] <- auc.f(e.loo.pred,e.loo.truth)
    
    # Adjustment for RF to calculate likelihood (removes perfect 1s and 0s)
    if(shortmod=="rf"){
      e.loo.pred[e.loo.pred==0] <- 1/n.trees
      e.loo.pred[e.loo.pred==1] <- 1-1/n.trees
    } else {
      e.loo.pred[e.loo.pred==0] <- min(e.loo.pred[e.loo.pred!=0])
      e.loo.pred[e.loo.pred==1] <- min(e.loo.pred[e.loo.pred!=1])
    }
    # Log likelihood (overall)
    e.loo[e.loo$Model==longmod,c("LogLik")] <- sum(dbinom(e.loo.truth, size=1, prob=e.loo.pred, log=T))
    # Log likelihood (per data point, to compare uneven folds)
    e.loo[e.loo$Model==longmod,c("LogLik_n")] <- e.loo[e.loo$Model==longmod,"LogLik"]/nrow(tr)
  
  } # LOOP for each model
  
  # Save the error table
  write.csv(e.loo, paste0("CV/Validation Table - LOO (n=",n.loo,", b=",buffer/1000,"km).csv"), row.names=F)
  #e.loo <- read.csv(paste0("CV/Validation Table - LOO (n=",n.loo,", b=",buffer/1000,"km).csv"))

} # LOOP for each buffer size



## PART 5 - Compile CV results
## ---------------------------

# Key table of CV and CV types for summaries
CV.type.table <- data.frame(
  CVlong=c("Resubstitution",
           paste0("Random-",ran.fold.list),
           paste0("Spatial Grid ",apply(block.arrange,1,function(z){paste0(z[1],"x",z[2])})),
           paste0("Cluster-",cluster.list),
           paste0("PCAcut-",2^pca.cut.list),
           paste0("LOO (n=",n.loo,", b=",buffer.list,"km)")),
  CV=c("resub",
       paste0("ran",ran.fold.list),
       apply(block.arrange,1,function(z){paste0(z[1],"x",z[2])}),
       paste0("clus",cluster.list),
       paste0("cut",2^pca.cut.list),
       paste0("loo",buffer.list)),
  CVtype=c("Resub",
           rep("Random",length(ran.fold.list)),
           rep("Spatial",nrow(block.arrange)),
           rep("Enviro",length(cluster.list)),
           rep("Enviro",length(2^pca.cut.list)),
           rep("LOO",length(buffer.list)))
)

# List of files to process
file.list <- list.files("CV",pattern="Validation Table")
bcv.list <- substr(file.list,20,nchar(file.list)-4)

# LOOP for all cross-validations
for(b in bcv.list){

  # Compile error tables
  b.err <- read.csv(paste0("CV/Validation Table - ",b,".csv"))
  
  # Compile the results tables
  if(b==bcv.list[1]){
    bcv.err  <- data.frame(b.err["Model"], CVlong=b, b.err[-1])
  } else {
    bcv.err  <- rbind(bcv.err,  data.frame(b.err["Model"], CVlong=b, b.err[-1]))
  }

} # LOOP for all cross-validations

# Add CV types
bcv.err$CVtype <- CV.type.table[match(bcv.err$CV,CV.type.table$CV),"CVtype"]

# Save combined results table
write.csv(bcv.err, "BCV Combined validation table.csv", row.names=F)

# Plot the results
# RMSE
png(file="Plots/CV Results by RMSE.png")
ggplot(bcv.err, aes(x=CVlong, y=RMSE, col=Model)) + 
  geom_point() + 
  facet_grid(~CVtype, drop=T, scales="free", space="free") +
  labs(x="", y="RMSE") +
  theme(legend.position="right", axis.text.x = element_text(angle = 90))
dev.off()
# AUC
png(file="Plots/CV Results by AUC.png")
ggplot(bcv.err, aes(x=CVlong, y=AUC, col=Model)) + 
  geom_point() + 
  facet_grid(~CVtype, drop=T, scales="free", space="free") +
  labs(x="", y="AUC") +
  theme(legend.position="right", axis.text.x = element_text(angle = 90))
dev.off()



## PART 6a - Measure geographic distances between training and testing points
## --------------------------------------------------------------------------

# Make an output directory
dir.create("Distances")

# subsample of 5000 points (for computational efficiency)
sub <- sample(1:nrow(tr),5000)
sub <- sub[order(sub)]
tr.sub <- tr[sub,]

# LOOP for all cross-validations
# Different method for LOO, resub, and other CVs
for(b in bcv.list){
  
  # LOO
  if(substr(b,1,3)=="LOO"){
    
    # Set the buffer size
    buffer <- as.numeric(substr(b,15,nchar(b)-3))
    
    # Calculate a distance matrix (caution: 5000x5000 matrix is large)
    loo.m <- as.matrix(dist(tr.sub[c("x","y")]/1000))
    near <- apply(loo.m, 1, function(x){min(x[x>buffer])})
    head(near)
    dist.tab <- data.frame(ID=tr.sub$ID, Fold=1, n=length(near), GeogDist=near)
    
  }

  # Resubstitution (all distances should = 0)!!!
  if(substr(b,1,3)=="Res"){
    
    # Calculate distances
    near <- sqrt(ann(as.matrix(tr.sub[c("y","x")]/1000), as.matrix(tr.sub[c("y","x")]/1000), k=1, verbose=F)$knnIndexDist[,2])

    # Make distance table
    dist.tab <- data.frame(ID=tr.sub$ID, Fold=1, GeogDist=near)

  }
  
  # All other CV methods
  if(!substr(b,1,3) %in% c("LOO","Res")){
    
    # Make a list of folds in the CV
    folds <- read.csv(paste0("CV/Fold List - ",b,".csv"))[,2]
    fold.list <- unique(folds[sub])[order(unique(folds[sub]))]
    folds.sub <- folds[sub]
    
    # LOOP for each fold
    for(i in fold.list){
      
      # Make output table (first fold loop only)
      if(i==fold.list[1]){dist.tab <- data.frame(tr.sub["ID"], Fold=folds.sub, GeogDist=NA)}
      
      # Model training data (XY in km)
      near.train <- tr.sub[folds.sub!=i,c("y","x")]/1000
      # Model testing data (XY in km)
      near.test  <- tr.sub[folds.sub==i,c("y","x")]/1000
      # Number of training points in subsample
      #n.p <- nrow(near.train)
      
      # Calculate minimum distances (km) using approximate nearest neighbour
      near <- sqrt(ann(as.matrix(near.train), as.matrix(near.test), k=1, verbose=F)$knnIndexDist[,2])
     
      # Populate distance table
      dist.tab[dist.tab$Fold==i, "GeogDist"] <- round(near,5)

    } # LOOP for each fold
    
    # Save the distance table
    write.csv(dist.tab, paste0("Distances/Min Geog Dist Table - ",b,".csv"), row.names=F)
    
  }
  
  # Make (first method loop only) or populate the summary table (with average minimum distances)
  if(b==bcv.list[1]){
    dist.fold.sum <- data.frame(CVlong=b, ddply(dist.tab, .(Fold), summarise, n=length(ID), avg_GeogDist=mean(GeogDist)))
  } else {
    dist.fold.sum <- rbind(dist.fold.sum, data.frame(CVlong=b, ddply(dist.tab, .(Fold), summarise, n=length(ID), avg_GeogDist=mean(GeogDist))))
  }
  
} # LOOP for each CV method

# Add CV types
dist.fold.sum$CV <- CV.type.table[match(dist.fold.sum$CVlong,CV.type.table$CVlong),"CV"]
dist.fold.sum$CVtype <- CV.type.table[match(dist.fold.sum$CVlong,CV.type.table$CVlong),"CVtype"]

# Save the summary table (by folds)
write.csv(dist.fold.sum, "Distances/Average Min Geog Dist (by Fold).csv", row.names=F)

# Make a summary table by method (weighted by n in each fold)
dist.cv.sum <- ddply(dist.fold.sum, .(CVlong), summarise, nFolds=length(Fold), avg_GeogDist=weighted.mean(avg_GeogDist,n))

# Save the summary table (by method)
write.csv(dist.cv.sum, "Distances/Average Min Geog Dist (by CV).csv", row.names=F)



## PART 6b - Measure environmental distances between training and testing points
## -----------------------------------------------------------------------------

# NOTE: this calculation is based on:
#   1) the same subsample as the spatial distance measurements (tr.sub) in 5a above, and
#   2) the same PCA data weighted by variance explained (pc.adj) in 3d above

# Make the subsampled data with adjusted PCA values
pca.sub <- data.frame(tr[sub,c("ID","x","y")], pc.adj[sub,])

# LOOP for all cross-validations
# Different method for LOO, resub, and other CVs
for(b in bcv.list){
  
  # LOO
  if(substr(b,1,3)=="LOO"){
    
    # Set the buffer size
    buffer <- as.numeric(substr(b,15,nchar(b)-3))

    # Calculate a distance matrix (caution: 5000x5000 matrix is large)
    loo.xy.m <- as.matrix(dist(pca.sub[c("x","y")]/1000))
    loo.env.m <- as.matrix(dist(pca.sub[modvars]))
    for(i in 1:nrow(loo.xy.m)){
      if(i==1){near <- c()}
      near[i] <- min(loo.env.m[which(loo.xy.m[,i]>buffer),i] )
    }
    env.tab <- data.frame(ID=tr.sub$ID, Fold=1, n=length(near), EnvDist=near)
    
  }

  # Resubstitution (all distances should = 0)!!!
  if(substr(b,1,3)=="Res"){
    
    # Calculate distances
    near <- sqrt(ann(as.matrix(pca.sub[modvars]), as.matrix(pca.sub[modvars]), k=1, verbose=F)$knnIndexDist[,2])

    # Make distance table
    env.tab <- data.frame(ID=tr.sub$ID, Fold=1, EnvDist=near)
    
  }
  
  # All other CV methods
  if(!substr(b,1,3) %in% c("LOO","Res")){
    
    # Make a list of folds in the CV
    folds <- read.csv(paste0("CV/Fold List - ",b,".csv"))[,2]
    fold.list <- unique(folds[sub])[order(unique(folds[sub]))]
    folds.sub <- folds[sub]
    
    # LOOP for each fold
    for(i in fold.list){
      
      # Make output table (first fold loop only)
      if(i==fold.list[1]){env.tab <- data.frame(pca.sub["ID"], Fold=folds.sub, EnvDist=NA)}
      
      # Model training data (XY in km)
      near.train <- pca.sub[folds.sub!=i,modvars]
      # Model testing data (XY in km)
      near.test  <- pca.sub[folds.sub==i,modvars]
      # Number of training points in subsample
      #n.p <- nrow(near.train)
      
      # Calculate minimum distances (km) using approximate nearest neighbour
      near <- sqrt(ann(as.matrix(near.train), as.matrix(near.test), k=1, verbose=F)$knnIndexDist[,2])
      
      # Populate distance table
      env.tab[env.tab$Fold==i, "EnvDist"] <- round(near,5)
      
    } # LOOP for each fold
    
    # Save the distance table
    write.csv(env.tab, paste0("Distances/Min Enviro Dist Table - ",b,".csv"), row.names=F)
    
  }
  
  # Make (first method loop only) or populate the summary table (with average minimum distances)
  if(b==bcv.list[1]){
    env.fold.sum <- data.frame(CVlong=b, ddply(env.tab, .(Fold), summarise, n=length(ID), avg_EnvDist=mean(EnvDist)))
  } else {
    env.fold.sum <- rbind(env.fold.sum, data.frame(CVlong=b, ddply(env.tab, .(Fold), summarise, n=length(ID), avg_EnvDist=mean(EnvDist))))
  }
  
} # LOOP for each CV method

# Add CV types
env.fold.sum$CV <- CV.type.table[match(env.fold.sum$CVlong,CV.type.table$CVlong),"CV"]
env.fold.sum$CVtype <- CV.type.table[match(env.fold.sum$CVlong,CV.type.table$CVlong),"CVtype"]

# Save the summary table (by folds)
write.csv(env.fold.sum, "Distances/Average Min Enviro Dist (by Fold).csv", row.names=F)

# Make a summary table by method (weighted by n in each fold)
env.cv.sum <- ddply(env.fold.sum, .(CVlong), summarise, nFolds=length(Fold), avg_EnvDist=weighted.mean(avg_EnvDist,n))

# Save the summary table (by method)
write.csv(dist.cv.sum, "Distances/Average Min Enviro Dist (by CV).csv", row.names=F)



## PART 7 - Compare distances (geog and env) to error estimates
## ------------------------------------------------------------

# Make the plots for AUC
plotdat <- merge(bcv.err,dist.cv.sum[c("CVlong","avg_GeogDist")], by="CVlong", all.x=T, all.y=F, sort=F)
plotdat <- merge(plotdat,env.cv.sum[c("CVlong","avg_EnvDist")], by="CVlong", all.x=T, all.y=F, sort=F)
# AUC Plot
png(file="Plots/AUC vs Geog and Enviro Distance.png", width=640, height=320)
p.auc.geog <- ggplot(plotdat[plotdat$Model=="Random Forest",], aes(x=avg_GeogDist,y=AUC,col=CVtype)) +
  geom_point() +
  stat_smooth(method="lm", level=0, lwd=0.5) +
  labs(x="Average minimum geographic distance") + 
  ggtitle("Random Forest") + theme(legend.position="bottom")
p.auc.env <- ggplot(plotdat[plotdat$Model=="Random Forest",], aes(x=avg_EnvDist,y=AUC,col=CVtype)) +
  geom_point() +
  stat_smooth(method="lm", level=0, lwd=0.5) +
  labs(x="Average minimum environmental distance") + 
  ggtitle("Random Forest") + theme(legend.position="bottom")
grid.arrange(p.auc.geog,p.auc.env, ncol=2)
dev.off()















