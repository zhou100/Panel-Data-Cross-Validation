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
####  BOX 3 - Blocking for genetic relatedness
####  ----------------------------------------
####
####  CONTENTS:
####  ---------
####  1. Set up the simulation and build data
####  2. Set up the cross-validation (CV)
####  3. Run the cross validations
####       a. Resubstitution CV
####       b. Random k-fold CV
####       c. Blocked k-fold CV
####       d. Buffered leave-one-out CV
####       e. Generalised Least Squares (GLS)
####   4. Summary tables and plots
####   5. Measuring autocorrelation
####       a. Semivariograms
####       b. Correlograms



####  PACKAGES AND FUNCTIONS

# Phylogenetics
#install.packages('phytools')
library(phytools)

# Stepwise model selection
#install.packages('MASS')
library(MASS)

# Model dredging
#install.packages('MuMIn')
library(MuMIn)

# Cluster analysis for blocked folds
#install.packages('cluster')
library(cluster)

# GLS models
#install.packages('nlme')
library(nlme)

# NMDS (for autocorrelation)
#install.packages('MASS')
library(MASS)

# Variograms
#install.packages('geoR')
library(geoR)

# Correlograms
#install.packages('ncf')
library(ncf)

# Data management
#install.packages('plyr')
library(plyr)

# Plotting
#install.packages('ggplot2')
library(ggplot2)
#install.packages('RColorBrewer')
library(RColorBrewer)



## PART 1 - Set up the simulation and build data
## ---------------------------------------------

# Set the number of nodes (i.e. hypothetical species)
n.nodes <- 50

# Set the number of simulation runs
n.sim <- 100

# Create a directory to save the simulation data
dir.create("Sample Data", showWarnings=F)

# Make the phylogenetic tree simulations
tree <- pbtree(n=n.nodes, tip.label=paste0("N",1:n.nodes), scale=1, nsim=n.sim)

# Individual trees can be visualised by plotting
# See the phytools documentation for more info
plot(tree[1])

# Tree data can be saved for reproducibility
save(tree, file=paste0("Sample Data/Trees (",n.nodes," nodes, ",n.sim," trees).RData"))
#load(paste0("Sample Data/Trees (",n.nodes," nodes, ",n.sim," trees).RData"))

# LOOP to generate the trait values and environmental variable values for each simulation
for(b in 1:n.sim){
  
  # Brownian motion (BM) simulation to generate the error terms based on phylogenetic structure
  P <- fastBM(tree[[b]], mu=0, sig2=1)
  
  # Generate random environmental values for each species
  # Values are ordered to simulate structure in phylogenetic tree
  E <- runif(n.nodes, min=0, max=90)
  E <- E[order(E)]
  
  # Generate 'true model' trait values
  R0 <- 25 + (1.5 * E) -  (0.01 * E^2)
  
  # Generate trait values with the BM autocorrelated error added
  R <- 25 + (1.5 * E) -  (0.01 * E^2) + P*10
  
  # Make lists of values to use in the analysis
  if(b==1){
    P.list <- list(Tree1=P)
    E.list <- list(Tree1=E)
    R.list <- list(Tree1=R)
    R0.list <- list(Tree1=R0)
  } else {
    P.list[[b]] <- P
    E.list[[b]] <- E
    R.list[[b]] <- R
    R0.list[[b]] <- R0
  }#if
  
}#for

# Lists of values can be saved for reproducibility (our data is provided)
save(P.list, file=paste0("Sample Data/BM value list (",n.nodes," nodes, ",n.sim," trees).RData"))
save(E.list, file=paste0("Sample Data/Enviro value list (",n.nodes," nodes, ",n.sim," trees).RData"))
save(R.list, file=paste0("Sample Data/Trait value list (",n.nodes," nodes, ",n.sim," trees).RData"))
save(R0.list, file=paste0("Sample Data/True model value list (",n.nodes," nodes, ",n.sim," trees).RData"))

# To use the publication data, un-comment and run this section of code.
# Uses 100 simulated trees of 50 nodes each.
# load("Appendix_6_Box_3_DATA_Trees.RData")
# load("Appendix_6_Box_3_DATA_BM_value_list.RData")         # P.list
# load("Appendix_6_Box_3_DATA_Enviro_value_list.RData")     # E.list
# load("Appendix_6_Box_3_DATA_Trait_value_list.RData")      # R.list
# load("Appendix_6_Box_3_DATA_True_model_value_list.RData") # R0.list
# P.list <- P.list[1:100]
# E.list <- E.list[1:100]
# R.list <- R.list[1:100]
# R0.list <- R0.list[1:100]

# The simulated trait values can be visualised for any simulation sample (e.g. 25)
b <- 1
plot(E.list[[b]],R.list[[b]],pch=20, xlab="Trait value (R)", ylab="Environmental variable (E)")
# The true model can be added as a fitted curve
lines(E.list[[b]],R0.list[[b]],col="red")



## PART 2 - Set up the cross-validation (CV)
## -----------------------------------------

# List the number of folds (k) to run in the CV.
# These k values will be used for the random and blocked CV.
fold.list <- c(2,5)

# List the buffered leave-one-out (LOO) radii to run in the CV.
# These are based on using scale=1 in the phylogenetic tree generation.
loo.list <- c(0,0.25,0.5,0.75,1)

# Create the output directory for the cross-validation results
dir.create("Results", showWarnings=F)
dir.create(paste0("Results/",n.nodes," nodes, ",n.sim," trees/"), showWarnings=F)
results.path <- paste0("Results/",n.nodes," nodes, ",n.sim," trees/")

# Check the expected OVERALL ERROR
TrueRes.list <- mapply('-', R.list,R0.list, SIMPLIFY=FALSE)
for(i in 1:n.sim){
  if(i==1){TrueRMSE <- data.frame(TREE=1:n.sim, TrueRMSE=NA)}
  TrueRMSE$TrueRMSE[i] <- sqrt(mean(TrueRes.list[[i]]^2)) 
}
write.csv(TrueRMSE, paste0(results.path,"True Model Residuals (",n.nodes," nodes, ",n.sim," trees).csv"), row.names=F)
#TrueRMSE <- read.csv(paste0(results.path,"True Model Residuals (",n.nodes," nodes, ",n.sim," trees).csv"))

# Plot the expected errors (RMSE)
hist.rmse <- hist(TrueRMSE$TrueRMSE, col="grey75", border="grey75", axes=T, xlab="True model RMSE", main="")
abline(v=median(TrueRMSE$TrueRMSE), lty=2); abline(v=mean(TrueRMSE$TrueRMSE), lty=2)
text(max(hist.rmse$breaks),max(hist.rmse$counts), labels=paste0("Mean = ",round(mean(TrueRMSE$TrueRMSE),2)), adj=1)
text(max(hist.rmse$breaks),max(hist.rmse$counts)*0.9, labels=paste0("Median = ",round(median(TrueRMSE$TrueRMSE),2)), adj=1)

# make the list of alternative model formulae (for model selection)
form.list <- list(
  form1 <- "R ~ E",
  form2 <- "R ~ E + I(E^2)",
  form3 <- "R ~ E + I(E^2) + I(E^3)",
  form4 <- "R ~ E + I(E^2) + I(E^3) + I(E^4)",
  form5 <- "R ~ E + I(E^2) + I(E^3) + I(E^4) + I(E^5)",
  form6 <- "R ~ E + I(E^2) + I(E^3) + I(E^4) + I(E^5) + I(E^6)",
  form7 <- "R ~ E + I(E^2) + I(E^3) + I(E^4) + I(E^5) + I(E^6) + I(E^7)",
  form8 <- "R ~ E + I(E^2) + I(E^3) + I(E^4) + I(E^5) + I(E^6) + I(E^7) + I(E^8)"
)



## PART 3a - Resubstitution CV
## ---------------------------

# LOOP for all simulations
for(b in 1:n.sim){
  
  # Pull the data for the simulation from the lists of values 
  P <- P.list[[b]]
  E <- E.list[[b]]
  R <- abs(R.list[[b]])
  
  # Create the output table (first simulation only)
  if(b==1){resub.tab <- data.frame(expand.grid(CV="Resub", MODEL=1:10, FoldBuf=1, TREE=1:n.sim, AIC=NA, RMSE=NA, COEF.0=NA, COEF.1=NA, COEF.2=NA, COEF.3=NA, COEF.4=NA, COEF.5=NA, COEF.6=NA, COEF.7=NA, COEF.8=NA), nParam=c(2:9,NA,NA))}
  
  # Fit a GLM for each of the 8 model structures
  # Record the AIC, RMSE, and coefficients for each model in the output table
  # LOOP for each model structure
  for(i in 1:8){
    resub.mod <- glm(as.character(form.list[i]))
    resub.tab[resub.tab$TREE==b,c("AIC")][i] <- round(AIC(resub.mod),3)
    resub.tab[resub.tab$TREE==b,c("RMSE")][i] <- round(sqrt(mean((R-predict(resub.mod))^2)),3)
    resub.tab[resub.tab$TREE==b & resub.tab$MODEL==i,paste0("COEF.",0:i)] <- resub.mod$coefficients
  }
  
  # Run the stepwise model selection
  resub.mod <- stepAIC(lm(as.character(form.list[8])), trace=0, direction="both")
  # Record stepwise-selected model in the output table 
  resub.tab[resub.tab$TREE==b,c("AIC")][9] <- round(AIC(resub.mod),3)
  resub.tab[resub.tab$TREE==b,c("RMSE")][9] <- round(sqrt(mean((R-predict(resub.mod))^2)),3)
  mc <- names(resub.mod$coefficients)[-1]
  if(sum(mc=="E")>0){mc[which(mc=="E")] <- "(E1)"}
  mc <- c(0, as.numeric(substr(mc,nchar(mc)-1,nchar(mc)-1)))
  resub.tab[resub.tab$TREE==b & resub.tab$MODEL==9,paste0("COEF.",mc)] <- resub.mod$coefficients
  resub.tab[resub.tab$TREE==b,c("nParam")][9] <- length(resub.mod$coefficients)
  
  # Run the dredging model selection
  resub.mod <- get.models(dredge(lm(as.character(form.list[8]),na.action="na.fail")), subset=1)[[1]]
  # Record dredging-selected model in the output table 
  resub.tab[resub.tab$TREE==b,c("AIC")][10] <- round(AIC(resub.mod),3)
  resub.tab[resub.tab$TREE==b,c("RMSE")][10] <- round(sqrt(mean((R-predict(resub.mod))^2)),3)
  mc <- names(resub.mod$coefficients)[-1]
  if(sum(mc=="E")>0){mc[which(mc=="E")] <- "(E1)"}
  mc <- c(0, as.numeric(substr(mc,nchar(mc)-1,nchar(mc)-1)))
  resub.tab[resub.tab$TREE==b & resub.tab$MODEL==10,paste0("COEF.",mc)] <- resub.mod$coefficients
  resub.tab[resub.tab$TREE==b,c("nParam")][10] <- length(resub.mod$coefficients)
  
} # FOR all simulations

# Save the Resubstitution output table
write.csv(resub.tab, paste0(results.path,"Full output data (",n.nodes," nodes, ",n.sim," trees) - Resub.csv"), row.names=F)

# Summarise the resubstitution output data by model
resub.by.mod <- ddply(resub.tab, .(MODEL), summarise, nTree=max(TREE), CV="Resub",
                      p50.nParam=median(nParam),  # Median number of model parameters
                      x.nParam=mean(nParam),      # Average number of model parameters
                      p50.AIC=median(AIC),        # Median AIC
                      x.AIC=mean(AIC),            # Average AIC
                      sd.AIC=sd(AIC),             # Standard deviation AIC
                      p50.RMSE=median(RMSE),      # Median RMSE
                      x.RMSE=mean(RMSE),          # Average RMSE
                      sd.RMSE=mean(RMSE))        # Standard deviation RMSE

# Save the summary-by-model data
write.csv(resub.by.mod, paste0(results.path,"Summary by Model (",n.nodes," nodes, ",n.sim," trees) - Resub.csv"), row.names=F)

# Summarise the resubstitution output data by simulation tree
resub.by.tree <- ddply(resub.tab, .(TREE), summarise,  CV="Resub",
                       best.model.AIC=which.min(AIC),    # Best model by AIC
                       best.model.RMSE=which.min(RMSE))  # Best model by RMSE

# Number of model parameters (incl intercept) for each model selection
resub.by.tree$nStep <- rowSums(!is.na(resub.tab[resub.tab$MODEL==9, paste0("COEF.",1:8)]))
resub.by.tree$nDredge <- rowSums(!is.na(resub.tab[resub.tab$MODEL==10, paste0("COEF.",1:8)]))
resub.by.tree$nAIC <- with(resub.by.tree, 
                           ifelse(best.model.AIC<9, best.model.AIC+1,
                                  ifelse(best.model.AIC==9, nStep, nDredge)))
resub.by.tree$nRMSE <- with(resub.by.tree, 
                            ifelse(best.model.RMSE<9, best.model.RMSE+1, 
                                   ifelse(best.model.RMSE==9, nStep, nDredge)))

# Save the summary-by-tree data
write.csv(resub.by.tree, paste0(results.path,"Summary by Tree (",n.nodes," nodes, ",n.sim," trees) - Resub.csv"), row.names=F)



## PART 3b - Random k-fold CV
## --------------------------

# LOOP for all simulations
for(b in 1:n.sim){
  
  # Pull the data for the simulation from the lists of values 
  P <- P.list[[b]]
  E <- E.list[[b]]
  R <- abs(R.list[[b]])
  
  # Create the output table (first simulation only)
  if(b==1){random.tab <- data.frame(expand.grid(CV=c("Random"), MODEL=1:8, FoldBuf=fold.list, TREE=1:n.sim, AIC=NA, sd.AIC=NA, RMSE=NA), nParam=c(2:9), min.DIST=NA, p50.DIST=NA, avg.DIST=NA)}
  
  # LOOP for all k-fold structures (i.e. number of folds)
  for(n.folds in fold.list){
    
    # Make fold assignments (while statement to ensure sample selects at least one node from each fold)
    RandomFold <- c()
    while(length(unique(RandomFold))<n.folds){RandomFold <- sample(1:n.folds,n.nodes,replace=T)}
    
    # Measure the phylogenetic distance (average and median) between folds
    tree.dist <- round(cophenetic.phylo(tree[[b]]),3)
    for(i in 1:n.folds){
      if(i==1){ran.dist <- data.frame(FOLD=1:n.folds, min.DIST=NA, p50.DIST=NA, avg.DIST=NA)}
      ran.dist[i,c("min.DIST","p50.DIST")] <- quantile(tree.dist[RandomFold==i,RandomFold!=i],c(0,0.5))
      ran.dist[i,c("avg.DIST")] <- mean(tree.dist[RandomFold==i,RandomFold!=i])
    }
    
    # Populate the random CV results table with the phylogenetic distances between folds
    random.tab[random.tab$TREE==b & random.tab$FoldBuf==n.folds,c("min.DIST","p50.DIST","avg.DIST")] <- matrix(rep(apply(ran.dist[-1],2,mean),each=8),nrow=8)
    
    # Fit a GLM for each of the 8 model structures
    # Record the AIC, RMSE, and coefficients for each model in the output table
    # LOOP for each model structure
    for(i in 1:8){
      
      # LOOP for each randomly assigned fold
      for(f in 1:n.folds){
        
        # Build the output table (first fold only)
        if(f==1){
          ran.cv <- data.frame(NODE=1:n.nodes, FOLD=RandomFold, PRED=NA)
          ran.aic <- data.frame(FOLD=1:n.folds, N=as.numeric(table(RandomFold)), AIC=NA)
        }
        
        # Make the data, create the GLM, calculate AIC and RMSE
        ran.train <- data.frame(E=E[RandomFold!=f], R=R[RandomFold!=f])
        ran.test  <- data.frame(E=E[RandomFold==f], R=R[RandomFold==f])
        ran.mod <- glm(as.character(form.list[i]), data=ran.train)
        ran.cv[ran.cv$FOLD==f,"PRED"] <- predict(ran.mod, ran.test, type="response")
        ran.aic$AIC[f] <- AIC(ran.mod)
        
      } # FOR each fold
      
      # Populate the output table
      random.tab[random.tab$TREE==b & random.tab$FoldBuf==f,c("AIC")][i] <- round(mean(ran.aic$AIC),2)
      random.tab[random.tab$TREE==b & random.tab$FoldBuf==f,c("sd.AIC")][i] <- round(sd(ran.aic$AIC),5)
      random.tab[random.tab$TREE==b & random.tab$FoldBuf==f,c("RMSE")][i] <- round(sqrt(mean((R-ran.cv$PRED)^2)),3)
      
    } # FOR each model structure
    
  } # FOR each fold structure
  
} # FOR each simulation

# Save the random k-fold output table
write.csv(random.tab, paste0(results.path,"Full output data (",n.nodes," nodes, ",n.sim," trees) - Random.csv"), row.names=F)

# Summarise the random k-fold output data by model
random.by.mod <- ddply(random.tab[random.tab$MODEL<=8,], .(FoldBuf,MODEL), summarise,
                       nTree=max(TREE,na.rm=T), CV="Random",   
                       p50.nParam=median(nParam,na.rm=T),       # Median number of model parameters
                       x.nParam=mean(nParam,na.rm=T),           # Average number of model parameters
                       p50.AIC=median(AIC[AIC!=-Inf],na.rm=T),  # Median AIC
                       x.AIC=mean(AIC[AIC!=-Inf],na.rm=T),      # Average AIC
                       sd.AIC=sd(AIC[AIC!=-Inf],na.rm=T),       # Standard deviation AIC
                       p50.RMSE=median(RMSE,na.rm=T),           # Median RMSE
                       x.RMSE=mean(RMSE,na.rm=T),               # Average RMSE
                       sd.RMSE=sd(RMSE,na.rm=T))                # Standard deviation RMSE

# Save the summary-by-model data
write.csv(random.by.mod, paste0(results.path,"Summary by Model (",n.nodes," nodes, ",n.sim," trees) - Random.csv"), row.names=F)

# Summarise the random k-fold output data by simulation tree
random.by.tree <- ddply(random.tab[random.tab$MODEL<=8,], .(FoldBuf,TREE), summarise, CV="Random",
                        min.DIST=round(mean(min.DIST,na.rm=T),3),  # Avg minimum distance between folds
                        p50.DIST=round(mean(p50.DIST,na.rm=T),3),  # Avg median distance between folds
                        avg.DIST=round(mean(avg.DIST,na.rm=T),3),  # Avg average distance between folds
                        best.model.AIC=which.min(AIC),             # Best model by AIC
                        best.model.RMSE=which.min(RMSE))           # Best model by RMSE

# Number of model parameters (incl intercept) based on AIC or RMSE
random.by.tree$nAIC <- random.by.tree$best.model.AIC+1
random.by.tree$nRMSE <- random.by.tree$best.model.RMSE+1

# Save the summary-by-tree data
write.csv(random.by.tree, paste0(results.path,"Summary by Tree (",n.nodes," nodes, ",n.sim," trees) - Random.csv"), row.names=F)



## PART 3c - Blocked k-fold CV
## ---------------------------

# LOOP for all simulations
for(b in 1:n.sim){
  
  # Pull the data for the simulation from the lists of values 
  P <- P.list[[b]]
  E <- E.list[[b]]
  R <- abs(R.list[[b]])
  
  # Create the output table (first simulation only)
  if(b==1){block.tab <- data.frame(expand.grid(CV="Blocked", MODEL=1:8, FoldBuf=fold.list, TREE=1:n.sim, AIC=NA, sd.AIC=NA, RMSE=NA), nParam=c(2:9), min.DIST=NA, p50.DIST=NA, avg.DIST=NA)}
  
  # LOOP for all k-fold structures (i.e. number of folds)
  for(n.folds in fold.list){
    
    # Measure phylogenetic distance between all nodes
    tree.dist <- round(cophenetic.phylo(tree[[b]]),3)
    
    # Use clustering to define distance-based folds (see 'cluster' package documentation for more info)
    BlockFold <- pam(tree.dist, k=n.folds)$clustering
    
    # Measure phylogenetic distances between all folds
    block.dist <- data.frame(FOLD=1:n.folds, min.DIST=NA, p50.DIST=NA, avg.DIST=NA)
    for(i in 1:n.folds){
      block.dist[i,c("min.DIST","p50.DIST")] <- quantile(tree.dist[BlockFold==i,BlockFold!=i],c(0,0.5))
      block.dist[i,c("avg.DIST")] <- mean(tree.dist[BlockFold==i,BlockFold!=i])
    }
    
    # Populate the blocked output data with phylogenetic distances between folds
    block.tab[block.tab$TREE==b & block.tab$FoldBuf==n.folds,c("min.DIST","p50.DIST","avg.DIST")] <- matrix(rep(apply(block.dist[-1],2,mean),each=8),nrow=8)
    
    # Fit a GLM for each of the 8 model structures
    # Record the AIC, RMSE, and coefficients for each model in the output table
    # LOOP for each model structure
    for(i in 1:8){
      
      # LOOP for each blocked fold     
      for(f in 1:n.folds){
        
        # Build the output table (first fold only)
        if(f==1){
          block.cv <- data.frame(NODE=1:n.nodes,FOLD=BlockFold,PRED=NA)
          block.aic <- data.frame(FOLD=1:n.folds, N=as.numeric(table(BlockFold)), AIC=NA)
        }
        
        # Make the data, create the GLM, calculate AIC and RMSE
        block.train <- data.frame(E=E[BlockFold!=f], R=R[BlockFold!=f])
        block.test  <- data.frame(E=E[BlockFold==f], R=R[BlockFold==f])
        block.mod <- glm(as.character(form.list[i]), data=block.train)
        block.cv[block.cv$FOLD==f,"PRED"] <- predict(block.mod, block.test, type="response")
        block.aic$AIC[f] <- AIC(block.mod)
        
      } # FOR each fold
      
      # Populate the output table
      block.tab[block.tab$TREE==b & block.tab$FoldBuf==f,c("AIC")][i] <- round(mean(block.aic$AIC),2)
      block.tab[block.tab$TREE==b & block.tab$FoldBuf==f,c("sd.AIC")][i] <- round(sd(block.aic$AIC),5)
      block.tab[block.tab$TREE==b & block.tab$FoldBuf==f,c("RMSE")][i] <- round(sqrt(mean((R-block.cv$PRED)^2)),3)
      
    } # FOR each model structure
    
  } # FOR each fold structure
  
} # FOR each simulation

# Save the blocked k-fold output table
write.csv(block.tab, paste0(results.path,"Full output data (",n.nodes," nodes, ",n.sim," trees) - Blocked.csv"), row.names=F)

# Summarise the blocked k-fold output data by model
block.by.mod <- ddply(block.tab[block.tab$MODEL<=8,], .(FoldBuf,MODEL), summarise, 
                      nTree=max(TREE,na.rm=T), CV="Blocked",   
                      p50.nParam=median(nParam,na.rm=T),       # Median number of model parameters
                      x.nParam=mean(nParam,na.rm=T),           # Average number of model parameters
                      p50.AIC=median(AIC[AIC!=-Inf],na.rm=T),  # Median AIC
                      x.AIC=mean(AIC[AIC!=-Inf],na.rm=T),      # Average AIC
                      sd.AIC=sd(AIC[AIC!=-Inf],na.rm=T),       # Standard deviation AIC
                      p50.RMSE=median(RMSE,na.rm=T),           # Median RMSE
                      x.RMSE=mean(RMSE,na.rm=T),               # Average RMSE
                      sd.RMSE=sd(RMSE,na.rm=T))                # Standard deviation RMSE

# Save the summary-by-model data
write.csv(block.by.mod, paste0(results.path,"Summary by Model (",n.nodes," nodes, ",n.sim," trees) - Blocked.csv"), row.names=F)

# Summarise the blocked k-fold output data by simulation tree
block.by.tree <- ddply(block.tab[block.tab$MODEL<=8,], .(FoldBuf,TREE), summarise, CV="Blocked",
                       min.DIST=round(mean(min.DIST,na.rm=T),3),  # Avg minimum distance between folds
                       p50.DIST=round(mean(p50.DIST,na.rm=T),3),  # Avg median distance between folds
                       avg.DIST=round(mean(avg.DIST,na.rm=T),3),  # Avg average distance between folds
                       best.model.AIC=which.min(AIC),             # Best model by AIC
                       best.model.RMSE=which.min(RMSE))           # Best model by RMSE

# Number of model parameters (incl intercept) based on AIC or RMSE
block.by.tree$nAIC <- block.by.tree$best.model.AIC+1
block.by.tree$nRMSE <- block.by.tree$best.model.RMSE+1

# Save the summary-by-tree data
write.csv(block.by.tree, paste0(results.path,"Summary by Tree (",n.nodes," nodes, ",n.sim," trees) - Blocked.csv"), row.names=F)



## PART 3d - Buffered leave-one-out CV
## -----------------------------------

# LOOP for all simulations
for(b in 1:n.sim){
  
  # Pull the data for the simulation from the lists of values 
  P <- P.list[[b]]
  E <- E.list[[b]]
  R <- abs(R.list[[b]])
  
  # Create the output table (first simulation only)
  if(b==1){loo.tab <- data.frame(expand.grid(CV=c("LOO"), MODEL=1:8, FoldBuf=loo.list, TREE=1:n.sim, nTrain=NA, AIC=NA, sd.AIC=NA, RMSE=NA), nParam=c(2:9), min.DIST=NA, p50.DIST=NA, avg.DIST=NA)}
  
  # LOOP for all LOO buffer sizes
  for(buffer in loo.list){
    
    # Measure phylogenetic distance between all nodes
    tree.dist <- round(cophenetic.phylo(tree[[b]]),3)
    
    # Measure phylogenetic distances between validation and training points 
    loo.dist <- data.frame(NODE=1:n.nodes, min.DIST=NA, p50.DIST=NA, avg.DIST=NA)
    for(i in 1:n.nodes){
      loo.dist[i,c("min.DIST","p50.DIST")] <- quantile(tree.dist[i,tree.dist[i,]>buffer],c(0,0.5))
      loo.dist[i,c("avg.DIST")] <- mean(tree.dist[i,tree.dist[i,]>buffer])
    }
    
    # Populate the LOO output data with phylogenetic distances
    loo.tab[loo.tab$TREE==b & loo.tab$FoldBuf==buffer,c("min.DIST","p50.DIST","avg.DIST")] <- matrix(rep(apply(loo.dist[-1],2,mean),each=8),nrow=8)
    
    # Fit a GLM for each of the 8 model structures
    # Record the AIC, RMSE, and coefficients for each model in the output table
    # LOOP for each model structure
    for(i in 1:8){
      
      # Loop for each validation point (each node)
      for(f in 1:n.nodes){
        
        # Build the output table (first point only)
        if(f==1){
          loo.cv <- data.frame(NODE=1:n.nodes,PRED=NA)
          loo.aic <- data.frame(NODE=1:n.nodes, N=NA, AIC=NA)
        }
        
        # Make the data, create the GLM, calculate AIC and RMSE
        loo.train <- data.frame(E=E[which(tree.dist[f,]>buffer)], R=R[tree.dist[f,]>buffer])
        loo.test  <- data.frame(E=E[f], R=R[f])
        loo.mod <- glm(as.character(form.list[i]), data=loo.train)
        loo.cv[f,"PRED"] <- predict(loo.mod, loo.test, type="response")
        loo.aic[f,"N"] <- nrow(loo.train)
        loo.aic[f,"AIC"] <- AIC(loo.mod)
        
      } # FOR each fold
      
      # Populate the output table
      loo.tab[loo.tab$TREE==b & loo.tab$FoldBuf==buffer,"AIC"][i] <- round(mean(loo.aic$AIC),2)
      loo.tab[loo.tab$TREE==b & loo.tab$FoldBuf==buffer,"sd.AIC"][i] <- round(sd(loo.aic$AIC),5)
      loo.tab[loo.tab$TREE==b & loo.tab$FoldBuf==buffer,"RMSE"][i] <- round(sqrt(mean((R-loo.cv$PRED)^2)),3)
      loo.tab[loo.tab$TREE==b & loo.tab$FoldBuf==buffer,"nTrain"][i] <- round(mean(loo.aic$N),2)
      
    } # FOR each model structure
    
  } # FOR each buffer size
  
} # FOR each simulation


# Save the LOO output table
write.csv(loo.tab, paste0(results.path,"Full output data (",n.nodes," nodes, ",n.sim," trees) - LOO.csv"), row.names=F)

# Summarise the LOO output data by model
loo.by.mod <- ddply(loo.tab[loo.tab$MODEL<=8,], .(FoldBuf,MODEL), summarise,
                    nTree=max(TREE,na.rm=T), CV="LOO",
                    p50.nParam=median(nParam,na.rm=T),       # Median number of model parameters
                    x.nParam=mean(nParam,na.rm=T),           # Average number of model parameters
                    p50.AIC=median(AIC[AIC!=-Inf],na.rm=T),  # Median AIC
                    x.AIC=mean(AIC[AIC!=-Inf],na.rm=T),      # Average AIC
                    sd.AIC=sd(AIC[AIC!=-Inf],na.rm=T),       # Standard deviation AIC
                    p50.RMSE=median(RMSE,na.rm=T),           # Median RMSE
                    x.RMSE=mean(RMSE,na.rm=T),               # Average RMSE
                    sd.RMSE=sd(RMSE,na.rm=T))                # Standard deviation RMSE

# Save the summary-by-model data
write.csv(loo.by.mod, paste0(results.path,"Summary by Model (",n.nodes," nodes, ",n.sim," trees) - LOO.csv"), row.names=F)

# Summarise the LOO output data by simulation tree
loo.by.tree <- ddply(loo.tab[loo.tab$MODEL<=8,], .(FoldBuf,TREE), summarise,  CV="LOO",
                     min.DIST=round(mean(min.DIST,na.rm=T),3),  # Avg minimum distance between folds
                     p50.DIST=round(mean(p50.DIST,na.rm=T),3),  # Avg median distance between folds
                     avg.DIST=round(mean(avg.DIST,na.rm=T),3),  # Avg average distance between folds
                     best.model.AIC=which.min(AIC),             # Best model by AIC
                     best.model.RMSE=which.min(RMSE))           # Best model by RMSE

# Number of model parameters (incl intercept) based on AIC or RMSE
loo.by.tree$nAIC <- loo.by.tree$best.model.AIC+1
loo.by.tree$nRMSE <- loo.by.tree$best.model.RMSE+1

# Save the summary-by-tree data
write.csv(loo.by.tree, paste0(results.path,"Summary by Tree (",n.nodes," nodes, ",n.sim," trees) - LOO.csv"), row.names=F)



## PART 3e - Generalised Least Squares (GLS)
## -----------------------------------------

# LOOP for all simulations
for(b in 1:n.sim){
  
  # Pull the data for the simulation from the lists of values 
  P <- P.list[[b]]
  E <- E.list[[b]]
  R <- abs(R.list[[b]])
  
  # Create the output table (first simulation only)
  if(b==1){gls.tab <- data.frame(expand.grid(CV="GLS", MODEL=1:10, FoldBuf=1, TREE=1:n.sim, AIC=NA, RMSE=NA, COEF.0=NA, COEF.1=NA, COEF.2=NA, COEF.3=NA, COEF.4=NA, COEF.5=NA, COEF.6=NA, COEF.7=NA, COEF.8=NA), nParam=c(2:9,NA,NA))}
  
  # Fit a GLM for each of the 8 model structures
  # Must use Method="ML" to use stepAIC.
  # Record the AIC, RMSE, and coefficients for each model in the output table
  # LOOP for each model structure
  for(i in 1:8){
    gls.mod <- gls(as.formula(as.character(form.list[i])), correlation=corBrownian(phy=tree[[b]]), method="ML")
    gls.tab[gls.tab$TREE==b,c("AIC")][i] <- round(AIC(gls.mod),3)
    gls.tab[gls.tab$TREE==b,c("RMSE")][i] <- round(sqrt(mean((R-predict(gls.mod))^2)),3)
    gls.tab[gls.tab$TREE==b & gls.tab$MODEL==i,paste0("COEF.",0:i)] <- gls.mod$coefficients
  }
  
  # Run the stepwise model selection
  gls.mod <- stepAIC(gls(as.formula(as.character(form.list[8])), correlation=corBrownian(phy=tree[[b]]), method="ML"), trace=0, direction="both")
  # Record stepwise-selected model in the output table 
  gls.tab[gls.tab$TREE==b,c("AIC")][9] <- round(AIC(gls.mod),3)
  gls.tab[gls.tab$TREE==b,c("RMSE")][9] <- round(sqrt(mean((R-predict(gls.mod))^2)),3)
  mc <- names(gls.mod$coefficients)[-1]
  if(sum(mc=="E")>0){mc[which(mc=="E")] <- "(E1)"}
  mc <- c(0, as.numeric(substr(mc,nchar(mc)-1,nchar(mc)-1)))
  gls.tab[gls.tab$TREE==b & gls.tab$MODEL==9,paste0("COEF.",mc)] <- gls.mod$coefficients
  gls.tab[gls.tab$TREE==b,c("nParam")][9] <- length(gls.mod$coefficients)
  
  # Run the dredging model selection
  gls.mod <- get.models(dredge(gls(as.formula(as.character(form.list[8])), correlation=corBrownian(phy=tree[[b]]), method="ML", na.action="na.fail")), subset=1)[[1]]
  # Record dredging-selected model in the output table 
  gls.tab[gls.tab$TREE==b,c("AIC")][10] <- round(AIC(gls.mod),3)
  gls.tab[gls.tab$TREE==b,c("RMSE")][10] <- round(sqrt(mean((R-predict(gls.mod))^2)),3)
  mc <- names(gls.mod$coefficients)[-1]
  if(sum(mc=="E")>0){mc[which(mc=="E")] <- "(E1)"}
  mc <- c(0, as.numeric(substr(mc,nchar(mc)-1,nchar(mc)-1)))
  gls.tab[gls.tab$TREE==b & gls.tab$MODEL==10,paste0("COEF.",mc)] <- gls.mod$coefficients
  gls.tab[gls.tab$TREE==b,c("nParam")][10] <- length(gls.mod$coefficients)
  
} # FOR all simulations

# Save the GLS output table
write.csv(gls.tab, paste0(results.path,"Full output data (",n.nodes," nodes, ",n.sim," trees) - GLS.csv"), row.names=F)

# Summarise the GLS output data by model
gls.by.mod <- ddply(gls.tab, .(MODEL), summarise, nTree=max(TREE), CV="Resub",
                    p50.nParam=median(nParam),  # Median number of model parameters
                    x.nParam=mean(nParam),      # Average number of model parameters
                    p50.AIC=median(AIC),        # Median AIC
                    x.AIC=mean(AIC),            # Average AIC
                    sd.AIC=sd(AIC),             # Standard deviation AIC
                    p50.RMSE=median(RMSE),      # Median RMSE
                    x.RMSE=mean(RMSE),          # Average RMSE
                    sd.RMSE=mean(RMSE))         # Standard deviation RMSE

# Save the summary-by-model data
write.csv(gls.by.mod, paste0(results.path,"Summary by Model (",n.nodes," nodes, ",n.sim," trees) - GLS.csv"), row.names=F)

# Summarise the GLS output data by simulation tree
gls.by.tree <- ddply(gls.tab, .(TREE), summarise,  CV="GLS",
                     best.model.AIC=which.min(AIC),    # Best model by AIC
                     best.model.RMSE=which.min(RMSE))  # Best model by RMSE

# Number of model parameters (incl intercept) for each model selection
gls.by.tree$nStep <- rowSums(!is.na(gls.tab[gls.tab$MODEL==9,paste0("COEF.",1:8)]))
gls.by.tree$nDredge <- rowSums(!is.na(gls.tab[gls.tab$MODEL==10,paste0("COEF.",1:8)]))
gls.by.tree$nAIC <- with(gls.by.tree, 
                         ifelse(best.model.AIC<9, best.model.AIC+1,
                                ifelse(best.model.AIC==9, nStep, nDredge)))
gls.by.tree$nRMSE <- with(gls.by.tree,
                          ifelse(best.model.RMSE<9, best.model.RMSE+1,
                                 ifelse(best.model.RMSE==9, nStep, nDredge)))

# Save the summary-by-tree data
write.csv(gls.by.tree, paste0(results.path,"Summary by Tree (",n.nodes," nodes, ",n.sim," trees) - GLS.csv"), row.names=F)



## Part 4 - Summary Tables and Plots
## ---------------------------------

# Combine the output tables
clean.tab <- rbind.fill(resub.tab, gls.tab, random.tab, block.tab, loo.tab)[c("CV","MODEL","FoldBuf","TREE","AIC","RMSE","nParam","min.DIST","p50.DIST","avg.DIST","nTrain")]
# Set resubstitution and GLS phylo distances to zero
clean.tab[clean.tab$CV=="Resub" | clean.tab$CV=="GLS",c("min.DIST","p50.DIST","avg.DIST")] <- 0
write.csv(clean.tab, paste0(results.path,"^ Clean output data (",n.nodes," nodes, ",n.sim," trees).csv"), row.names=F)
#clean.tab <- read.csv(paste0(results.path,"^ Clean output data (",n.nodes," nodes, ",n.sim," trees).csv"))

# Tabulate the selected models based on AIC and RMSE criteria
best.mod <- ddply(clean.tab, .(CV,FoldBuf,TREE), summarise, 
                  BestAIC=which.min(AIC),             # AIC-selected model
                  BestAIC.n=nParam[which.min(AIC)],   # N parameters in AIC-selected model
                  BestAIC.AIC=min(AIC),               # AIC of AIC-selected model
                  BestAIC.RMSE=RMSE[which.min(AIC)],  # RMSE of AIC-selected model
                  BestRMSE=which.min(RMSE),           # RMSE-selected model
                  BestRMSE.n=nParam[which.min(RMSE)], # N parameters in RMSE-selected model
                  BestRMSE.AIC=AIC[which.min(RMSE)],  # AIC of RMSE-selected model
                  BestRMSE.RMSE=min(RMSE))            # RMSE of RMSE-selected model
# Rename models #9 and #10 to stepwise- or dredge-selected 
best.mod$BestAIC <- ifelse(best.mod$BestAIC<9, best.mod$BestAIC,
                           ifelse(best.mod$BestAIC==9, "Step", "Dredge"))
best.mod$BestRMSE <- ifelse(best.mod$BestRMSE<9, best.mod$BestRMSE,
                            ifelse(best.mod$BestRMSE==9, "Step", "Dredge"))

# Summarise the selected models with the appropriate criteria
# LM and GLS models should be AIC-selected
# Cross-Validated models (Resubstitution, Random, Blocked, and LOO) should be RMSE-selected
select.mod.aic <- with(best.mod[best.mod$CV %in% c("LM","GLS"),], data.frame(CV=CV, FoldBuf=FoldBuf, Group=paste(CV,FoldBuf), TREE=TREE, BestParam=BestAIC.n, BestRMSE=BestAIC.RMSE))
select.mod.rmse <- with(best.mod[!best.mod$CV %in% c("LM","GLS"),], data.frame(CV=CV, FoldBuf=FoldBuf, Group=paste(CV,FoldBuf), TREE=TREE, BestParam=BestRMSE.n, BestRMSE=BestRMSE.RMSE))
select.mod <- rbind(select.mod.aic,select.mod.rmse)

# Save the selected model summary table
write.csv(select.mod, paste0(results.path,"Summary table of selected models by AIC and RMSE.csv"), row.names=F)

# Plot the RMSE of the selected models
# True RMSE (median shown as vertical line) is included for comparison
rmse.plot <- rbind(select.mod, data.frame(CV="Truth", FoldBuf=0, Group="Truth", TREE=TrueRMSE$TREE, BestParam=3, BestRMSE=TrueRMSE$TrueRMSE))
rmse.plot$Group <- factor(rmse.plot$Group, levels=c("Truth","GLS 1","Resub 1","Random 2","Random 5","Blocked 2","Blocked 5","LOO 0","LOO 0.25","LOO 0.5","LOO 0.75","LOO 1"))
rmse.plot$CV <- factor(rmse.plot$CV, levels=c("Truth","GLS","Resub","Random","Blocked","LOO"))
# Colour palette for plotting
group.cols <- c("black","#358A3F","#9DD352",                       # Truth, GLS, Resub
                colorRampPalette(brewer.pal(7,"PuRd"))(5)[4:3],    # Random
                colorRampPalette(brewer.pal(7,"YlOrBr"))(5)[5:4],  # Blocked
                colorRampPalette(brewer.pal(5,"YlGnBu"))(6)[2:6])  # SLOO
ggplot(rmse.plot, aes(x=BestRMSE, col=Group)) +
  geom_vline(xintercept=median(TrueRMSE$TrueRMSE), lty=2) +
  geom_density() +
  facet_grid(CV~.) +
  scale_x_continuous("RMSE of selected model", limits=c(0,15)) +
  scale_color_manual(values=group.cols)

# Summarise the selected models for each tree and CV approach
# LM and GLS select best model by AIC
# Cross-validations select best model by RMSE
best.n.param <- rbind(
  ddply(clean.tab[clean.tab$CV %in% c("Resub","GLS"),], .(CV,FoldBuf,TREE), summarise, 
        BestCriteria="AIC",
        nParam=nParam[which.min(AIC)]),
  ddply(clean.tab[clean.tab$CV %in% c("Random","Blocked","LOO"),], .(CV,FoldBuf,TREE), summarise,
        BestCriteria="RMSE",
        nParam=nParam[which.min(RMSE)]))
best.n.param$Group <- paste(best.n.param$CV,best.n.param$FoldBuf,sep="_")
best.n.param[best.n.param$Group=="Resub_1","Group"] <- "LM"
best.n.param[best.n.param$Group=="GLS_1","Group"] <- "GLS"
best.n.param$CV <- factor(best.n.param$CV, levels=c("Resub","GLS","Random","Blocked","LOO"))
best.n.param$Group <- factor(best.n.param$Group, levels=c("LM","GLS",paste0("Random_",rev(fold.list)),paste0("Blocked_",rev(fold.list)),paste0("LOO_",loo.list)))

# Boxplots of selected model parameter counts
boxplot(best.n.param$nParam~best.n.param$Group, las=2, col="grey90")
abline(h=3, col="blue") # True model n parameters = 3



## Part 5a - Measuring autocorrelation with semivariograms
## -------------------------------------------------------

# LOOP for all simulations
for(b in 1:n.sim){
  
  # Extract the data values for the tree
  E <- E.list[[b]]
  R <- R.list[[b]]
  
  #  Phylogenetic distance matrix of the nodes
  tree.dist <- round(cophenetic.phylo(tree[[b]]),3)
  
  # NMDS to calculate hypothetical coordinates of the nodes
  # http://stackoverflow.com/questions/17271509/coordinates-from-distance-matrix-in-r
  mds <- cmdscale(tree.dist, k=2)
  x <- mds[,1]
  y <- mds[,2]
  
  # Coordinates should be stretched across known max distance (=2)
  dist.factor <- 2/max(dist(cbind(x,y)))
  x <- x * dist.factor
  y <- y * dist.factor
  
  # Explain with the true (data generating) model
  mod.form <- as.formula(R ~ E + I(E^2))
  full.mod <- glm(mod.form)
  
  # Calculate model residuals
  res <- R - predict(full.mod, type="response")
  
  # Create the variogram & the summary table (for plotting)
  v <- variog(coords=cbind(x,y), data=res, breaks=seq(0.1,2,0.1), messages=F)
  v.sum <- data.frame(Dist=v$u, SV=v$v, N=v$n)
  
  # Compile all variogram data into a single table
  if(b==1){v.out <- v.sum} else {v.out <- rbind(v.out,v.sum)}
}

# Save the semivariogram data
write.csv(v.out, paste0(results.path,"Autocorrelation by semivariogram (",n.sim," trees, ",n.nodes," nodes).csv"), row.names=F)

# Make plotting data (combine simulations)
v.plot <- ddply(v.out, .(Dist), summarise,
                n.rows=length(SV), # Number of distance intervals
                x.SV=mean(SV),     # Average semivariance within each interval
                x.n=mean(N),       # Average number of data points within the interval
                SD=sd(SV))         # SD of the semivariance within each interval
v.plot$SE <- v.plot$SD/sqrt(v.plot$n.rows) # Standard error

# Plot the semivariogram
v.plot$err.high <- v.plot$x.SV+v.plot$SE
v.plot$err.low <- v.plot$x.SV-v.plot$SE
plot(x.SV~Dist, data=v.plot, type="l", col="grey75",
     ylim=range(c(v.plot$err.high,v.plot$err.low)),
     xlab="Phylogenetic distance", ylab="Semivariance",
     main="Phylogenetic autocorrelation")
with(v.plot, arrows(Dist,err.low,Dist,err.high, length=0.05, angle=90, code=3, col="grey75"))
points(x.SV~Dist, data=v.plot, pch=20)



## Part 5b - Measuring autocorrelation using correlograms with Moran's I
## ---------------------------------------------------------------------

# Set the number of data resamplings (see ncf package documentation)
n.resample <- 1

# Set the increment for the distance classes
increment <- 0.05

# LOOP for all simulations
for(b in 1:n.sim){
  
  # Extract the data values for the tree
  E <- E.list[[b]]
  R <- R.list[[b]]
  
  #  Phylogenetic distance matrix of the nodes
  tree.dist <- round(cophenetic.phylo(tree[[b]]),3)
  
  # NMDS to calculate hypothetical coordinates of the nodes
  # http://stackoverflow.com/questions/17271509/coordinates-from-distance-matrix-in-r
  mds <- cmdscale(tree.dist, k=2)
  x <- mds[,1]
  y <- mds[,2]


  # Explain with the true (data generating) model
  mod.form <- as.formula(R ~ E + I(E^2))
  full.mod <- glm(mod.form)
  
  # Calculate residuals
  res <- R - predict(full.mod, type="response")


  # Build the correlogram data
  ncf.t <- correlog(x, y, res, increment=increment, resamp=n.resample, quiet=T)
  ncf.t <- data.frame(Tree=b, ncf.t$mean.of.class, ncf.t$correlation, ncf.t$p, ncf.t$n)
  names(ncf.t) <- c("Tree","Dist","Moran","p","n")
  if(b==1){ncf.out <- ncf.t} else {ncf.out <- rbind(ncf.out,ncf.t)}
}
  
# Save the correlogram data
write.csv(ncf.out, paste0(results.path,"Autocorrelation by correlogram (",n.sim," trees, ",n.nodes," nodes, ",n.resample," resamples, ",increment," increment).csv"), row.names=F)

# Make plotting data (combine simulations)
# Compiles statistics across all simulations
# Rounding makes distance classes
ncf.out$DistRound <- round(ncf.out$Dist,1)
ncf.plot <- ddply(ncf.out[-1], .(DistRound), summarise,
                  n.rows=length(Moran),  # Number of rows within the distance class
                  x.Moran=mean(Moran),   # Average Moran's I within the distance class
                  x.p=mean(p),           # Average p-value for distance class
                  x.n=mean(n),           # Average number of values within the distance class
                  SD=sd(Moran))          # SD of Moran's I within the distance class
ncf.plot$SE <- ncf.plot$SD/sqrt(ncf.plot$n.rows)  # Standard Error

# Plot the correlogram
# Note that correlations are unreliable with small samples (n.rows)
ncf.plot$err.high <- ncf.plot$x.Moran+ncf.plot$SE
ncf.plot$err.low <- ncf.plot$x.Moran-ncf.plot$SE
plot(x.Moran~DistRound, data=ncf.plot, type="l", col="grey75",
     ylim=range(c(ncf.plot$err.high,ncf.plot$err.low)),
     xlab="Phylogenetic distance", ylab="Moran's I",
     main="Phylogenetic autocorrelation")
with(ncf.plot, arrows(DistRound,err.low,DistRound,err.high, length=0.05, angle=90, code=3, col="grey75"))
points(x.Moran~DistRound, data=ncf.plot, pch=20)
