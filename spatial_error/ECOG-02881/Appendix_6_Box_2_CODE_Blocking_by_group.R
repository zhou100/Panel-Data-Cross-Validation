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
####
####  BOX 2 - Blocking for individuals and groups
####  -------------------------------------------
####
####  CONTENTS:
####  ---------
####  1. Import the elk data.
####  2. Collinearity and multi-collinearity check for chosen Resource Selection Functions' predictors.
####  3. Visualization of results of the sensitivity analys.
####  4. GLMM fit, normality of random effect check, and spatial residual autocorrelation analysis.
####  5. Model evaluation sensu Boyce et al. 2002 - Resubstitution (model trained and tested on the same data, no data splitting).  
####  6. Model evaluation sensu Boyce et al. 2002 - Random cross-validation.     
####  7. Model evaluation sensu Boyce et al. 2002 - Individual block cross-validation.     
####  8. Model evaluation sensu Boyce et al. 2002 - Spatially independent individual block cross-validation.     
####  9. Details on how individuals were allocated to folds: cluster analysis fit and visualization.     
####  10. Visualization of main results     

 

####  PACKAGES

# Data management 
#install.packages("reshape2")
library(reshape2)

# Generalized additive models
#install.packages("mgcv")
library(mgcv)

# GLMMs (needed to build resource selection functions)
#install.packages("lme4")
library(lme4)

# Plotting 
#install.packages("lattice")
library(lattice)

# Spatial residual autocorrelation
#install.packages("gstat")
library(gstat)
#install.packages("sp")
library(sp)

#### FUNCTIONS

# Variance Inflation Factor (VIF) functions to check for multi-collinearity among model predictors.  
# Citation: Mixed effects models and extensions in ecology with R. (2009).
#           Zuur, AF, Ieno, EN, Walker, N, Saveliev, AA, and Smith, GM.
#           Springer. http://www.highstat.com/book2.htm
corvif <- function(dataz) {
  dataz <- as.data.frame(dataz)
  form    <- formula(paste("fooy ~ ",paste(strsplit(names(dataz)," "),collapse=" + ")))
  dataz   <- data.frame(fooy=1,dataz)
  lm_mod  <- lm(form,dataz)
  cat("\n\nVariance inflation factors\n\n")
  print(myvif(lm_mod))
}
myvif <- function(mod) {
  v <- vcov(mod)
  assign <- attributes(model.matrix(mod))$assign
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  } else warning("No intercept: vifs may not be sensible.")
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("The model contains fewer than 2 terms")
  if (length(assign) > dim(v)[1] ) {
    diag(tmp_cor)<-0
    if (any(tmp_cor==1.0)){
      return("Sample size is too small, 100% collinearity is present")
    } else {
      return("Sample size is too small")
    }
  }
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) {
    result <- data.frame(GVIF=result[, 1])
  } else {
    result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  }
  invisible(result)
}

# Handy Variogram function to check for spatial residual autocorrelation.
# Citation: Mixed effects models and extensions in ecology with R. (2009).
#           Zuur, AF, Ieno, EN, Walker, N, Saveliev, AA, and Smith, GM.
#           Springer. http://www.highstat.com/book2.htm
MyVariogram <- function(x,y,z, MyDistance) {
  mydata      <- data.frame(z, x, y)
  coordinates(mydata)    <- c("x", "y")  
  Var <- variogram(z ~ 1, mydata, cutoff = MyDistance)
  data.frame(Var$np, Var$dist, Var$gamma)
}

# Make output directory where main plots are stored as png files (compatible with all operating systems).
dir.create("Results")


## PART 1 - Import the elk data 
## ----------------------------

# Import the data
load("Appendix_6_Box_2_DATA_db_2008.RData")

## These data were compiled and provided by:
##
##    Simone Ciuti, the Montane Elk Research Program, http://montaneelk.com/
##    Individuals: 43 adult female elk
##    Sample size: 27148 satellite elk relocations combined with random available locations.
##    Study period: July 1st to Aug 31st 2008.
##
## Data details
## [1] "ElkID"            Elk identity        
## [2] "Easting"          UTM coordinates NAD1983 UTM ZONE 11N      
## [3] "Northing"         UTM coordinates NAD1983 UTM ZONE 11N
## [4] "dem"              Digital elevation model
## [5] "cc_abmt"          Canopy cover model
## [6] "lc_30"            Land cover 30 m resolution (categorical, 8 levels)
## [7] "use"              Elk relocations (1) and random available points (0)
## [8] "lc_30f"           Reclassified land cover map 30 m resolution (4 levels)          ***
## [9] "new_herd"         Elk herd
## [10] "LogRugg"         Terrain Ruggedeness index (log-transformed)                     ***
## [11] "SqrtDistroad"    Distance to the closest road (sqrt-transformed)                 ***
## [12] "NDVI2"           Normalized Difference Vegetation Index (square-transformed)     ***
## 
##  *** Environmental predictors chosen for the resource selection models. See text for more details.



## PART 2 - Collinearity and multi-collinearity check for chosen Resource Selection Functions' predictors  
## ------------------------------------------------------------------------------------------------------

db_2008$lc <- as.numeric(db_2008$lc_30f)   # lc_30f converted to numeric for the collinearity check. 
db_2008 <- db_2008[!is.na(db_2008$lc), ]   # removed NAs (i.e., unclassified habitat types in the land cover map).
cor(db_2008[,10:13])                       # all correlations |rp| << 0.5.
corvif(db_2008[,10:13])                    # no issue with multicollinearity (VIF << 3)



## PART 3 - Visualization of results of the sensitivity analysis
## -------------------------------------------------------------

# Load beta estimates from the sensitivity analysis.
load(file = "Appendix_6_Box_2_DATA_all_betas.RData")
colnames(all_betas) <- c("simulation", 2:101)

## NOTE: Beta estimates were obtained by running a simple for-loop (which requires roughly 2 days to run with the elk dataset). The for-loop assigns a given number of random locations to each elk used location: ratios of used:available points were 1:0.05, 1:0.1, 1:0.5, 1:1, 1:2, . . . 1:13. Parameters for each use:available ratio have been estimated 100 times using the GLMM's structure specified in the main manuscript, each run with a new selection of random available points.

# Visualise the beta estimates calculated in the sensitivity analysis.
# Here we reported the estimates for the quadratic term for canopy cover as an example.
newdata <- melt(all_betas, id.vars="simulation")
png("Results/Beta estimate sensitivity.png")
plot(value ~ simulation, data=newdata, 
     main="", bty="n", col="black", bg=1, ylim=c(-0.40, 0), xaxt="n", las=2, pch=1, cex=0.8,
     ylab="Beta estimates (crown cover^2)", xlab="Number of available points per 'used' point") 
axis(1, at=c(0.10,1:13), labels=c(0.10, 1:13), las=1)
data <- data.frame(simulation=seq(0,13,0.5))
model <- gam(value ~ s(simulation), data=newdata, family=gaussian)
predictions <- predict(model, data, type="response", se=T)
lines(data$simulation, predictions$fit, col="grey", lty=2)
lines(data$simulation, predictions$fit + 1.96 * predictions$se.fit, col = "grey", lty=1)
lines(data$simulation, predictions$fit - 1.96 * predictions$se.fit, col = "grey", lty=1)
dev.off() # go to the 'Results' directory to visualized the plot.


## PART 4 - GLMM fit, normality of random effect check, and spatial residual autocorrelation analysis   
## --------------------------------------------------------------------------------------------------

# Define model formula
formula3 <- use ~ scale(LogRugg) + I(scale(LogRugg)^2) + scale(SqrtDistroad) + I(scale(SqrtDistroad)^2) + scale(NDVI2) + I(scale(NDVI2)^2) + lc_30f + scale(LogRugg) * scale(SqrtDistroad) + scale(NDVI2) * scale(SqrtDistroad) + (1|ElkID) 

#################################  SLOW STEP  ##################################
model3 <- glmer(formula3, db_2008, family=binomial)
save(model3, file = "model3.RData")
#load("model3.RData")

# Check the normal distribution of random effect
png("Results/random_effect_distribution.png")
r_int <- ranef(model3)$ElkID$`(Intercept)`
qqnorm(r_int, main = "Normal Q-Q plot of random intercepts (elk individuals)",  bty="n")
qqline(r_int, col = 2, lwd = 3, lty = 3)
(foo <- shapiro.test(r_int))
text(0,min(r_int), paste0("Shapiro-Wilk normality test: W = ",round(foo$statistic,3),", p = ",round(foo$p.value,3)), col = 2)
dev.off()

# Plot model intercept for each individual animal
png("Results/individual_intercepts.png")
dotplot(ranef(model3, condVar=TRUE))
dev.off()

# Setting the scene for the variograms - spatial autocorrelation of residuals
model3.proj <- residuals(model3)
res <- data.frame(db_2008[c("Easting","Northing", "use", "ElkID")],  RES_inter=model3.proj)

# Variogram for presence-available data
res <- data.frame(db_2008[c("Easting","Northing", "use", "ElkID")],  RES_inter=model3.proj)
set.seed(5)
res <- res[sample(1:nrow(res), 27000, replace=FALSE),] #select a subsample to speed up the analysis
newdata2a <- MyVariogram(res$Easting, res$Northing, res$RES_inter, MyDistance = 8000)

# Variogram for presence-only data
res <- data.frame(db_2008[c("Easting","Northing", "use", "ElkID")],  RES_inter=model3.proj)
res <- res[res$use == 1,]               
newdata2b <- MyVariogram(res$Easting, res$Northing, res$RES_inter, MyDistance = 8000)

# Plot variogram (presence-only)
# Plot variogram (presence-availability)
png("Results/Variograms for pres-avl and pres-only.png")
par(mfrow=c(2,1), mar=c(5,5,2,1)); layout.show(2)
# Presence-availability
plot(newdata2a$Var.dist, newdata2a$Var.gamma, bty="n", type="l", col="red", lwd=2, las=1,
     xlim=c(0,8000), ylim=c(0,0.6), xlab="Distance (m)", ylab="Semivariance", main="Variogram (presence-availability)")
# Presence-only
plot(newdata2b$Var.dist, newdata2b$Var.gamma, bty="n", type="l", col="red", lwd=2, las=1,
     xlim=c(0,8000), ylim=c(0,0.04), xlab="Distance (m)", ylab="Semivariance", main="Variogram (presence-only)")
dev.off()



## PART 5 - Model evaluation sensu Boyce et al. 2002 - Resubstitution
##          (Model trained and tested on the same data, no data splitting.)
## ------------------------------------------------------------------------

# Calculate the contribution to RSF scores by the levels of the categorical predictor (reference category: coniferous forest)
db_2008$lc_30_contribute <- as.numeric(ifelse(db_2008$lc_30f == "Coniferous forest", 0,
                                             ifelse(db_2008$lc_30f == "Deciduous forest", model3@beta[8],
                                                    ifelse(db_2008$lc_30f == "Mixed forest", model3@beta[9],
                                                           ifelse(db_2008$lc_30f == "Open areas", model3@beta[10],"ERROR")))))

# Calculate RSF scores assuming the exponential form.
db_2008$RSFscores1 <- exp(model3@beta[2] * scale(db_2008$LogRugg) +
                           model3@beta[3] * I(scale(db_2008$LogRugg)^2) +
                           model3@beta[4] * scale(db_2008$SqrtDistroad) +
                           model3@beta[5] * I(scale(db_2008$SqrtDistroad)^2) +
                           model3@beta[6] * scale(db_2008$NDVI2) +
                           model3@beta[7] * I(scale(db_2008$NDVI2)^2) +
                           db_2008$lc_30_contribute +
                           model3@beta[11] * scale(db_2008$LogRugg) * scale(db_2008$SqrtDistroad) +
                           model3@beta[12] * scale(db_2008$SqrtDistroad) * scale(db_2008$NDVI2))

# Run evaluation sensu Boyce et al. 2002
newdata <- db_2008[complete.cases(db_2008[,"RSFscores1"]),]
q.pp1 <- quantile(newdata$RSFscores1, probs=seq(0,1,.1))   # computing quantiles of RSF scores
bin1 <- rep(NA,length(newdata$RSFscores1))                 # binning RSF scores (10 bins)
for (i in 1:10){
  bin1[newdata$RSFscores1 >= q.pp1[i] & newdata$RSFscores1 < q.pp1[i+1]] = i
}
used1 <- newdata$use
table(used1, bin1)         # Occurrence of presence and available data by bin
a1 <- table(used1, bin1) 
a1 <- t(a1) 
a1 <- as.data.frame.matrix(a1)      # the next few lines compute area-adjusted frequency of categories (bins) of RSF scores 
a1$areaadjusted <- rep(NA,length(10))
sum0 <- sum(a1[,1])
sum1 <- sum(a1[,2])
a1$areaadjusted <- (a1[,2] / sum1 ) / (a1[,1] / sum0)
a1$bins <- seq(1, 10 ,by=1) ; a1
with(a1, cor.test(bins, areaadjusted, method="spearman"))    # Spearman correlation coefficient between RSF bin ranks and area-adjusted frequencies
png("Results/Binned RSF scores for the resubstitution.png")
with(a1, plot(bins, areaadjusted,  type="b", col = "red", cex=1.4, pch = 20))  
abline(h = 1, lty = 3)
dev.off()

rm(newdata)
plot_full = a1[,3] # storing plot data for final visualization (see below)



## PART 6 - Model evaluation sensu Boyce et al. 2002 - Random cross-validation
## ---------------------------------------------------------------------------
##
## 5-fold cross-validation with all data (GPS location fixes) assigned to folds randomly (i.e. all elk contribute with telemetry fixes to all folds), resulting in 5 folds with 20% of telemetry relocations each.     

# Split the data into 5 folds by simply selecting fixes randomly
set.seed(5)
db_2008$rand.vec <- sample(1:5,nrow(db_2008),replace=TRUE) 

# Fit the model in all folds but one.
#################################  SLOW STEP  ##################################
mod_inter_randomf1 <- glmer(formula3, db_2008[db_2008$rand.vec != 1,], family = binomial)
mod_inter_randomf2 <- glmer(formula3, db_2008[db_2008$rand.vec != 2,], family = binomial)
mod_inter_randomf3 <- glmer(formula3, db_2008[db_2008$rand.vec != 3,], family = binomial)
mod_inter_randomf4 <- glmer(formula3, db_2008[db_2008$rand.vec != 4,], family = binomial)
mod_inter_randomf5 <- glmer(formula3, db_2008[db_2008$rand.vec != 5,], family = binomial)
save(mod_inter_randomf1, file = "mod_inter_randomf1.RData")
save(mod_inter_randomf2, file = "mod_inter_randomf2.RData")
save(mod_inter_randomf3, file = "mod_inter_randomf3.RData")
save(mod_inter_randomf4, file = "mod_inter_randomf4.RData")
save(mod_inter_randomf5, file = "mod_inter_randomf5.RData")
# load(file = "mod_inter_randomf1.RData")
# load(file = "mod_inter_randomf2.RData")
# load(file = "mod_inter_randomf3.RData")
# load(file = "mod_inter_randomf4.RData")
# load(file = "mod_inter_randomf5.RData")

# Remove columns used for evaluations fit earlier
db_2008$lc_30_contribute <- NULL
db_2008$RSFscores <- NULL

# Calculate the contribution to RSF scores by the levels of the categorical predictor (reference category: coniferous forest) and calculate RSF scores assuming the exponential form. Repeating it for each withhold fold.

# NOTE THE ATTACH COMMAND HERE - DETACH BELOW
attach(db_2008)
db_2008$lc_30_contribute[rand.vec == 1] <- as.numeric(ifelse(db_2008$lc_30f[rand.vec == 1] == "Coniferous forest", 0,
                                                            ifelse(db_2008$lc_30f[rand.vec == 1] == "Deciduous forest", mod_inter_randomf1@beta[8],
                                                                   ifelse(db_2008$lc_30f[rand.vec == 1] == "Mixed forest", mod_inter_randomf1@beta[9],
                                                                          ifelse(db_2008$lc_30f[rand.vec == 1] == "Open areas", mod_inter_randomf1@beta[10], "ERROR")))))
db_2008$RSFscores[rand.vec == 1] <- exp(mod_inter_randomf1@beta[2] * scale(db_2008$LogRugg[rand.vec == 1]) +
                                         mod_inter_randomf1@beta[3] * I(scale(db_2008$LogRugg[rand.vec == 1])^2) +
                                         mod_inter_randomf1@beta[4] * scale(db_2008$SqrtDistroad[rand.vec == 1]) +
                                         mod_inter_randomf1@beta[5] * I(scale(db_2008$SqrtDistroad[rand.vec == 1])^2) +
                                         mod_inter_randomf1@beta[6] * scale(db_2008$NDVI2[rand.vec == 1]) +
                                         mod_inter_randomf1@beta[7] * I(scale(db_2008$NDVI2[rand.vec == 1])^2) +
                                         db_2008$lc_30_contribute[rand.vec == 1] +
                                         mod_inter_randomf1@beta[11] * scale(db_2008$LogRugg[rand.vec == 1]) * scale(db_2008$SqrtDistroad[rand.vec == 1]) +
                                         mod_inter_randomf1@beta[12] * scale(db_2008$SqrtDistroad[rand.vec == 1]) * scale(db_2008$NDVI2[rand.vec == 1]))

db_2008$lc_30_contribute[rand.vec == 2] <- as.numeric(ifelse(db_2008$lc_30f[rand.vec == 2] == "Coniferous forest", 0,
                                                            ifelse(db_2008$lc_30f[rand.vec == 2] == "Deciduous forest", mod_inter_randomf2@beta[8],
                                                                   ifelse(db_2008$lc_30f[rand.vec == 2] == "Mixed forest", mod_inter_randomf2@beta[9],
                                                                          ifelse(db_2008$lc_30f[rand.vec == 2] == "Open areas", mod_inter_randomf2@beta[10], "ERROR")))))
db_2008$RSFscores[rand.vec == 2] <- exp(mod_inter_randomf2@beta[2] * scale(db_2008$LogRugg[rand.vec == 2]) +
                                         mod_inter_randomf2@beta[3] * I(scale(db_2008$LogRugg[rand.vec == 2])^2) +
                                         mod_inter_randomf2@beta[4] * scale(db_2008$SqrtDistroad[rand.vec == 2]) +
                                         mod_inter_randomf2@beta[5] * I(scale(db_2008$SqrtDistroad[rand.vec == 2])^2) +
                                         mod_inter_randomf2@beta[6] * scale(db_2008$NDVI2[rand.vec == 2]) +
                                         mod_inter_randomf2@beta[7] * I(scale(db_2008$NDVI2[rand.vec == 2])^2) +
                                         db_2008$lc_30_contribute[rand.vec == 2] +
                                         mod_inter_randomf2@beta[11] * scale(db_2008$LogRugg[rand.vec == 2]) * scale(db_2008$SqrtDistroad[rand.vec == 2]) +
                                         mod_inter_randomf2@beta[12] * scale(db_2008$SqrtDistroad[rand.vec == 2]) * scale(db_2008$NDVI2[rand.vec == 2]))

db_2008$lc_30_contribute[rand.vec == 3] <- as.numeric(ifelse(db_2008$lc_30f[rand.vec == 3] == "Coniferous forest", 0,
                                                            ifelse(db_2008$lc_30f[rand.vec == 3] == "Deciduous forest", mod_inter_randomf3@beta[8],
                                                                   ifelse(db_2008$lc_30f[rand.vec == 3] == "Mixed forest", mod_inter_randomf3@beta[9],
                                                                          ifelse(db_2008$lc_30f[rand.vec == 3] == "Open areas", mod_inter_randomf3@beta[10], "ERROR")))))
db_2008$RSFscores[rand.vec == 3] <- exp(mod_inter_randomf3@beta[2] * scale(db_2008$LogRugg[rand.vec == 3]) +
                                         mod_inter_randomf3@beta[3] * I(scale(db_2008$LogRugg[rand.vec == 3])^2) +
                                         mod_inter_randomf3@beta[4] * scale(db_2008$SqrtDistroad[rand.vec == 3]) +
                                         mod_inter_randomf3@beta[5] * I(scale(db_2008$SqrtDistroad[rand.vec == 3])^2) +
                                         mod_inter_randomf3@beta[6] * scale(db_2008$NDVI2[rand.vec == 3]) +
                                         mod_inter_randomf3@beta[7] * I(scale(db_2008$NDVI2[rand.vec == 3])^2) +
                                         db_2008$lc_30_contribute[rand.vec == 3] +
                                         mod_inter_randomf3@beta[11] * scale(db_2008$LogRugg[rand.vec == 3]) * scale(db_2008$SqrtDistroad[rand.vec == 3]) +
                                         mod_inter_randomf3@beta[12] * scale(db_2008$SqrtDistroad[rand.vec == 3]) * scale(db_2008$NDVI2[rand.vec == 3]))

db_2008$lc_30_contribute[rand.vec == 4] <- as.numeric(ifelse(db_2008$lc_30f[rand.vec == 4] == "Coniferous forest", 0,
                                                            ifelse(db_2008$lc_30f[rand.vec == 4] == "Deciduous forest", mod_inter_randomf4@beta[8],
                                                                   ifelse(db_2008$lc_30f[rand.vec == 4] == "Mixed forest", mod_inter_randomf4@beta[9],
                                                                          ifelse(db_2008$lc_30f[rand.vec == 4] == "Open areas", mod_inter_randomf4@beta[10], "ERROR")))))
db_2008$RSFscores[rand.vec == 4] <- exp(mod_inter_randomf4@beta[2] * scale(db_2008$LogRugg[rand.vec == 4]) +
                                         mod_inter_randomf4@beta[3] * I(scale(db_2008$LogRugg[rand.vec == 4])^2) +
                                         mod_inter_randomf4@beta[4] * scale(db_2008$SqrtDistroad[rand.vec == 4]) +
                                         mod_inter_randomf4@beta[5] * I(scale(db_2008$SqrtDistroad[rand.vec == 4])^2) +
                                         mod_inter_randomf4@beta[6] * scale(db_2008$NDVI2[rand.vec == 4]) +
                                         mod_inter_randomf4@beta[7] * I(scale(db_2008$NDVI2[rand.vec == 4])^2) +
                                         db_2008$lc_30_contribute[rand.vec == 4] +
                                         mod_inter_randomf4@beta[11] * scale(db_2008$LogRugg[rand.vec == 4]) * scale(db_2008$SqrtDistroad[rand.vec == 4]) +
                                         mod_inter_randomf4@beta[12] * scale(db_2008$SqrtDistroad[rand.vec == 4]) * scale(db_2008$NDVI2[rand.vec == 4]))

db_2008$lc_30_contribute[rand.vec == 5] <- as.numeric(ifelse(db_2008$lc_30f[rand.vec == 5] == "Coniferous forest", 0,
                                                            ifelse(db_2008$lc_30f[rand.vec == 5] == "Deciduous forest", mod_inter_randomf5@beta[8],
                                                                   ifelse(db_2008$lc_30f[rand.vec == 5] == "Mixed forest", mod_inter_randomf5@beta[9],
                                                                          ifelse(db_2008$lc_30f[rand.vec == 5] == "Open areas", mod_inter_randomf5@beta[10], "ERROR")))))
db_2008$RSFscores[rand.vec == 5] <- exp(mod_inter_randomf5@beta[2] * scale(db_2008$LogRugg[rand.vec == 5]) +
                                         mod_inter_randomf5@beta[3] * I(scale(db_2008$LogRugg[rand.vec == 5])^2) +
                                         mod_inter_randomf5@beta[4] * scale(db_2008$SqrtDistroad[rand.vec == 5]) +
                                         mod_inter_randomf5@beta[5] * I(scale(db_2008$SqrtDistroad[rand.vec == 5])^2) +
                                         mod_inter_randomf5@beta[6] * scale(db_2008$NDVI2[rand.vec == 5]) +
                                         mod_inter_randomf5@beta[7] * I(scale(db_2008$NDVI2[rand.vec == 5])^2) +
                                         db_2008$lc_30_contribute[rand.vec == 5] +
                                         mod_inter_randomf5@beta[11] * scale(db_2008$LogRugg[rand.vec == 5]) * scale(db_2008$SqrtDistroad[rand.vec == 5]) +
                                         mod_inter_randomf5@beta[12] * scale(db_2008$SqrtDistroad[rand.vec == 5]) * scale(db_2008$NDVI2[rand.vec == 5]))
detach(db_2008)


# Run the k-fold CV evaluation sensu Boyce et al. 2002
dataset <- db_2008[complete.cases(db_2008[,"RSFscores"]),]
rho_model <- numeric(5) ## it will store Spearman's coefficients
fold <- subset(dataset,rand.vec==unique(dataset$rand.vec)[1]) # run the procedure for the 1st fold - the for-loop below will work on folds 2 to 5
q.pp <- quantile(fold$RSFscores,probs=seq(0,1,.1))   ## computing quantiles of RSF scores
# --------------------------------------------------------
bin <- rep(NA,length(fold$RSFscores))                ## binning RSF scores (10 bins)
for (j in 1:10){
  bin[fold$RSFscores>=q.pp[j]& fold$RSFscores<q.pp[j+1]] = j
}
used <- fold$use
# --------------------------------------------------------
a <- table(used,bin) ## Occurrence of presence and available data by bin
a <- t(a) #transpose the table
a <- as.data.frame.matrix(a) ## the next few lines compute area-adjusted frequency of categories (bins) of RSF scores 
a$areaadjusted <- rep(NA,length(10))
sum0 <- sum(a[,1])
sum1 <- sum(a[,2])
a$areaadjusted <- (a[,2] / sum1 ) / (a[,1] / sum0)
a$bins <- seq(1,10,by=1);a

# --------------------------------------------------------
rho_model[1] <- with(a,cor.test(bins,areaadjusted,method="spearm"))$estimate ## Spearman correlation coefficient between RSF bin ranks and area-adjusted frequencies
# --------------------------------------------------------

# Run the procedure for the other folds and plot the binned RSF scores from the cross-validation
png("Results/Binned RSF scores from the random CV.png", height=360, width=600)
par(oma=c(1,2,1,1)) 
par(mar=c(4.2,4.2,2,2))
with(a,plot(bins,areaadjusted, ylab="Area adjusted frequency", xlab="Binned RSF scores",xlim=c(1,10),ylim=c(0,2),type="b",cex=1.4,cex.lab=1.5,las=1,main="K-fold validation"))
abline(h = 1, lty = 3) ## Spearman correlation coefficient between RSF bin ranks and area-adjusted frequencies
plot_random <- data.frame(a[,3], b = NA, c = NA, d = NA, e = NA)
for (i in 2:5){
  fold <- subset(dataset,rand.vec==unique(dataset$rand.vec)[i]) ## run the procedure on folds 2 to 5
  # --------------------------------------------------------
  q.pp <- quantile(fold$RSFscores,probs=seq(0,1,.1))  ## computing quantiles of RSF scores
  # --------------------------------------------------------
  bin <- rep(NA,length(fold$RSFscores)) ## binning RSF scores (10 bins)
  for (j in 1:10){
    bin[fold$RSFscores>=q.pp[j]& fold$RSFscores<q.pp[j+1]] = j
  }
  used = fold$use
  # --------------------------------------------------------
  a <- table(used,bin) ## area adjusted freq in used/available for each bin
  a <- t(a) #transpose the table
  a <- as.data.frame.matrix(a) ## the next few lines compute area-adjusted frequency of categories (bins) of RSF scores 
  a$areaadjusted <- rep(NA,length(10))
  sum0 <- sum(a[,1])
  sum1 <- sum(a[,2])
  a$areaadjusted <- (a[,2] / sum1 ) / (a[,1] / sum0)
  a$bins <- seq(1,10,by=1);a
  # --------------------------------------------------------
  rho_model[i] <- with(a,cor.test(bins,areaadjusted,method="spearm"))$estimate  ## store Spearman correlation coefficients between RSF bin ranks and area-adjusted frequencies
  # --------------------------------------------------------
  par(oma=c(1,2,1,1)) 
  par(mar=c(4.2,4.2,2,2))
  with(a,points(bins,areaadjusted,col=i,type="b")) 
  plot_random[,i] = a[,3]
}
dev.off()

## store Spearman correlation coefficients that will be used for final plots below ##
Rho_random_fixes <- rho_model 



## PART 7 - Model evaluation sensu Boyce et al. 2002 - Individual block cross-validation
## -------------------------------------------------------------------------------------
##
## 5-fold cross-validation with data split by randomly assigning all GPS fixes from a single individual to a given fold, resulting in 3 folds with 9 individuals each and 2 folds with 8 individuals each. Home ranges of individuals assigned to different folds may overlap.    

# Remove columns used for evaluations fit earlier
db_2008$rand.vec <- NULL   

# split individuals randomly
newdata <- data.frame(ElkID = unique(db_2008$ElkID))
set.seed(5)
random_sample <- data.frame(ElkID = sample(newdata$ElkID, 43))
random_sample$rand.vec <- 0
random_sample$rand.vec[1:9] <- 1
random_sample$rand.vec[10:17] <- 2
random_sample$rand.vec[18:26] <- 3
random_sample$rand.vec[27:34] <- 4
random_sample$rand.vec[35:43] <- 5
with(random_sample, tapply(ElkID, rand.vec, length))
db_2008 <- merge(db_2008, random_sample, by = "ElkID", all.x = T )
with(db_2008, tapply(ElkID, rand.vec, unique))

# Fit the model in all folds but one.
#################################  SLOW STEP  ##################################
mod_inter_randomIND1 <- glmer(formula3, db_2008[db_2008$rand.vec != 1,], family = binomial)
mod_inter_randomIND2 <- glmer(formula3, db_2008[db_2008$rand.vec != 2,], family = binomial)
mod_inter_randomIND3 <- glmer(formula3, db_2008[db_2008$rand.vec != 3,], family = binomial)
mod_inter_randomIND4 <- glmer(formula3, db_2008[db_2008$rand.vec != 4,], family = binomial)
mod_inter_randomIND5 <- glmer(formula3, db_2008[db_2008$rand.vec != 5,], family = binomial)
save(mod_inter_randomIND1, file = "mod_inter_randomIND1.RData")
save(mod_inter_randomIND2, file = "mod_inter_randomIND2.RData")
save(mod_inter_randomIND3, file = "mod_inter_randomIND3.RData")
save(mod_inter_randomIND4, file = "mod_inter_randomIND4.RData")
save(mod_inter_randomIND5, file = "mod_inter_randomIND5.RData")
# load(file = "mod_inter_randomIND1.RData")
# load(file = "mod_inter_randomIND2.RData")
# load(file = "mod_inter_randomIND3.RData")
# load(file = "mod_inter_randomIND4.RData")
# load(file = "mod_inter_randomIND5.RData")

# Remove columns used for evaluations fit earlier
db_2008$lc_30_contribute <- NULL
db_2008$RSFscores <- NULL

# Calculation of the contribution to RSF scores by the levels of the categorical predictor (reference category: coniferous forest) and calculation of RSF scores assuming the exponential form. Repeating it for each withhold fold.

# NOTE THE ATTACH COMMAND HERE - DETACH BELOW
attach(db_2008)
db_2008$lc_30_contribute[rand.vec == 1] <- as.numeric(ifelse(db_2008$lc_30f[rand.vec == 1] == "Coniferous forest", 0,
                                                            ifelse(db_2008$lc_30f[rand.vec == 1] == "Deciduous forest", mod_inter_randomIND1@beta[8],
                                                                   ifelse(db_2008$lc_30f[rand.vec == 1] == "Mixed forest", mod_inter_randomIND1@beta[9],
                                                                          ifelse(db_2008$lc_30f[rand.vec == 1] == "Open areas", mod_inter_randomIND1@beta[10], "ERROR")))))
db_2008$RSFscores[rand.vec == 1] = exp(mod_inter_randomIND1@beta[2] * scale(db_2008$LogRugg[rand.vec == 1]) +
                                         mod_inter_randomIND1@beta[3] * I(scale(db_2008$LogRugg[rand.vec == 1])^2) +
                                         mod_inter_randomIND1@beta[4] * scale(db_2008$SqrtDistroad[rand.vec == 1]) +
                                         mod_inter_randomIND1@beta[5] * I(scale(db_2008$SqrtDistroad[rand.vec == 1])^2) +
                                         mod_inter_randomIND1@beta[6] * scale(db_2008$NDVI2[rand.vec == 1]) +
                                         mod_inter_randomIND1@beta[7] * I(scale(db_2008$NDVI2[rand.vec == 1])^2) +
                                         db_2008$lc_30_contribute[rand.vec == 1] +
                                         mod_inter_randomIND1@beta[11] * scale(db_2008$LogRugg[rand.vec == 1]) * scale(db_2008$SqrtDistroad[rand.vec == 1]) +
                                         mod_inter_randomIND1@beta[12] * scale(db_2008$SqrtDistroad[rand.vec == 1]) * scale(db_2008$NDVI2[rand.vec == 1]))

db_2008$lc_30_contribute[rand.vec == 2] <- as.numeric(ifelse(db_2008$lc_30f[rand.vec == 2] == "Coniferous forest", 0,
                                                            ifelse(db_2008$lc_30f[rand.vec == 2] == "Deciduous forest", mod_inter_randomIND2@beta[8],
                                                                   ifelse(db_2008$lc_30f[rand.vec == 2] == "Mixed forest", mod_inter_randomIND2@beta[9],
                                                                          ifelse(db_2008$lc_30f[rand.vec == 2] == "Open areas", mod_inter_randomIND2@beta[10], "ERROR")))))
db_2008$RSFscores[rand.vec == 2] <- exp(mod_inter_randomIND2@beta[2] * scale(db_2008$LogRugg[rand.vec == 2]) +
                                         mod_inter_randomIND2@beta[3] * I(scale(db_2008$LogRugg[rand.vec == 2])^2) +
                                         mod_inter_randomIND2@beta[4] * scale(db_2008$SqrtDistroad[rand.vec == 2]) +
                                         mod_inter_randomIND2@beta[5] * I(scale(db_2008$SqrtDistroad[rand.vec == 2])^2) +
                                         mod_inter_randomIND2@beta[6] * scale(db_2008$NDVI2[rand.vec == 2]) +
                                         mod_inter_randomIND2@beta[7] * I(scale(db_2008$NDVI2[rand.vec == 2])^2) +
                                         db_2008$lc_30_contribute[rand.vec == 2] +
                                         mod_inter_randomIND2@beta[11] * scale(db_2008$LogRugg[rand.vec == 2]) * scale(db_2008$SqrtDistroad[rand.vec == 2]) +
                                         mod_inter_randomIND2@beta[12] * scale(db_2008$SqrtDistroad[rand.vec == 2]) * scale(db_2008$NDVI2[rand.vec == 2]))

db_2008$lc_30_contribute[rand.vec == 3] <- as.numeric(ifelse(db_2008$lc_30f[rand.vec == 3] == "Coniferous forest", 0,
                                                            ifelse(db_2008$lc_30f[rand.vec == 3] == "Deciduous forest", mod_inter_randomIND3@beta[8],
                                                                   ifelse(db_2008$lc_30f[rand.vec == 3] == "Mixed forest", mod_inter_randomIND3@beta[9],
                                                                          ifelse(db_2008$lc_30f[rand.vec == 3] == "Open areas", mod_inter_randomIND3@beta[10], "ERROR")))))
db_2008$RSFscores[rand.vec == 3] <- exp(mod_inter_randomIND3@beta[2] * scale(db_2008$LogRugg[rand.vec == 3]) +
                                         mod_inter_randomIND3@beta[3] * I(scale(db_2008$LogRugg[rand.vec == 3])^2) +
                                         mod_inter_randomIND3@beta[4] * scale(db_2008$SqrtDistroad[rand.vec == 3]) +
                                         mod_inter_randomIND3@beta[5] * I(scale(db_2008$SqrtDistroad[rand.vec == 3])^2) +
                                         mod_inter_randomIND3@beta[6] * scale(db_2008$NDVI2[rand.vec == 3]) +
                                         mod_inter_randomIND3@beta[7] * I(scale(db_2008$NDVI2[rand.vec == 3])^2) +
                                         db_2008$lc_30_contribute[rand.vec == 3] +
                                         mod_inter_randomIND3@beta[11] * scale(db_2008$LogRugg[rand.vec == 3]) * scale(db_2008$SqrtDistroad[rand.vec == 3]) +
                                         mod_inter_randomIND3@beta[12] * scale(db_2008$SqrtDistroad[rand.vec == 3]) * scale(db_2008$NDVI2[rand.vec == 3]))

db_2008$lc_30_contribute[rand.vec == 4] <- as.numeric(ifelse(db_2008$lc_30f[rand.vec == 4] == "Coniferous forest", 0,
                                                            ifelse(db_2008$lc_30f[rand.vec == 4] == "Deciduous forest", mod_inter_randomIND4@beta[8],
                                                                   ifelse(db_2008$lc_30f[rand.vec == 4] == "Mixed forest", mod_inter_randomIND4@beta[9],
                                                                          ifelse(db_2008$lc_30f[rand.vec == 4] == "Open areas", mod_inter_randomIND4@beta[10], "ERROR")))))
db_2008$RSFscores[rand.vec == 4] <- exp(mod_inter_randomIND4@beta[2] * scale(db_2008$LogRugg[rand.vec == 4]) +
                                         mod_inter_randomIND4@beta[3] * I(scale(db_2008$LogRugg[rand.vec == 4])^2) +
                                         mod_inter_randomIND4@beta[4] * scale(db_2008$SqrtDistroad[rand.vec == 4]) +
                                         mod_inter_randomIND4@beta[5] * I(scale(db_2008$SqrtDistroad[rand.vec == 4])^2) +
                                         mod_inter_randomIND4@beta[6] * scale(db_2008$NDVI2[rand.vec == 4]) +
                                         mod_inter_randomIND4@beta[7] * I(scale(db_2008$NDVI2[rand.vec == 4])^2) +
                                         db_2008$lc_30_contribute[rand.vec == 4] +
                                         mod_inter_randomIND4@beta[11] * scale(db_2008$LogRugg[rand.vec == 4]) * scale(db_2008$SqrtDistroad[rand.vec == 4]) +
                                         mod_inter_randomIND4@beta[12] * scale(db_2008$SqrtDistroad[rand.vec == 4]) * scale(db_2008$NDVI2[rand.vec == 4]))

db_2008$lc_30_contribute[rand.vec == 5] <- as.numeric(ifelse(db_2008$lc_30f[rand.vec == 5] == "Coniferous forest", 0,
                                                            ifelse(db_2008$lc_30f[rand.vec == 5] == "Deciduous forest", mod_inter_randomIND5@beta[8],
                                                                   ifelse(db_2008$lc_30f[rand.vec == 5] == "Mixed forest", mod_inter_randomIND5@beta[9],
                                                                          ifelse(db_2008$lc_30f[rand.vec == 5] == "Open areas", mod_inter_randomIND5@beta[10], "ERROR")))))
db_2008$RSFscores[rand.vec == 5] <- exp(mod_inter_randomIND5@beta[2] * scale(db_2008$LogRugg[rand.vec == 5]) +
                                         mod_inter_randomIND5@beta[3] * I(scale(db_2008$LogRugg[rand.vec == 5])^2) +
                                         mod_inter_randomIND5@beta[4] * scale(db_2008$SqrtDistroad[rand.vec == 5]) +
                                         mod_inter_randomIND5@beta[5] * I(scale(db_2008$SqrtDistroad[rand.vec == 5])^2) +
                                         mod_inter_randomIND5@beta[6] * scale(db_2008$NDVI2[rand.vec == 5]) +
                                         mod_inter_randomIND5@beta[7] * I(scale(db_2008$NDVI2[rand.vec == 5])^2) +
                                         db_2008$lc_30_contribute[rand.vec == 5] +
                                         mod_inter_randomIND5@beta[11] * scale(db_2008$LogRugg[rand.vec == 5]) * scale(db_2008$SqrtDistroad[rand.vec == 5]) +
                                         mod_inter_randomIND5@beta[12] * scale(db_2008$SqrtDistroad[rand.vec == 5]) * scale(db_2008$NDVI2[rand.vec == 5]))
detach(db_2008)

# Run the k-fold CV evaluation sensu Boyce et al. 2002
dataset <- db_2008[complete.cases(db_2008[,"RSFscores"]),]
rho_model <- numeric(5) ## it will store Spearman's coefficients

fold <- subset(dataset,rand.vec==unique(dataset$rand.vec)[1]) # run the procedure for the 1st fold - the for-loop below will work on folds 2 to 5
q.pp <- quantile(fold$RSFscores,probs=seq(0,1,.1)) ## computing quantiles of RSF scores
# --------------------------------------------------------
bin <- rep(NA,length(fold$RSFscores))
for (j in 1:10){
  bin[fold$RSFscores>=q.pp[j]& fold$RSFscores<q.pp[j+1]] = j ## binning RSF scores (10 bins)
}
used <- fold$use
# --------------------------------------------------------
a <- table(used,bin) ## Occurrence of presence and available data by bin
a <- t(a) #transpose the table
a <- as.data.frame.matrix(a) ## the next few lines compute area-adjusted frequency of categories (bins) of RSF scores 
a$areaadjusted <- rep(NA,length(10))
sum0 <- sum(a[,1])
sum1 <- sum(a[,2])
a$areaadjusted <- (a[,2] / sum1 ) / (a[,1] / sum0)
a$bins <- seq(1,10,by=1);a

# --------------------------------------------------------
rho_model[1] <- with(a,cor.test(bins,areaadjusted,method="spearm"))$estimate ## Spearman correlation coefficient between RSF bin ranks and area-adjusted frequencies
# --------------------------------------------------------

# Run the procedure for the other folds and plot the binned RSF scores from the cross-validation
png("Results/Binned RSF scores from the individual CV.png", height=360, width=600)
par(oma=c(1,2,1,1)) 
par(mar=c(4.2,4.2,2,2))
with(a,plot(bins,areaadjusted, ylab="Area adjusted frequency", xlab="Binned RSF scores",xlim=c(1,10),ylim=c(0,2.5),type="b",cex=1.4,cex.lab=1.5,las=1,main="K-fold validation"))
abline(h = 1, lty = 3) ## Spearman correlation coefficient between RSF bin ranks and area-adjusted frequencies
plot_random_ind <- data.frame(a[,3], b = NA, c = NA, d = NA, e = NA)
for (i in 2:5){
  fold <- subset(dataset,rand.vec==unique(dataset$rand.vec)[i]) ## run the procedure on folds 2 to 5
  # --------------------------------------------------------
  q.pp <- quantile(fold$RSFscores,probs=seq(0,1,.1)) ## computing quantiles of RSF scores
  # --------------------------------------------------------
  bin <- rep(NA,length(fold$RSFscores))
  for (j in 1:10){
    bin[fold$RSFscores>=q.pp[j]& fold$RSFscores<q.pp[j+1]] = j  ## binning RSF scores (10 bins)
  }
  used <- fold$use
  # --------------------------------------------------------
  a <- table(used,bin) ## area adjusted freq in used/available for each bin
  a <- t(a) #transpose the table
  a <- as.data.frame.matrix(a) ## the next few lines compute area-adjusted frequency of categories (bins) of RSF scores 
  a$areaadjusted <- rep(NA,length(10))
  sum0 <- sum(a[,1])
  sum1 <- sum(a[,2])
  a$areaadjusted <- (a[,2] / sum1 ) / (a[,1] / sum0)
  a$bins <- seq(1,10,by=1);a
  # --------------------------------------------------------
  rho_model[i] <- with(a,cor.test(bins,areaadjusted,method="spearm"))$estimate ## store Spearman correlation coefficients between RSF bin ranks and area-adjusted frequencies
  # --------------------------------------------------------
  par(oma=c(1,2,1,1)) 
  par(mar=c(4.2,4.2,2,2))
  with(a,points(bins,areaadjusted,col=i,type="b"))  
  plot_random_ind[,i] = a[,3]
}
dev.off()

## store Spearman correlation coefficients that will be used for final plots below ##
Rho_random_individuals <- rho_model



## PART 8 - Model evaluation sensu Boyce et al. 2002 - Spatially independent individual block cross-validation
## -----------------------------------------------------------------------------------------------------------
## 
## 5-fold cross-validation with data split by spatially independent individuals. Data from each individual contribute only to one fold and individuals closer than 20 km are never allocated to the same fold. This results in 3 folds with 8 individuals each, 1 fold with 7 individuals, and 1 fold with 12 individuals. Home ranges of individuals assigned to different folds do not overlap

# Remove columns used for evaluations fit earlier
db_2008$rand.vec = NULL   

# Select spatially independent individuals (see PART 9 for how these animals have been selected)
db_2008$rand.vec = ifelse(db_2008$ElkID == "E019" | db_2008$ElkID == "E043"| db_2008$ElkID == "E080"| db_2008$ElkID == "E087"| db_2008$ElkID == "E070"| db_2008$ElkID == "E073"| db_2008$ElkID == "E038"| db_2008$ElkID == "E040", 1,ifelse(db_2008$ElkID == "E085"| db_2008$ElkID == "E027"| db_2008$ElkID == "E044"| db_2008$ElkID == "E021"| db_2008$ElkID == "E023"| db_2008$ElkID == "E026"| db_2008$ElkID == "E083", 2, ifelse(db_2008$ElkID == "E062"| db_2008$ElkID == "E018"| db_2008$ElkID == "E054"| db_2008$ElkID == "E058"| db_2008$ElkID == "E056"| db_2008$ElkID == "E090"| db_2008$ElkID == "E061"| db_2008$ElkID == "E079", 3, ifelse(db_2008$ElkID == "E057"| db_2008$ElkID == "E074"| db_2008$ElkID == "E077"| db_2008$ElkID == "E004"| db_2008$ElkID == "E053"| db_2008$ElkID == "E075"| db_2008$ElkID == "E086"| db_2008$ElkID == "E009"| db_2008$ElkID == "E016"| db_2008$ElkID == "E001"| db_2008$ElkID == "E003"| db_2008$ElkID == "E067", 4, ifelse(db_2008$ElkID == "E055"| db_2008$ElkID == "E066"| db_2008$ElkID == "E017"| db_2008$ElkID == "E059"| db_2008$ElkID == "E060"| db_2008$ElkID == "E005"| db_2008$ElkID == "E051"| db_2008$ElkID == "E052", 5, "ERROR")))))

## visualize sample size split by folds ##
with(db_2008, tapply(ElkID, rand.vec, unique))

# Fit the model in all folds but one.
#################################  SLOW STEP  ##################################
mod_inter_indIND1 <- glmer(formula3, db_2008[db_2008$rand.vec != 1,], family = binomial)
mod_inter_indIND2 <- glmer(formula3, db_2008[db_2008$rand.vec != 2,], family = binomial)
mod_inter_indIND3 <- glmer(formula3, db_2008[db_2008$rand.vec != 3,], family = binomial)
mod_inter_indIND4 <- glmer(formula3, db_2008[db_2008$rand.vec != 4,], family = binomial)
mod_inter_indIND5 <- glmer(formula3, db_2008[db_2008$rand.vec != 5,], family = binomial)
save(mod_inter_indIND1, file = "mod_inter_indIND1.RData")
save(mod_inter_indIND2, file = "mod_inter_indIND2.RData")
save(mod_inter_indIND3, file = "mod_inter_indIND3.RData")
save(mod_inter_indIND4, file = "mod_inter_indIND4.RData")
save(mod_inter_indIND5, file = "mod_inter_indIND5.RData")
# load(file = "mod_inter_indIND1.RData")
# load(file = "mod_inter_indIND2.RData")
# load(file = "mod_inter_indIND3.RData")
# load(file = "mod_inter_indIND4.RData")
# load(file = "mod_inter_indIND5.RData")

# Remove columns used for evaluations fit earlier
db_2008$lc_30_contribute <- NULL
db_2008$RSFscores <- NULL

# Calculation of the contribution to RSF scores by the levels of the categorical predictor (reference category: coniferous forest) and calculation of RSF scores assuming the exponential form. Repeating it for each withhold fold.
# NOTE THE ATTACH COMMAND HERE - DETACH BELOW
attach(db_2008)
db_2008$lc_30_contribute[rand.vec == 1] <- as.numeric(ifelse(db_2008$lc_30f[rand.vec == 1] == "Coniferous forest", 0,
                                                            ifelse(db_2008$lc_30f[rand.vec == 1] == "Deciduous forest", mod_inter_indIND1@beta[8],
                                                                   ifelse(db_2008$lc_30f[rand.vec == 1] == "Mixed forest", mod_inter_indIND1@beta[9],
                                                                          ifelse(db_2008$lc_30f[rand.vec == 1] == "Open areas", mod_inter_indIND1@beta[10], "ERROR")))))
db_2008$RSFscores[rand.vec == 1] <- exp(mod_inter_indIND1@beta[2] * scale(db_2008$LogRugg[rand.vec == 1]) +
                                         mod_inter_indIND1@beta[3] * I(scale(db_2008$LogRugg[rand.vec == 1])^2) +
                                         mod_inter_indIND1@beta[4] * scale(db_2008$SqrtDistroad[rand.vec == 1]) +
                                         mod_inter_indIND1@beta[5] * I(scale(db_2008$SqrtDistroad[rand.vec == 1])^2) +
                                         mod_inter_indIND1@beta[6] * scale(db_2008$NDVI2[rand.vec == 1]) +
                                         mod_inter_indIND1@beta[7] * I(scale(db_2008$NDVI2[rand.vec == 1])^2) +
                                         db_2008$lc_30_contribute[rand.vec == 1] +
                                         mod_inter_indIND1@beta[11] * scale(db_2008$LogRugg[rand.vec == 1]) * scale(db_2008$SqrtDistroad[rand.vec == 1]) +
                                         mod_inter_indIND1@beta[12] * scale(db_2008$SqrtDistroad[rand.vec == 1]) * scale(db_2008$NDVI2[rand.vec == 1]))

db_2008$lc_30_contribute[rand.vec == 2] <- as.numeric(ifelse(db_2008$lc_30f[rand.vec == 2] == "Coniferous forest", 0,
                                                            ifelse(db_2008$lc_30f[rand.vec == 2] == "Deciduous forest", mod_inter_indIND2@beta[8],
                                                                   ifelse(db_2008$lc_30f[rand.vec == 2] == "Mixed forest", mod_inter_indIND2@beta[9],
                                                                          ifelse(db_2008$lc_30f[rand.vec == 2] == "Open areas", mod_inter_indIND2@beta[10], "ERROR")))))
db_2008$RSFscores[rand.vec == 2] <- exp(mod_inter_indIND2@beta[2] * scale(db_2008$LogRugg[rand.vec == 2]) +
                                         mod_inter_indIND2@beta[3] * I(scale(db_2008$LogRugg[rand.vec == 2])^2) +
                                         mod_inter_indIND2@beta[4] * scale(db_2008$SqrtDistroad[rand.vec == 2]) +
                                         mod_inter_indIND2@beta[5] * I(scale(db_2008$SqrtDistroad[rand.vec == 2])^2) +
                                         mod_inter_indIND2@beta[6] * scale(db_2008$NDVI2[rand.vec == 2]) +
                                         mod_inter_indIND2@beta[7] * I(scale(db_2008$NDVI2[rand.vec == 2])^2) +
                                         db_2008$lc_30_contribute[rand.vec == 2] +
                                         mod_inter_indIND2@beta[11] * scale(db_2008$LogRugg[rand.vec == 2]) * scale(db_2008$SqrtDistroad[rand.vec == 2]) +
                                         mod_inter_indIND2@beta[12] * scale(db_2008$SqrtDistroad[rand.vec == 2]) * scale(db_2008$NDVI2[rand.vec == 2]))

db_2008$lc_30_contribute[rand.vec == 3] <- as.numeric(ifelse(db_2008$lc_30f[rand.vec == 3] == "Coniferous forest", 0,
                                                            ifelse(db_2008$lc_30f[rand.vec == 3] == "Deciduous forest", mod_inter_indIND3@beta[8],
                                                                   ifelse(db_2008$lc_30f[rand.vec == 3] == "Mixed forest", mod_inter_indIND3@beta[9],
                                                                          ifelse(db_2008$lc_30f[rand.vec == 3] == "Open areas", mod_inter_indIND3@beta[10], "ERROR")))))
db_2008$RSFscores[rand.vec == 3] <- exp(mod_inter_indIND3@beta[2] * scale(db_2008$LogRugg[rand.vec == 3]) +
                                         mod_inter_indIND3@beta[3] * I(scale(db_2008$LogRugg[rand.vec == 3])^2) +
                                         mod_inter_indIND3@beta[4] * scale(db_2008$SqrtDistroad[rand.vec == 3]) +
                                         mod_inter_indIND3@beta[5] * I(scale(db_2008$SqrtDistroad[rand.vec == 3])^2) +
                                         mod_inter_indIND3@beta[6] * scale(db_2008$NDVI2[rand.vec == 3]) +
                                         mod_inter_indIND3@beta[7] * I(scale(db_2008$NDVI2[rand.vec == 3])^2) +
                                         db_2008$lc_30_contribute[rand.vec == 3] +
                                         mod_inter_indIND3@beta[11] * scale(db_2008$LogRugg[rand.vec == 3]) * scale(db_2008$SqrtDistroad[rand.vec == 3]) +
                                         mod_inter_indIND3@beta[12] * scale(db_2008$SqrtDistroad[rand.vec == 3]) * scale(db_2008$NDVI2[rand.vec == 3]))

db_2008$lc_30_contribute[rand.vec == 4] <- as.numeric(ifelse(db_2008$lc_30f[rand.vec == 4] == "Coniferous forest", 0,
                                                            ifelse(db_2008$lc_30f[rand.vec == 4] == "Deciduous forest", mod_inter_indIND4@beta[8],
                                                                   ifelse(db_2008$lc_30f[rand.vec == 4] == "Mixed forest", mod_inter_indIND4@beta[9],
                                                                          ifelse(db_2008$lc_30f[rand.vec == 4] == "Open areas", mod_inter_indIND4@beta[10], "ERROR")))))
db_2008$RSFscores[rand.vec == 4] <- exp(mod_inter_indIND4@beta[2] * scale(db_2008$LogRugg[rand.vec == 4]) +
                                         mod_inter_indIND4@beta[3] * I(scale(db_2008$LogRugg[rand.vec == 4])^2) +
                                         mod_inter_indIND4@beta[4] * scale(db_2008$SqrtDistroad[rand.vec == 4]) +
                                         mod_inter_indIND4@beta[5] * I(scale(db_2008$SqrtDistroad[rand.vec == 4])^2) +
                                         mod_inter_indIND4@beta[6] * scale(db_2008$NDVI2[rand.vec == 4]) +
                                         mod_inter_indIND4@beta[7] * I(scale(db_2008$NDVI2[rand.vec == 4])^2) +
                                         db_2008$lc_30_contribute[rand.vec == 4] +
                                         mod_inter_indIND4@beta[11] * scale(db_2008$LogRugg[rand.vec == 4]) * scale(db_2008$SqrtDistroad[rand.vec == 4]) +
                                         mod_inter_indIND4@beta[12] * scale(db_2008$SqrtDistroad[rand.vec == 4]) * scale(db_2008$NDVI2[rand.vec == 4]))

db_2008$lc_30_contribute[rand.vec == 5] <- as.numeric(ifelse(db_2008$lc_30f[rand.vec == 5] == "Coniferous forest", 0,
                                                            ifelse(db_2008$lc_30f[rand.vec == 5] == "Deciduous forest", mod_inter_indIND5@beta[8],
                                                                   ifelse(db_2008$lc_30f[rand.vec == 5] == "Mixed forest", mod_inter_indIND5@beta[9],
                                                                          ifelse(db_2008$lc_30f[rand.vec == 5] == "Open areas", mod_inter_indIND5@beta[10], "ERROR")))))
db_2008$RSFscores[rand.vec == 5] <- exp(mod_inter_indIND5@beta[2] * scale(db_2008$LogRugg[rand.vec == 5]) +
                                         mod_inter_indIND5@beta[3] * I(scale(db_2008$LogRugg[rand.vec == 5])^2) +
                                         mod_inter_indIND5@beta[4] * scale(db_2008$SqrtDistroad[rand.vec == 5]) +
                                         mod_inter_indIND5@beta[5] * I(scale(db_2008$SqrtDistroad[rand.vec == 5])^2) +
                                         mod_inter_indIND5@beta[6] * scale(db_2008$NDVI2[rand.vec == 5]) +
                                         mod_inter_indIND5@beta[7] * I(scale(db_2008$NDVI2[rand.vec == 5])^2) +
                                         db_2008$lc_30_contribute[rand.vec == 5] +
                                         mod_inter_indIND5@beta[11] * scale(db_2008$LogRugg[rand.vec == 5]) * scale(db_2008$SqrtDistroad[rand.vec == 5]) +
                                         mod_inter_indIND5@beta[12] * scale(db_2008$SqrtDistroad[rand.vec == 5]) * scale(db_2008$NDVI2[rand.vec == 5]))
detach(db_2008)

# Run the k-fold CV evaluation sensu Boyce et al. 2002
dataset <- db_2008[complete.cases(db_2008[,"RSFscores"]),]
rho_model <- numeric(5) ## it will store Spearman's coefficients

fold <- subset(dataset,rand.vec==unique(dataset$rand.vec)[1]) # run the procedure for the 1st fold - the for-loop below will work on folds 2 to 5
q.pp <- quantile(fold$RSFscores,probs=seq(0,1,.1))  ## computing quantiles of RSF scores
# --------------------------------------------------------
bin <- rep(NA,length(fold$RSFscores)) ## binning RSF scores (10 bins)
for (j in 1:10){
  bin[fold$RSFscores>=q.pp[j]& fold$RSFscores<q.pp[j+1]] = j
}
used <- fold$use
# --------------------------------------------------------
a <- table(used,bin) ## Occurrence of presence and available data by bin
a <- t(a) #transpose the table
a <- as.data.frame.matrix(a) ## the next few lines compute area-adjusted frequency of categories (bins) of RSF scores 
a$areaadjusted <- rep(NA,length(10))
sum0 <- sum(a[,1])
sum1 <- sum(a[,2])
a$areaadjusted <- (a[,2] / sum1 ) / (a[,1] / sum0)
a$bins <- seq(1,10,by=1);a

# --------------------------------------------------------
rho_model[1] <- with(a,cor.test(bins,areaadjusted,method="spearm"))$estimate ## Spearman correlation coefficient between RSF bin ranks and area-adjusted frequencies
# --------------------------------------------------------


# Run the procedure for the other folds and plot the binned RSF scores from the cross-validation
png("Results/Binned RSF scores from the spatially independent CV.png", height=360, width=600)
par(oma=c(1,2,1,1)) 
par(mar=c(4.2,4.2,2,2))
with(a,plot(bins,areaadjusted, ylab="Area adjusted frequency", xlab="Binned RSF scores",xlim=c(1,10),ylim=c(0,3),type="b",cex=1.4,cex.lab=1.5,las=1,main="K-fold validation"))
abline(h = 1, lty = 3) ## Spearman correlation coefficient between RSF bin ranks and area-adjusted frequencies
plot_ind_ind <- data.frame(a[,3], b = NA, c = NA, d = NA, e = NA)
for (i in 2:5){
  fold <- subset(dataset,rand.vec==unique(dataset$rand.vec)[i]) ## run the procedure on folds 2 to 5
  # --------------------------------------------------------
  q.pp <- quantile(fold$RSFscores,probs=seq(0,1,.1)) ## computing quantiles of RSF scores
  # --------------------------------------------------------
  bin <- rep(NA,length(fold$RSFscores))
  for (j in 1:10){
    bin[fold$RSFscores>=q.pp[j]& fold$RSFscores<q.pp[j+1]] = j ## binning RSF scores (10 bins)
  }
  used = fold$use
  # --------------------------------------------------------
  a <- table(used,bin) ## area adjusted freq in used/available for each bin
  a <- t(a) #transpose the table
  a <- as.data.frame.matrix(a) ## the next few lines compute area-adjusted frequency of categories (bins) of RSF scores
  a$areaadjusted <- rep(NA,length(10))
  sum0 <- sum(a[,1])
  sum1 <- sum(a[,2])
  a$areaadjusted <- (a[,2] / sum1 ) / (a[,1] / sum0)
  a$bins <- seq(1,10,by=1);a
  # --------------------------------------------------------
  rho_model[i] <- with(a,cor.test(bins,areaadjusted,method="spearm"))$estimate  ## store Spearman correlation coefficients between RSF bin ranks and area-adjusted frequencies
  # --------------------------------------------------------
  par(oma=c(1,2,1,1)) 
  par(mar=c(4.2,4.2,2,2))
  with(a,points(bins,areaadjusted,col=i,type="b"))  
  plot_ind_ind[,i] = a[,3]
}
dev.off()

## store Spearman correlation coefficients that will be used for final plots below ##
Rho_independent_individuals <- rho_model

# Make combined table of all Spearman coefficients (RHOs) previously stored after each cross-validation (it will be used below for final visualization)
RHO_coeff <- data.frame(Rho_full = c(1,NA,NA,NA,NA), Rho_random_fixes = Rho_random_fixes, Rho_random_individuals = Rho_random_individuals, Rho_independent_individuals = Rho_independent_individuals)



## PART 9 - Details on how individuals were allocated to folds: cluster analysis fit and visualization.
## ----------------------------------------------------------------------------------------------------

# Load the file centres2008 (geometric centres of elk Home Ranges)
load("Appendix_6_Box_2_DATA_centres2008.RData")

# Some data visualization
png("Results/Individual home range centres.png")
plot(data$Easting, data$Northing, xlab = "Easting", ylab = "Northing")
text(data$Easting, data$Northing, labels = data$ElkID, cex = 0.6, pos = 3)
dev.off()

# Compute distance matrix and fit cluster analysis
data$newcol <- as.factor(do.call(paste, c(data[c(1,3)], sep = " ")))
HRcentres <- as.matrix(data[,c("Easting", "Northing")]) 
rownames(HRcentres) <- data[,"newcol"]
dist_matrix <- dist(HRcentres, method = "euclidean")
fit <- hclust(dist_matrix, method = "complete")

# Plot cluster tree with spatially independent individuals (20 km threshold)
png("Results/Dendrogram of spatially independent individuals.png", height=800, width=800)
par(mar=c(9.5,5,2,1))
plot(fit) # display dendogram
abline(h = 20000, lwd = 3, col = "green", lty = 3)
### grouping to achieve spatial independency
rect(35.5, -60000, 43.5, -55000, col = "red", border = "red"); text(39, -57600, "fold 5")
rect(23.5, -60000, 35.5, -55000, col = "blue", border = "blue") ; text(29, -57600, "fold 4", col = "white")
rect(15.5, -60000, 23.5, -55000, col = "green", border = "green"); text(19.5, -57600, "fold 3")
rect(8.5, -60000, 15.5, -55000, col = "orange", border = "orange") ; text(12.5, -57600, "fold 2")
rect(0.5, -60000, 8.5, -55000, col = "brown", border = "brown"); text(4.5, -57600, "fold 1")
dev.off()

#plot cluster tree with visualization of randomly selected individuals 
data$newcol <- as.factor(do.call(paste, c(data[c(1,3)], sep = " ")))
HRcentres <- as.matrix(data[,c("Easting", "Northing")]) 
rownames(HRcentres) <- data[,"newcol"]
dist_matrix <- dist(HRcentres, method = "euclidean")
fit <- hclust(dist_matrix, method = "complete")

# Plot the binned RSF scores from the random cross-validation
png("Results/Dendrogram of random individuals.png", height=800, width=800)
par(mar=c(7,5,2,1))
plot(fit) # display dendogram
text(0, -55000, "fold:", font = 2, cex = 0.7)
# Fold 1
for(i in c(1,3,14,16,25,28,31,37,40))
{text(i, -55000, "1", col = "brown", font = 4)}
# Fold 2
for(i in c(2,6,11,12,13,17,32,41))
{text(i, -55000, "2", col = "orange", font = 2)}
# Fold 3
for(i in c(5,7,15,18,20,21,24,27,39))
{text(i, -55000, "3", col = "green", font = 3)}
# Fold 4
for(i in c(4,9,19,22,29,33,36,43))
{text(i, -55000, "4", col = "blue")}
# Fold 5
for(i in c(8,10,23,26,30,34,35,38,42))
{text(i, -55000, "5", col = "red")}
dev.off()



## PART 10 - Visualization of main results
## ---------------------------------------

# Plot all results together
png("Results/Results plot (Rhos, betas, scores).png", height=750, width=1050)

# Spearman rank correlations by CV method
layout(matrix(c(1,1,2,2,1,1,2,2,3,4,5,6), 3, 4, byrow = TRUE))
attach(RHO_coeff)
plot(1:5, Rho_full, las=1, xaxt="n", col="black", bg="black", pch=21, cex=2,  ylim=c(0.75, 1), xlim=(c(0.5,4.5)),
     ylab="Spearman-rank correlations (rs)", xlab="Blocking design", cex.lab=1.2, bty="n")
x1 = c(1:4)
axis(1, at=x1, labels=c("full model","CV","BCV (ri)","BCV (sbi)"), cex.axis=0.9)
x = c(1.9, 1.95, 2, 2.05, 2.1)
points(x,   Rho_random_fixes,            col="black", bg="gray60", pch=21, cex = 2)
points(x+1, Rho_random_individuals,      col="black", bg="gray90", pch=21, cex = 2)
points(x+2, Rho_independent_individuals, col="black", bg="white",  pch=21, cex = 2)
detach(RHO_coeff)
text(4, 0.75, "Spearman's Rho", col="red", cex=1.5, font=2, pos=2)
legend(1, 0.8, pch=c(21,21,21,21), 
       pt.bg=c("black","gray60","gray90","white"),
       legend=c("Resubstitution","Random fixes","Randomly blocked individuals","Spatially blocked individuals"),
       col=rep("black",4), cex=1.1, pt.cex=2, bty="n")

# Beta estimates by CV method
x = 1:13
y = seq(1:13)
plot(x,y, pch="", bty="n", las=1, xaxt="n", ylim = c(-0.25, 0.9),
     ylab = "Beta estimates (GLMMs)", xlab="", cex.lab=1.2)
axis(1, at=c(2:12), cex.axis=0.9, las=2,
     labels=c("tri","tri^2","d","d^2","NDVI","NDVI^2","Decid","Mixed","Open","tri * d","d * NDVI"))
x=2:12; i=1.8; a=(-0.50); b=(-0.25); c=0; d=0.25
points(x + a, model3@beta[2:12],               col="black", bg="black",  pch=21, cex=i)
points(x + b, mod_inter_randomf1@beta[2:12],   col="black", bg="gray60", pch=21, cex=i)
points(x + b, mod_inter_randomf2@beta[2:12],   col="black", bg="gray60", pch=21, cex=i)
points(x + b, mod_inter_randomf3@beta[2:12],   col="black", bg="gray60", pch=21, cex=i)
points(x + b, mod_inter_randomf4@beta[2:12],   col="black", bg="gray60", pch=21, cex=i)
points(x + b, mod_inter_randomf5@beta[2:12],   col="black", bg="gray60", pch=21, cex=i)
points(x + c, mod_inter_randomIND1@beta[2:12], col="black", bg="gray90", pch=21, cex=i)
points(x + c, mod_inter_randomIND2@beta[2:12], col="black", bg="gray90", pch=21, cex=i)
points(x + c, mod_inter_randomIND3@beta[2:12], col="black", bg="gray90", pch=21, cex=i)
points(x + c, mod_inter_randomIND4@beta[2:12], col="black", bg="gray90", pch=21, cex=i)
points(x + c, mod_inter_randomIND5@beta[2:12], col="black", bg="gray90", pch=21, cex=i)
points(x + d, mod_inter_indIND1@beta[2:12],    col="black", bg="white",  pch=21, cex=i)
points(x + d, mod_inter_indIND2@beta[2:12],    col="black", bg="white",  pch=21, cex=i)
points(x + d, mod_inter_indIND3@beta[2:12],    col="black", bg="white",  pch=21, cex=i)
points(x + d, mod_inter_indIND4@beta[2:12],    col="black", bg="white",  pch=21, cex=i)
points(x + d, mod_inter_indIND5@beta[2:12],    col="black", bg="white",  pch=21, cex=i)
text(12, -0.25, "Beta estimates", col="red", cex=1.5, font=2, pos=2)

# Four plots of Binned RSF scores for each CV approach

y.limits <- c(0,max(max(plot_full),max(plot_random),max(plot_random_ind),max(plot_ind_ind)))

# Resubstitution (no blocking)
plot(1:10, plot_full, las=1, bty="n", col="black", bg="black", pch=21, cex=1.5, type="b", ylim=y.limits,
     ylab="Area adjusted frequency", xlab="Binned RSF scores", cex.lab=1.2)
abline(h=1, lty=3, lwd=0.1, col="gray")
text(10, 0.2, "Resubstitution", col="red", cex=1, font=2, pos=2)

# Random
plot(1:10, plot_random[,1], las=1, bty="n", col="black", bg="gray60", pch=21, cex=1.5, type="b", ylim=y.limits,
     ylab="Area adjusted frequency", xlab="Binned RSF scores", cex.lab=1.2)
points(1:10, plot_random[,2], bty="n", col="black", bg="gray60", pch=21, cex=1.5, type="b")
points(1:10, plot_random[,3], bty="n", col="black", bg="gray60", pch=21, cex=1.5, type="b")
points(1:10, plot_random[,4], bty="n", col="black", bg="gray60", pch=21, cex=1.5, type="b")
points(1:10, plot_random[,5], bty="n", col="black", bg="gray60", pch=21, cex=1.5, type="b")
abline(h = 1, lty = 3, lwd = 0.1, col = "gray")
text(10, 0.2, "Random fixes", col="red", cex=1, font=2, pos=2)

# Random individuals
plot(1:10, plot_random_ind[,1], las = 1, bty="n", col="black", bg="gray90", pch=21, cex=1.5, type="b", ylim=y.limits, 
     ylab="Area adjusted frequency", xlab="Binned RSF scores", cex.lab=1.2)
points(1:10, plot_random_ind[,2], bty="n", col="black", bg="gray90", pch=21, cex=1.5, type="b")
points(1:10, plot_random_ind[,3], bty="n", col="black", bg="gray90", pch=21, cex=1.5, type="b")
points(1:10, plot_random_ind[,4], bty="n", col="black", bg="gray90", pch=21, cex=1.5, type="b")
points(1:10, plot_random_ind[,5], bty="n", col="black", bg="gray90", pch=21, cex=1.5, type="b")
abline(h = 1, lty = 3, lwd = 0.1, col = "gray")
text(10, 0.2, "Randomly blocked individuals", col="red", cex=1, font=2, pos=2)

# Spatially independent individuals
plot(1:10, plot_ind_ind[,1], las=1, bty="n", col="black", bg="white", pch=21, cex=1.5, type="b", ylim=y.limits, 
     ylab="Area adjusted frequency", xlab="Binned RSF scores", cex.lab=1.2)
points(1:10, plot_ind_ind[,2], bty="n", col="black", bg="white", pch=21, cex=1.5, type="b")
points(1:10, plot_ind_ind[,3], bty="n", col="black", bg="white", pch=21, cex=1.5, type="b")
points(1:10, plot_ind_ind[,4], bty="n", col="black", bg="white", pch=21, cex=1.5, type="b")
points(1:10, plot_ind_ind[,5], bty="n", col="black", bg="white", pch=21, cex=1.5, type="b")
abline(h = 1, lty = 3, lwd = 0.1, col = "gray")
text(10, 0.2, "Spatially blocked individuals", col="red", cex=1, font=2, pos=2)

dev.off()



