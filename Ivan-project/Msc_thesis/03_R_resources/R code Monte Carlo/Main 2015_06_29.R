setwd("C:\\Users\\andre\\Documents\\New folder")
source("SimulationFunctions.r")
library(fun)

aggList = c(1, 2, 4, 5, 6, 10, 12)


################################################################################
### SPECIAL RUNS

nRep = 2

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.9"
OLSEst(c(1, 1, 0.5, 0.9), randomTreat = 2, aggList, disSD = 1, nRep, fullRun = FALSE, disSteps = 2, chartTitle, xType = "wwx")


chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.9"
STSLSEst(c(1, 1, 0.5, 0.9), 2, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.9"
STSLSEst(c(1, 1, 0.5, 0.9), randomTreat = 2, aggList, disSD = 1, nRep, fullRun = FALSE, disSteps = 6, chartTitle, xType = "wwx")

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.9"
OLSEst(c(1, 1, 0.5, 0.9), randomTreat = 3, aggList, disSD = 1, nRep, fullRun = FALSE, disSteps = 2, chartTitle, xType = "wwx")

chartTitle = "Y random, X random, T random"
OLSEst(c(1, 1, 0.5, 0), randomTreat = 0, aggList, disSD = 1, nRep, fullRun = TRUE, chartTitle, xType = "random")

chartTitle = "Y random, X random, T random"
OLSEst(c(1, 1, 0.5, 0), randomTreat = 0, aggList, disSD = 1, nRep, fullRun = FALSE, disSteps = 5, chartTitle, xType = "random", useBinaryAgg = FALSE)


# Scenario 5: T: X dependent | X: correlated (X=W^2 K) | Y: Random, corr = 0.3
chartTitle = "Spatial Autocorrelation X, T, Y"
OLSEst(trueCoeff = c(1, 1, 0.5, 0), randomTreat = 2, aggList, disSD = 1, nRep, fullRun = TRUE, chartTitle, xType = "wwx", useBinaryAgg = TRUE, corrXT = 0.3)

chartTitle = "Spatial Autocorrelation X, T, Y"
OLSEst(trueCoeff = c(1, 1, 0.5, 0), randomTreat = 2, aggList, disSD = 1, nRep, fullRun = FALSE, disSteps = c(10,12), chartTitle = chartTitle, xType = "wwx", useBinaryAgg = TRUE, corrXT = 0.3)


# Simple Matching
chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.9"
NonSpatialMatching(trueCoeff = c(1, 1, 0.5, 0), randomTreat = 2, aggList, disSD = 1, nRep, fullRun = TRUE, chartTitle, xType = "wwx", useBinaryAgg = TRUE, corrXT = 0.5)


### Study of Bump in SE

chartTitle = "Traditional Literature Case"
SEMEst(trueCoeff = c(1, 1, 0.5, 0), randomTreat = 1, aggList, disSD = 1, nRep, fullRun = FALSE, disSteps = 6, chartTitle, xType = "random", useBinaryAgg = TRUE, corrXT = 0)

SEMEst2(trueCoeff = c(1, 1, 0.5, 0), randomTreat = 1, aggList, disSD = 1, nRep, fullRun = FALSE, disSteps = 6, chartTitle, xType = "random", useBinaryAgg = TRUE, corrXT = 0)


OLSEst(trueCoeff = c(1, 1, 0.5, 0), randomTreat = 1, aggList, disSD = 1, nRep, fullRun = FALSE, disSteps = 6, chartTitle, xType = "random", useBinaryAgg = TRUE, corrXT = 0)


################################################################################
### OLS

nRep = 2

chartTitle = "Traditional Literature Case"
OLSEst(c(1, 1, 0.5, 0), 1, aggList, 1, nRep, chartTitle, fullRun = FALSE, disSteps = c(1, 2, 4, 5, 6, 10, 12), xType = "random")

chartTitle = "Spatial Lag Y 0.5"
OLSEst(c(1, 1, 0.5, 0.5), 0, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Lag Y 0.9"
OLSEst(c(1, 1, 0.5, 0.9), 0, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
OLSEst(c(1, 1, 0.5, 0), 1, aggList, 1, nRep, chartTitle, fullRun = FALSE, disSteps = 6, xType = "lagged")

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
OLSEst(c(1, 1, 0.5, 0.5), 1, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.9"
OLSEst(c(1, 1, 0.5, 0.9), 1, aggList, 1, nRep, chartTitle, xType = "wwx")

chartTitle = "Spatial Autocorrelation X, T, Y"
OLSEst(c(1, 1, 0.5, 0), 2, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
OLSEst(c(1, 1, 0.5, 0.5), 2, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.9"
OLSEst(c(1, 1, 0.5, 0.9), 2, aggList, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
OLSEst(c(1, 1, 0.5, 0), 3, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
OLSEst(c(1, 1, 0.5, 0.5), 3, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.9"
OLSEst(c(1, 1, 0.5, 0.9), 3, aggList, 1, nRep, chartTitle)


### Study of Spatial Heterogeneity

chartTitle = "Spatial Lag Y 0.5"
OLSEst(c(1, 1, 0.5, 0.5), 0, aggList, 0, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
OLSEst(c(1, 1, 0.5, 0), 1, aggList, 0, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
OLSEst(c(1, 1, 0.5, 0.5), 1, aggList, 0, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
OLSEst(c(1, 1, 0.5, 0), 2, aggList, 0, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
OLSEst(c(1, 1, 0.5, 0.5), 2, aggList, 0, nRep, chartTitle)

chartTitle = "Spatial Lag Y 0.5"
OLSEst(c(1, 1, 0.5, 0.5), 0, aggList, 5, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
OLSEst(c(1, 1, 0.5, 0), 1, aggList, 5, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5" #HERE
OLSEst(c(1, 1, 0.5, 0.5), 1, aggList, 5, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
OLSEst(c(1, 1, 0.5, 0), 2, aggList, 5, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
OLSEst(c(1, 1, 0.5, 0.5), 2, aggList, 5, nRep, chartTitle)




################################################################################
### STSLS

nRep = 1000

chartTitle = "No Spatial Process"
STSLSEst(c(1, 1, 0.5, 0), 0, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Lag Y 0.5"
STSLSEst(c(1, 1, 0.5, 0.5), 0, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Lag Y 0.9"
STSLSEst(c(1, 1, 0.5, 0.9), 0, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
STSLSEst(c(1, 1, 0.5, 0), 1, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
STSLSEst(c(1, 1, 0.5, 0.5), 1, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.9"
STSLSEst(c(1, 1, 0.5, 0.9), 1, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
STSLSEst(c(1, 1, 0.5, 0), 2, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
STSLSEst(c(1, 1, 0.5, 0.5), 2, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.9"
STSLSEst(c(1, 1, 0.5, 0.9), 2, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
STSLSEst(c(1, 1, 0.5, 0), 3, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
STSLSEst(c(1, 1, 0.5, 0.5), 3, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.9"
STSLSEst(c(1, 1, 0.5, 0.9), 3, aggList, 1, nRep, chartTitle)


### Study of Spatial Heterogeneity

chartTitle = "Spatial Lag Y 0.5"
STSLSEst(c(1, 1, 0.5, 0.5), 0, aggList, 0, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
STSLSEst(c(1, 1, 0.5, 0), 1, aggList, 0, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
STSLSEst(c(1, 1, 0.5, 0.5), 1, aggList, 0, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
STSLSEst(c(1, 1, 0.5, 0), 2, aggList, 0, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
STSLSEst(c(1, 1, 0.5, 0.5), 2, aggList, 0, nRep, chartTitle)

chartTitle = "Spatial Lag Y 0.5"
STSLSEst(c(1, 1, 0.5, 0.5), 0, aggList, 5, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
STSLSEst(c(1, 1, 0.5, 0), 1, aggList, 5, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
STSLSEst(c(1, 1, 0.5, 0.5), 1, aggList, 5, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
STSLSEst(c(1, 1, 0.5, 0), 2, aggList, 5, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
STSLSEst(c(1, 1, 0.5, 0.5), 2, aggList, 5, nRep, chartTitle)




################################################################################
### Non-Spatial DID

nRep = 1000

chartTitle = "No Spatial Process"
NonSpatialDID(c(1, 1, 0.5, 0), 0, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Lag Y 0.5"
NonSpatialDID(c(1, 1, 0.5, 0.5), 0, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Lag Y 0.9"
NonSpatialDID(c(1, 1, 0.5, 0.9), 0, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
NonSpatialDID(c(1, 1, 0.5, 0), 1, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
NonSpatialDID(c(1, 1, 0.5, 0.5), 1, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.9"
NonSpatialDID(c(1, 1, 0.5, 0.9), 1, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
NonSpatialDID(c(1, 1, 0.5, 0), 2, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
NonSpatialDID(c(1, 1, 0.5, 0.5), 2, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.9"
NonSpatialDID(c(1, 1, 0.5, 0.9), 2, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
NonSpatialDID(c(1, 1, 0.5, 0), 3, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
NonSpatialDID(c(1, 1, 0.5, 0.5), 3, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.9"
NonSpatialDID(c(1, 1, 0.5, 0.9), 3, aggList, 1, nRep, chartTitle)
shutdown()

### Study of Spatial Heterogeneity

chartTitle = "Spatial Lag Y 0.5"
NonSpatialDID(c(1, 1, 0.5, 0.5), 0, aggList, 0, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
NonSpatialDID(c(1, 1, 0.5, 0), 1, aggList, 0, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
NonSpatialDID(c(1, 1, 0.5, 0.5), 1, aggList, 0, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
NonSpatialDID(c(1, 1, 0.5, 0), 2, aggList, 0, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
NonSpatialDID(c(1, 1, 0.5, 0.5), 2, aggList, 0, nRep, chartTitle)

chartTitle = "Spatial Lag Y 0.5"
NonSpatialDID(c(1, 1, 0.5, 0.5), 0, aggList, 5, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
NonSpatialDID(c(1, 1, 0.5, 0), 1, aggList, 5, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
NonSpatialDID(c(1, 1, 0.5, 0.5), 1, aggList, 5, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
NonSpatialDID(c(1, 1, 0.5, 0), 2, aggList, 5, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
NonSpatialDID(c(1, 1, 0.5, 0.5), 2, aggList, 5, nRep, chartTitle)




### Spatial DID

nRep = 50

chartTitle = "No Spatial Process"
SpatialDID(c(1, 1, 0.5, 0), 0, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Lag Y 0.5"
SpatialDID(c(1, 1, 0.5, 0.5), 0, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Lag Y 0.9"
SpatialDID(c(1, 1, 0.5, 0.9), 0, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
SpatialDID(c(1, 1, 0.5, 0), 1, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
SpatialDID(c(1, 1, 0.5, 0.5), 1, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.9"
SpatialDID(c(1, 1, 0.5, 0.9), 1, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
SpatialDID(c(1, 1, 0.5, 0), 2, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
SpatialDID(c(1, 1, 0.5, 0.5), 2, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.9"
SpatialDID(c(1, 1, 0.5, 0.9), 2, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
SpatialDID(c(1, 1, 0.5, 0), 3, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
SpatialDID(c(1, 1, 0.5, 0.5), 3, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.9"
SpatialDID(c(1, 1, 0.5, 0.9), 3, aggList, 1, nRep, chartTitle)


### Study of Spatial Heterogeneity

chartTitle = "Spatial Lag Y 0.5"
SpatialDID(c(1, 1, 0.5, 0.5), 0, aggList, 0, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
SpatialDID(c(1, 1, 0.5, 0), 1, aggList, 0, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
SpatialDID(c(1, 1, 0.5, 0.5), 1, aggList, 0, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
SpatialDID(c(1, 1, 0.5, 0), 2, aggList, 0, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
SpatialDID(c(1, 1, 0.5, 0.5), 2, aggList, 0, nRep, chartTitle)

chartTitle = "Spatial Lag Y 0.5"
SpatialDID(c(1, 1, 0.5, 0.5), 0, aggList, 5, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
SpatialDID(c(1, 1, 0.5, 0), 1, aggList, 5, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
SpatialDID(c(1, 1, 0.5, 0.5), 1, aggList, 5, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
SpatialDID(c(1, 1, 0.5, 0), 2, aggList, 5, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
SpatialDID(c(1, 1, 0.5, 0.5), 2, aggList, 5, nRep, chartTitle)




################################################################################
### Non-Spatial Matching

nRep = 1000

chartTitle = "No Spatial Process"
NonSpatialMatching(c(1, 1, 0.5, 0), 0, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Lag Y 0.5"
NonSpatialMatching(c(1, 1, 0.5, 0.5), 0, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Lag Y 0.9"
NonSpatialMatching(c(1, 1, 0.5, 0.9), 0, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
NonSpatialMatching(c(1, 1, 0.5, 0), 1, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
NonSpatialMatching(c(1, 1, 0.5, 0.5), 1, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.9"
NonSpatialMatching(c(1, 1, 0.5, 0.9), 1, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
NonSpatialMatching(c(1, 1, 0.5, 0), 2, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5" #MISSING
NonSpatialMatching(c(1, 1, 0.5, 0.5), 2, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.9"
NonSpatialMatching(c(1, 1, 0.5, 0.9), 2, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y" #HERE
NonSpatialMatching(c(1, 1, 0.5, 0), 3, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
NonSpatialMatching(c(1, 1, 0.5, 0.5), 3, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.9"
NonSpatialMatching(c(1, 1, 0.5, 0.9), 3, aggList, 1, nRep, chartTitle)
shutdown()

### Study of Spatial Heterogeneity

chartTitle = "Spatial Lag Y 0.5"
NonSpatialMatching(c(1, 1, 0.5, 0.5), 0, aggList, 0, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
NonSpatialMatching(c(1, 1, 0.5, 0), 1, aggList, 0, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
NonSpatialMatching(c(1, 1, 0.5, 0.5), 1, aggList, 0, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
NonSpatialMatching(c(1, 1, 0.5, 0), 2, aggList, 0, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
NonSpatialMatching(c(1, 1, 0.5, 0.5), 2, aggList, 0, nRep, chartTitle)

chartTitle = "Spatial Lag Y 0.5"
NonSpatialMatching(c(1, 1, 0.5, 0.5), 0, aggList, 5, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
NonSpatialMatching(c(1, 1, 0.5, 0), 1, aggList, 5, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
NonSpatialMatching(c(1, 1, 0.5, 0.5), 1, aggList, 5, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
NonSpatialMatching(c(1, 1, 0.5, 0), 2, aggList, 5, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
NonSpatialMatching(c(1, 1, 0.5, 0.5), 2, aggList, 5, nRep, chartTitle)




################################################################################
### Spatial Matching

nRep = 10

chartTitle = "No Spatial Process"
SpatialMatching(c(1, 1, 0.5, 0), 0, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Lag Y 0.5"
SpatialMatching(c(1, 1, 0.5, 0.5), 0, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Lag Y 0.9"
SpatialMatching(c(1, 1, 0.5, 0.9), 0, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
SpatialMatching(c(1, 1, 0.5, 0), 1, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
SpatialMatching(c(1, 1, 0.5, 0.5), 1, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.9"
SpatialMatching(c(1, 1, 0.5, 0.9), 1, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
SpatialMatching(c(1, 1, 0.5, 0), 2, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
SpatialMatching(c(1, 1, 0.5, 0.5), 2, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.9"
SpatialMatching(c(1, 1, 0.5, 0.9), 2, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
SpatialMatching(c(1, 1, 0.5, 0), 3, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
SpatialMatching(c(1, 1, 0.5, 0.5), 3, aggList, 1, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.9"
SpatialMatching(c(1, 1, 0.5, 0.9), 3, aggList, 1, nRep, chartTitle)


### Study of Spatial Heterogeneity

chartTitle = "Spatial Lag Y 0.5"
SpatialMatching(c(1, 1, 0.5, 0.5), 0, aggList, 0, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
SpatialMatching(c(1, 1, 0.5, 0), 1, aggList, 0, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
SpatialMatching(c(1, 1, 0.5, 0.5), 1, aggList, 0, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
SpatialMatching(c(1, 1, 0.5, 0), 2, aggList, 0, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
SpatialMatching(c(1, 1, 0.5, 0.5), 2, aggList, 0, nRep, chartTitle)

chartTitle = "Spatial Lag Y 0.5"
SpatialMatching(c(1, 1, 0.5, 0.5), 0, aggList, 5, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
SpatialMatching(c(1, 1, 0.5, 0), 1, aggList, 5, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
SpatialMatching(c(1, 1, 0.5, 0.5), 1, aggList, 5, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y"
SpatialMatching(c(1, 1, 0.5, 0), 2, aggList, 5, nRep, chartTitle)

chartTitle = "Spatial Autocorrelation X, T, Y + Lag 0.5"
SpatialMatching(c(1, 1, 0.5, 0.5), 2, aggList, 5, nRep, chartTitle)
