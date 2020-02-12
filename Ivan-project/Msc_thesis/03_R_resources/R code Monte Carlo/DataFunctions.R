################################################################################
###                          MAUP Package                                    ###
###                                                                          ###
###                Andre F. T. Avelino (ACE/UIUC) 2014                       ###
################################################################################

library(Matrix)
library(spdep)


# <summary>
# Square-aggregates variables from a spatial grid. Requires mod(matSize/aggLevel) = 0.
# </summary>
# <param name="X"> Original variable </param>
# <param name="aggLevel"> Desired aggregation level (for 2x2 type 2, 3x3 type 3, etc.) </param>
# <param name="binaryAgg"> Binary transformation for treatment, if true returns binary </param>
# <returns> Aggregated variable </returns>
Aggregate <- function(X, aggLevel, binaryAgg = FALSE){

     if (aggLevel <= 1) return(X)
     
     matSize = length(X)^0.5
     m = (ceiling(matSize/aggLevel)-1)
     Xnew = rep(0, ceiling(matSize/aggLevel)^2)

     for (j in 0:m){                        # Row Position (aggregated matrix)
         for (i in 1:(m+1)){                # Column Position (aggregated matrix)
             for(dj in 0:(aggLevel-1)){     # Row Position (original matrix)
                    for(di in 1:aggLevel){  # Column Position (original matrix)
                           Xnew[(j*(m+1)) + i] = Xnew[(j*(m+1)) + i] + X[((j*aggLevel)+dj)*matSize + ((aggLevel*(i-1))+di)]
                    }
             }         
             Xnew[(j*(m+1)) + i] = Xnew[(j*(m+1)) + i] / (aggLevel^2)
             if(binaryAgg){if(Xnew[(j*(m+1)) + i] > 0.5){Xnew[(j*(m+1)) + i] = 1} else {Xnew[(j*(m+1)) + i] = 0}}
         }
     }

     return (Xnew)
     
}


# <summary>
# Square-disaggregates variables from a spatial grid using a normal distribution.
# </summary>
# <param name="X"> Original variable </param>
# <param name="disLevel"> Desired disaggregation level (for 2x2 type 2, 3x3 type 3, etc.) </param>
# <param name="stdDev"> Standard deviation </param>
# <returns> Aggregated variable </returns>
Disaggregate <- function(X, disLevel, stdDev){
     
     if (disLevel <= 1) return(X)
     
     newDim = (length(X)^0.5)*disLevel
     Xnew = rep(0, length(X)*disLevel^2)
     j = 0
     i = 1
      
     for (k in 1:length(X)){
         Xaux = rnorm(disLevel^2, mean = X[k], sd = stdDev)
         n = 1
         for(dj in 0:(disLevel-1)){
                for(di in 1:disLevel){
                       Xnew[((j*disLevel)+dj)*newDim + ((disLevel*(i-1))+di)] = Xaux[n]
                       n = n + 1
                }
         }
         if ((k - (j*(length(X)^0.5)) + 1) / (length(X)^0.5) > 1){
              i = 1
              j = j + 1
         }else{
              i = i + 1
         }
     }
     
     return (Xnew)

}


# <summary>
# Creates the treatment dummy for a polygon of cells (Xij to Xlm).
# </summary>
# <param name="Xini"> Column position of the upper left corner of the polygon </param>
# <param name="Yini"> Row position of the upper left corner of the polygon </param>
# <param name="Xend"> Column position of the lower right corner of the polygon </param>
# <param name="Yend"> Row position of the lower right corner of the polygon </param>
# <param name="matSize"> Row dimension of the spatial grid </param>
# <returns> Treatment variable </returns>
GenerateTreatment <- function(Xini, Yini, Xend, Yend, matSize){

        Xnew = rep(0, matSize^2)
        
        for (j in Xini:Xend){
            for (i in Yini:Yend){
                Xnew[(j-1)*matSize + i] = 1
            }
        }

        return (Xnew)

}


# <summary>
# Randomly creates the treatment dummy for a given number of cells.
# </summary>
# <param name="nTreatment"> Number of treatment units </param>
# <param name="matSize"> Row dimension of the spatial grid </param>
# <returns> Treatment variable </returns>
GenerateRandomTreatment <- function(nTreatment, matSize){
  
  Xnew = rep(0, matSize^2)
  
  selectedCells = sample(seq(1, by = 1, len = (matSize^2)), nTreatment, replace = FALSE)
  
  for (i in 1:length(selectedCells)){
    Xnew[selectedCells[i]] = 1
  }
  
  return (Xnew)
  
}


# <summary>
# Creates the treatment dummy for a polygon of cells (Xij to Xlm).
# </summary>
# <param name="Xini"> Column position of the upper left corner of the polygon </param>
# <param name="Yini"> Row position of the upper left corner of the polygon </param>
# <param name="Xend"> Column position of the lower right corner of the polygon </param>
# <param name="Yend"> Row position of the lower right corner of the polygon </param>
# <param name="matSize"> Row dimension of the spatial grid </param>
# <returns> Treatment variable </returns>
GenerateInteractiveTreatment <- function(Xint, nTreat = 0, customCorr = 0, mySeed = 1){
  
  if (nTreat == 0){nTreat = length(Xint)/2}
  if (nTreat > length(Xint)){nTreat = length(Xint)}
  
  dT = rep(0, length(Xint))
  
  data = cbind(seq(1, length(Xint), 1), Xint)
  
  if (customCorr == 0){
    
    sortedData = data[order(data[,2], decreasing = TRUE),]
    for (j in 1:nTreat){dT[sortedData[j, 1]] = 1}
  
  }else{
    
    set.seed(mySeed)
    
    dT = rep(0, length(Xint))
    sortedData = data[order(data[,2], decreasing = TRUE),]
    for (j in 1:nTreat){dT[sortedData[j, 1]] = 1}
    removeList = sample(1:nTreat, nTreat, replace = F)
    insertList = sample(nTreat:(nTreat*2), nTreat, replace = F)
    
    numSwaps = 0
    
    repeat{
      numSwaps = numSwaps + 1
      dT[sortedData[removeList[numSwaps], 1]] = 0
      dT[sortedData[insertList[numSwaps], 1]] = 1
      if((cor(dT, Xint) <= (customCorr + 0.02)) | (numSwaps == nTreat)){break}
      
    }
    
  }
  
  return (dT)
  
}


# <summary>
# Generates accummulation clusters for a given variable.
# </summary>
# <param name="X"> Variable to create hot spots </param>
# <param name="nHotSpots"> Number of hot spots </param>
# <param name="spread"> Number of lags in the accumulation points (default = 2) </param>
# <param name="mType"> Spatial weight matrix type: "rook" or "queen" </param>
# <param name="minUnif"> Minimum of the uniform distribution interval </param>
# <param name="maxUnif"> Maximum of the uniform distribution interval </param>
# <returns> Modified variable </returns>
GenerateHotSpots <- function(X, nHotSpots, spread = 2, mType = "rook", minUnif = 22, maxUnif = 30){

    xSize = length(X)^.5

    sHotSpots = sample(seq(1, by=1, len = (xSize^2)), nHotSpots, replace = FALSE)
    
    listHS = rep(0, xSize^2)
    for (i in 1:length(sHotSpots)){
        listHS[sHotSpots[i]] = 1
    }
    
    Wl = nb2listw(cell2nb(xSize, xSize, type = mType))
    HS = lag.listw(Wl, lag.listw(Wl, listHS))
    
    for (i in 1:spread){
        HS = lag.listw(Wl, HS)
    }
    
    for (i in 1:length(X)){
        X[i] = ifelse(HS[i] > 0, runif(1, min = minUnif, max = maxUnif), X[i])
    }
    
    image(t(matrix(X, size, size)))
    
    return(X)
       
}


# <summary>
# Aggregation Variance: calculates the variance of the aggregated mean
# </summary>
# <param name="X"> Original variable </param>
# <param name="aggLevel"> Desired aggregation level (for 2x2 type 2, 3x3 type 3, etc.) </param>
# <param name="average"> Mean average or matrix </param>
# <returns> Aggregation Variance </returns>
AggVariance <- function(X, aggLevel, average = TRUE){
  

  if (aggLevel <= 1) {if (average) return(0) else return(rep(0, ceiling(matSize/aggLevel)^2))}
  
  matSize = length(X)^0.5
  m = (ceiling(matSize/aggLevel)-1)
  Xnew = rep(0, ceiling(matSize/aggLevel)^2)
  Xvar = matrix(0, aggLevel, aggLevel)
  
  
  for (j in 0:m){                        # Row Position (aggregated matrix)
    for (i in 1:(m+1)){                # Column Position (aggregated matrix)
      for(dj in 0:(aggLevel-1)){     # Row Position (original matrix)
        for(di in 1:aggLevel){  # Column Position (original matrix)
          Xvar[dj+1,di] = X[((j*aggLevel)+dj)*matSize + ((aggLevel*(i-1))+di)]
        }
      }
      Xnew[(j*(m+1)) + i] = var(as.vector(Xvar))
    }
  }
  
  if (average) return(mean(as.vector(Xnew))) else return(Xnew)
  
}


# <summary>
# Calculates the p-value of the ATE or ATT coefficient
# </summary>
# <param name="result"> Output from Match </param>
# <returns> P-Value </returns>
PValueMatching <- function(result){
  
  xbar = result$est
  se = result$se
  n = result$wnobs
  t = xbar/se
  
  return(2*pt(-abs(t), df = n-2))
  
}