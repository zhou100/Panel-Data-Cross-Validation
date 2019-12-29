source("DataFunctions.R")
library(splm)
library(Matching)
library(MatchIt)

# randomTreat Index:
#   0 : purely random allocation (50% treatment)
#   1 : purely clustered allocation (left half of the grid - 50%)
#   2 : dependent on X allocation (50%)
#   3 : purely clustered allocation (upper left quarter of the grid - 25%)


# xType Options:
#   "lagged" : lagged X using 0.9
#   "wwx"    : W*W*X
#   "random" : purely random


################################################################################
###       OLS                                                                ###
################################################################################

OLSEst <- function(trueCoeff, randomTreat, aggList, disSD, nRep, chartTitle, fullRun = TRUE, disSteps = 1, xType = "random", useBinaryAgg = TRUE, corrXT = 0){

    cat("----- OLS ----- \n\n")
    cat("Specification: \n")
    cat(paste("    Coefficients: const =", trueCoeff[1] ,"| X =", trueCoeff[2], "| T =", trueCoeff[3], ifelse(trueCoeff[4] > 0, paste("| Spatial Lag Y =", trueCoeff[4]), "| No Spatial Lag Y"),"\n"))
    cat(paste("    Treatment Type: ", ifelse((randomTreat == 0),"Random Allocation", ifelse((randomTreat == 1), "Contiguous (50%)", ifelse((randomTreat == 2), paste("X Dependent", ifelse(corrXT == 0, "0.75", corrXT)), "Contiguous (25%)"))), "\n"))
    cat(paste("    X Type: ", xType, "\n"))
    cat(paste("    Disaggregation SD: ", disSD, "\n\n"))
    
    cat("    --> Creating Spatial Weight Matrices: \n")
    
    if(fullRun){disSteps = aggList}
    set.seed(1987)

    pb2 = txtProgressBar(0, length(disSteps), style = 3)
    
    for (p in 1:length(disSteps)){
      
      size = (120 / disSteps[p])
      
      eval(parse(text = paste("Wn", p, "= as(as_dgRMatrix_listw(nb2listw(cell2nb(size, size, type = \"queen\"))), \"CsparseMatrix\")", sep = "")))
      eval(parse(text = paste("Wl", p, "= nb2listw(cell2nb(size, size, type = \"queen\"))", sep = "")))
      
      setTxtProgressBar(pb2, p)
      
    }
    
    cat("\n\n")
    cat("    --> Running Simulations: \n")
    
    estBETA0Tab = matrix(0, nrow = nRep, ncol = length(aggList)*length(disSteps))
    estSE0Tab = matrix(0, nrow = nRep, ncol = length(aggList)*length(disSteps))
    estCIL0Tab = matrix(0, nrow = nRep, ncol = length(aggList)*length(disSteps))
    estCIU0Tab = matrix(0, nrow = nRep, ncol = length(aggList)*length(disSteps))
    estBETA1Tab = matrix(0, nrow = nRep, ncol = length(aggList)*length(disSteps))
    estSE1Tab = matrix(0, nrow = nRep, ncol = length(aggList)*length(disSteps))
    estCIL1Tab = matrix(0, nrow = nRep, ncol = length(aggList)*length(disSteps))
    estCIU1Tab = matrix(0, nrow = nRep, ncol = length(aggList)*length(disSteps))
      
    for (d in 1:length(disSteps)){
      
      size = (120 / disSteps[d])
      
      combTab = matrix(0, nrow = length(aggList), ncol = 10)
      
      
      eval(parse(text = paste("Wn = Wn", d, sep = "")))   
      if (randomTreat > -1){eval(parse(text = paste("Wl = Wl", d, sep = "")))}
           
      cat("\n\n")
      cat("    Run: ", d, "- Spatial Process at", size, "x", size, "level\n")
      
      pb = txtProgressBar(0, nRep, style = 3)
      
      
      for (k in 1:nRep){
        
        X = runif(size^2, min = 0, max = 5)   
        
        if (xType != "random"){
          if (xType == "lagged"){
            X = as(powerWeights(Wn, 0.9, X = as.matrix(X)), "matrix")[,1]
          }else{
            X = as(lag.listw(Wl, lag.listw(Wl, X)), "matrix")[,1]
          }
        }

        if (randomTreat == 0){dT = GenerateRandomTreatment((size*size)/2, size)}
        if (randomTreat == 1){dT = GenerateTreatment(1, 1, size, size/2, size)}
        if (randomTreat == 2){dT = GenerateInteractiveTreatment(X, customCorr = corrXT, mySeed = k)}
        if (randomTreat == 3){dT = GenerateTreatment(1, 1, size/2, size/2, size)}
        
        e = rnorm(size^2, mean = 0, sd = 1)
        
        Y = powerWeights(Wn, trueCoeff[4], X = as.matrix(trueCoeff[1] + trueCoeff[2]*X + trueCoeff[3]*dT + e))
        
        # Disaggregation Process
        
        Yd = Disaggregate(as(Y, "matrix")[,1], disSteps[d], disSD)
        Xd = Disaggregate(X, disSteps[d], disSD)
        Td = Disaggregate(dT, disSteps[d], 0)
        aSize = size * disSteps[d]
        
        resultTab = matrix(0, nrow = length(aggList), ncol = 10)
        
        for (i in 1:length(aggList)){
          SPGX = Aggregate(Xd, aggList[i])
          SPGY = Aggregate(Yd, aggList[i])
          SPGT = Aggregate(Td, aggList[i], binaryAgg = useBinaryAgg)
          OLS = lm(SPGY ~ SPGX + SPGT)
          CI = confint.default(OLS) # confidence interval for 95%
          resultTab[i,1] = AggVariance(Xd, aggList[i])
          resultTab[i,2] = summary(OLS)$coef[2,1]
          resultTab[i,3] = summary(OLS)$coef[2,2]
          resultTab[i,4] = summary(OLS)$coef[3,1]
          resultTab[i,5] = summary(OLS)$coef[3,2]
          resultTab[i,7] = CI[2,1] # X ci lower bound
          resultTab[i,8] = CI[2,2] # X ci upper bound
          resultTab[i,9] = CI[3,1] # T ci lower bound
          resultTab[i,10] = CI[3,2] # T ci upper bound       
          if(summary(OLS)$coef[3,4] < 0.05){resultTab[i,6] = resultTab[i,6] + 1}
          estBETA0Tab[k,i+((d-1)*length(disSteps))] = summary(OLS)$coef[2,1]
          estSE0Tab[k,i+((d-1)*length(disSteps))] = summary(OLS)$coef[2,2]
          estCIL0Tab[k,i+((d-1)*length(disSteps))] = CI[2,1]
          estCIU0Tab[k,i+((d-1)*length(disSteps))] = CI[2,2]      
          estBETA1Tab[k,i+((d-1)*length(disSteps))] = summary(OLS)$coef[3,1]
          estSE1Tab[k,i+((d-1)*length(disSteps))] = summary(OLS)$coef[3,2]
          estCIL1Tab[k,i+((d-1)*length(disSteps))] = CI[3,1]
          estCIU1Tab[k,i+((d-1)*length(disSteps))] = CI[3,2] 
        }
        
        combTab = combTab + resultTab
        
        setTxtProgressBar(pb, k)
        
      }
      
      if (d == 1){
        mainTab = (combTab / nRep)
      }else{
        mainTab = rbind(mainTab, (combTab / nRep))
      }
    
    
    colnames(mainTab) = c("Avg.Variance", "Coeff (b)", "s.e.(b)", "Coeff (t)", "s.e.(t)", "Perc Reject Null", "CI(X) LB", "CI(X) UB", "CI(T) LB", "CI(T) UB")
    if (d == 1){
      myRowNames = aggList
    }else{
      if (length(disSteps) > 1){
        myRowNames = c(myRowNames, aggList)
      }
    }
    #if (length(disSteps) > 1){for(d in 2:length(aggList)){myRowNames = c(myRowNames, aggList)}}
    
      
    rownames(mainTab) = myRowNames
    as.table(mainTab)
    
    myCSVName = paste("OLS_SD",disSD, "_Rho", (trueCoeff[4]*10)%%10, "_T-", ifelse((randomTreat == 0),"Rand", ifelse((randomTreat == 1), "Half", ifelse((randomTreat == 2), "Xdep", "Quarter"))),"_X-",  xType, ifelse(useBinaryAgg, "", "_ContAgg"), "_Rep", nRep,".csv", sep = "")
    write.table(mainTab, file = myCSVName, row.names = TRUE, col.names = TRUE, sep = ",")
    
    myCSVName = paste("OLS_SD",disSD, "_Rho", (trueCoeff[4]*10)%%10, "_T-", ifelse((randomTreat == 0),"Rand", ifelse((randomTreat == 1), "Half", ifelse((randomTreat == 2), "Xdep", "Quarter"))),"_X-",  xType, ifelse(useBinaryAgg, "", "_ContAgg"), "_Rep", nRep,"_BETA0.csv", sep = "")
    write.table(estBETA0Tab, file = myCSVName, row.names = FALSE, col.names = FALSE, sep = ",")
    
    myCSVName = paste("OLS_SD",disSD, "_Rho", (trueCoeff[4]*10)%%10, "_T-", ifelse((randomTreat == 0),"Rand", ifelse((randomTreat == 1), "Half", ifelse((randomTreat == 2), "Xdep", "Quarter"))),"_X-",  xType, ifelse(useBinaryAgg, "", "_ContAgg"), "_Rep", nRep,"_SE0.csv", sep = "")
    write.table(estSE0Tab, file = myCSVName, row.names = FALSE, col.names = FALSE, sep = ",")
    
    myCSVName = paste("OLS_SD",disSD, "_Rho", (trueCoeff[4]*10)%%10, "_T-", ifelse((randomTreat == 0),"Rand", ifelse((randomTreat == 1), "Half", ifelse((randomTreat == 2), "Xdep", "Quarter"))),"_X-",  xType, ifelse(useBinaryAgg, "", "_ContAgg"), "_Rep", nRep,"_CIL0.csv", sep = "")
    write.table(estCIL0Tab, file = myCSVName, row.names = FALSE, col.names = FALSE, sep = ",")
    
    myCSVName = paste("OLS_SD",disSD, "_Rho", (trueCoeff[4]*10)%%10, "_T-", ifelse((randomTreat == 0),"Rand", ifelse((randomTreat == 1), "Half", ifelse((randomTreat == 2), "Xdep", "Quarter"))),"_X-",  xType, ifelse(useBinaryAgg, "", "_ContAgg"), "_Rep", nRep,"_CIU0.csv", sep = "")
    write.table(estCIU0Tab, file = myCSVName, row.names = FALSE, col.names = FALSE, sep = ",")
    
    myCSVName = paste("OLS_SD",disSD, "_Rho", (trueCoeff[4]*10)%%10, "_T-", ifelse((randomTreat == 0),"Rand", ifelse((randomTreat == 1), "Half", ifelse((randomTreat == 2), "Xdep", "Quarter"))),"_X-",  xType, ifelse(useBinaryAgg, "", "_ContAgg"), "_Rep", nRep,"_BETA1.csv", sep = "")
    write.table(estBETA1Tab, file = myCSVName, row.names = FALSE, col.names = FALSE, sep = ",")
    
    myCSVName = paste("OLS_SD",disSD, "_Rho", (trueCoeff[4]*10)%%10, "_T-", ifelse((randomTreat == 0),"Rand", ifelse((randomTreat == 1), "Half", ifelse((randomTreat == 2), "Xdep", "Quarter"))),"_X-",  xType, ifelse(useBinaryAgg, "", "_ContAgg"), "_Rep", nRep,"_SE1.csv", sep = "")
    write.table(estSE1Tab, file = myCSVName, row.names = FALSE, col.names = FALSE, sep = ",")
    
    myCSVName = paste("OLS_SD",disSD, "_Rho", (trueCoeff[4]*10)%%10, "_T-", ifelse((randomTreat == 0),"Rand", ifelse((randomTreat == 1), "Half", ifelse((randomTreat == 2), "Xdep", "Quarter"))),"_X-",  xType, ifelse(useBinaryAgg, "", "_ContAgg"), "_Rep", nRep,"_CIL1.csv", sep = "")
    write.table(estCIL1Tab, file = myCSVName, row.names = FALSE, col.names = FALSE, sep = ",")
    
    myCSVName = paste("OLS_SD",disSD, "_Rho", (trueCoeff[4]*10)%%10, "_T-", ifelse((randomTreat == 0),"Rand", ifelse((randomTreat == 1), "Half", ifelse((randomTreat == 2), "Xdep", "Quarter"))),"_X-",  xType, ifelse(useBinaryAgg, "", "_ContAgg"), "_Rep", nRep,"_CIU1.csv", sep = "")
    write.table(estCIU1Tab, file = myCSVName, row.names = FALSE, col.names = FALSE, sep = ",")
    
    myChartName = paste("OLS_SD",disSD, "_Rho", (trueCoeff[4]*10)%%10, "_T-", ifelse((randomTreat == 0),"Rand", ifelse((randomTreat == 1), "Half", ifelse((randomTreat == 2), "Xdep", "Quarter"))),"_X-",  xType, ifelse(useBinaryAgg, "", "_ContAgg"), "_Rep", nRep,"_", sep = "")
    
    # Charts
      png(filename = paste(myChartName, "Level",(120 / disSteps[d]),".png", sep = ""))
      listRows = seq(length(aggList)*(d-1) + 1,length(aggList)*d)
      xAxis = cbind(1, 4, 16, 25, 36, 100, 144)
      polyX = c(rev(xAxis), xAxis)  
      polyYX = c(rev(mainTab[listRows,2] - mainTab[listRows,3]),(mainTab[listRows,2] + mainTab[listRows,3]))
      polyYT = c(rev(mainTab[listRows,4] - mainTab[listRows,5]),(mainTab[listRows,4] + mainTab[listRows,5]))
      
      plot(c(1,144), c(min(mainTab[listRows,2], mainTab[listRows,4], polyYX, polyYT), max(mainTab[listRows,2], mainTab[listRows,4], polyYX, polyYT)), type = "n", main = paste(chartTitle, "(OLS) \n (Spatial Process at the", (disSteps[d]), "x", (disSteps[d]), "Agg Level)"), xlab = "Number of Units in Each Agg Group", ylab = "")
      
      polygon(polyX, polyYX, col = 'grey80', border = NA)
      lines(xAxis, (mainTab[listRows,2] - mainTab[listRows,3]), lty = 'dashed', col = 'red')
      lines(xAxis, (mainTab[listRows,2] + mainTab[listRows,3]), lty = 'dashed', col = 'red')
      
      polygon(polyX, polyYT, col = 'grey80', border = NA)
      lines(xAxis, (mainTab[listRows,4] - mainTab[listRows,5]), lty = 'dashed', col = 'red')
      lines(xAxis, (mainTab[listRows,4] + mainTab[listRows,5]), lty = 'dashed', col = 'red')
      
      #lines(xAxis, mainTab[listRows,1], col = "red", lty = 1, lwd = 2.5)
      lines(xAxis, mainTab[listRows,2], col = "green", lty = 1, lwd = 2.5)
      lines(xAxis, mainTab[listRows,4], col = "blue", lty = 1, lwd = 2.5)
      legend(0, 3.5, c("X", "Treatment Dummy"), cex = 0.8, col = c("green", "blue"), lty = c(1, 1), lwd = c(2.5, 2.5), bty = "n")
      abline(h = trueCoeff[2], col = "black", lty = 2)
      text(140, (trueCoeff[2] + 0.05), "True X coeff", cex = .6)
      abline(h = trueCoeff[3], col = "black", lty = 2)
      text(140, (trueCoeff[3] + 0.05), "True T coeff", cex = .6)
      abline(v = (disSteps[d]^2), col = "black", lty = 3)
      dev.off()
      
    }
    
    cat("\n\n** Monte Carlo Completed \n\n")

}




################################################################################
###       STSLS                                                              ###
################################################################################

STSLSEst<- function(trueCoeff, randomTreat, aggList, disSD, nRep, chartTitle, fullRun = TRUE, disSteps = 1, xType = "random"){
  
    cat("----- STSLS ----- \n\n")
    cat("    --> Creating Spatial Weight Matrices: \n")
    
    if(fullRun){disSteps = aggList}
    set.seed(1987)
    
    pb2 = txtProgressBar(0, length(disSteps), style = 3)
    
    for (p in 1:length(aggList)){
      
      size = (120 / aggList[p])
      
      eval(parse(text = paste("Wn", aggList[p], "= as(as_dgRMatrix_listw(nb2listw(cell2nb(size, size, type = \"queen\"))), \"CsparseMatrix\")", sep = "")))
      
      eval(parse(text = paste("trMat", aggList[p], "= trW(as(as_dgRMatrix_listw(nb2listw(cell2nb(size, size, type = \"queen\"))), \"CsparseMatrix\"), type = \"MC\")", sep = "")))
      
      eval(parse(text = paste("Wl", aggList[p], "= nb2listw(cell2nb(size, size, type = \"queen\"))", sep = "")))
      
      setTxtProgressBar(pb2, p)
      
    }
    
    cat("\n\n    --> Running Simulations: \n")

    for (d in 1:length(disSteps)){
      
      size = (120 / disSteps[d])  
      combTab = matrix(0, nrow = length(aggList), ncol = 20)
      
      eval(parse(text = paste("Wn = Wn", disSteps[d], sep = "")))   
      if (randomTreat > 0){eval(parse(text = paste("Wl = Wl", disSteps[d], sep = "")))}
      
      cat("    Run: ", d, "- Spatial Process at", size, "x", size, "level\n")
      
      pb = txtProgressBar(0, nRep, style = 3)   
      
      for (k in 1:nRep){
        
        X = runif(size^2, min = 0, max = 5)
        
        if (xType != "random"){
          if (xType == "lagged"){
            X = as(powerWeights(Wn, 0.9, X = as.matrix(X)), "matrix")[,1]
          }else{
            X = as(lag.listw(Wl, lag.listw(Wl, X)), "matrix")[,1]
          }
        }
        
        if (randomTreat == 0){dT = GenerateRandomTreatment((size*size)/2, size)}
        if (randomTreat == 1){dT = GenerateTreatment(1, 1, size, size/2, size)}
        if (randomTreat == 2){dT = GenerateInteractiveTreatment(X)}
        if (randomTreat == 3){dT = GenerateTreatment(1, 1, size/2, size/2, size)}
        
        e = rnorm(size^2, mean = 0, sd = 1)
        
        Y = powerWeights(Wn, trueCoeff[4], X = as.matrix(trueCoeff[1] + trueCoeff[2]*X + trueCoeff[3]*dT + e))
        
        # Disaggregation Process
        
        Yd = Disaggregate(as(Y, "matrix")[,1], disSteps[d], disSD)
        Xd = Disaggregate(X, disSteps[d], disSD)
        Td = Disaggregate(dT, disSteps[d], 0)
        aSize = size * disSteps[d]
        
        resultTab = matrix(0, nrow = length(aggList), ncol = 20)
        
        for (i in 1:length(aggList)){
          SPGX = Aggregate(Xd, aggList[i])
          SPGY = Aggregate(Yd, aggList[i])
          SPGT = Aggregate(Td, aggList[i], RT = TRUE)
          eval(parse(text = paste("trMat = trMat", aggList[i], sep = ""))) 
          eval(parse(text = paste("Wls = Wl", aggList[i], sep = ""))) 
          SAR = stsls(SPGY ~ SPGX + SPGT, listw = Wls)
          eval(parse(text = paste("trMat = trMat", aggList[i], sep = ""))) 
          impactTab = summary(impacts(SAR, tr = trMat, R = 2000), zstats = TRUE, short = FALSE)
          resultTab[i,1] = AggVariance(Xd, aggList[i]) # Avg. Agg. Variance
          resultTab[i,2] = summary(SAR)$Coef[3,1]      # Coefficient X
          resultTab[i,3] = summary(SAR)$Coef[3,2]      # S.E. X
          resultTab[i,4] = summary(SAR)$Coef[4,1]      # Coefficient T
          resultTab[i,5] = summary(SAR)$Coef[4,2]      # S.E. T
          resultTab[i,6] = summary(SAR)$Coef[1,1]      # Rho
          resultTab[i,7] = summary(SAR)$Coef[1,2]      # S.E. Rho
          resultTab[i,8] = impactTab$direct_sum$statistics[1,1]       # X Direct Impacts
          resultTab[i,9] = impactTab$indirect_sum$statistics[1,1]     # X Indirect Impacts
          resultTab[i,10] = impactTab$total_sum$statistics[1,1]       # X Total Impacts
          resultTab[i,11] = impactTab$direct_sum$statistics[2,1]      # T Direct Impacts
          resultTab[i,12] = impactTab$indirect_sum$statistics[2,1]    # T Indirect Impacts
          resultTab[i,13] = impactTab$total_sum$statistics[2,1]       # T Total Impacts
          resultTab[i,14] = impactTab$direct_sum$statistics[1,4]      # X Direct Impacts (Time-Series SE)
          resultTab[i,15] = impactTab$indirect_sum$statistics[1,4]    # X Indirect Impacts (Time-Series SE)
          resultTab[i,16] = impactTab$total_sum$statistics[1,4]       # X Total Impacts (Time-Series SE)
          resultTab[i,17] = impactTab$direct_sum$statistics[2,4]      # T Direct Impacts (Time-Series SE)
          resultTab[i,18] = impactTab$indirect_sum$statistics[2,4]    # T Indirect Impacts (Time-Series SE)
          resultTab[i,19] = impactTab$total_sum$statistics[2,4]       # T Total Impacts (Time-Series SE)
          if(summary(SAR)$Coef[4,4] < 0.05){resultTab[i,20] = resultTab[i,20] + 1}
        }
        
        combTab = combTab + resultTab
        
        setTxtProgressBar(pb, k)
        
      }
      
      if (d == 1){
        mainTab = combTab
      }else{
        mainTab = rbind(mainTab, combTab)
      }
    
    mainTab = mainTab / nRep
    
    colnames(mainTab) <- c("Corr(r)", "Coeff (X)", "s.e.(X)", "Coeff (T)", "s.e.(T)", "Rho", "s.e.(rho)", "DI (X)", "II (X)", "TI (X)", "DI (T)", "II (T)", "TI (T)", "s.e. DI (X)", "s.e. II (X)", "s.e. TI (X)", "s.e. DI (T)", "s.e. II (T)", "s.e. TI (T)", "Prob.Rej.Null")
    if (d == 1){
      myRowNames = aggList
    }else{
      if (length(disSteps) > 1){
        myRowNames = c(myRowNames, aggList)
      }
    }
    rownames(mainTab) = myRowNames
    as.table(mainTab)
      
    myCSVName = paste("STSLS_SD",disSD, "_Rho", (trueCoeff[4]*10)%%10, "_T-", ifelse((randomTreat == 0),"Rand", ifelse((randomTreat == 1), "Half", ifelse((randomTreat == 2), "Xdep", "Quarter"))),"_X-", xType,"_Rep", nRep,".csv", sep = "")
    write.table(mainTab, file = myCSVName, row.names = TRUE, col.names = TRUE, sep = ",")
    
    myChartName = paste("STSLS_SD",disSD, "_Rho", (trueCoeff[4]*10)%%10, "_T-", ifelse((randomTreat == 0),"Rand", ifelse((randomTreat == 1), "Half", ifelse((randomTreat == 2), "Xdep", "Quarter"))),"_X-", xType,"_Rep", nRep,"_", sep = "")
    
    
    #Charts  
      png(filename = paste(myChartName, d,".png", sep = ""))
      listRows = seq(length(aggList)*(d-1) + 1,length(aggList)*d)
      xAxis = cbind(1, 4, 16, 25, 36, 100, 144)
      polyX = c(rev(xAxis), xAxis)  
      polyYX = c(rev(mainTab[listRows,2] - mainTab[listRows,3]),(mainTab[listRows,2] + mainTab[listRows,3]))
      polyYT = c(rev(mainTab[listRows,4] - mainTab[listRows,5]),(mainTab[listRows,4] + mainTab[listRows,5]))
      polyYR = c(rev(mainTab[listRows,6] - mainTab[listRows,7]),(mainTab[listRows,6] + mainTab[listRows,7]))
      
      plot(c(1,144), c(min(mainTab[listRows,2], mainTab[listRows,4], polyYX, polyYT, polyYR), max(mainTab[listRows,2], mainTab[listRows,4], polyYX, polyYT, polyYR)), type = "n", main = paste(chartTitle, "(STSLS) \n (Spatial Process at the", (disSteps[d]), "x", (disSteps[d]), "Agg Level)"), xlab = "Number of Units in Each Agg Group", ylab = "")
      
      polygon(polyX, polyYX, col = 'grey80', border = NA)
      lines(xAxis, (mainTab[listRows,2] - mainTab[listRows,3]), lty = 'dashed', col = 'red')
      lines(xAxis, (mainTab[listRows,2] + mainTab[listRows,3]), lty = 'dashed', col = 'red')
      
      polygon(polyX, polyYT, col = 'grey80', border = NA)
      lines(xAxis, (mainTab[listRows,4] - mainTab[listRows,5]), lty = 'dashed', col = 'red')
      lines(xAxis, (mainTab[listRows,4] + mainTab[listRows,5]), lty = 'dashed', col = 'red')
      
      polygon(polyX, polyYR, col = 'grey80', border = NA)
      lines(xAxis, (mainTab[listRows,6] - mainTab[listRows,7]), lty = 'dashed', col = 'red')
      lines(xAxis, (mainTab[listRows,6] + mainTab[listRows,7]), lty = 'dashed', col = 'red')
      
      #lines(xAxis, mainTab[listRows,1], col = "red", lty = 1, lwd = 2.5)
      lines(xAxis, mainTab[listRows,2], col = "green", lty = 1, lwd = 2.5)
      lines(xAxis, mainTab[listRows,4], col = "blue", lty = 1, lwd = 2.5)
      lines(xAxis, mainTab[listRows,6], col = "yellow", lty = 1, lwd = 2.5)
      #legend(0, 3.5, c("Correlation (Y,X)", "X", "Treatment Dummy", "Rho"), cex = 0.8, col = c("red", "green", "blue", "yellow"), lty = c(1, 1, 1, 1), lwd = c(2.5, 2.5, 2.5, 2.5), bty = "n")
      abline(h = trueCoeff[2], col = "black", lty = 2)
      text(140, (trueCoeff[2]), "True X coeff", cex = .6)
      abline(h = trueCoeff[3], col = "black", lty = 2)
      text(140, (trueCoeff[3]), "True T coeff", cex = .6)
      abline(h = trueCoeff[4], col = "black", lty = 2)
      text(140, (trueCoeff[4]), "True rho coeff", cex = .6)
      abline(v = (disSteps[d]^2), col = "black", lty = 3)
      dev.off()
      
    }
    
    # Marginal Effects Charts
    
    myChartName = paste("STSLS_ME_SD",disSD, "_Rho", (trueCoeff[4]*10)%%10, "_T-", ifelse((randomTreat == 0),"Rand", ifelse((randomTreat == 1), "Half", ifelse((randomTreat == 2), "Xdep", "Quarter"))),"_X-", xType,"_Rep", nRep,"_", sep = "")
    
    for (d in 1:length(disSteps)){
      
      png(filename = paste(myChartName, d,".png", sep = ""))
      listRows = seq(length(aggList)*(d-1) + 1,length(aggList)*d)
      xAxis = cbind(1, 4, 16, 25, 36, 100, 144)
      
      plot(c(1,144), c(min(mainTab[listRows,8], mainTab[listRows,10], mainTab[listRows,11], mainTab[listRows,13]), max(mainTab[listRows,8], mainTab[listRows,10], mainTab[listRows,11], mainTab[listRows,13])), type = "n", main = paste(chartTitle, "(STSLS) - Marginal Effects \n (Spatial Process at the", (disSteps[d]), "x", (disSteps[d]), "Agg Level)"), xlab = "Number of Units in Each Agg Group", ylab = "")
      
      lines(xAxis, mainTab[listRows,8], col = "red", lty = 1, lwd = 2.5)
      lines(xAxis, mainTab[listRows,10], col = "red", lty = 2, lwd = 2.5)
      lines(xAxis, mainTab[listRows,11], col = "blue", lty = 1, lwd = 2.5)
      lines(xAxis, mainTab[listRows,13], col = "blue", lty = 2, lwd = 2.5)
      #legend(100, 10, c("Direct Impact (X)", "Total Impact (X)", "Direct Impact (T)", "Total Impact (T)"), cex = 0.8, col = c("red", "red", "blue", "blue"), lty = c(1, 2, 1, 2), lwd = c(2.5, 2.5, 2.5, 2.5), bty = "n")
      abline(v = (disSteps[d]^2), col = "black", lty = 3)
      dev.off()
      
    }
    
    cat("\n\n** Monte Carlo Completed \n\n")

}




################################################################################
###       Diff-in-Diff                                                       ###
################################################################################


#####  Simple DID  #####

NonSpatialDID <- function(trueCoeff, randomTreat, aggList, disSD, nRep, chartTitle, fullRun = TRUE, disSteps = 1, xType = "lagged"){
  
    cat("----- NonSpatial DID via OLS ----- \n\n")
    cat("    --> Creating Spatial Weight Matrices: \n")
    
    if(fullRun){disSteps = aggList}
    set.seed(1987)
    
    pb2 = txtProgressBar(0, length(disSteps), style = 3)
    
    for (p in 1:length(disSteps)){
      
      size = (120 / disSteps[p])
      
      eval(parse(text = paste("Wn", p, "= as(as_dgRMatrix_listw(nb2listw(cell2nb(size, size, type = \"queen\"))), \"CsparseMatrix\")", sep = "")))
      eval(parse(text = paste("Wl", p, "= nb2listw(cell2nb(size, size, type = \"queen\"))", sep = "")))
      
      setTxtProgressBar(pb2, p)
      
    }
    
    cat("\n\n    --> Running Simulations: \n")
      
    for (d in 1:length(disSteps)){
         
      size = (120 / disSteps[d])
      
      combTab = matrix(0, nrow = length(aggList), ncol = 7)
      
      eval(parse(text = paste("Wn = Wn", d, sep = "")))   
      if (randomTreat > 0){eval(parse(text = paste("Wl = Wl", d, sep = "")))}
      
      cat("    Run: ", d, "- Spatial Process at", size, "x", size, "level\n")
      
      pb = txtProgressBar(0, nRep, style = 3)
     
      for (k in 1:nRep){
        
        X1 = runif(size^2, min = 0, max = 5)   
        
        if (randomTreat > 0){
          if (xType == "lagged"){
            X1 = as(powerWeights(Wn, 0.9, X = as.matrix(X1)), "matrix")[,1]
          }else{
            X1 = as(lag.listw(Wl, lag.listw(Wl, X1)), "matrix")[,1]
          }
        }
                
        if (randomTreat == 0){dT = GenerateRandomTreatment((size*size)/2, size)}
        if (randomTreat == 1){dT = GenerateTreatment(1, 1, size, size/2, size)}
        if (randomTreat == 2){dT = GenerateInteractiveTreatment(X1)}
        if (randomTreat == 3){dT = GenerateTreatment(1, 1, size/2, size/2, size)}
        
        X2 = ifelse(dT == 1, X1 + trueCoeff[3], X1)
        
        ei = rnorm(size^2, mean = 0, sd = 1)
        et1 = rep(rnorm(1, mean = 0, sd = 1), size^2)
        et2 = rep(rnorm(1, mean = 1, sd = .5), size^2)
        
        
        Y1 = powerWeights(Wn, trueCoeff[4], X = as.matrix(trueCoeff[1] + trueCoeff[2]*X1 + ei + et1))   
        Y2 = powerWeights(Wn, trueCoeff[4], X = as.matrix(trueCoeff[1] + trueCoeff[2]*X2 + ei + et2))
        
        
        # Disaggregation Process
        
        Yd1 = Disaggregate(as(Y1, "matrix")[,1], disSteps[d], disSD)
        Yd2 = Disaggregate(as(Y2, "matrix")[,1], disSteps[d], disSD)
        Td = Disaggregate(dT, disSteps[d], 0)
        aSize = size * disSteps[d]
        
        resultTab = matrix(0, nrow = length(aggList), ncol = 7)
        
        for (i in 1:length(aggList)){
          SPGY1 = Aggregate(Yd1, aggList[i])
          SPGY2 = Aggregate(Yd2, aggList[i])
          SPGT = Aggregate(Td, aggList[i], RT = TRUE)
          Y = c(SPGY1, SPGY2) 
          yearDummy = c(rep(0, (aSize/aggList[i])^2), rep(1, (aSize/aggList[i])^2))   
          treatDummy = c(SPGT, SPGT)    
          interaction = yearDummy*treatDummy
          linearDID = lm(Y ~ yearDummy + treatDummy + interaction)
          resultTab[i,1] = summary(linearDID)$coef[2,1]
          resultTab[i,2] = summary(linearDID)$coef[2,2]
          resultTab[i,3] = summary(linearDID)$coef[3,1]
          resultTab[i,4] = summary(linearDID)$coef[3,2]
          resultTab[i,5] = summary(linearDID)$coef[4,1]
          resultTab[i,6] = summary(linearDID)$coef[4,2]
          if(summary(linearDID)$coef[4,4] < 0.05){resultTab[i,7] = resultTab[i,7] + 1}
        }
        
        combTab = combTab + resultTab
        
        setTxtProgressBar(pb, k)
        
      }
      
      if (d == 1){
        mainTab = combTab
      }else{
        mainTab = rbind(mainTab, combTab)
      }
      
    }
    
    mainTab = mainTab / nRep
    
    colnames(mainTab) <- c("Coeff (year)", "s.e.(year)", "Coeff (T)", "s.e.(T)", "Coeff (TE)", "s.e.(TE)", "Prob.Rej.Null")
    myRowNames = aggList
    if (length(disSteps) > 1){for(d in 2:length(aggList)){myRowNames = c(myRowNames, aggList)}}
    rownames(mainTab) <- myRowNames
    as.table(mainTab)
    
    myCSVName = paste("nsDID_SD",disSD, "_Rho", (trueCoeff[4]*10)%%10, "_", ifelse((randomTreat == 0),"R", ifelse((randomTreat == 1), "H", ifelse((randomTreat == 2), "X", "Q"))),"_Rep", nRep,".csv", sep = "")
    write.table(mainTab, file = myCSVName, row.names = TRUE, col.names = TRUE, sep = ",")
    
    myChartName = paste("nsDID_SD",disSD, "_Rho", (trueCoeff[4]*10)%%10, "_", ifelse((randomTreat == 0),"R", ifelse((randomTreat == 1), "H", ifelse((randomTreat == 2), "X", "Q"))),"_Rep", nRep,"_", sep = "")
    
    for (d in 1:length(disSteps)){
      
      png(filename = paste(myChartName, d,".png", sep = ""))
      listRows = seq(length(aggList)*(d-1) + 1,length(aggList)*d)
      xAxis = cbind(1, 4, 16, 25, 36, 100, 144)
      polyX = c(rev(xAxis), xAxis)
      polyYT = c(rev(mainTab[listRows,5] - mainTab[listRows,6]),(mainTab[listRows,5] + mainTab[listRows,6]))
      
      plot(c(1,144), c(min(mainTab[listRows,5], polyYT), max(mainTab[listRows,5], polyYT)), type = "n", main = paste(chartTitle, "(Non-Spatial DID) \n (Spatial Process at the", (disSteps[d]), "x", (disSteps[d]), "Agg Level)"), xlab = "Number of Units in Each Agg Group", ylab = "")
      
      polygon(polyX, polyYT, col = 'grey80', border = NA)
      lines(xAxis, (mainTab[listRows,5] - mainTab[listRows,6]), lty = 'dashed', col = 'red')
      lines(xAxis, (mainTab[listRows,5] + mainTab[listRows,6]), lty = 'dashed', col = 'red')
      
      lines(xAxis, mainTab[listRows,5], col = "blue", lty = 1, lwd = 2.5)
      #legend(0, 3.5, c("Treatment Effect"), cex = 0.8, col = c("blue"), lty = c(1), lwd = c(2.5), bty = "n")
      abline(h = (trueCoeff[3]), col = "black", lty = 2)
      text(140, (trueCoeff[3] + 0.02), "True T coeff", cex = .6)
      abline(v = (disSteps[d]^2), col = "black", lty = 3)
      dev.off()
      
    }
    
    cat("\n\n** Monte Carlo Completed \n\n")

}



#####  Spatial DID  #####

SpatialDID <- function(trueCoeff, randomTreat, aggList, disSD, nRep, chartTitle){
  
    cat("----- Spatial DID via GMM ----- \n\n")
    cat("    --> Creating Spatial Weight Matrices: \n")
    
    disSteps = aggList
    set.seed(1987)
    
    pb2 = txtProgressBar(1, length(disSteps), style = 3)
    
    for (p in 1:length(disSteps)){
      
      size = (120 / disSteps[p])
      
      eval(parse(text = paste("Wn", p, "= as(as_dgRMatrix_listw(nb2listw(cell2nb(size, size, type = \"queen\"))), \"CsparseMatrix\")", sep = "")))
      eval(parse(text = paste("Wl", p, "= nb2listw(cell2nb(size, size, type = \"queen\"))", sep = "")))
      
      setTxtProgressBar(pb2, p)
      
    }
    
    cat("\n\n    --> Running Simulations: \n")

    for (d in 1:length(disSteps)){
      
      size = (120 / disSteps[d])
      
      combTab = matrix(0, nrow = length(aggList), ncol = 9)
      
      eval(parse(text = paste("Wn = Wn", d, sep = "")))   
      if (randomTreat > 0){eval(parse(text = paste("Wl = Wl", d, sep = "")))}
      
      cat("    Run: ", d, "\n")
      
      pb = txtProgressBar(0, nRep, style = 3)
      
      for (k in 1:nRep){
        
        X1 = runif(size^2, min = 0, max = 5)
        if (randomTreat > 0){X1 = as(powerWeights(Wn, 0.9, X = as.matrix(X1)), "matrix")[,1]}
        
        if (randomTreat == 0){dT = GenerateRandomTreatment((size*size)/2, size)}
        if (randomTreat == 1){dT = GenerateTreatment(1, 1, size, size/2, size)}
        if (randomTreat == 2){dT = GenerateInteractiveTreatment(X1)}
        if (randomTreat == 3){dT = GenerateTreatment(1, 1, size/2, size/2, size)}
        
        X2 = ifelse(dT == 1, X1 + trueCoeff[3], X1)
        
        ei = rnorm(size^2, mean = 0, sd = 1)
        et1 = rep(rnorm(1, mean = 0, sd = 1), size^2)
        et2 = rep(rnorm(1, mean = 1, sd = .5), size^2)
        
        
        Y1 = powerWeights(Wn, trueCoeff[4], X = as.matrix(trueCoeff[1] + trueCoeff[2]*X1 + ei + et1))   
        Y2 = powerWeights(Wn, trueCoeff[4], X = as.matrix(trueCoeff[1] + trueCoeff[2]*X2 + ei + et2))
        
        
        # Disaggregation Process
        
        Yd1 = Disaggregate(as(Y1, "matrix")[,1], disSteps[d], disSD)
        Yd2 = Disaggregate(as(Y2, "matrix")[,1], disSteps[d], disSD)
        Td = Disaggregate(dT, disSteps[d], 0)
        aSize = size * disSteps[d]
        
        resultTab = matrix(0, nrow = length(aggList), ncol = 9)
        
        for (i in 1:length(aggList)){
          SPGY1 = Aggregate(Yd1, aggList[i])
          SPGY2 = Aggregate(Yd2, aggList[i])
          SPGT = Aggregate(Td, aggList[i], RT = TRUE)
          eval(parse(text = paste("Wls = Wl", i, sep = ""))) 
          Y = c(SPGY1, SPGY2)    
          yearDummy = c(rep(0, (aSize/aggList[i])^2), rep(1, (aSize/aggList[i])^2))   
          treatDummy = c(SPGT, SPGT)    
          interaction = yearDummy*treatDummy
          id = c(seq(1, (aSize/aggList[i])^2, by = 1), seq(1, (aSize/aggList[i])^2, by = 1))
          time = c(rep(2000,(aSize/aggList[i])^2), rep(2001, (aSize/aggList[i])^2))
          myData = as.data.frame(cbind(id, time, Y, yearDummy, treatDummy, interaction))   
          spatialDID = spgm((Y ~ yearDummy + treatDummy + interaction), data = myData, listw = Wls, model = "within", lag = TRUE, moments = "fullweights", method = "w2sls")
          resultTab[i,1] = summary(spatialDID)$Coef[1,1]
          resultTab[i,2] = summary(spatialDID)$Coef[1,2]
          resultTab[i,3] = summary(spatialDID)$Coef[2,1]
          resultTab[i,4] = summary(spatialDID)$Coef[2,2]
          resultTab[i,5] = summary(spatialDID)$Coef[3,1]
          resultTab[i,6] = summary(spatialDID)$Coef[3,2]
          resultTab[i,7] = summary(spatialDID)$Coef[4,1]
          resultTab[i,8] = summary(spatialDID)$Coef[4,2]
          if(summary(spatialDID)$Coef[4,4] < 0.05){resultTab[i,9] = resultTab[i,9] + 1}
        }
        
        combTab = combTab + resultTab
        
        setTxtProgressBar(pb, k)
        
      }
      
      if (d == 1){
        mainTab = combTab
      }else{
        mainTab = rbind(mainTab, combTab)
      }
      
    }
    
    mainTab = mainTab / nRep
    
    colnames(mainTab) <- c("Coeff (lambda)", "s.e.(lambda)", "Coeff (year)", "s.e.(year)", "Coeff (T)", "s.e.(T)", "Coeff (TE)", "s.e.(TE)", "Prob.Rej.Null")
    myRowNames = aggList
    for(d in 2:length(disSteps)){myRowNames = c(myRowNames, aggList)}
    rownames(mainTab) <- myRowNames
    as.table(mainTab)
    
    myCSVName = paste("sDID_SD",disSD, "_Rho", (trueCoeff[4]*10)%%10, "_", ifelse((randomTreat == 0),"R", ifelse((randomTreat == 1), "H", ifelse((randomTreat == 2), "X", "Q"))),"_Rep", nRep,".csv", sep = "")
    write.table(mainTab, file = myCSVName, row.names = TRUE, col.names = TRUE, sep = ",")
    
    myChartName = paste("sDID_SD",disSD, "_Rho", (trueCoeff[4]*10)%%10, "_", ifelse((randomTreat == 0),"R", ifelse((randomTreat == 1), "H", ifelse((randomTreat == 2), "X", "Q"))),"_Rep", nRep,"_", sep = "")
    
    for (d in 1:length(disSteps)){
      
      png(filename = paste(myChartName, d,".png", sep = ""))
      listRows = seq(length(aggList)*(d-1) + 1,length(aggList)*d)
      xAxis = cbind(1, 4, 16, 25, 36, 100, 144)
      polyX = c(rev(xAxis), xAxis)  
      polyYT = c(rev(mainTab[listRows,7] - mainTab[listRows,8]),(mainTab[listRows,7] + mainTab[listRows,8]))
      
      plot(c(1,144), c(min(mainTab[listRows,7], polyYT), max(mainTab[listRows,7], polyYT)), type = "n", main = paste(chartTitle, "(Spatial DID) \n (Spatial Process at the", (disSteps[d]), "x", (disSteps[d]), "Agg Level)"), xlab = "Number of Units in Each Agg Group", ylab = "")
      
      polygon(polyX, polyYT, col = 'grey80', border = NA)
      lines(xAxis, (mainTab[listRows,7] - mainTab[listRows,8]), lty = 'dashed', col = 'red')
      lines(xAxis, (mainTab[listRows,7] + mainTab[listRows,8]), lty = 'dashed', col = 'red')
      
      lines(xAxis, mainTab[listRows,7], col = "blue", lty = 1, lwd = 2.5)
      #legend(0, 3.5, c("Average Treatment Effect"), cex = 0.8, col = c("blue"), lty = c(1), lwd = c(2.5), bty = "n")
      abline(h = trueCoeff[3], col = "black", lty = 2)
      #text(140, 5.1, "True T coeff", cex = .6)
      abline(v = (disSteps[d]^2), col = "black", lty = 3)
      dev.off()
      
    }
    
    cat("\n\n** Monte Carlo Completed \n\n")
    
}




################################################################################
###       Matching Estimator                                                 ###
################################################################################


#####  Simple Matching  #####

NonSpatialMatching <- function(trueCoeff, randomTreat, aggList, disSD, nRep, chartTitle, fullRun = TRUE, disSteps = 1, xType = "random", useBinaryAgg = TRUE, corrXT = 0){
  
    cat("----- Non-Spatial Matching ----- \n\n")
    cat("Specification: \n")
    cat(paste("    Coefficients: const =", trueCoeff[1] ,"| X =", trueCoeff[2], "| T =", trueCoeff[3], ifelse(trueCoeff[4] > 0, paste("| Spatial Lag Y =", trueCoeff[4]), "| No Spatial Lag Y"),"\n"))
    cat(paste("    Treatment Type: ", ifelse((randomTreat == 0),"Random Allocation", ifelse((randomTreat == 1), "Contiguous (50%)", ifelse((randomTreat == 2), paste("X Dependent", ifelse(corrXT == 0, "0.75", corrXT)), "Contiguous (25%)"))), "\n"))
    cat(paste("    X Type: ", xType, "\n"))
    cat(paste("    Disaggregation SD: ", disSD, "\n\n"))
     
    cat("    --> Creating Spatial Weight Matrices: \n")
    
    if(fullRun){disSteps = aggList}
    set.seed(1987)
    
    pb2 = txtProgressBar(0, length(disSteps), style = 3)
    
    for (p in 1:length(disSteps)){
      
      size = (120 / disSteps[p])
      
      eval(parse(text = paste("Wn", p, "= as(as_dgRMatrix_listw(nb2listw(cell2nb(size, size, type = \"queen\"))), \"CsparseMatrix\")", sep = "")))
      eval(parse(text = paste("Wl", p, "= nb2listw(cell2nb(size, size, type = \"queen\"))", sep = "")))
      
      setTxtProgressBar(pb2, p)
      
    }
    
    cat("\n\n")
    cat("    --> Running Simulations: \n")

    estATETab = matrix(0, nrow = nRep, ncol = length(aggList)*length(disSteps))
    estATEseTab = matrix(0, nrow = nRep, ncol = length(aggList)*length(disSteps))
    estATErejTab = matrix(0, nrow = nRep, ncol = length(aggList)*length(disSteps))
    estATTTab = matrix(0, nrow = nRep, ncol = length(aggList)*length(disSteps))
    estATTseTab = matrix(0, nrow = nRep, ncol = length(aggList)*length(disSteps))
    estATTrejTab = matrix(0, nrow = nRep, ncol = length(aggList)*length(disSteps))
     
    
    
    for (d in 1:length(disSteps)){
      
      size = (120 / disSteps[d])
      
      combTab = matrix(0, nrow = length(aggList), ncol = 6) # Avg. Treat. Effect (Est + se + Prob. Rej. Null), Avg. Treat. Effect for the Treated (Est + se + Prob. Rej. Null)
      
      eval(parse(text = paste("Wn = Wn", d, sep = "")))   
      if (randomTreat > -1){eval(parse(text = paste("Wl = Wl", d, sep = "")))}
      
      cat("\n\n")
      cat("    Run: ", d, "- Spatial Process at", size, "x", size, "level\n")
      
      pb = txtProgressBar(0, nRep, style = 3)

      for (k in 1:nRep){
        
        X = runif(size^2, min = 0, max = 5)
        
        if (xType != "random"){
          if (xType == "lagged"){
            X = as(powerWeights(Wn, 0.9, X = as.matrix(X)), "matrix")[,1]
          }else{
            X = as(lag.listw(Wl, lag.listw(Wl, X)), "matrix")[,1]
          }
        }
        
        if (randomTreat == 0){dT = GenerateRandomTreatment((size*size)/2, size)}
        if (randomTreat == 1){dT = GenerateTreatment(1, 1, size, size/2, size)}
        if (randomTreat == 2){dT = GenerateInteractiveTreatment(X, customCorr = corrXT, mySeed = k)}
        if (randomTreat == 3){dT = GenerateTreatment(1, 1, size/2, size/2, size)}
        
        e = rnorm(size^2, mean = 0, sd = 1)
        
        Y = powerWeights(Wn, trueCoeff[4], X = as.matrix(trueCoeff[1] + trueCoeff[2]*X + trueCoeff[3]*dT + e))
        
        # Disaggregation Process
        
        Yd = Disaggregate(as(Y, "matrix")[,1], disSteps[d], disSD)
        Xd = Disaggregate(X, disSteps[d], disSD)
        Td = Disaggregate(dT, disSteps[d], 0)
        aSize = size * disSteps[d]
        
        resultTab = matrix(0, nrow = length(aggList), ncol = 6)
        
        for (i in 1:length(aggList)){
          SPGX = Aggregate(Xd, aggList[i])
          SPGY = Aggregate(Yd, aggList[i])
          SPGT = Aggregate(Td, aggList[i], binaryAgg = useBinaryAgg)
          glm1 = glm(SPGT ~ SPGX, family = binomial)
          
          mout = Match(Y = SPGY, Tr = SPGT, X = glm1$fitted, estimand = "ATE")   
          resultTab[i,1] = mout$est
          resultTab[i,2] = mout$se
          if(PValueMatching(mout) < 0.05){resultTab[i,3] = resultTab[i,3] + 1}
          estATETab[k,i+((d-1)*length(disSteps))] = mout$est
          estATEseTab[k,i+((d-1)*length(disSteps))] = mout$se
          estATErejTab[k,i+((d-1)*length(disSteps))] = PValueMatching(mout)
          
          mout = Match(Y = SPGY, Tr = SPGT, X = glm1$fitted, estimand = "ATT")
          resultTab[i,4] = mout$est
          resultTab[i,5] = mout$se
          if(PValueMatching(mout) < 0.05){resultTab[i,6] = resultTab[i,6] + 1} 
          estATTTab[k,i+((d-1)*length(disSteps))] = mout$est     
          estATTseTab[k,i+((d-1)*length(disSteps))] = mout$se
          estATTrejTab[k,i+((d-1)*length(disSteps))] = PValueMatching(mout)

        }
        
        combTab = combTab + resultTab
        
        setTxtProgressBar(pb, k)
        
      }
      
      if (d == 1){
        mainTab = (combTab / nRep)
      }else{
        mainTab = rbind(mainTab, (combTab / nRep))
      }
      
    colnames(mainTab) <- c("ATE Treat Est", "ATE s.e.", "ATE Prob.Rej.Null", "ATT Treat Est", "ATT s.e.", "ATT Prob.Rej.Null")
    if (d == 1){
      myRowNames = aggList
    }else{
      if (length(disSteps) > 1){
        myRowNames = c(myRowNames, aggList)
      }
    }

    rownames(mainTab) = myRowNames
    as.table(mainTab)
    
    myCSVName = paste("nsMatch_SD",disSD, "_Rho", (trueCoeff[4]*10)%%10, "_T-", ifelse((randomTreat == 0),"Rand", ifelse((randomTreat == 1), "Half", ifelse((randomTreat == 2), "Xdep", "Quarter"))),"_X-",  xType, ifelse(useBinaryAgg, "", "_ContAgg"), "_Rep", nRep,".csv", sep = "")
    write.table(mainTab, file = myCSVName, row.names = TRUE, col.names = TRUE, sep = ",")

    myCSVName = paste("nsMatch_SD",disSD, "_Rho", (trueCoeff[4]*10)%%10, "_T-", ifelse((randomTreat == 0),"Rand", ifelse((randomTreat == 1), "Half", ifelse((randomTreat == 2), "Xdep", "Quarter"))),"_X-",  xType, ifelse(useBinaryAgg, "", "_ContAgg"), "_Rep", nRep,"_ATE.csv", sep = "")
    write.table(estATETab, file = myCSVName, row.names = FALSE, col.names = FALSE, sep = ",")

    myCSVName = paste("nsMatch_SD",disSD, "_Rho", (trueCoeff[4]*10)%%10, "_T-", ifelse((randomTreat == 0),"Rand", ifelse((randomTreat == 1), "Half", ifelse((randomTreat == 2), "Xdep", "Quarter"))),"_X-",  xType, ifelse(useBinaryAgg, "", "_ContAgg"), "_Rep", nRep,"_ATEse.csv", sep = "")
    write.table(estATEseTab, file = myCSVName, row.names = FALSE, col.names = FALSE, sep = ",")
    
    myCSVName = paste("nsMatch_SD",disSD, "_Rho", (trueCoeff[4]*10)%%10, "_T-", ifelse((randomTreat == 0),"Rand", ifelse((randomTreat == 1), "Half", ifelse((randomTreat == 2), "Xdep", "Quarter"))),"_X-",  xType, ifelse(useBinaryAgg, "", "_ContAgg"), "_Rep", nRep,"_ATErej.csv", sep = "")
    write.table(estATErejTab, file = myCSVName, row.names = FALSE, col.names = FALSE, sep = ",")

    myCSVName = paste("nsMatch_SD",disSD, "_Rho", (trueCoeff[4]*10)%%10, "_T-", ifelse((randomTreat == 0),"Rand", ifelse((randomTreat == 1), "Half", ifelse((randomTreat == 2), "Xdep", "Quarter"))),"_X-",  xType, ifelse(useBinaryAgg, "", "_ContAgg"), "_Rep", nRep,"_ATT.csv", sep = "")
    write.table(estATTTab, file = myCSVName, row.names = FALSE, col.names = FALSE, sep = ",")
    
    myCSVName = paste("nsMatch_SD",disSD, "_Rho", (trueCoeff[4]*10)%%10, "_T-", ifelse((randomTreat == 0),"Rand", ifelse((randomTreat == 1), "Half", ifelse((randomTreat == 2), "Xdep", "Quarter"))),"_X-",  xType, ifelse(useBinaryAgg, "", "_ContAgg"), "_Rep", nRep,"_ATTse.csv", sep = "")
    write.table(estATTseTab, file = myCSVName, row.names = FALSE, col.names = FALSE, sep = ",")
    
    myCSVName = paste("nsMatch_SD",disSD, "_Rho", (trueCoeff[4]*10)%%10, "_T-", ifelse((randomTreat == 0),"Rand", ifelse((randomTreat == 1), "Half", ifelse((randomTreat == 2), "Xdep", "Quarter"))),"_X-",  xType, ifelse(useBinaryAgg, "", "_ContAgg"), "_Rep", nRep,"_ATTrej.csv", sep = "")
    write.table(estATTrejTab, file = myCSVName, row.names = FALSE, col.names = FALSE, sep = ",")

    
    myChartName = paste("nsMatch_SD",disSD, "_Rho", (trueCoeff[4]*10)%%10, "_T-", ifelse((randomTreat == 0),"Rand", ifelse((randomTreat == 1), "Half", ifelse((randomTreat == 2), "Xdep", "Quarter"))),"_X-",  xType, ifelse(useBinaryAgg, "", "_ContAgg"), "_Rep", nRep,"_", sep = "")
    
    # Charts  
    png(filename = paste(myChartName, "Level",(120 / disSteps[d]),".png", sep = ""))
    listRows = seq(length(aggList)*(d-1) + 1,length(aggList)*d)
    xAxis = cbind(1, 4, 16, 25, 36, 100, 144)
    polyX = c(rev(xAxis), xAxis) 
    polyATE = c(rev(mainTab[listRows,1] - mainTab[listRows,2]),(mainTab[listRows,1] + mainTab[listRows,2]))
    polyATT = c(rev(mainTab[listRows,4] - mainTab[listRows,5]),(mainTab[listRows,4] + mainTab[listRows,5]))
      
    plot(c(1,144), c(min(mainTab[listRows,1], polyATE, polyATT), max(mainTab[listRows,1], polyATE, polyATT)), type = "n", main = paste(chartTitle, "(Simple Matching) \n (Spatial Process at the", (disSteps[d]), "x", (disSteps[d]), "Agg Level)"), xlab = "Number of Units in Each Agg Group", ylab = "")
      
    polygon(polyX, polyATE, col = 'grey80', border = NA)
    lines(xAxis, (mainTab[listRows,1] - mainTab[listRows,2]), lty = 'dashed', col = 'red')
    lines(xAxis, (mainTab[listRows,1] + mainTab[listRows,2]), lty = 'dashed', col = 'red')
      
    polygon(polyX, polyATT, col = 'grey80', border = NA)
    lines(xAxis, (mainTab[listRows,4] - mainTab[listRows,5]), lty = 'dashed', col = 'red')
    lines(xAxis, (mainTab[listRows,4] + mainTab[listRows,5]), lty = 'dashed', col = 'red')
      
    lines(xAxis, mainTab[listRows,1], col = "blue", lty = 1, lwd = 2.5)
    lines(xAxis, mainTab[listRows,4], col = "green", lty = 1, lwd = 2.5)
    #legend(0, 3.5, c("ATE", "ATT"), cex = 0.8, col = c("blue", "green"), lty = c(1, 1), lwd = c(2.5, 2.5), bty = "n")
    abline(h = trueCoeff[3], col = "black", lty = 2)
    text(140, (trueCoeff[3] + 0.1), "True T coeff", cex = .6)
    abline(v = (disSteps[d]^2), col = "black", lty = 3)
    dev.off()
      
  }
    
    cat("\n\n** Monte Carlo Completed \n\n")

}



#####  Spatial Matching (Abadie and Imbens (2006)) #####

SpatialMatching <- function(trueCoeff, randomTreat, aggList, disSD, nRep, chartTitle){
  
    cat("----- Spatial Matching ----- \n\n")
    cat("    --> Creating Spatial Weight Matrices: \n")
    
    disSteps = aggList
    set.seed(1987)
    
    pb2 = txtProgressBar(1, length(disSteps), style = 3)
    
    for (p in 1:length(disSteps)){
      
      size = (120 / disSteps[p])
      
      eval(parse(text = paste("Wn", p, "= as(as_dgRMatrix_listw(nb2listw(cell2nb(size, size, type = \"queen\"))), \"CsparseMatrix\")", sep = "")))
      eval(parse(text = paste("Wl", p, "= nb2listw(cell2nb(size, size, type = \"queen\"))", sep = "")))
      
      setTxtProgressBar(pb2, p)
      
    }
    
    cat("\n\n    --> Running Simulations: \n")
    
    for (d in 1:length(disSteps)){
      
      size = (120 / disSteps[d])
      
      combTab = matrix(0, nrow = length(aggList), ncol = 6) # Avg. Treat. Effect (Est + se + Prob. Rej. Null), Avg. Treat. Effect for the Treated (Est + se + Prob. Rej. Null)
      
      eval(parse(text = paste("Wn = Wn", d, sep = "")))   
      if (randomTreat > 0){eval(parse(text = paste("Wl = Wl", d, sep = "")))}
      
      cat("    Run: ", d, "\n")
      
      pb = txtProgressBar(0, nRep, style = 3)
 
      for (k in 1:nRep){
        
        X = runif(size^2, min = 0, max = 5)
        if (randomTreat > 0){X = as(powerWeights(Wn, 0.9, X = as.matrix(X)), "matrix")[,1]}
        
        if (randomTreat == 0){dT = GenerateRandomTreatment((size*size)/2, size)}
        if (randomTreat == 1){dT = GenerateTreatment(1, 1, size, size/2, size)}
        if (randomTreat == 2){dT = GenerateInteractiveTreatment(X)}
        if (randomTreat == 3){dT = GenerateTreatment(1, 1, size/2, size/2, size)}
        
        e = rnorm(size^2, mean = 0, sd = 1)
        
        Y = powerWeights(Wn, trueCoeff[4], X = as.matrix(trueCoeff[1] + trueCoeff[2]*X + trueCoeff[3]*dT + e))
        
        
        # Disaggregation Process
        
        Yd = Disaggregate(as(Y, "matrix")[,1], disSteps[d], disSD)
        Xd = Disaggregate(X, disSteps[d], disSD)
        Td = Disaggregate(dT, disSteps[d], 0)
        aSize = size * disSteps[d]
        
        resultTab = matrix(0, nrow = length(aggList), ncol = 6)
        
        for (i in 1:length(aggList)){
          SPGX = Aggregate(Xd, aggList[i])
          SPGY = Aggregate(Yd, aggList[i])
          SPGT = Aggregate(Td, aggList[i], RT = TRUE)
          eval(parse(text = paste("Wls = Wl", i, sep = "")))
          glm1 = glm(SPGT ~ SPGX + lag.listw(Wls, SPGX), family = binomial)
          mout = Match(Y = SPGY, Tr = SPGT, X = glm1$fitted, BiasAdj = TRUE, estimand = "ATE")   
          resultTab[i,1] = mout$est
          resultTab[i,2] = mout$se
          if(PValueMatching(mout) < 0.05){resultTab[i,3] = resultTab[i,3] + 1}
          mout = Match(Y = SPGY, Tr = SPGT, X = glm1$fitted, estimand = "ATT", BiasAdj = TRUE)
          resultTab[i,4] = mout$est
          resultTab[i,5] = mout$se
          if(PValueMatching(mout) < 0.05){resultTab[i,6] = resultTab[i,6] + 1}
        }
        
        combTab = combTab + resultTab
        
        setTxtProgressBar(pb, k)
        
      }
      
      if (d == 1){
        mainTab = combTab
      }else{
        mainTab = rbind(mainTab, combTab)
      }
      
    }
    
    mainTab = mainTab / nRep
    
    colnames(mainTab) <- c("ATE Treat Est", "ATE s.e.", "ATE Prob.Rej.Null", "ATT Treat Est", "ATT s.e.", "ATT Prob.Rej.Null")
    myRowNames = aggList
    for(d in 2:length(disSteps)){myRowNames = c(myRowNames, aggList)}
    rownames(mainTab) <- myRowNames
    as.table(mainTab)
    
    myCSVName = paste("sMatch_SD",disSD, "_Rho", (trueCoeff[4]*10)%%10, "_", ifelse((randomTreat == 0),"R", ifelse((randomTreat == 1), "H", ifelse((randomTreat == 2), "X", "Q"))),"_Rep", nRep,".csv", sep = "")
    write.table(mainTab, file = myCSVName, row.names = TRUE, col.names = TRUE, sep = ",")
    
    myChartName = paste("sMatch_SD",disSD, "_Rho", (trueCoeff[4]*10)%%10, "_", ifelse((randomTreat == 0),"R", ifelse((randomTreat == 1), "H", ifelse((randomTreat == 2), "X", "Q"))),"_Rep", nRep,"_", sep = "")
    
    for (d in 1:length(disSteps)){
      
      png(filename = paste(myChartName, d,".png", sep = ""))
      listRows = seq(length(aggList)*(d-1) + 1,length(aggList)*d)
      xAxis = cbind(1, 4, 16, 25, 36, 100, 144)
      polyX = c(rev(xAxis), xAxis) 
      polyATE = c(rev(mainTab[listRows,1] - mainTab[listRows,2]),(mainTab[listRows,1] + mainTab[listRows,2]))
      polyATT = c(rev(mainTab[listRows,4] - mainTab[listRows,5]),(mainTab[listRows,4] + mainTab[listRows,5]))
      
      plot(c(1,144), c(min(mainTab[listRows,1], polyATE, polyATT), max(mainTab[listRows,1], polyATE, polyATT)), type = "n", main = paste(chartTitle, "(Spatial Matching) \n (Spatial Process at the", (disSteps[d]), "x", (disSteps[d]), "Agg Level)"), xlab = "Number of Units in Each Agg Group", ylab = "")
      
      polygon(polyX, polyATE, col = 'grey80', border = NA)
      lines(xAxis, (mainTab[listRows,1] - mainTab[listRows,2]), lty = 'dashed', col = 'red')
      lines(xAxis, (mainTab[listRows,1] + mainTab[listRows,2]), lty = 'dashed', col = 'red')
      
      polygon(polyX, polyATT, col = 'grey80', border = NA)
      lines(xAxis, (mainTab[listRows,4] - mainTab[listRows,5]), lty = 'dashed', col = 'red')
      lines(xAxis, (mainTab[listRows,4] + mainTab[listRows,5]), lty = 'dashed', col = 'red')
      
      lines(xAxis, mainTab[listRows,1], col = "blue", lty = 1, lwd = 2.5)
      lines(xAxis, mainTab[listRows,4], col = "green", lty = 1, lwd = 2.5)
      #legend(0, 3.5, c("ATE", "ATT"), cex = 0.8, col = c("blue", "green"), lty = c(1, 1), lwd = c(2.5, 2.5), bty = "n")
      abline(h = trueCoeff[3], col = "black", lty = 2)
      text(140, (trueCoeff[3] + 0.1), "True T coeff", cex = .6)
      abline(v = (disSteps[d]^2), col = "black", lty = 3)
      dev.off()
      
    }
    
    cat("\n\n** Monte Carlo Completed \n\n")
  
}
