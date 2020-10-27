
probBin = function(m, minEvents = 500, nbin = 11, ...) 
{
  lenVec = length(m)
  if (lenVec <= nbin) 
    stop("length of m needs to greater than nbin")
  
  if (floor(lenVec/nbin) == ceiling(lenVec/nbin)) {
    idx = seq(1, lenVec, by = floor(lenVec/nbin))
  } 
  else {
    idx = c(seq(1, lenVec, by = floor(lenVec/nbin)), lenVec)
  }
  sm = sort(m)
  binEdges = unique(sm[idx])
  histm = hist(m, binEdges)
  ret = list(binEdges = binEdges, counts = histm$counts)
  return(ret)
} 


Analyze2 <- function(WT_NSG, WT_C57, KO_NSG, KO_C57, rawWT_NSG, rawWT_C57, rawKO_NSG, rawKO_C57, timevector){
  
  #Thinning all of the MCMC Analysis chains to 1000 data points beginning from point 1001 to account for burn in
  
  Thin_WT_NSG <- WT_NSG[seq(1001,10000,9),]
  Thin_WT_C57 <- WT_C57[seq(1001,10000,9),]
  
  Thin_KO_NSG <- KO_NSG[seq(1001,10000,9),]
  Thin_KO_C57 <- KO_C57[seq(1001,10000,9),]
  
  #Comparing NSG to C57 in WT
  
  Kp_WT <- na.omit(log10(Thin_WT_C57[,1]/Thin_WT_NSG[,1]))
  
  To_WT <- na.omit(log10(c(Thin_WT_C57[,2], Thin_WT_C57[,3], Thin_WT_C57[,4], Thin_WT_C57[,5], Thin_WT_C57[,6]))-
                log10(c(Thin_WT_NSG[,2], Thin_WT_NSG[,3], Thin_WT_NSG[,4], Thin_WT_NSG[,5], Thin_WT_NSG[,6])))
  
  #Comparing NSG to C57 in KO
  
  Kp_KO <- na.omit(log10(Thin_KO_C57[,1]/Thin_KO_NSG[,1]))
  
  To_KO <- na.omit(log10(c(Thin_KO_C57[,2], Thin_KO_C57[,3], Thin_KO_C57[,4], Thin_KO_C57[,5], Thin_KO_C57[,6]))-
                 log10(c(Thin_KO_NSG[,2], Thin_KO_NSG[,3], Thin_KO_NSG[,4], Thin_KO_NSG[,5], Thin_KO_NSG[,6])))
  
  #Organizing the values of the parameters from smallest to largest to determine 95% enclosing points and average
  
  #Slope
  S_WT_NSG <- sort(Thin_WT_NSG[,1])
  S_WT_C57 <- sort(Thin_WT_C57[,1])
  
  S_KO_NSG <- sort(Thin_KO_NSG[,1])
  S_KO_C57 <- sort(Thin_KO_C57[,1])
  
  #Intercept is based on the intercept values for all 5 individual mice in the cohort
  I_WT_NSG <- sort(c(Thin_WT_NSG[,2], Thin_WT_NSG[,3], Thin_WT_NSG[,4], Thin_WT_NSG[,5], Thin_WT_NSG[,6]))
  I_WT_C57 <- sort(c(Thin_WT_C57[,2], Thin_WT_C57[,3], Thin_WT_C57[,4], Thin_WT_C57[,5], Thin_WT_C57[,6]))
  
  I_KO_NSG <- sort(c(Thin_KO_NSG[,2], Thin_KO_NSG[,3], Thin_KO_NSG[,4], Thin_KO_NSG[,5], Thin_KO_NSG[,6]))
  I_KO_C57 <- sort(c(Thin_KO_C57[,2], Thin_KO_C57[,3], Thin_KO_C57[,4], Thin_KO_C57[,5], Thin_KO_C57[,6]))
  
  #Lower Trendlines are based on the 2.5% data point for the slope and intercept of the dataset (25th point for slopes, 125th point for intercepts)
  
  Low_WT_NSG <- S_WT_NSG[25]*timevector + I_WT_NSG[125]
  Low_WT_C57 <- S_WT_C57[25]*timevector + I_WT_C57[125]
  
  Low_KO_NSG <- S_KO_NSG[25]*timevector + I_KO_NSG[125]
  Low_KO_C57 <- S_KO_C57[25]*timevector + I_KO_C57[125]
  
  #Upper Trendlines are based on the 97.5% data point for the slope and intercept of the dataset (975th point for slopes, 4875th point for intercepts)
  
  High_WT_NSG <- S_WT_NSG[975]*timevector + I_WT_NSG[4875]
  High_WT_C57 <- S_WT_C57[975]*timevector + I_WT_C57[4875]
  
  High_KO_NSG <- S_KO_NSG[975]*timevector + I_KO_NSG[4875]
  High_KO_C57 <- S_KO_C57[975]*timevector + I_KO_C57[4875]
  
  #Average Trendlines are based on the median data point for the slope and intercept of the dataset (500th point for slopes, 2500th point for intercepts)
  
  Avg_WT_NSG <- S_WT_NSG[500]*timevector + I_WT_NSG[2500]
  Avg_WT_C57 <- S_WT_C57[500]*timevector + I_WT_C57[2500]
  
  Avg_KO_NSG <- S_KO_NSG[500]*timevector + I_KO_NSG[2500]
  Avg_KO_C57 <- S_KO_C57[500]*timevector + I_KO_C57[2500]
  
  #Plotting the Results
  oldPar <- par(mfrow = c(1, 4), pty = "s", mar = c(3,4,3,2))

  #NSG Plot
  plot(rawWT_NSG[,1], rawWT_NSG[,2], type = "p", col = "black", pch = 15, cex = 1.2, cex.lab = 1.3, cex.axis = 1.2,
       #main = "NSG", 
       ylab = " Tumor size (10^mm3)", xlab = "Time (Days)",
       xlim = range(timevector), ylim = c(0,4))
  points(rawKO_NSG[,1], rawKO_NSG[,2], col = "red", pch = 19)

  lines(timevector, Avg_WT_NSG, col = "black", lty = "solid")
  lines(timevector, Low_WT_NSG, col = "black", lty = "dashed")
  lines(timevector, High_WT_NSG, col = "black", lty = "dashed")

  lines(timevector, Avg_KO_NSG, col = "red", lty = "solid")
  lines(timevector, Low_KO_NSG, col = "red", lty = "dashed")
  lines(timevector, High_KO_NSG, col = "red", lty = "dashed")

  #C57 Plot
  plot(rawWT_C57[,1], rawWT_C57[,2], type = "p", col = "black", pch = 15, cex = 1.2, cex.lab = 1.3, cex.axis = 1.2,
       #main = "C57", 
       ylab = "Tumor size (10^mm3)", xlab = "Time (Days)",
       xlim = range(timevector), ylim = c(0,4))
  points(rawKO_C57[,1], rawKO_C57[,2], col = "red", pch = 19)

  lines(timevector, Avg_WT_C57, col = "black", lty = "solid")
  lines(timevector, Low_WT_C57, col = "black", lty = "dashed")
  lines(timevector, High_WT_C57, col = "black", lty = "dashed")

  lines(timevector, Avg_KO_C57, col = "red", lty = "solid")
  lines(timevector, Low_KO_C57, col = "red", lty = "dashed")
  lines(timevector, High_KO_C57, col = "red", lty = "dashed")

  #Slope Plot
  h1 <- density(Kp_WT)
  h2 <- density(Kp_KO)
  maxh1h2 <- max(c(h1$y, h2$y))
  h1$y <- h1$y/maxh1h2
  h2$y <- h2$y/maxh1h2
  ymax <- max(c(h1$y, h2$y)) * 1.05
  plot(h1$x, h1$y, type = "l", col = "black", cex = 1.2, cex.lab = 1.3, cex.axis = 1.2,
       xlim = c(-.5,.5),
       ylim = c(0, ymax), #main = "Log10 Ratio of C57 to NSG \nKp in WT and KO",
       xlab = "Value of Kp", ylab = "Normalized Density"
       )
  lines(h2$x, h2$y, col = "red")

  #Initial Tumor Size Plot
  h1 <- density(To_WT)
  h2 <- density(To_KO)
  maxh1h2 <- max(c(h1$y, h2$y))
  h1$y <- h1$y/maxh1h2
  h2$y <- h2$y/maxh1h2
  ymax <- max(c(h1$y, h2$y)) * 1.05
  plot(h1$x, h1$y, type = "l", col = "black", cex = 1.2, cex.lab = 1.3, cex.axis = 1.2,
       xlim = c(-3,3),
       ylim = c(0, ymax), #main = "Log10 Ratio of C57 to NSG \nTo in WT and KO",
       xlab = "Value of To", ylab = "Normalized Density"
       )
  lines(h2$x, h2$y, col = "red")

  par(oldPar)
  
  HWT <- probBin(Kp_WT)
  TbinEdges <- HWT$binEdges
  TbinEdges[length(TbinEdges)] <- max(max(Kp_WT), max(Kp_KO))
  HWT <- hist(Kp_WT, breaks = TbinEdges)
  HKO <- hist(Kp_KO, breaks = TbinEdges)


  CTable <- cbind(HWT$counts, HKO$counts)
  chisq.test(CTable, simulate.p.value = TRUE, B = 100000)

  Xsq <- chisq.test(CTable, simulate.p.value = TRUE, B = 100000)
 #______________________________________
  HWT <- probBin(To_WT)
  TbinEdges <- HWT$binEdges
  TbinEdges[length(TbinEdges)] <- max(max(To_WT), max(To_KO))
  HWT <- hist(To_WT, breaks = TbinEdges)
  HKO <- hist(To_KO, breaks = TbinEdges)


  CTable <- cbind(HWT$counts, HKO$counts)
  chisq.test(CTable, simulate.p.value = TRUE, B = 100000)

  Xsq2 <- chisq.test(CTable, simulate.p.value = TRUE, B = 100000)


  return(data.frame(Xsq$statistic, Xsq$p.value, Xsq2$statistic, Xsq2$p.value))
}

#Creating Figures for all experiments except WISP1

S1_Shift <- Analyze2(T_SWT_NSG1, T_SWT_C571, T_S1C6_NSG1, T_S1C6_C571, SWT_NSG, SWT_C57, S1C6_NSG, S1C6_C57, timevector = 5:26)
S2_Shift <- Analyze2(T_SWT_NSG1, T_SWT_C571, T_S2C2_NSG1, T_S2C2_C571, SWT_NSG, SWT_C57, S2C2_NSG, S2C2_C57, timevector = 5:26)

#NOTE: If you get an error, swap To_WT and To_KO everywhere in 161-165
D_shift <- Analyze2(T_DPWT_NSG1, T_DPWT_C571, T_D_NSG1, T_D_C571, DPWT_NSG, DPWT_C57, D_NSG, D_C57, timevector = 5:26)
#NOTE: Both Kp and To need swapping, so swap To_WT and To_KO for analysis in 149-153
P_Shift <- Analyze2(T_DPWT_NSG1, T_DPWT_C571, T_P_NSG1, T_P_C571, DPWT_NSG, DPWT_C57, P_NSG, P_C57, timevector = 5:26)


#The Analyze 2 function needs minor changes for WISP1 because n = 8 instead of n = 5 like the other experiments

Analyze2 <- function(WT_NSG, WT_C57, KO_NSG, KO_C57, rawWT_NSG, rawWT_C57, rawKO_NSG, rawKO_C57, timevector){
  
  #Thinning all of the MCMC Analysis chains to 1000 data points beginning from point 1001 to account for burn in
  
  Thin_WT_NSG <- WT_NSG[seq(1001,10000,9),]
  Thin_WT_C57 <- WT_C57[seq(1001,10000,9),]
  
  Thin_KO_NSG <- KO_NSG[seq(1001,10000,9),]
  Thin_KO_C57 <- KO_C57[seq(1001,10000,9),]
  
  #Comparing NSG to C57 in WT
  
  Kp_WT <- log10(Thin_WT_C57[,1]/Thin_WT_NSG[,1])
  
  To_WT <- na.omit(sample(log10(c(Thin_WT_C57[,2], Thin_WT_C57[,3], Thin_WT_C57[,4], Thin_WT_C57[,5], Thin_WT_C57[,6],Thin_WT_C57[,7],Thin_WT_C57[,8],Thin_WT_C57[,9])), size = 5000)-
                     log10(c(Thin_WT_NSG[,2], Thin_WT_NSG[,3], Thin_WT_NSG[,4], Thin_WT_NSG[,5], Thin_WT_NSG[,6])))
  
  #Comparing NSG to C57 in KO
  
  Kp_KO <- log10(Thin_KO_C57[,1]/Thin_KO_NSG[,1])
  
  To_KO <- na.omit(sample(log10(c(Thin_KO_C57[,2], Thin_KO_C57[,3], Thin_KO_C57[,4], Thin_KO_C57[,5], Thin_KO_C57[,6],Thin_KO_C57[,7],Thin_KO_C57[,8],Thin_KO_C57[,9])), size = 5000)-
                     log10(c(Thin_KO_NSG[,2], Thin_KO_NSG[,3], Thin_KO_NSG[,4], Thin_KO_NSG[,5], Thin_KO_NSG[,6])))
  
  
  #Organizing the values of the parameters from smallest to largest to determine 95% enclosing points and average
  
  #Slope
  S_WT_NSG <- sort(Thin_WT_NSG[,1])
  S_WT_C57 <- sort(Thin_WT_C57[,1])
  
  S_KO_NSG <- sort(Thin_KO_NSG[,1])
  S_KO_C57 <- sort(Thin_KO_C57[,1])
  
  #Intercept is based on the intercept values for all individual mice in the cohort
  I_WT_NSG <- sort(c(Thin_WT_NSG[,2], Thin_WT_NSG[,3], Thin_WT_NSG[,4], Thin_WT_NSG[,5], Thin_WT_NSG[,6]))
  I_WT_C57 <- sort(c(Thin_WT_C57[,2], Thin_WT_C57[,3], Thin_WT_C57[,4], Thin_WT_C57[,5], Thin_WT_C57[,6],Thin_WT_C57[,7],Thin_WT_C57[,8],Thin_WT_C57[,9]))
  
  I_KO_NSG <- sort(c(Thin_KO_NSG[,2], Thin_KO_NSG[,3], Thin_KO_NSG[,4], Thin_KO_NSG[,5], Thin_KO_NSG[,6]))
  I_KO_C57 <- sort(c(Thin_KO_C57[,2], Thin_KO_C57[,3], Thin_KO_C57[,4], Thin_KO_C57[,5], Thin_KO_C57[,6],Thin_KO_C57[,7],Thin_KO_C57[,8],Thin_KO_C57[,9]))
  
  #Lower Trendlines are based on the 2.5% data point for the slope and intercept of the dataset (25th point for slopes, 125th (NSG) or 200th (C57) point for intercepts)
  
  Low_WT_NSG <- S_WT_NSG[25]*timevector + I_WT_NSG[125]
  Low_WT_C57 <- S_WT_C57[25]*timevector + I_WT_C57[200]
  
  Low_KO_NSG <- S_KO_NSG[25]*timevector + I_KO_NSG[125]
  Low_KO_C57 <- S_KO_C57[25]*timevector + I_KO_C57[200]
  
  #Upper Trendlines are based on the 97.5% data point for the slope and intercept of the dataset (975th point for slopes,4875 (NSG) or 7800th (C57) point for intercepts)
  
  High_WT_NSG <- S_WT_NSG[975]*timevector + I_WT_NSG[4875]
  High_WT_C57 <- S_WT_C57[975]*timevector + I_WT_C57[7800]
  
  High_KO_NSG <- S_KO_NSG[975]*timevector + I_KO_NSG[4875]
  High_KO_C57 <- S_KO_C57[975]*timevector + I_KO_C57[7800]
  
  #Average Trendlines are based on the median data point for the slope and intercept of the dataset (500th point for slopes, 2500th (NSG) or 4000th (C57) point for intercepts)
  
  Avg_WT_NSG <- S_WT_NSG[500]*timevector + I_WT_NSG[2500]
  Avg_WT_C57 <- S_WT_C57[500]*timevector + I_WT_C57[4000]
  
  Avg_KO_NSG <- S_KO_NSG[500]*timevector + I_KO_NSG[2500]
  Avg_KO_C57 <- S_KO_C57[500]*timevector + I_KO_C57[4000]
  
  #Plotting the Results
  oldPar <- par(mfrow = c(1, 4), pty = "s", mar = c(3,4,3,2))
  
  #NSG Plot
  plot(rawWT_NSG[,1], rawWT_NSG[,2], type = "p", col = "black", pch = 15, cex = 1.2, cex.lab = 1.3, cex.axis = 1.2, 
      # main = "NSG", 
       ylab = " Tumor size (10^mm3)", xlab = "Time (Days)", 
       xlim = range(timevector), ylim = c(0,4))
  points(rawKO_NSG[,1], rawKO_NSG[,2], col = "red", pch = 19)
  
  lines(timevector, Avg_WT_NSG, col = "black", lty = "solid")
  lines(timevector, Low_WT_NSG, col = "black", lty = "dashed")
  lines(timevector, High_WT_NSG, col = "black", lty = "dashed")
  
  lines(timevector, Avg_KO_NSG, col = "red", lty = "solid")
  lines(timevector, Low_KO_NSG, col = "red", lty = "dashed")
  lines(timevector, High_KO_NSG, col = "red", lty = "dashed")
  
  #C57 Plot 
  plot(rawWT_C57[,1], rawWT_C57[,2], type = "p", col = "black", pch = 15, cex = 1.2, cex.lab = 1.3, cex.axis = 1.2, 
      # main = "C57", 
       ylab = "Tumor size (10^mm3)", xlab = "Time (Days)",
       xlim = range(timevector), ylim = c(0,4))
  points(rawKO_C57[,1], rawKO_C57[,2], col = "red", pch = 19)
  
  lines(timevector, Avg_WT_C57, col = "black", lty = "solid")
  lines(timevector, Low_WT_C57, col = "black", lty = "dashed")
  lines(timevector, High_WT_C57, col = "black", lty = "dashed")
  
  lines(timevector, Avg_KO_C57, col = "red", lty = "solid")
  lines(timevector, Low_KO_C57, col = "red", lty = "dashed")
  lines(timevector, High_KO_C57, col = "red", lty = "dashed")
  
  #Slope Plot
  h1 <- density(Kp_WT)
  h2 <- density(Kp_KO)
  maxh1h2 <- max(c(h1$y, h2$y))
  h1$y <- h1$y/maxh1h2
  h2$y <- h2$y/maxh1h2
  ymax <- max(c(h1$y, h2$y)) * 1.05
  plot(h1$x, h1$y, type = "l", col = "black", cex = 1.2, cex.lab = 1.3, cex.axis = 1.2,
       xlim = c(-.5,.5), 
       ylim = c(0, ymax), #main = "Natural Log Ratio of C57 to NSG \nKp in WT and KO", 
       xlab = "Value of Kp", ylab = "Density")
  lines(h2$x, h2$y, col = "red")
  
  #Initial Tumor Size Plot
  h1 <- density(To_WT)
  h2 <- density(To_KO)
  maxh1h2 <- max(c(h1$y, h2$y))
  h1$y <- h1$y/maxh1h2
  h2$y <- h2$y/maxh1h2
  ymax <- max(c(h1$y, h2$y)) * 1.05
  plot(h1$x, h1$y, type = "l", col = "black",cex = 1.2, cex.lab = 1.3, cex.axis = 1.2, 
       xlim = c(-3,3),
       ylim = c(0, ymax), #main = "Natural Log Ratio of C57 to NSG \nTo in WT and KO",
       xlab = "Value of To", ylab = "Density")
  lines(h2$x, h2$y, col = "red")
  
  par(oldPar)
   HWT <- probBin(Kp_KO)
  TbinEdges <- HWT$binEdges
  TbinEdges[length(TbinEdges)] <- max(max(Kp_WT), max(Kp_KO))
  HWT <- hist(Kp_KO, breaks = TbinEdges)
  HKO <- hist(Kp_WT, breaks = TbinEdges)


  CTable <- cbind(HWT$counts, HKO$counts)
  chisq.test(CTable, simulate.p.value = TRUE, B = 100000)

  Xsq <- chisq.test(CTable, simulate.p.value = TRUE, B = 100000)

  To_WT <- To_WT[!is.infinite(To_WT)]
  To_KO <- To_KO[!is.infinite(To_KO)]
 # ______________________________________
  HWT <- probBin(To_KO)
  TbinEdges <- HWT$binEdges
  TbinEdges[length(TbinEdges)] <- max(max(To_WT), max(To_KO))
  HWT <- hist(To_KO, breaks = TbinEdges)
  HKO <- hist(To_WT, breaks = TbinEdges)


  CTable <- cbind(HWT$counts, HKO$counts)
  chisq.test(CTable, simulate.p.value = TRUE, B = 100000)

  Xsq2 <- chisq.test(CTable, simulate.p.value = TRUE, B = 100000)


  return(data.frame(Xsq$statistic, Xsq$p.value, Xsq2$statistic, Xsq2$p.value))
}

W_Shift <- Analyze2(T_WWT_NSG1, T_WWT_C571, T_W_NSG1, T_W_C571, WWT_NSG, WWT_C57, W_NSG, W_C57, timevector = 5:26)
