

##NOTE: This script is designed to further analyze the results of the Markov Chain Monte Carlo 
#analysis presented in the B16_InVivoScreen_MCMC.R file in combination with the raw tumor growth data. 
#The function Analyze ultimately creates a figure with 4 plots representing the raw tumor growth data in NSG (plot 1) 
#and C57 mice (plot 2) with an average trendline (solid line) and 95% confidence interval (dashed lines) along with density plots for the tumor
#growth rate Kp (plot 3) and the initial tumor injection size To (plot 4) for both WT and KO cell lines (black and red respectively)
#in both NSG and C57 (solid and dashed lines respectively)

Analyze <- function(WT_NSG, WT_C57, KO_NSG, KO_C57, rawWT_NSG, rawWT_C57, rawKO_NSG, rawKO_C57, timevector){

  #Thinning all of the MCMC Analysis chains to 1000 data points beginning from point 1001 to account for burn in

  Thin_WT_NSG <- WT_NSG[seq(1001,10000,9),]
  Thin_WT_C57 <- WT_C57[seq(1001,10000,9),]
  
  Thin_KO_NSG <- KO_NSG[seq(1001,10000,9),]
  Thin_KO_C57 <- KO_C57[seq(1001,10000,9),]

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
  oldPar <- par(mfrow = c(1, 4), pty = "s")
  
  #NSG Plot
  plot(rawWT_NSG[,1], rawWT_NSG[,2], type = "p", col = "black", pch = 15, 
       main = "NSG", ylab = " Log Tumor size (mm3)", xlab = "Days Post Tumor Injection", 
       xlim = range(timevector), ylim = c(0,max(High_WT_NSG, High_KO_NSG)))
  points(rawKO_NSG[,1], rawKO_NSG[,2], col = "red", pch = 19)
  
  lines(timevector, Avg_WT_NSG, col = "black", lty = "solid")
  lines(timevector, Low_WT_NSG, col = "black", lty = "dashed")
  lines(timevector, High_WT_NSG, col = "black", lty = "dashed")
  
  lines(timevector, Avg_KO_NSG, col = "red", lty = "solid")
  lines(timevector, Low_KO_NSG, col = "red", lty = "dashed")
  lines(timevector, High_KO_NSG, col = "red", lty = "dashed")
  
  #C57 Plot 
  plot(rawWT_C57[,1], rawWT_C57[,2], type = "p", col = "black", pch = 15, 
       main = "C57", ylab = "Log Tumor size (mm3)", xlab = "Days Post Tumor Injection",
       xlim = range(timevector), ylim = c(0,max(High_WT_NSG, High_KO_NSG)))
  points(rawKO_C57[,1], rawKO_C57[,2], col = "red", pch = 19)
  
  lines(timevector, Avg_WT_C57, col = "black", lty = "solid")
  lines(timevector, Low_WT_C57, col = "black", lty = "dashed")
  lines(timevector, High_WT_C57, col = "black", lty = "dashed")
  
  lines(timevector, Avg_KO_C57, col = "red", lty = "solid")
  lines(timevector, Low_KO_C57, col = "red", lty = "dashed")
  lines(timevector, High_KO_C57, col = "red", lty = "dashed")
    
 #Slope Plot
  h1 <- density(S_WT_NSG)
  h2 <- density(S_WT_C57)
  h3 <- density(S_KO_NSG)
  h4 <- density(S_KO_C57)
  ymax <- max(c(h1$y, h2$y, h3$y, h4$y )) * 1.05
  plot(h1$x, h1$y, type = "l", col = "black", 
       xlim = c(min(S_WT_C57, S_KO_C57, S_WT_NSG, S_KO_NSG),max(S_WT_C57, S_KO_C57, S_WT_NSG, S_KO_NSG)), 
       ylim = c(0, ymax), main = "Kp", 
       xlab = "Value of Kp", ylab = "Normalized Density")
  lines(h2$x, h2$y, col = "black", lty = "dashed")
  lines(h3$x, h3$y, col = "red")
  lines(h4$x, h4$y, col = "red", lty = "dashed")
  
  #Intercept Plot
  h1 <- density(I_WT_NSG)
  h2 <- density(I_WT_C57)
  h3 <- density(I_KO_NSG)
  h4 <- density(I_KO_C57)
  ymax <- max(c(h1$y, h2$y, h3$y, h4$y )) * 1.05
  plot(h1$x, h1$y, type = "l", col = "black",
       xlim = c(min(I_WT_C57, I_KO_C57, I_WT_NSG, I_KO_NSG), max(I_WT_C57, I_KO_C57, I_WT_NSG, I_KO_NSG)), 
       ylim = c(0, ymax), main = "To", 
       xlab = "Value of To", ylab = "Normalized Density")
  lines(h2$x, h2$y, col = "black", lty = "dashed")
  lines(h3$x, h3$y, col = "red")
  lines(h4$x, h4$y, col = "red", lty = "dashed")
  
  par(oldPar)
  return(data.frame(S_WT_NSG, S_WT_C57, S_KO_NSG, S_KO_C57))
}

#Creating Figures for all experiments
DNMT3A_Kp <- Analyze(T_DPWT_NSG1, T_DPWT_C571, T_D_NSG1, T_D_C571, DPWT_NSG, DPWT_C57, D_NSG, D_C57, timevector = 0:24)
PTPN11_Kp <- Analyze(T_DPWT_NSG1, T_DPWT_C571, T_P_NSG1, T_P_C571, DPWT_NSG, DPWT_C57, P_NSG, P_C57, timevector = 0:26)
SPARC1C6_Kp <- Analyze(T_SWT_NSG1, T_SWT_C571, T_S1C6_NSG1, T_S1C6_C571, SWT_NSG, SWT_C57, S1C6_NSG, S1C6_C57, timevector = 0:24)
SPARC2C2_Kp <- Analyze(T_SWT_NSG1, T_SWT_C571, T_S2C2_NSG1, T_S2C2_C571, SWT_NSG, SWT_C57, S2C2_NSG, S2C2_C57, timevector = 0:24)
WISP1_Kp <- Analyze(T_WWT_NSG1, T_WWT_C571, T_W_NSG1, T_W_C571, WWT_NSG, WWT_C57, W_NSG, W_C57, timevector = 0:26)

#Natural log of average C57 growth rate divided by average NSG growth rate. 
#Positive indicates C57 grows faster. Negative indicates NSG grows faster

#WT DNMT3A and PTPN11 (small difference)
log(DNMT3A_Kp[500,2]/DNMT3A_Kp[500,1])
#KO DNMT3A (C57 grows faster)
log(DNMT3A_Kp[500,4]/DNMT3A_Kp[500,3])
#KO PTPN11 (C57 grows faster)
log(PTPN11_Kp[500,4]/PTPN11_Kp[500,3])


#WT SPARC 1C6 (NSG grows faster)
log(SPARC1C6_Kp[500,2]/SPARC1C6_Kp[500,1])
#KO SPARC 1C6 (NSG grows faster)
log(SPARC1C6_Kp[500,4]/SPARC1C6_Kp[500,3])


#WT SPARC 2C2 (NSG grows faster)
log(SPARC2C2_Kp[500,2]/SPARC2C2_Kp[500,1])
#KO SPARC 2C2 (NSG grows faster)
log(SPARC2C2_Kp[500,4]/SPARC2C2_Kp[500,3])


#WT WISP1 (NSG grows faster)
log(WISP1_Kp[500,2]/WISP1_Kp[500,1])
#KO WISP1 (NSG grows faster)
log(WISP1_Kp[500,4]/WISP1_Kp[500,3])