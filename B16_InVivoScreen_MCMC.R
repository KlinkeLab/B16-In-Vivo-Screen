##NOTE: This script is the second of 2 to complete a Markov Chain analysis on tumor volume growth curves.
#This script contains the Markov Chain function and the plotting function used to create the results printouts. 

#Basic Clean-up
graphics.off()  #closes all open figures
cat("\014")     #clears the console

#MARKOV CHAINS
#___________________________________________________________________________________________________________
#___________________________________________________________________________________________________________
require(stats)
library(colorspace)
library(MASS)
library(EBImage)
library(flowViz)
library(coda)

setwd("E:/Experiment Files/B16 In Vivo Screen/MCMC")

#Define Metropolis-Hastings algorithm

MHmcmc <- function(sigma, posterior, data, steps = 1000, target = 0.2, randomSeed = NULL, startValue = NULL) 
{
  if (steps < 100) {
    warning("Function should take at least 100 steps")
  }
  #determine number of parameter dimensions
  np <- length(sigma)
  if (any(sigma <= 0)) 
    stop("All standard deviations must be strictly non-zero and positive")
  
  targetSample <- matrix(rep(0, np*steps), nrow = steps, byrow = TRUE) #empty matrix we'll populate later
  
  if (!is.null(randomSeed)) 
    set.seed(randomSeed)
  z <- rnorm(steps, 0, sigma[1]) #normal distribution centered at 0 with stdev of sigma [1] to create array for 1st parameter
  for (n in 2:np){
    z <- cbind(z, rnorm(steps, 0, sigma[n])) #repeat for all parameters
  }
  u <- runif(steps) #proposal distribution
  if (is.null(startValue)) 
    startValue <- z[1,] #specifying where you're gonna start the chain
  targetSample[1,] <- startValue
  g <- rep(0, steps)
  proposal <- matrix(rep(0, np*steps), nrow = steps, byrow = TRUE)
  alpha <- rep(0, steps)
  naf <- rep(0, steps)
  nz <- rep(0, steps)
  g[1] <- posterior(targetSample[1,], data)
  af <- 1
  sigma1 <- sigma[1]
  i1 <- 1
  nstep = 1
  accept = 1
  
  for (n in 2:steps) {
    proposal[n,] <- targetSample[i1,] + z[n,] #proposing the random step you want to take
    g[n] <- posterior(proposal[n,], data) #calculate posterior
    k3 <- g[n] #proposed point
    k4 <- g[i1] #where you currently are
    alpha[n] <- ifelse(k3/k4 > 1, 1, k3/k4)
    if (u[n] >= alpha[n]) {
      targetSample[n,] <- targetSample[i1,] #reject the step, copy the current point into the next element
    }
    else {
      targetSample[n,] <- proposal[n,] #accept the step, put the new poin into the next element
      i1 <- n
      accept <- accept + 1
    }
    if (nstep >= 200){
      af <- accept/nstep
      if (af > target){
        sigma1 <- sigma1 *1.1
        z <- z * 1.1 #take greater risks. Expand the area you're going to explore
      } else if (af < target){
        sigma1 <- sigma1 * 0.9
        z <- z * 0.9 #take fewer risks. Reduce the area
      }
      nstep = 0
      accept = 0
    } else {
      nstep = nstep + 1
    }
    nz[n] <- sigma1
    naf[n] <- af
  } #end of the Markov Chain Monte Carlo
  
  return(data.frame(MC_param = targetSample, Accept_Fraction = naf, Posterior_Eval = g, Sigma_n = nz))
}

Multiplot <- function(chain1, chain2, chain3){
  
  oldPar <- par(mfrow = c(2, 3), pty = "s")
  
  h1 <- density(chain1$MC_param.1)
  h2 <- density(chain2$MC_param.1)
  h3 <- density(chain3$MC_param.1)
  ymax <- max(c(h1$y, h2$y, h3$y)) * 1.05
  plot(h1$x, h1$y, type = "l", col = "black", xlim = range(chain1$MC_param.1), 
       ylim = c(0, ymax), main = "Posterior Parameter 1", 
       xlab = "x", ylab = "Density")
  lines(h2$x, h2$y, col = "red")
  lines(h3$x, h3$y, col = "blue")
  
  h1 <- density(chain1$MC_param.2)
  h2 <- density(chain2$MC_param.2)
  h3 <- density(chain3$MC_param.2)
  ymax <- max(c(h1$y, h2$y, h3$y)) * 1.05
  plot(h1$x, h1$y, type = "l", col = "black", xlim = range(chain1$MC_param.2), 
       ylim = c(0, ymax), main = "Posterior Parameter 2", 
       xlab = "x", ylab = "Density")
  lines(h2$x, h2$y, col = "red")
  lines(h3$x, h3$y, col = "blue")
  
  h1 <- density(chain1$MC_param.3)
  h2 <- density(chain2$MC_param.3)
  h3 <- density(chain3$MC_param.3)
  ymax <- max(c(h1$y, h2$y, h3$y)) * 1.05
  plot(h1$x, h1$y, type = "l", col = "black", xlim = range(chain1$MC_param.3), 
       ylim = c(0, ymax), main = "Posterior Parameter 3", 
       xlab = "x", ylab = "Density")
  lines(h2$x, h2$y, col = "red")
  lines(h3$x, h3$y, col = "blue")
  
  plot(chain1$MC_param.1, type = "l", col = "black", main = "", ylab = "Target 1 Sample")
  lines(chain2$MC_param.1, col = "red")
  lines(chain3$MC_param.1, col = "blue")
  
  plot(chain1$MC_param.2, type = "l", main = "", ylab = "Target 2 Sample")
  lines(chain2$MC_param.2, col = "red")
  lines(chain3$MC_param.2, col = "blue")
  
  plot(chain1$MC_param.3, type = "l", main = "", ylab = "Target 3 Sample")
  lines(chain2$MC_param.3, col = "red")
  lines(chain3$MC_param.3, col = "blue")
  
  h1 <- density(chain1$MC_param.4)
  h2 <- density(chain2$MC_param.4)
  h3 <- density(chain3$MC_param.4)
  ymax <- max(c(h1$y, h2$y, h3$y)) * 1.05
  plot(h1$x, h1$y, type = "l", col = "black", xlim = range(chain1$MC_param.4), 
       ylim = c(0, ymax), main = "Posterior Parameter 4", 
       xlab = "x", ylab = "Density")
  lines(h2$x, h2$y, col = "red")
  lines(h3$x, h3$y, col = "blue")
  
  h1 <- density(chain1$MC_param.5)
  h2 <- density(chain2$MC_param.5)
  h3 <- density(chain3$MC_param.5)
  ymax <- max(c(h1$y, h2$y, h3$y)) * 1.05
  plot(h1$x, h1$y, type = "l", col = "black", xlim = range(chain1$MC_param.5), 
       ylim = c(0, ymax), main = "Posterior Parameter 5", 
       xlab = "x", ylab = "Density")
  lines(h2$x, h2$y, col = "red")
  lines(h3$x, h3$y, col = "blue")
  
  h1 <- density(chain1$MC_param.6)
  h2 <- density(chain2$MC_param.6)
  h3 <- density(chain3$MC_param.6)
  ymax <- max(c(h1$y, h2$y, h3$y)) * 1.05
  plot(h1$x, h1$y, type = "l", col = "black", xlim = range(chain1$MC_param.6), 
       ylim = c(0, ymax), main = "Posterior Parameter 6", 
       xlab = "x", ylab = "Density")
  lines(h2$x, h2$y, col = "red")
  lines(h3$x, h3$y, col = "blue")

  plot(chain1$MC_param.4, type = "l", main = "", ylab = "Target 4 Sample")
  lines(chain2$MC_param.4, col = "red")
  lines(chain3$MC_param.4, col = "blue")
  
  plot(chain1$MC_param.5, type = "l", main = "", ylab = "Target 5 Sample")
  lines(chain2$MC_param.5, col = "red")
  lines(chain3$MC_param.5, col = "blue")
  
  plot(chain1$MC_param.6, type = "l", main = "", ylab = "Target 6 Sample")
  lines(chain2$MC_param.6, col = "red")
  lines(chain3$MC_param.6, col = "blue")
  
  plot(chain1$Accept_Fraction, type = "p", main = "Acceptance Fraction", xlab = "index", col = "black")
  points(chain2$Accept_Fraction, col = "red")
  points(chain3$Accept_Fraction, col = "blue")
  
  par(oldPar)
}

#Some common likelihood functions that will be used multiple times. Likelihood functions are typically modified for each data
#set based on the number of data points you have available for each mouse. Each Yhat analyzes the data from one mouse.

#Likelihood evaluation for 5 time points
LLH5 <- function(theta, data){
  Yhat1 <- theta[1] * data[1:5,1] + theta[2] 
  SSE1 <- sum((data[1:5,2] - Yhat1)^2) 
  
  Yhat2 <- theta[1] * data[6:10,1] + theta[3] 
  SSE2 <- sum((data[6:10,2] - Yhat2)^2) 
  
  Yhat3 <- theta[1] * data[11:15,1] + theta[4] 
  SSE3 <- sum((data[11:15,2] - Yhat3)^2) 
  
  Yhat4 <- theta[1] * data[16:20,1] + theta[5] 
  SSE4 <- sum((data[16:20,2] - Yhat4)^2)
  
  Yhat5 <- theta[1] * data[21:25,1] + theta[6] 
  SSE5 <- sum((data[21:25,2] - Yhat5)^2) 
  
  post <- (SSE1+SSE2+SSE3+SSE4+SSE5)^(-nrow(data)/2) 
  
  return(post)
}

#Some commonly used Likelihood functions for when no data points are missing in a set.

#Likelihood evaluation for 7 time points
LLH7 <- function(theta, data){
  Yhat1 <- theta[1] * data[1:7,1] + theta[2] 
  SSE1 <- sum((data[1:7,2] - Yhat1)^2) 
  
  Yhat2 <- theta[1] * data[8:14,1] + theta[3] 
  SSE2 <- sum((data[8:14,2] - Yhat2)^2) 
  
  Yhat3 <- theta[1] * data[15:21,1] + theta[4] 
  SSE3 <- sum((data[9:21,2] - Yhat3)^2) 
  
  Yhat4 <- theta[1] * data[22:28,1] + theta[5] 
  SSE4 <- sum((data[22:28,2] - Yhat4)^2)
  
  Yhat5 <- theta[1] * data[29:35,1] + theta[6] 
  SSE5 <- sum((data[29:35,2] - Yhat5)^2) 
  
  post <- (SSE1+SSE2+SSE3+SSE4+SSE5)^(-nrow(data)/2) 
  
  return(post)
}


#Analysis and PDF creation

#_______________________________________________________________
Sig <- c(.1, 1, 1, 1, 1, 1)

#Analysis of NSG AND C57B/6 mice injected with B16F0 WT cells as control group for both DNMT3A and PTPM11 KO groups
T_DPWT_NSG1 <- MHmcmc(Sig, LLH7, DPWT_NSG, steps = 100000, target = 0.2, randomSeed = 1, startValue = c(.3, 3, 3, 3, 3, 3))
T_DPWT_NSG2 <- MHmcmc(Sig, LLH7, DPWT_NSG, steps = 100000, target = 0.2, randomSeed = 2, startValue = c(.2, 2.75, 2.5, 2.75, 2.25, 2.5))
T_DPWT_NSG3 <- MHmcmc(Sig, LLH7, DPWT_NSG, steps = 100000, target = 0.2, randomSeed = 3, startValue = c(.5, 1, 1, 1, 1, 1))

DPWT_NSG_List <- mcmc.list(DPWT_NSG1_MCMC, DPWT_NSG2_MCMC, DPWT_NSG3_MCMC)

#Gelman Rubin Convergence Diagnostic to determine if the three Markov Chains for each dataset are 
#converging to the same point. A multivariate psrf value close to 1 indicates convergence.
#Verify convergence before generating a pdf to save time

#Convergence
DPWT_NSG1_MCMC <- mcmc(data = T_DPWT_NSG1[,1:6])
DPWT_NSG2_MCMC <- mcmc(data = T_DPWT_NSG2[,1:6])
DPWT_NSG3_MCMC <- mcmc(data = T_DPWT_NSG3[,1:6])

DPWT_NSG_List <- mcmc.list(DPWT_NSG1_MCMC, DPWT_NSG2_MCMC, DPWT_NSG3_MCMC)

gelman.diag( DPWT_NSG_List, confidence = 0.95, transform = FALSE, autoburnin = TRUE, multivariate = TRUE)

# pdf("DNMT3A AND PTPN11 WT Control in NSG MCMC.pdf", width = 11, height = 8)
# Multiplot(T_DPWT_NSG1, T_DPWT_NSG2, T_DPWT_NSG3)
# dev.off() # End plot 


#Likelihood evaluation for DNMT3A and PTPN11 Control in C57B/6 Mice
LLH_DPWT_C57 <- function(theta, data){
  Yhat1 <- theta[1] * data[1:3,1] + theta[2] 
  SSE1 <- sum((data[1:3,2] - Yhat1)^2) 
  
  Yhat2 <- theta[1] * data[4:7,1] + theta[3] 
  SSE2 <- sum((data[4:7,2] - Yhat2)^2) 
  
  Yhat3 <- theta[1] * data[8:13,1] + theta[4] 
  SSE3 <- sum((data[8:13,2] - Yhat3)^2) 
  
  Yhat4 <- theta[1] * data[14:21,1] + theta[5] 
  SSE4 <- sum((data[14:21,2] - Yhat4)^2)
  
  Yhat5 <- theta[1] * data[22:25,1] + theta[6] 
  SSE5 <- sum((data[22:25,2] - Yhat5)^2) 
  
  post <- (SSE1+SSE2+SSE3+SSE4+SSE5)^(-nrow(data)/2) 
  
  return(post)
}

T_DPWT_C571 <- MHmcmc(Sig, LLH_DPWT_C57, DPWT_C57, steps = 100000, target = 0.2, randomSeed = 4, startValue = c(.3, 3, 3, 3, 3, 3))
T_DPWT_C572 <- MHmcmc(Sig, LLH_DPWT_C57, DPWT_C57, steps = 100000, target = 0.2, randomSeed = 5, startValue = c(.2, 2.75, 2.5, 2.75, 2.25, 2.5))
T_DPWT_C573 <- MHmcmc(Sig, LLH_DPWT_C57, DPWT_C57, steps = 100000, target = 0.2, randomSeed = 6, startValue = c(.5, 1, 1, 1, 1, 1))

#Convergence
DPWT_C571_MCMC <- mcmc(data = T_DPWT_C571[,1:6])
DPWT_C572_MCMC <- mcmc(data = T_DPWT_C572[,1:6])
DPWT_C573_MCMC <- mcmc(data = T_DPWT_C573[,1:6])

DPWT_C57_List <- mcmc.list(DPWT_C571_MCMC, DPWT_C572_MCMC, DPWT_C573_MCMC)

gelman.diag( DPWT_C57_List, confidence = 0.95, transform = FALSE, autoburnin = TRUE, multivariate = TRUE)

# pdf("DNMT3A AND PTPN11 WT Control in C57 MCMC.pdf", width = 11, height = 8)
# Multiplot(T_DPWT_C571, T_DPWT_C572, T_DPWT_C573)
# dev.off() # End plot 


#_______________________________________________________________
#Analysis of NSG and C57 mice injected with B16F0 DNMT3A KO cells
Sig <- c(.1, 1, 1, 1, 1, 1)

T_D_NSG1 <- MHmcmc(Sig, LLH7, D_NSG, steps = 100000, target = 0.2, randomSeed = 7, startValue = c(.3, 3, 3, 3, 3, 3))
T_D_NSG2 <- MHmcmc(Sig, LLH7, D_NSG, steps = 100000, target = 0.2, randomSeed = 8, startValue = c(.2, 2.75, 2.5, 2.75, 2.25, 2.5))
T_D_NSG3 <- MHmcmc(Sig, LLH7, D_NSG, steps = 100000, target = 0.2, randomSeed = 9, startValue = c(.5, 1, 1, 1, 1, 1))

#Convergence
D_NSG1_MCMC <- mcmc(data = T_D_NSG1[,1:6])
D_NSG2_MCMC <- mcmc(data = T_D_NSG2[,1:6])
D_NSG3_MCMC <- mcmc(data = T_D_NSG3[,1:6])

D_NSG_List <- mcmc.list(D_NSG1_MCMC, D_NSG2_MCMC, D_NSG3_MCMC)

gelman.diag( D_NSG_List, confidence = 0.95, transform = FALSE, autoburnin = TRUE, multivariate = TRUE)

# pdf("DNMT3A KO in NSG MCMC.pdf", width = 11, height = 8)
# Multiplot(T_D_NSG1, T_D_NSG2, T_D_NSG3)
# dev.off() # End plot 

#Likelihood evaluation for DNMT3A KO in C57 Mice
LLH_D_C57 <- function(theta, data){
  Yhat1 <- theta[1] * data[1:6,1] + theta[2] 
  SSE1 <- sum((data[1:6,2] - Yhat1)^2) 
  
  Yhat2 <- theta[1] * data[7:10,1] + theta[3] 
  SSE2 <- sum((data[7:10,2] - Yhat2)^2) 
  
  Yhat3 <- theta[1] * data[11:14,1] + theta[4] 
  SSE3 <- sum((data[11:14,2] - Yhat3)^2) 
  
  Yhat4 <- theta[1] * data[15:20,1] + theta[5] 
  SSE4 <- sum((data[15:20,2] - Yhat4)^2)
  
  Yhat5 <- theta[1] * data[21:24,1] + theta[6] 
  SSE5 <- sum((data[21:24,2] - Yhat5)^2) 
  
  post <- (SSE1+SSE2+SSE3+SSE4+SSE5)^(-nrow(data)/2) 
  
  return(post)
}

T_D_C571 <- MHmcmc(Sig, LLH_D_C57, D_C57, steps = 100000, target = 0.2, randomSeed = 10, startValue = c(.3, 3, 3, 3, 3, 3))
T_D_C572 <- MHmcmc(Sig, LLH_D_C57, D_C57, steps = 100000, target = 0.2, randomSeed = 11, startValue = c(.2, 2.75, 2.5, 2.75, 2.25, 2.5))
T_D_C573 <- MHmcmc(Sig, LLH_D_C57, D_C57, steps = 100000, target = 0.2, randomSeed = 12, startValue = c(.5, 1, 1, 1, 1, 1))

#Convergence
D_C571_MCMC <- mcmc(data = T_D_C571[,1:6])
D_C572_MCMC <- mcmc(data = T_D_C572[,1:6])
D_C573_MCMC <- mcmc(data = T_D_C573[,1:6])

D_C57_List <- mcmc.list(D_C571_MCMC, D_C572_MCMC, D_C573_MCMC)

gelman.diag( D_C57_List, confidence = 0.95, transform = FALSE, autoburnin = TRUE, multivariate = TRUE)

# pdf("DNMT3A KO in C57 MCMC.pdf", width = 11, height = 8)
# Multiplot(T_D_C571, T_D_C572, T_D_C573)
# dev.off() # End plot 


#_______________________________________________________________
#Analysis of NSG and C57 mice injected with B16F0 PTPN11 KO cells  
Sig <- c(.1, 1, 1, 1, 1, 1)

T_P_NSG1 <- MHmcmc(Sig, LLH7, P_NSG, steps = 100000, target = 0.2, randomSeed = 13, startValue = c(.3, 3, 3, 3, 3, 3))
T_P_NSG2 <- MHmcmc(Sig, LLH7, P_NSG, steps = 100000, target = 0.2, randomSeed = 14, startValue = c(.2, 2.75, 2.5, 2.75, 2.25, 2.5))
T_P_NSG3 <- MHmcmc(Sig, LLH7, P_NSG, steps = 100000, target = 0.2, randomSeed = 15, startValue = c(.5, 1, 1, 1, 1, 1))

#Convergence
P_NSG1_MCMC <- mcmc(data = T_P_NSG1[,1:6])
P_NSG2_MCMC <- mcmc(data = T_P_NSG2[,1:6])
P_NSG3_MCMC <- mcmc(data = T_P_NSG3[,1:6])

P_NSG_List <- mcmc.list(P_NSG1_MCMC, P_NSG2_MCMC, P_NSG3_MCMC)

gelman.diag( P_NSG_List, confidence = 0.95, transform = FALSE, autoburnin = TRUE, multivariate = TRUE)

# pdf("PTPN11 KO in NSG MCMC.pdf", width = 11, height = 8)
# Multiplot(T_P_NSG1, T_P_NSG2, T_P_NSG3)
# dev.off() # End plot 

#Likelihood evaluation for PTPN11 KO in C57 Mice
LLH_P_C57 <- function(theta, data){
  Yhat1 <- theta[1] * data[1:7,1] + theta[2] 
  SSE1 <- sum((data[1:7,2] - Yhat1)^2) 
  
  Yhat2 <- theta[1] * data[8:11,1] + theta[3] 
  SSE2 <- sum((data[8:11,2] - Yhat2)^2) 
  
  Yhat3 <- theta[1] * data[12:14,1] + theta[4] 
  SSE3 <- sum((data[12:14,2] - Yhat3)^2) 
  
  Yhat4 <- theta[1] * data[15:21,1] + theta[5] 
  SSE4 <- sum((data[15:21,2] - Yhat4)^2)
  
  Yhat5 <- theta[1] * data[22:26,1] + theta[6] 
  SSE5 <- sum((data[22:26,2] - Yhat5)^2) 
  
  post <- (SSE1+SSE2+SSE3+SSE4+SSE5)^(-nrow(data)/2) 
  
  return(post)
}

Sig <- c(.1, .5, 1, 1, 1, 1)

T_P_C571 <- MHmcmc(Sig, LLH_P_C57, P_C57, steps = 100000, target = 0.2, randomSeed = 16, startValue = c(.4, 2, 2, 2, 2, 2))
T_P_C572 <- MHmcmc(Sig, LLH_P_C57, P_C57, steps = 100000, target = 0.2, randomSeed = 17, startValue = c(.3, 1, 1, 1, 1, 1))
T_P_C573 <- MHmcmc(Sig, LLH_P_C57, P_C57, steps = 100000, target = 0.2, randomSeed = 18, startValue = c(.6, 0, 0, 0, 0, 0))

#Convergence
P_C571_MCMC <- mcmc(data = T_P_C571[,1:6])
P_C572_MCMC <- mcmc(data = T_P_C572[,1:6])
P_C573_MCMC <- mcmc(data = T_P_C573[,1:6])

P_C57_List <- mcmc.list(P_C571_MCMC, P_C572_MCMC, P_C573_MCMC)

gelman.diag( P_C57_List, confidence = 0.95, transform = FALSE, autoburnin = TRUE, multivariate = TRUE)

# pdf("PTPN11 KO in C57 MCMC.pdf", width = 11, height = 8)
# Multiplot(T_P_C571, T_P_C572, T_P_C573)
# dev.off() # End plot 

#_______________________________________________________________
#Analysis of NSG and C57B/6 mice injected with B16F0 WT cells as control group for both the SPARC KO group
Sig <- c(.1, .5, 1, 1, 1, 1)

T_SWT_NSG1 <- MHmcmc(Sig, LLH5, SWT_NSG, steps = 100000, target = 0.2, randomSeed = 19, startValue = c(.4, 2, 2, 2, 2, 2))
T_SWT_NSG2 <- MHmcmc(Sig, LLH5, SWT_NSG, steps = 100000, target = 0.2, randomSeed = 20, startValue = c(.3, 1, 1, 1, 1, 1))
T_SWT_NSG3 <- MHmcmc(Sig, LLH5, SWT_NSG, steps = 100000, target = 0.2, randomSeed = 21, startValue = c(.6, 0, 0, 0, 0, 0))

#Convergence
SWT_NSG1_MCMC <- mcmc(data = T_SWT_NSG1[,1:6])
SWT_NSG2_MCMC <- mcmc(data = T_SWT_NSG2[,1:6])
SWT_NSG3_MCMC <- mcmc(data = T_SWT_NSG3[,1:6])

SWT_NSG_List <- mcmc.list(SWT_NSG1_MCMC, SWT_NSG2_MCMC, SWT_NSG3_MCMC)

gelman.diag( SWT_NSG_List, confidence = 0.95, transform = FALSE, autoburnin = TRUE, multivariate = TRUE)

# pdf("SPARC WT Control in NSG MCMC.pdf", width = 11, height = 8)
# Multiplot(T_SWT_NSG1, T_SWT_NSG2, T_SWT_NSG3)
# dev.off() # End plot 

#Likelihood evaluation for SPARC Control in C57B/6 Mice
LLH_SWT_C57 <- function(theta, data){
  Yhat1 <- theta[1] * data[1:8,1] + theta[2] 
  SSE1 <- sum((data[1:8,2] - Yhat1)^2) 
  
  Yhat2 <- theta[1] * data[9:15,1] + theta[3] 
  SSE2 <- sum((data[9:15,2] - Yhat2)^2) 
  
  Yhat3 <- theta[1] * data[16:23,1] + theta[4] 
  SSE3 <- sum((data[16:23,2] - Yhat3)^2) 
  
  Yhat4 <- theta[1] * data[24:30,1] + theta[5] 
  SSE4 <- sum((data[24:30,2] - Yhat4)^2)
  
  Yhat5 <- theta[1] * data[31:38,1] + theta[6] 
  SSE5 <- sum((data[31:38,2] - Yhat5)^2) 
  
  post <- (SSE1+SSE2+SSE3+SSE4+SSE5)^(-nrow(data)/2) 
  
  return(post)
}

Sig <- c(.1, .5, 1, 1, 1, 1)

T_SWT_C571 <- MHmcmc(Sig, LLH_SWT_C57, SWT_C57, steps = 100000, target = 0.2, randomSeed = 22, startValue = c(.3, 3, 3, 3, 3, 3))
T_SWT_C572 <- MHmcmc(Sig, LLH_SWT_C57, SWT_C57, steps = 100000, target = 0.2, randomSeed = 23, startValue = c(.2, 2.75, 2.5, 2.75, 2.25, 2.5))
T_SWT_C573 <- MHmcmc(Sig, LLH_SWT_C57, SWT_C57, steps = 100000, target = 0.2, randomSeed = 24, startValue = c(.5, 1, 1, 1, 1, 1))

#Convergence
SWT_C571_MCMC <- mcmc(data = T_SWT_C571[,1:6])
SWT_C572_MCMC <- mcmc(data = T_SWT_C572[,1:6])
SWT_C573_MCMC <- mcmc(data = T_SWT_C573[,1:6])

SWT_C57_List <- mcmc.list(SWT_C571_MCMC, SWT_C572_MCMC, SWT_C573_MCMC)

gelman.diag( SWT_C57_List, confidence = 0.95, transform = FALSE, autoburnin = TRUE, multivariate = TRUE)

# pdf("SPARC WT Control in C57 MCMC.pdf", width = 11, height = 8)
# Multiplot(T_SWT_C571, T_SWT_C572, T_SWT_C573)
# dev.off() # End plot 

#_______________________________________________________________
#Analysis of NSG and C57 mice injected with B16F0 SPARC 1C6 KO cells 

#Likelihood evaluation for SPARC 1C6 KO IN NSG Mice
LLH_S1C6_NSG <- function(theta, data){
  Yhat1 <- theta[1] * data[1:5,1] + theta[2] 
  SSE1 <- sum((data[1:5,2] - Yhat1)^2) 
  
  Yhat2 <- theta[1] * data[6:10,1] + theta[3] 
  SSE2 <- sum((data[6:10,2] - Yhat2)^2) 
  
  Yhat3 <- theta[1] * data[11:15,1] + theta[4] 
  SSE3 <- sum((data[11:15,2] - Yhat3)^2) 
  
  Yhat4 <- theta[1] * data[16:19,1] + theta[5] 
  SSE4 <- sum((data[16:19,2] - Yhat4)^2)
  
  Yhat5 <- theta[1] * data[20:24,1] + theta[6] 
  SSE5 <- sum((data[20:24,2] - Yhat5)^2) 
  
  post <- (SSE1+SSE2+SSE3+SSE4+SSE5)^(-nrow(data)/2) 
  
  return(post)
}

Sig <- c(.1, 1, 1, 1, 1, 1)

T_S1C6_NSG1 <- MHmcmc(Sig, LLH_S1C6_NSG, S1C6_NSG, steps = 100000, target = 0.2, randomSeed = 25, startValue = c(.4, 0, 0, 0, 0, 0))
T_S1C6_NSG2 <- MHmcmc(Sig, LLH_S1C6_NSG, S1C6_NSG, steps = 100000, target = 0.2, randomSeed = 26, startValue = c(.6, 2, 2, 2, 2, 2))
T_S1C6_NSG3 <- MHmcmc(Sig, LLH_S1C6_NSG, S1C6_NSG, steps = 100000, target = 0.2, randomSeed = 27, startValue = c(.5, 1, 1, 1, 1, 1))

#Convergence
S1C6_NSG1_MCMC <- mcmc(data = T_S1C6_NSG1[,1:6])
S1C6_NSG2_MCMC <- mcmc(data = T_S1C6_NSG2[,1:6])
S1C6_NSG3_MCMC <- mcmc(data = T_S1C6_NSG3[,1:6])

S1C6_NSG_List <- mcmc.list(S1C6_NSG1_MCMC, S1C6_NSG2_MCMC, S1C6_NSG3_MCMC)

gelman.diag( S1C6_NSG_List, confidence = 0.95, transform = FALSE, autoburnin = TRUE, multivariate = TRUE)

# pdf("SPARC 1C6 KO in NSG MCMC.pdf", width = 11, height = 8)
# Multiplot(T_S1C6_NSG1, T_S1C6_NSG2, T_S1C6_NSG3)
# dev.off() # End plot 

#Likelihood function for SPARC 1C6 KO in C57B/6 Mice
LLH_S1C6_C57 <- function(theta, data){
  Yhat1 <- theta[1] * data[1:8,1] + theta[2] 
  SSE1 <- sum((data[1:8,2] - Yhat1)^2) 
  
  Yhat2 <- theta[1] * data[9:14,1] + theta[3] 
  SSE2 <- sum((data[9:14,2] - Yhat2)^2) 
  
  Yhat3 <- theta[1] * data[15:22,1] + theta[4] 
  SSE3 <- sum((data[15:22,2] - Yhat3)^2) 
  
  Yhat4 <- theta[1] * data[23:30,1] + theta[5] 
  SSE4 <- sum((data[23:30,2] - Yhat4)^2)
  
  Yhat5 <- theta[1] * data[31:38,1] + theta[6] 
  SSE5 <- sum((data[31:38,2] - Yhat5)^2) 
  
  post <- (SSE1+SSE2+SSE3+SSE4+SSE5)^(-nrow(data)/2) 
  
  return(post)
}
Sig <- c(.1, 1, 1, 1, 1, 1)

T_S1C6_C571 <- MHmcmc(Sig, LLH_S1C6_C57, S1C6_C57, steps = 100000, target = 0.2, randomSeed = 28, startValue = c(.3, 3, 3, 3, 3, 3))
T_S1C6_C572 <- MHmcmc(Sig, LLH_S1C6_C57, S1C6_C57, steps = 100000, target = 0.2, randomSeed = 29, startValue = c(.2, 2.75, 2.5, 2.75, 2.25, 2.5))
T_S1C6_C573 <- MHmcmc(Sig, LLH_S1C6_C57, S1C6_C57, steps = 100000, target = 0.2, randomSeed = 30, startValue = c(.5, 1, 1, 1, 1, 1))

#Convergence
S1C6_C571_MCMC <- mcmc(data = T_S1C6_C571[,1:6])
S1C6_C572_MCMC <- mcmc(data = T_S1C6_C572[,1:6])
S1C6_C573_MCMC <- mcmc(data = T_S1C6_C573[,1:6])

S1C6_C57_List <- mcmc.list(S1C6_C571_MCMC, S1C6_C572_MCMC, S1C6_C573_MCMC)

gelman.diag( S1C6_C57_List, confidence = 0.95, transform = FALSE, autoburnin = TRUE, multivariate = TRUE)

# pdf("SPARC 1C6 KO in C57 MCMC.pdf", width = 11, height = 8)
# Multiplot(T_S1C6_C571, T_S1C6_C572, T_S1C6_C573)
# dev.off() # End plot 


#_______________________________________________________________
#Analysis of NSG and C57 mice injected with B16F0 SPARC 2C2 KO cells 

#Likelihood evaluation for SPARC 2C2 KO IN NSG Mice
LLH_S2C2_NSG <- function(theta, data){
  Yhat1 <- theta[1] * data[1:5,1] + theta[2] 
  SSE1 <- sum((data[1:5,2] - Yhat1)^2) 
  
  Yhat2 <- theta[1] * data[6:9,1] + theta[3] 
  SSE2 <- sum((data[6:9,2] - Yhat2)^2) 
  
  Yhat3 <- theta[1] * data[10:14,1] + theta[4] 
  SSE3 <- sum((data[10:14,2] - Yhat3)^2) 
  
  Yhat4 <- theta[1] * data[15:19,1] + theta[5] 
  SSE4 <- sum((data[15:19,2] - Yhat4)^2)
  
  Yhat5 <- theta[1] * data[20:24,1] + theta[6] 
  SSE5 <- sum((data[20:24,2] - Yhat5)^2) 
  
  post <- (SSE1+SSE2+SSE3+SSE4+SSE5)^(-nrow(data)/2) 
  
  return(post)
}

Sig <- c(.1, 1, 1, 1, 1, 1)

T_S2C2_NSG1 <- MHmcmc(Sig, LLH_S2C2_NSG, S2C2_NSG, steps = 100000, target = 0.2, randomSeed = 31, startValue = c(.3, 3, 3, 3, 3, 3))
T_S2C2_NSG2 <- MHmcmc(Sig, LLH_S2C2_NSG, S2C2_NSG, steps = 100000, target = 0.2, randomSeed = 32, startValue = c(.2, 2.75, 2.5, 2.75, 2.25, 2.5))
T_S2C2_NSG3 <- MHmcmc(Sig, LLH_S2C2_NSG, S2C2_NSG, steps = 100000, target = 0.2, randomSeed = 33, startValue = c(.5, 1, 1, 1, 1, 1))

S2C2_NSG1_MCMC <- mcmc(data = T_S2C2_NSG1[,1:6])
S2C2_NSG2_MCMC <- mcmc(data = T_S2C2_NSG2[,1:6])
S2C2_NSG3_MCMC <- mcmc(data = T_S2C2_NSG3[,1:6])

S2C2_NSG_List <- mcmc.list(S2C2_NSG1_MCMC, S2C2_NSG2_MCMC, S2C2_NSG3_MCMC)

gelman.diag( S2C2_NSG_List, confidence = 0.95, transform = FALSE, autoburnin = TRUE, multivariate = TRUE)

# pdf("SPARC 2C2 KO in NSG MCMC.pdf", width = 11, height = 8)
# Multiplot(T_S2C2_NSG1, T_S2C2_NSG2, T_S2C2_NSG3)
# dev.off() # End plot 

#Likelihood function for SPARC 2C2 KO in C57B/6 Mice
LLH_S2C2_C57 <- function(theta, data){
  Yhat1 <- theta[1] * data[1:7,1] + theta[2] 
  SSE1 <- sum((data[1:7,2] - Yhat1)^2) 
  
  Yhat2 <- theta[1] * data[8:14,1] + theta[3] 
  SSE2 <- sum((data[8:14,2] - Yhat2)^2) 
  
  Yhat3 <- theta[1] * data[15:21,1] + theta[4] 
  SSE3 <- sum((data[15:21,2] - Yhat3)^2) 
  
  Yhat4 <- theta[1] * data[22:26,1] + theta[5] 
  SSE4 <- sum((data[22:26,2] - Yhat4)^2)
  
  Yhat5 <- theta[1] * data[27:31,1] + theta[6] 
  SSE5 <- sum((data[27:31,2] - Yhat5)^2) 
  
  post <- (SSE1+SSE2+SSE3+SSE4+SSE5)^(-nrow(data)/2) 
  
  return(post)
}

Sig <- c(.1, 1, 1, 1, 1, 1)

T_S2C2_C571 <- MHmcmc(Sig, LLH_S2C2_C57, S2C2_C57, steps = 100000, target = 0.2, randomSeed = 34, startValue = c(.3, 0, 0, 0, 0, 0))
T_S2C2_C572 <- MHmcmc(Sig, LLH_S2C2_C57, S2C2_C57, steps = 100000, target = 0.2, randomSeed = 35, startValue = c(.4, 2, 2, 2, 2, 2))
T_S2C2_C573 <- MHmcmc(Sig, LLH_S2C2_C57, S2C2_C57, steps = 100000, target = 0.2, randomSeed = 36, startValue = c(.5, 1, 1, 1, 1, 1))

#Convergence
S2C2_C571_MCMC <- mcmc(data = T_S2C2_C571[,1:6])
S2C2_C572_MCMC <- mcmc(data = T_S2C2_C572[,1:6])
S2C2_C573_MCMC <- mcmc(data = T_S2C2_C573[,1:6])

S2C2_C57_List <- mcmc.list(S2C2_C571_MCMC, S2C2_C572_MCMC, S2C2_C573_MCMC)

gelman.diag( S2C2_C57_List, confidence = 0.95, transform = FALSE, autoburnin = TRUE, multivariate = TRUE)

# pdf("SPARC 2C2 KO in C57 MCMC.pdf", width = 11, height = 8)
# Multiplot(T_S2C2_C571, T_S2C2_C572, T_S2C2_C573)
# dev.off() # End plot 

#_______________________________________________________________
#Analysis of NSG and C57B/6 mice injected with B16F0 WT cells as the control group for the WISP1 KO group
Sig <- c(.1, 1, 1, 1, 1, 1)

T_WWT_NSG1 <- MHmcmc(Sig, LLH7, WWT_NSG, steps = 100000, target = 0.2, randomSeed = 37, startValue = c(.3, 3, 3, 3, 3, 3))
T_WWT_NSG2 <- MHmcmc(Sig, LLH7, WWT_NSG, steps = 100000, target = 0.2, randomSeed = 38, startValue = c(.2, 2.75, 2.5, 2.75, 2.25, 2.5))
T_WWT_NSG3 <- MHmcmc(Sig, LLH7, WWT_NSG, steps = 100000, target = 0.2, randomSeed = 39, startValue = c(.5, 1, 1, 1, 1, 1))

#Convergence
WWT_NSG1_MCMC <- mcmc(data = T_WWT_NSG1[,1:6])
WWT_NSG2_MCMC <- mcmc(data = T_WWT_NSG2[,1:6])
WWT_NSG3_MCMC <- mcmc(data = T_WWT_NSG3[,1:6])

WWT_NSG_List <- mcmc.list(WWT_NSG1_MCMC, WWT_NSG2_MCMC, WWT_NSG3_MCMC)

gelman.diag( WWT_NSG_List, confidence = 0.95, transform = FALSE, autoburnin = TRUE, multivariate = TRUE)

# pdf("WISP1 WT Control in NSG MCMC.pdf", width = 11, height = 8)
# Multiplot(T_WWT_NSG1, T_WWT_NSG2, T_WWT_NSG3)
# dev.off() # End plot 

#Likelihood evaluation for WISP1 WT Control in C57 Mice
LLH_WWT_C57 <- function(theta, data){
  Yhat1 <- theta[1] * data[1:4,1] + theta[2] 
  SSE1 <- sum((data[1:4,2] - Yhat1)^2) 
  
  Yhat2 <- theta[1] * data[5:7,1] + theta[3] 
  SSE2 <- sum((data[5:7,2] - Yhat2)^2) 
  
  Yhat3 <- theta[1] * data[8:15,1] + theta[4] 
  SSE3 <- sum((data[8:15,2] - Yhat3)^2) 
  
  Yhat4 <- theta[1] * data[16:18,1] + theta[5] 
  SSE4 <- sum((data[16:18,2] - Yhat4)^2)
  
  Yhat5 <- theta[1] * data[19:26,1] + theta[6] 
  SSE5 <- sum((data[19:26,2] - Yhat5)^2) 
  
  post <- (SSE1+SSE2+SSE3+SSE4+SSE5)^(-nrow(data)/2) 
  
  return(post)
}

Sig <- c(.1, 1, 1, 1, 1, 1)

T_WWT_C571 <- MHmcmc(Sig, LLH_WWT_C57, WWT_C57, steps = 100000, target = 0.2, randomSeed = 40, startValue = c(.3, 3, 3, 3, 3, 3))
T_WWT_C572 <- MHmcmc(Sig, LLH_WWT_C57, WWT_C57, steps = 100000, target = 0.2, randomSeed = 41, startValue = c(.2, 2.75, 2.5, 2.75, 2.25, 2.5))
T_WWT_C573 <- MHmcmc(Sig, LLH_WWT_C57, WWT_C57, steps = 100000, target = 0.2, randomSeed = 42, startValue = c(0, 4, 4, 4, 4, 4))

#Convergence
WWT_C571_MCMC <- mcmc(data = T_WWT_C571[,1:6])
WWT_C572_MCMC <- mcmc(data = T_WWT_C572[,1:6])
WWT_C573_MCMC <- mcmc(data = T_WWT_C573[,1:6])

WWT_C57_List <- mcmc.list(WWT_C571_MCMC, WWT_C572_MCMC, WWT_C573_MCMC)

gelman.diag( WWT_C57_List, confidence = 0.95, transform = FALSE, autoburnin = TRUE, multivariate = TRUE)

# pdf("WISP1 WT Control in C57 MCMC.pdf", width = 11, height = 8)
# Multiplot(T_WWT_C571, T_WWT_C572, T_WWT_C573)
# dev.off() # End plot 

#_______________________________________________________________
#Analysis of NSG and C57B/6 mice injected with B16F0 WISP1 KO cells

T_W_NSG1 <- MHmcmc(Sig, LLH7, W_NSG, steps = 100000, target = 0.2, randomSeed = 43, startValue = c(.3, 3, 3, 3, 3, 3))
T_W_NSG2 <- MHmcmc(Sig, LLH7, W_NSG, steps = 100000, target = 0.2, randomSeed = 44, startValue = c(.2, 2.75, 2.5, 2.75, 2.25, 2.5))
T_W_NSG3 <- MHmcmc(Sig, LLH7, W_NSG, steps = 100000, target = 0.2, randomSeed = 45, startValue = c(.5, 1, 1, 1, 1, 1))

#Convergence
W_NSG1_MCMC <- mcmc(data = T_W_NSG1[,1:6])
W_NSG2_MCMC <- mcmc(data = T_W_NSG2[,1:6])
W_NSG3_MCMC <- mcmc(data = T_W_NSG3[,1:6])

W_NSG_List <- mcmc.list(W_NSG1_MCMC, W_NSG2_MCMC, W_NSG3_MCMC)

gelman.diag( W_NSG_List, confidence = 0.95, transform = FALSE, autoburnin = TRUE, multivariate = TRUE)

# pdf("WISP1 KO in NSG MCMC.pdf", width = 11, height = 8)
# Multiplot(T_W_NSG1, T_W_NSG2, T_W_NSG3)
# dev.off() # End plot 

#Likelihood evaluation for WISP1 WT Control in C57 Mice
LLH_W_C57 <- function(theta, data){
  Yhat1 <- theta[1] * data[1:7,1] + theta[2] 
  SSE1 <- sum((data[1:7,2] - Yhat1)^2) 
  
  Yhat2 <- theta[1] * data[8:19,1] + theta[3] 
  SSE2 <- sum((data[8:19,2] - Yhat2)^2) 
  
  Yhat3 <- theta[1] * data[20:27,1] + theta[4] 
  SSE3 <- sum((data[20:27,2] - Yhat3)^2) 
  
  Yhat4 <- theta[1] * data[28:31,1] + theta[5] 
  SSE4 <- sum((data[28:31,2] - Yhat4)^2)
  
  Yhat5 <- theta[1] * data[32:43,1] + theta[6] 
  SSE5 <- sum((data[32:43,2] - Yhat5)^2) 
  
  post <- (SSE1+SSE2+SSE3+SSE4+SSE5)^(-nrow(data)/2) 
  
  return(post)
}

T_W_C571 <- MHmcmc(Sig, LLH_W_C57, W_C57, steps = 100000, target = 0.2, randomSeed = 46, startValue = c(.3, 3, 3, 3, 3, 3))
T_W_C572 <- MHmcmc(Sig, LLH_W_C57, W_C57, steps = 100000, target = 0.2, randomSeed = 47, startValue = c(.2, 2.75, 2.5, 2.75, 2.25, 2.5))
T_W_C573 <- MHmcmc(Sig, LLH_W_C57, W_C57, steps = 100000, target = 0.2, randomSeed = 48, startValue = c(.5, 1, 1, 1, 1, 1))

#Convergence
W_C571_MCMC <- mcmc(data = T_W_C571[,1:6])
W_C572_MCMC <- mcmc(data = T_W_C572[,1:6])
W_C573_MCMC <- mcmc(data = T_W_C573[,1:6])

W_C57_List <- mcmc.list(W_C571_MCMC, W_C572_MCMC, W_C573_MCMC)

gelman.diag( W_C57_List, confidence = 0.95, transform = FALSE, autoburnin = TRUE, multivariate = TRUE)

# pdf("WISP1 KO in C57 MCMC.pdf", width = 11, height = 8)
# Multiplot(T_W_C571, T_W_C572, T_W_C573)
# dev.off() # End plot 