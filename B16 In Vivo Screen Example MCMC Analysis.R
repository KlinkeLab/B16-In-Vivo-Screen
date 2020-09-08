#B16 In Vivo Screen MCMC Example Code
#Klinke Lab
#9/8/2020

#Basic Clean-up
graphics.off()  #closes all open figures
rm(list=ls())   #clears the variables
cat("\014")     #clears the console

# Defining Raw Growth Data

#First, the data must be organized into a 2 column matrix for input into the Markov Chain analysis.
#The tumor volume data from each mouse arranged into matrix with column 1 indicating the day after tumor challenge
#and column 2 indicating the tumor volume on that day. 
#After a matrix is created for each mouse, these are bound together using the abind function to create a single
#large 2 column matrix with the number of rows equal to the total number of data points collected.

require(stats)
library(colorspace)
library(MASS)
library(EBImage)
library(flowViz)
library(coda)

#Wild Type data for comparison with DNMT3A and PTPM11 KO data

#Defining X axis as Days 8 - 22 with data points every 2 days
C57_DP_Days  <- seq(8,22,2)

#5 C57BL/6 mice injected with B16FO as the control group for both DNMT3A and PTPN11 KO groups
#Data vector is PER MOUSE, each point is an individual mouse tumor volume in mm3
#Data from mice 1, 2, 8, 9, and 15
DPWT_C571 <- matrix(c(C57_DP_Days[1:3], log(c(86.4,	185.2,	547.1))), nrow = 3, byrow = FALSE)
DPWT_C572 <- matrix(c(C57_DP_Days[1:4], log(c(77.2,	151.5,	471.0,	965.9))), nrow = 4, byrow = FALSE)
DPWT_C573 <- matrix(c(C57_DP_Days[1:6], log(c(76.3,	150.8,	278.1,	665.9,	544.8,	1453.3))), nrow = 6, byrow = FALSE)
DPWT_C574 <- matrix(c(C57_DP_Days,      log(c(40.8,	167.7,	411.7,	304.0,	975.7,	883.5,	1669.2,	2417.9))), nrow = 8, byrow = FALSE)
DPWT_C575 <- matrix(c(C57_DP_Days[1:4], log(c(68.1,	165.6,	355.2,	597.8))), nrow = 4, byrow = FALSE)

DPWT_C57 <- abind(DPWT_C571, DPWT_C572, DPWT_C573, DPWT_C574, DPWT_C575, along = 1)
rm(list = 'DPWT_C571', 'DPWT_C572', 'DPWT_C573', 'DPWT_C574', 'DPWT_C575')

#________________________________________________________________________________________________
# Defining the Markov Chain Analysis Function

#This function performs the Markov Chain analysis for a single data set and set of starting parameters.
#Sigma refers to the step size between test points for each parameter.
#Posterior refers to the likelihood function you wish to use to evaluate the test point.
#Data refers to the dataset you wish to analyze. 

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
  
#________________________________________________________________________________________________
# Automating the Plotting Output for the Markov Chain Analysis
  
#This function is used to automate the graphing for the PDF outputs of the Markov Chain Analysis.
#It is designed to simultaneously plot three chains for each of the 6 parameters.
#The final plot is the acceptance fraction of the test points, which should hover at 20%.

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

#________________________________________________________________________________________________
# Creating Indivualized Likelihood Functions for Each Dataset

#Next, define the likelihood function for your analysis.
#In this case, we are looking at linear equations and hoping to define a common slope
#for all of the mice (representing the growth rate of the cell line in vivo)
#and a y intercept for each mouse (representing the initial injection size).
#A linear equation defining the tumor growth for each mouse is written, and the sum squared error
#for each equation is summed together to determine the overall error for the potential data point.


#Likelihood evaluation for DNMT3A and PTPN11 Control in C57B/6 Mice
#Note that each mouse has a different number of data points based on initial tumor visibility and time of sacrifice. Ensure
#each Yhat analyzed based on the appropriate number of data points for each mouse.
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

#________________________________________________________________________________________________
#Data Analysis

#Define sigma to represent the step size you want to start with between data points. In this case, the common slope has a much
#smaller step size than each of the individual y intercepts

Sig <- c(.1, 1, 1, 1, 1, 1)

#Use the MHmcmc function three times for each dataset. Vary the starting points to ensure the result isn't dependent.

#Analysis of C57B/6 mice injected with B16F0 WT cells as control group for both DNMT3A and PTPM11 KO groups

T_DPWT_C571 <- MHmcmc(Sig, LLH_DPWT_C57, DPWT_C57, steps = 100000, target = 0.2, randomSeed = 4, startValue = c(.3, 3, 3, 3, 3, 3))
T_DPWT_C572 <- MHmcmc(Sig, LLH_DPWT_C57, DPWT_C57, steps = 100000, target = 0.2, randomSeed = 5, startValue = c(.2, 2.75, 2.5, 2.75, 2.25, 2.5))
T_DPWT_C573 <- MHmcmc(Sig, LLH_DPWT_C57, DPWT_C57, steps = 100000, target = 0.2, randomSeed = 6, startValue = c(.5, 1, 1, 1, 1, 1))

#Gelman Rubin Convergence Diagnostic to determine if the three Markov Chains for each dataset are 
#converging to the same point. A multivariate psrf value close to 1 indicates convergence.
#Verify convergence before generating a pdf of your results to save time

#Convergence for C57B/6 Data
DPWT_C571_MCMC <- mcmc(data = T_DPWT_C571[,1:6])
DPWT_C572_MCMC <- mcmc(data = T_DPWT_C572[,1:6])
DPWT_C573_MCMC <- mcmc(data = T_DPWT_C573[,1:6])

DPWT_C57_List <- mcmc.list(DPWT_C571_MCMC, DPWT_C572_MCMC, DPWT_C573_MCMC)

gelman.diag( DPWT_C57_List, confidence = 0.95, transform = FALSE, autoburnin = TRUE, multivariate = TRUE)

#________________________________________________________________________________________________
#Generating a PDF of the resulting data

#A density plot is generated for each parameter. The Markov Chain is also plotted for each parameter to verify your
#analysis is adequately covering the entire sample space. Finally, the acceptance fraction is plotted to ensure you
#are reaching the desired 20% acceptance rate in adequate time

pdf("Example.pdf", width = 11, height = 8)
Multiplot(T_DPWT_C571, T_DPWT_C572, T_DPWT_C573)
dev.off() # End plot DPWT_C573)
dev.off() # End plot 