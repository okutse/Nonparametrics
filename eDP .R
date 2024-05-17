## enriched Dirichlet Process Prior
## Bayesian Nonparametrics
## Amos Okutse
## Adapted from https://github.com/jasonroy0/EDP_causal/blob/master/CppCode/cluster_continuous2.cpp



expit <- function(x) {
  # The expit function, also known as the logistic function,
  # maps any real-valued number to the (0, 1) interval, which
  # can be used to convert a linear regression output to a probability.
  return(1 / (1 + exp(-x)))
}


# Define the rmultinomF function in R
rmultinomF <- function(p) {
	#' that simulates a single draw from a multinomial distribution. It takes a vector of 
	#' probabilities p as input, calculates the cumulative sum of these probabilities, generates a random uniform number, 
	and then determines the index at which this random number exceeds the cumulative probability.
  # Normalize the probabilities so they sum to 1
  p <- p / sum(p)
  
  # Calculate the cumulative sum of probabilities
  csp <- cumsum(p)
  
  # Generate a random uniform number between 0 and 1
  rnd <- runif(1)
  
  # Initialize the result
  res <- 0
  
  # Get the size of the probability vector
  psize <- length(p)
  
  # Loop over the cumulative probabilities
  for(i in 1:psize) {
    # If the random number is greater than the current cumulative probability,
    # increment the result
    if(rnd > csp[i]) {
      res <- res + 1
    }
  }
  
  # The function returns the index of the interval that contains the random number,
  # incremented by 1 because R indices start at 1 instead of 0
  return(res + 1)
}


# Load the MASS library for the mvrnorm function
library(MASS) ## check out if functions are in LaplacesDemon package!!!
# Define the mvrnorm function in R
mvrnormR <- function(mu, Sigma) {
  # The mvrnorm function generates a random draw from a multivariate normal distribution.
  return(mvrnorm(n = 1, mu = mu, Sigma = Sigma))
}

# Define the rinvchisq function in R
rinvchisq <- function(df, scale) {
  # The rinvchisq function generates a random draw from an inverse chi-square distribution.
  return(1 / rchisq(n = 1, df = df, ncp = 0) / scale)
}

# Load the mvtnorm library for the dmvnorm function
library(mvtnorm)
# Define the dmvn function in R
dmvnR <- function(x, mu, Sigma, logt = FALSE) {
  # The dmvnorm function calculates the density of a multivariate normal distribution.
  # If logt is TRUE, it returns the log density; otherwise, it returns the density.
  return(dmvnorm(x = x, mean = mu, sigma = Sigma, log = logt))
}


# Define the updatevar function in R
updatevar <- function(x, nu0, tau0, c0, mu0) {
  # The function updates the variance in a Bayesian updating scheme.
  # that updates the variance in a Bayesian updating scheme. It takes as input a data point x, prior degrees of freedom nu0, prior scale tau0, prior sample size c0, and prior mean mu0.
  # Initialize the variance
  varx <- 0
  # Get the size of the data
  xsize <- length(x)
  # Update the degrees of freedom
  newdf <- nu0 + xsize
  
  # Calculate the variance if there is more than one data point
  if (xsize > 1) { 
    varx <- var(x) 
  }
  # Calculate the numerator of the scale parameter
  numer <- nu0 * tau0 + (xsize - 1) * varx + (c0 * xsize / (c0 + xsize)) * (mean(x) - mu0)^2
  # Generate a random draw from an inverse chi-square distribution
  newval <- 1 / rchisq(n = 1, df = newdf, ncp = 0) / (numer / newdf)
  # Return the updated variance
  return(newval)
}


# Define the updatemean function in R
updatemean <- function(x, tau, c0, mu0) {
  # The function updates the mean in a Bayesian updating scheme.
  # Get the size of the data
  xsize <- length(x)
  
  # Calculate the new variance
  newvar <- 1 / (c0 / tau + xsize / tau)
  
  # Calculate the new mean
  newmean <- (mu0 * c0 / tau + mean(x) * xsize / tau) * newvar
  
  # Generate a random draw from a normal distribution with the new mean and standard deviation
  newval <- rnorm(n = 1, mean = newmean, sd = sqrt(newvar))
  
  # Return the updated mean
  return(newval)
}


# Define the newbetafunction in R
newbetafunction <- function(curbet, x, y, betainit, diagbetacov0) {
  # The function performs a Metropolis-Hastings step in a Bayesian updating scheme.
  
  # Get the dimension of the regression coefficients
  dimbet <- length(curbet)
  
  # Define the proposal variance
  propvar <- diag(dimbet) * 0.01
  
  # Generate a proposed value for the regression coefficients
  proposed <- MASS::mvrnorm(n = 1, mu = curbet, Sigma = propvar)
  
  # Calculate the log likelihood of the proposed value
  loglikenew <- dbinom(y, size = 1, prob = plogis(sum(x * proposed)), log = TRUE) +
    mvtnorm::dmvnorm(x = proposed, mean = betainit, sigma = diagbetacov0, log = TRUE) +
    mvtnorm::dmvnorm(x = curbet, mean = proposed, sigma = propvar, log = TRUE)
  
  # Calculate the log likelihood of the current value
  loglikeold <- dbinom(y, size = 1, prob = plogis(sum(x * curbet)), log = TRUE) +
    mvtnorm::dmvnorm(x = curbet, mean = betainit, sigma = diagbetacov0, log = TRUE) +
    mvtnorm::dmvnorm(x = proposed, mean = curbet, sigma = propvar, log = TRUE)
  
  # Calculate the acceptance probability
  ans <- min(exp(loglikenew - loglikeold), 1)
  
  # If a uniform random number is less than the acceptance probability, update the regression coefficients
  if(runif(1) < ans)
    curbet <- proposed
  
  # Return the updated regression coefficients
  return(curbet)
}


### For Continuous Data
##################################

# Define the newregvar function in R
newregvar <- function(x, y, betaa0, betab0, beta0) {
  # Calculate the size of the vector x
  allp <- length(x)
  
  # Update the parameter 'an' for the gamma distribution
  an <- betaa0 + 0.5
  
  # Create a precision matrix with 'allp' dimensions, initialized to the identity matrix
  prec0 <- diag(allp)
  
  # Calculate the new precision matrix
  newprec <- tcrossprod(x) + prec0
  
  # Calculate the new beta vector
  betan <- solve(newprec, prec0 %*% beta0 + x * y)
  
  # Calculate the first part of the scale parameter for the gamma distribution
  part1 <- t(beta0) %*% prec0 %*% beta0
  
  # Calculate the second part of the scale parameter for the gamma distribution
  part2 <- t(betan) %*% newprec %*% betan
  
  # Check for dimension errors
  if(any(dim(part1) > 1) || any(dim(part2) > 1)) {
    stop("Something wrong in newregvar function")
  }
  
  # Update the parameter 'bn' for the gamma distribution
  bn <- betab0 + 0.5 * (y^2 + part1[1,1] - part2[1,1])
  
  # Calculate the new variance as the inverse of a gamma distribution draw
  outvar <- 1 / rgamma(1, shape = an, rate = bn)
  
  # Return the updated variance
  return(outvar)
}

# Define the newbet function in R
newbet <- function(x, y, sig, beta0) {
  # The function performs a Bayesian updating step for the regression coefficients.
  
  # Get the size of the vector x
  allp <- length(x)
  
  # Create a precision matrix with 'allp' dimensions, initialized to the identity matrix
  prec0 <- diag(allp)
  
  # Calculate the new precision matrix
  newprec <- tcrossprod(x) + prec0
  
  # Calculate the new beta vector
  betan <- solve(newprec, prec0 %*% beta0 + x * y)
  
  # Generate a proposed value for the regression coefficients
  bet <- MASS::mvrnorm(n = 1, mu = betan, Sigma = sig * solve(newprec))
  
  # Return the updated regression coefficients
  return(bet)
}




# Define the cluster function in R
cluster <- function(y, X, matX, Sy, Sx,
                    betaY, sig2, xPiPars, xMuPars, xSigPars,
                    alphapsi, alphatheta, h0y, h0i,
                    uniqueS, 
                    c0, mu0,
                    nu0, tau,
                    a0, b0,
                    betaa0, betab0,
                    betainit, diagbetacov0,
                    p1, ptx, p2) {
  # The function performs a clustering operation on continuous data in a Bayesian framework.
  
  # Get the number of observations and predictors
  n <- nrow(X)
  p <- ncol(X)
  
  # Initialize some variables
  indY <- integer(0)
  indX <- integer(0)
  numY <- 0
  uniqueY <- integer(0)
  numX <- integer(0)
  numTotalCluster <- 0
  ind_dummy <- integer(0)
  ind_dummy2 <- integer(0)
  num_dummy <- 0
  num_dummy2 <- 0
  count <- 0
  likeregy <- 0
  prodx <- 0
  prodx2 <- 0
  newCluster <- 0
  newpipars <- rep(0, p1 + ptx)
  newmupars <- rep(0, p2)
  newsigpars <- rep(0, p2)
  betadraw <- numeric(0)
  newbeta <- numeric(0)
  newsig <- 0
  
  # Loop through each person and change cluster memberships
  for(i in 1:n) {
    # Check if ith person is the lone person in his cluster
    ind_dummy <- which(Sy == Sy[i] & Sx == Sx[i])
    num_dummy <- length(ind_dummy)
    
    if(num_dummy == 1) { # If lone person in X-Y cluster
      # Delete associated coefficients in Y and X cluster
      ind_dummy2 <- which(Sy == Sy[i]) # Check if only person in Y cluster too
      num_dummy2 <- length(ind_dummy2)
      
      # Delete Y coef if only one in Y cluster
      if(num_dummy2 == 1) {
        betaY <- betaY[-Sy[i], ]
        sig2 <- sig2[-Sy[i]]
      }
      
      # Delete X coef
      # Should find row in uniqueS that corresponds to person i
      ind_dummy <- which(uniqueS[, 1] == Sy[i] & uniqueS[, 2] == Sx[i])
      
      xPiPars <- xPiPars[-ind_dummy, ]
      xMuPars <- xMuPars[-ind_dummy, ]
      xSigPars <- xSigPars[-ind_dummy, ]
      
      # Relabel X cluster
      for(j in 1:length(Sx)) {
        if(Sy[j] == Sy[i] && Sx[j] > Sx[i])
          Sx[j] <- Sx[j] - 1
      }
      
      for(j in 1:nrow(uniqueS)) {
        if(uniqueS[j, 1] == Sy[i] && uniqueS[j, 2] > Sx[i])
          uniqueS[j, 2] <- uniqueS[j, 2] - 1
      }
      
      # Relabel Y cluster (if needed)
      if(num_dummy2 == 1) {
        for(j in 1:length(Sy)) {
          if(Sy[j] > Sy[i])
            Sy[j] <- Sy[j] - 1
        }
        
        for(j in 1:nrow(uniqueS)) {
          if(uniqueS[j, 1] > Sy[i])
            uniqueS[j, 1] <- uniqueS[j, 1] - 1
        }
      }
      
      uniqueS <- uniqueS[-ind_dummy, ] # Get rid of row
    }
    
    # Need to delete row of Sy and Sx
    Sy <- Sy[-i]
    Sx <- Sx[-i]
    
    # Recalculate number of unique clusters
    numY <- nrow(betaY)
    numTotalCluster <- nrow(xMuPars)
    totalposs <- numY + numTotalCluster + 1
    probs <- numeric(totalposs)
    count <- 0
    
    # Counts for # in appropriate Y and X cluster, excluding the ith person
    njwoi <- 0
    nljwoi <- 0
    
    # # of X clusters within each Y cluster
    numXj <- 0
    
    for(j in 1:numY) {
      # Fill in probs for existing clusters
      
      # Get count of number of X clusters within jth Y cluster
      ind_dummy <- which(uniqueS[, 1] == (j + 1))
      numXj <- length(ind_dummy)
      
      # Get number of subjects within jth cluster
      ind_dummy <- which(Sy == (j + 1))
      njwoi <- length(ind_dummy)
      
      # Likelihood for each existing Y cluster
      likeregy <- dnorm(y[i], sum(matX[i, ] * betaY[j, ]), sqrt(sig2[j]), log = FALSE)
      
      for(k in 1:numXj) {
        prodx <- 1
        prodx2 <- 1
        
        ind_dummy <- which(Sy == (j + 1) & Sx == (k + 1))
        nljwoi <- length(ind_dummy)
        
        # Likelihood for binary covariates
        for(l in 1:(ptx + p1)) {
          prodx <- prodx * dbinom(X[i, l], size = 1, prob = xPiPars[count, l], log = FALSE)
        }
        
        # Likelihood for continuous covariates
        for(l in 1:p2) {
          prodx2 <- prodx2 * dnorm(X[i, ptx + p1 + l], mean = xMuPars[count, l], sd = sqrt(xSigPars[count, l]), log = FALSE)
        }
        
        probs[count] <- ((njwoi * nljwoi) / (njwoi + alphapsi)) * likeregy * prodx * prodx2
        count <- count + 1
      }
    }
    
    
    for(j in 1:numY) {
      # Fill in probs for new X clusters in existing Y clusters
      
      ind_dummy <- which(Sy == (j + 1))
      njwoi <- length(ind_dummy)
      
      # Likelihood for each existing Y cluster
      likeregy <- dnorm(y[i], sum(matX[i, ] * betaY[j, ]), sqrt(sig2[j]), log = FALSE)
      probs[numTotalCluster + j] <- ((njwoi * alphapsi) / (njwoi + alphapsi)) * likeregy * h0i[i]
    }
    
    probs[numY + numTotalCluster] <- alphatheta * h0y[i] * h0i[i] # Prob for new Y cluster
    
    # Use multinomial distribution to choose new cluster
    newCluster <- sample(1:length(probs), size = 1, prob = probs)
    probs <- rep(0, length(probs))
    
    # Need to map this integer to one of the clusters (or a new cluster)
    if(newCluster <= numTotalCluster) {
      Sy <- append(Sy, uniqueS[newCluster - 1, 1], after = i - 1)
      Sx <- append(Sx, uniqueS[newCluster - 1, 2], after = i - 1)
    } else {
      # Find out whether this is a new Y cluster or X cluster
      if(newCluster == (numTotalCluster + numY + 1)) {
        Sx <- append(Sx, 1, after = i - 1)
        
        # Need functions rm

# Update parameters
      newsig <- newregvar(t(matX[i, ]), y[i], betaa0, betab0, betainit)
      newbeta <- newbet(t(matX[i, ]), y[i], newsig, betainit)
      betaY <- rbind(betaY, t(newbeta))
      sig2 <- append(sig2, newsig)
      
      for(j in 1:(p1 + ptx)) {
        newpipars[j] <- rbeta(1, X[i, j] + a0, 1 - X[i, j] + b0)
      }
      
      for(j in 1:p2) {
        newsigpars[j] <- updatevar(X[i, (p1 + ptx + j)], nu0, tau, c0, mu0)
        newmupars[j] <- updatemean(X[i, (p1 + ptx + j)], tau, c0, mu0)
      }
      
      xPiPars <- rbind(xPiPars, newpipars)
      xSigPars <- rbind(xSigPars, newsigpars)
      xMuPars <- rbind(xMuPars, newmupars)
      
      Sy <- append(Sy, max(Sy) + 1, after = i)
      uniqueS <- rbind(uniqueS, c(max(Sy), Sx[i]))
    } else {
      # If new X cluster in existing Y cluster
      Sy <- append(Sy, newCluster - numTotalCluster, after = i)
      ind_dummy <- which(uniqueS[, 1] == Sy[i])
      Sx <- append(Sx, length(ind_dummy) + 1, after = i)
      
      # Update parameters
      for(j in 1:(p1 + ptx)) {
        newpipars[j] <- rbeta(1, X[i, j] + a0, 1 - X[i, j] + b0)
      }
      
      for(j in 1:p2) {
        newsigpars[j] <- updatevar(X[i, (p1 + ptx + j)], nu0, tau, c0, mu0)
        newmupars[j] <- updatemean(X[i, (p1 + ptx + j)], tau, c0, mu0)
      }
      
      ind_dummy <- which(uniqueS[, 1] <= Sy[i])
      num_dummy <- length(ind_dummy)
      
      xPiPars <- rbind(xPiPars[1:num_dummy, ], newpipars, xPiPars[(num_dummy + 1):nrow(xPiPars), ])
      xSigPars <- rbind(xSigPars[1:num_dummy, ], newsigpars, xSigPars[(num_dummy + 1):nrow(xSigPars), ])
      xMuPars <- rbind(xMuPars[1:num_dummy, ], newmupars, xMuPars[(num_dummy + 1):nrow(xMuPars), ])
      
      uniqueS <- rbind(uniqueS[1:num_dummy, ], c(Sy[i], Sx[i]), uniqueS[(num_dummy + 1):nrow(uniqueS), ])
    }
  }
  
  # Return a list containing the updated cluster memberships and parameters
  return(list(Sy = Sy, Sx = Sx, uniqueS = uniqueS, beta = betaY, sig2 = sig2, pipars = xPiPars, mupars = xMuPars, sigpars = xSigPars))
}













