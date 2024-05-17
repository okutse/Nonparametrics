## Bayesian Generative Machine Learning and Missing Outcome Data
## Amos Okutse
## Brown University School of Public Health
## Version last edited Sunday Apr 14th 2024
##-----------------------------------------------------------------------------

# NOTE: All analysis are under missing outcome data corresponding to 
# Case II settings 1 to 3 in the main paper.

## load all required packages
# function to install missing packages
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE, repos='http://cran.rstudio.com/')
  sapply(pkg, require, character.only = TRUE)
}
packages =c( "tidyverse","knitr", "kableExtra","skimr", "MatchIt", "RItools","optmatch", "ggplot2", "tufte", "tufterhandout", "plotly", "snowfall", "rstan", "gridExtra", "knitr", "gtsummary", "data.table", "GGally", "MASS", "broom", "boot", "foreach", "doParallel", "glmnet", "tidymodels" , "usemodels", "magrittr", "modelr")
ipak(packages)

## DGP
## -----------------------------------------------------------------------------
## data generation function for all simulation settings
## the saved data file should contain both the Z and the X for adjustment under correct and incorrect model specifications. 
## the linear model adjusting for the Z's is the correctly specified model and the efficiency gains given this model represent the optimal expected gains. [no misspecified model should outperform this model]
set.seed(123)
dgp2 <- function(n, ss, theta, eta0, eta1) {
  #' theta is the modifying effect of Z_1 on the treatment, A
  #' eta0 is the effect of A on the missing data mechanism
  #' ss is the residual standard deviation on the true outcome model
  #' n is the desired sample size
  #' eta1 is the effect modification of Z_1 on A in the missing data model
  
  # treatment variable
  A <- c()
  for (i in 1:n){if (i <= n/2){ A[i] = 1} else {if (i > n/2){A[i] = 0}}}
  # create the data frame with the true covariates
  mat <- data.frame(
    z1 <- rnorm(n = n, mean = 0, sd = 1), 
    z2 <- rnorm(n = n, mean = 0, sd = 1),
    z3 <- rnorm(n = n, mean = 0, sd = 1), 
    z4 <- rnorm(n = n, mean = 0, sd = 1))
  colnames(mat) <- c("z1", "z2", "z3", "z4")
  # generate the error term
  e <- rnorm(n = n, mean = 0, sd = ss)
  # generate the outcome variable based on the specified model
  y <- 210 + 50*A + theta*A*mat[, 1] + 27.4*mat[ ,1] + 13.7*mat[, 2] + 13.7*mat[, 3] + 13.7*mat[, 4] + e
  # generate the xi's actually observed by analyst
  x1 <- exp(mat[, 1]/2); x2 <- (mat[, 2]/ (1 + exp(mat[, 1]))) + 10; x3 <- (((mat[, 1]*mat[, 3])/25) + 0.6)^3; x4 <- (mat[,2] + mat[,4] + 20)^2
  # create the missing data variable based on the treatment
  pi <- locfit::expit(0 - mat[, 1] + eta0*A + eta1*A*mat[, 1] + 0.5*mat[, 2] - 0.25*mat[, 3] - 0.1*mat[, 4])
  R = rbinom(n = n, size = 1, prob = pi)
  # save variables as a data frame
  df <-data.frame(y, A, x1, x2, x3, x4, R, e, mat)
} 

##------------------------------------------------------------------------------
## Data sets for setting 1: set eta1 = theta = eta0 = 0
##------------------------------------------------------------------------------
set.seed(123)
# The organization of the files for each simulation setting hold as below
#sd = 1; n = 500
s_onea <- dgp2(n = 500, ss = 1, theta = 0, eta0 = 0, eta1 = 0)
# sd = 45; n = 500
s_oneb <- dgp2(n = 500, ss = 45, theta = 0, eta0 = 0, eta1 = 0)
# sd = 1; n = 2000
s_onec <- dgp2(n = 2000, ss = 1, theta = 0, eta0 = 0, eta1 = 0)
# sd = 45; n = 2000
s_oned <- dgp2(n = 2000, ss = 45, theta = 0, eta0 = 0, eta1 = 0)

## 1000 simulated data sets saved as a list
cores <- parallel::detectCores()
doParallel::registerDoParallel(cores - 1)
s1a = foreach(1:1000) %dopar% dgp2(n = 500, ss = 1, theta = 0, eta0 = 0, eta1 = 0)
s1b = foreach(1:1000) %dopar% dgp2(n = 500, ss = 45, theta = 0, eta0 = 0, eta1 = 0)
s1c = foreach(1:1000) %dopar% dgp2(n = 2000, ss = 1, theta = 0, eta0 = 0, eta1 = 0)
s1d = foreach(1:1000) %dopar% dgp2(n = 2000, ss = 45, theta = 0, eta0 = 0, eta1 = 0)

## Save all data sets as one for analyses under simulation setting one (one.RData)
save("s_onea", "s_oneb", "s_onec", "s_oned", "s1a", "s1b", "s1c", "s1d", file = "/Users/aokutse/Desktop/PhD/Spring 2024/Bayesian/FinalProject/one.RData")




##------------------------------------------------------------------------------
## Setting 2: theta1 = eta1 = 0 and eta0 = 1
##------------------------------------------------------------------------------

set.seed(123)
#sd = 1; n = 500
s_twoa <- dgp2(n = 500, ss = 1, theta = 0, eta0 = 1, eta1 = 0)
# sd = 45; n = 500
s_twob <- dgp2(n = 500, ss = 45, theta = 0, eta0 = 1, eta1 = 0)
# sd = 1; n = 2000
s_twoc <- dgp2(n = 2000, ss = 1, theta = 0, eta0 = 1, eta1 = 0)
# sd = 45; n = 2000
s_twod <- dgp2(n = 2000, ss = 45, theta = 0, eta0 = 1, eta1 = 0)


## Simulate 1000 data sets and save them as a list for further analysis
cores <- parallel::detectCores()
doParallel::registerDoParallel(cores - 1)
s2a = foreach(1:1000) %dopar% dgp2(n = 500, ss = 1, theta = 0, eta0 = 1, eta1 = 0)
s2b = foreach(1:1000) %dopar% dgp2(n = 500, ss = 45, theta = 0, eta0 = 1, eta1 = 0)
s2c = foreach(1:1000) %dopar% dgp2(n = 2000, ss = 1, theta = 0, eta0 = 1, eta1 = 0)
s2d = foreach(1:1000) %dopar% dgp2(n = 2000, ss = 45, theta = 0, eta0 = 1, eta1 = 0)

## save all data sets as one for analyses under simulation setting two data (two.RData)
save("s_twoa", "s_twob", "s_twoc", "s_twod", "s2a", "s2b", "s2c", "s2d", file = "/Users/aokutse/Desktop/PhD/Spring 2024/Bayesian/FinalProject/two.RData")




##------------------------------------------------------------------------------
## Setting 3: eta0 = eta1 = 1 and theta = 0
##------------------------------------------------------------------------------

set.seed(123)
## data sets for simulation setting three
#sd = 1; n = 500
s_threea <- dgp2(n = 500, ss = 1, theta = 0, eta0 = 1, eta1 = 1)
# sd = 45; n = 500
s_threeb <- dgp2(n = 500, ss = 45, theta = 0, eta0 = 1, eta1 = 1)
# sd = 1; n = 2000
s_threec <- dgp2(n = 2000, ss = 1, theta = 0, eta0 = 1, eta1 = 1)
# sd = 45; n = 2000
s_threed <- dgp2(n = 2000, ss = 45, theta = 0, eta0 = 1, eta1 = 1)


## Simulate 1000 data sets and save them as a list for further analysis
cores <- parallel::detectCores()
doParallel::registerDoParallel(cores - 1)
s3a = foreach(1:1000) %dopar% dgp2(n = 500, ss = 1, theta = 0, eta0 = 1, eta1 = 1)
s3b = foreach(1:1000) %dopar% dgp2(n = 500, ss = 45, theta = 0, eta0 = 1, eta1 = 1)
s3c = foreach(1:1000) %dopar% dgp2(n = 2000, ss = 1, theta = 0, eta0 = 1, eta1 = 1)
s3d = foreach(1:1000) %dopar% dgp2(n = 2000, ss = 45, theta = 0, eta0 = 1, eta1 = 1)

## Save all data sets as one for analyses corresponding to simulation setting three data (one.RData)
save("s_threea", "s_threeb", "s_threec", "s_threed", "s3a", "s3b", "s3c", "s3d", file = "/Users/aokutse/Desktop/PhD/Spring 2024/Bayesian/FinalProject/three.RData")





