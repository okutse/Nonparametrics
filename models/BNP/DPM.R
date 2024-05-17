## Dirichlet Process Mixtures

## clean the environment
rm(list = ls())
## load the data 
load("~/Desktop/PhD/Spring 2024/Bayesian/FinalProject/three.RData")

library(ChiRP)
library(tidyverse)
library(kableExtra)
library(latex2exp)
library(foreach)
library(doParallel)

## CONVERGENCE CHECKS
##-----------------------------------------------------------------------------
## n = 500
dpmod <- function(df, burnin, iter, init_k){
  df <- base::subset(df, R == 1, select = -c(e, z1, z2, z3, z4, R))
  N = dim(df)[1]
  d1 = df
  d1$A <- 1
  d0 <- df
  d0$A <- 0
  d_test <- rbind(d1, d0)
  res = ChiRP::fDPMix(d_train = df, formula = y ~ A + x1 + x2 + x3 + x4, burnin = burnin, iter = iter, d_test = d_test)
  et1 <- rowMeans(t(res$test[1:N,])) 
  et0 <- rowMeans(t(res$test[(N+1):(2*N),]))
  ATE = et1 - et0
  return(data.frame("et1" = et1, "et0" = et1, "ATE" = ATE))
}


## create MCMC chains for convergence and mixing analysis 
set.seed(123)
mod1 <- coda::mcmc(dpmod(df = s_threea, burnin = 2000, iter = 6000, init_k = 10))
mod2 <- coda::mcmc(dpmod(df = s_threea, burnin = 2000, iter = 6000, init_k = 100))
mod3 <- coda::mcmc(dpmod(df = s_threea, burnin = 2000, iter = 6000, init_k = 50))
mod4 <- coda::mcmc(dpmod(df = s_threea, burnin = 2000, iter = 6000, init_k = 200))


## merge all MCMC chains
allchains <- coda::mcmc.list(mod2, mod4)
allchains2 <- coda::mcmc.list(mod1, mod2, mod3, mod4)


## plot the chains to examine mixing
## plot(allchains)
plot(allchains2)

# ggplot()+
#   geom_line(data = data.frame(allchains2[[1]]), aes(x = c(1:length(ATE)), y = ATE), color = "blue")+
#   geom_line(data = data.frame(allchains2[[2]]), aes(x = c(1:length(ATE)), y = ATE), color = "red")+
#   geom_line(data = data.frame(allchains2[[3]]), aes(x = c(1:length(ATE)), y = ATE), color = "gray") +
#   geom_line(data = data.frame(allchains2[[4]]), aes(x = c(1:length(ATE)), y = ATE), color = "green", alpha = .5)+
#   labs(xlab = "Iteration, M", ylab = "Treatment Effect")+
#   scale_color_manual(name = "Chain", values = c("Chain 1"="blue", "red", "gray", "green"), labels = c("Chain 1", "Chain 2", "Chain 3", "Chain 4"))+
#   theme_bw()
  

## Mixing of the chains
par(mfrow = c(1, 2))
plot(x = data.frame(allchains2[[1]])$ATE, type = 'l', xlab = "Iterations, M", ylab = TeX("Treatment effect, $\\tau$"), col = 1)
lines(x = data.frame(allchains2[[2]])$ATE, type = 'l', xlab = "Iterations, M", ylab = TeX("Treatment effect, $\\tau$"), col = 2)
lines(x = data.frame(allchains2[[3]])$ATE, type = 'l', xlab = "Iterations, M", ylab = TeX("Treatment effect, $\\tau$"), col = 3)
lines(x = data.frame(allchains2[[4]])$ATE, type = 'l', xlab = "Iterations, M", ylab = TeX("Treatment effect, $\\tau$"), col = 4)


## Density of the estimated effect over chains
plot(density(x = data.frame(allchains2[[1]])$ATE, bw = 0.875), lwd = 2,
     col = 1, main = " ")


## Gelman-Rubin Statistic
gelmanRubin <- function(mcmc_chains) {
  n_chains <- length(mcmc_chains)
  n_iter <- nrow(mcmc_chains[[1]])
  n_params <- ncol(mcmc_chains[[1]])
  xnames <- names(mcmc_chains)
  # Calculate means and variances for each chain using lapply
  chain_means <- t(sapply(mcmc_chains, colMeans))
  within_variances <- t(sapply(mcmc_chains, function(x) apply(x, 2, var)))
  # Calculate between-chain variance and average within-chain variance
  b <- n_iter * apply(chain_means, 2, var)
  w <- rowMeans(within_variances)
  # Calculate weighted variance
  var_plus <- ((n_iter - 1) / n_iter) * w + mean(b) / n_iter
  # Calculate potential scale reduction factor
  psrf <- sqrt(var_plus / w)
  rownames(psrf) <- xnames
  return(psrf)
}
gelmanRubin(allchains2)





## Model 1: Observed Data Analysis with All Covariates 
##------------------------------------------------------------------------------
model_oneb <- function(df){
  df <- base::subset(df, R == 1, select = -c(e, z1, z2, z3, z4, R))
  N = dim(df)[1]
  
  d1 = df
  d1$A <- 1
  
  d0 <- df
  d0$A <- 0
  
  d_test <- rbind(d1, d0)
  res = ChiRP::fDPMix(d_train = df, formula = y ~ A + x1 + x2 + x3 + x4, burnin = 2000, iter = 6000, d_test = d_test)
  
  ATE_adjusted = mean(res$test[1:N,]) - mean(res$test[(N+1):(2*N),])
  bias_adjusted = ATE_adjusted - 50
  rslt = data.frame("ATE_adjusted" = ATE_adjusted, "bias_adjusted" = bias_adjusted)
  return(rslt)
}



## apply the function under model one to each data set in the list
cores <- detectCores()-1
registerDoParallel(cores)
cl = makeCluster(cores)
pkgs = clusterEvalQ(cl, c(library(tidyverse), library(tidymodels), library(magrittr), set.seed(123)))
# get the results under each  simulated data set

timea <- Sys.time()
oneaa = parLapply(cl, s3a, model_oneb)
timeb <- Sys.time()
time = timeb - timea
time ## 12 hours

onebb = parLapply(cl, s3b, model_oneb)
onecc = parLapply(cl, s3c, model_oneb)

t1 <- Sys.time()
onedd = parLapply(cl, s3d, model_oneb)
t2 <- Sys.time()
t2-t1
stopCluster(cl)


## convert the estimates to a data frame
oneaa <- oneaa %>% map_dfr(data.frame)
onebb <- onebb %>% map_dfr(data.frame)
onecc <- onecc %>% map_dfr(data.frame)
onedd <- onedd %>% map_dfr(data.frame)


#### Model 1 Results
##------------------------------------------------------------------------------
## case 1 [n = 500, SD = 1]
##------------------------------------------------------------------------------
obs <- c(n = nrow(subset(s_threea, R == 1)), ate = mean(oneaa$ATE_adjusted), sd = sd(oneaa$ATE_adjusted), bias = mean(oneaa$bias_adjusted))
obs
##------------------------------------------------------------------------------
## case 2 [n = 500, SD = 45]
##------------------------------------------------------------------------------
## observed
obs_2 <- c(n = nrow(subset(s_threeb, R == 1)), ate = mean(onebb$ATE_adjusted), sd = sd(onebb$ATE_adjusted), bias = mean(onebb$bias_adjusted))
obs_2
##------------------------------------------------------------------------------
## case 3 [n = 2000, SD = 1]
##------------------------------------------------------------------------------
## observed
obs_3 <- c(n = nrow(subset(s_threec, R == 1)), ate = mean(onecc$ATE_adjusted), sd = sd(onecc$ATE_adjusted), bias = mean(onecc$bias_adjusted))
obs_3
##------------------------------------------------------------------------------
## case 4 [n = 2000, SD = 45]
##------------------------------------------------------------------------------
## observed
obs_4 <- c(n = nrow(subset(s_threed, R == 1)), ate = mean(onedd$ATE_adjusted), sd = sd(onedd$ATE_adjusted), bias = mean(onedd$bias_adjusted))
obs_4
## write results to table
sim3_model_one = bind_rows(list("n = 500, SD = 1" = obs, 
                                "n = 500, SD = 45" = obs_2,
                                "n = 2000, SD = 1" = obs_3, 
                                "n = 2000, SD = 45"=  obs_4),
                           .id = "DGP") 


## Model 2: Observed Data Analysis With X1 Excluded
##---------------------------------------------------------------
# model 2 under simulation setting three where adjustment model excludes X1
model_twob <- function(df){
  df <- base::subset(df, R == 1, select = -c(e, z1, z2, z3, z4, R))
  N = dim(df)[1]
  
  d1 = df
  d1$A <- 1
  
  d0 <- df
  d0$A <- 0
  
  d_test <- rbind(d1, d0)
  res = ChiRP::fDPMix(d_train = df, formula = y ~ A + x2 + x3 + x4, burnin = 2000, iter = 6000, d_test = d_test)
  
  ATE_adjusted = mean(res$test[1:N,]) - mean(res$test[(N+1):(2*N),])
  bias_adjusted = ATE_adjusted - 50
  rslt = data.frame("ATE_adjusted" = ATE_adjusted, "bias_adjusted" = bias_adjusted)
  return(rslt)
}
  

## apply the function under model one to each data set in the list
cores <- detectCores()-1
registerDoParallel(cores)
cl = makeCluster(cores)
pkgs = clusterEvalQ(cl, c(library(tidyverse), library(tidymodels), library(magrittr), set.seed(123)))
# get the results under each  simulated data set
twoaa = parLapply(cl, s3a, model_twob) ## 12 hrs
twobb = parLapply(cl, s3b, model_twob) ## 2 hrs
tt1 = Sys.time()
twocc = parLapply(cl, s3c, model_twob) ## 12.9 hrs
tt2 = Sys.time()
tt2-tt1
tt2 = Sys.time()
twodd = parLapply(cl, s3d, model_twob) ##12.9 hrs
tt3 = Sys.time()
tt3 - tt2
stopCluster(cl)

## convert the estimates to a data frame
twoaa <- twoaa %>% map_dfr(data.frame)
twobb <- twobb %>% map_dfr(data.frame)
twocc <- twocc %>% map_dfr(data.frame)
twodd <- twodd %>% map_dfr(data.frame)


## Model 2 results
##------------------------------------------------------------------------------
## case 1 [n = 500, SD = 1]
##------------------------------------------------------------------------------
obs <- c(n = nrow(subset(s_threea, R == 1)), ate = mean(twoaa$ATE_adjusted), sd = sd(twoaa$ATE_adjusted), bias = mean(twoaa$bias_adjusted))
obs
##------------------------------------------------------------------------------
## case 2 [n = 500, SD = 45]
##------------------------------------------------------------------------------
## observed
obs_2 <- c(n = nrow(subset(s_threeb, R == 1)), ate = mean(twobb$ATE_adjusted), sd = sd(twobb$ATE_adjusted), bias = mean(twobb$bias_adjusted))
obs_2
##------------------------------------------------------------------------------
## case 3 [n = 2000, SD = 1]
##------------------------------------------------------------------------------
## observed
obs_3 <- c(n = nrow(subset(s_threec, R == 1)), ate = mean(twocc$ATE_adjusted), sd = sd(twocc$ATE_adjusted), bias = mean(twocc$bias_adjusted))
obs_3
##------------------------------------------------------------------------------
## case 4 [n = 2000, SD = 45]
##------------------------------------------------------------------------------
## observed
obs_4 <- c(n = nrow(subset(s_threed, R == 1)), ate = mean(twodd$ATE_adjusted), sd = sd(twodd$ATE_adjusted), bias = mean(twodd$bias_adjusted))
obs_4
## write results to table
sim3_model_two = bind_rows(list("n = 500, SD = 1" = obs, 
                                "n = 500, SD = 45" = obs_2,
                                "n = 2000, SD = 1" = obs_3, 
                                "n = 2000, SD = 45"=  obs_4),
                           .id = "DGP") 


## load saved model one results
sim3_model_one2 <- read.csv("/Users/aokutse/Desktop/PhD/Spring 2024/Bayesian/FinalProject/models/BNP/sim3_model_one.csv")
dpmm <- rbind(sim3_model_one2, sim3_model_two)

## save the CMLR model [first 4 rows are model 1 then rest are 2]
write.csv(dpmm, file = "/Users/aokutse/Desktop/PhD/Spring 2024/Bayesian/FinalProject/models/BNP/dpmm.csv", row.names = FALSE)
