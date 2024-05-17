## Dirichlet Process Mixtures Prelims

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
jpeg("conv_checks.jpeg", width = 6, height = 4, units = 'in', res = 300)
par(mfrow = c(1, 2))
cex = 1
plot(x = data.frame(allchains2[[1]])$ATE, type = 'l', xlab = "Iterations, M", ylab = TeX("Treatment effect, $\\hat{\\tau}$"), col = 1, cex.lab=cex, cex.main=cex, main = "(a)")
lines(x = data.frame(allchains2[[2]])$ATE, type = 'l', xlab = "Iterations, M", ylab = TeX("Treatment effect, $\\tau$"), col = 2, cex.lab=cex, cex.main=cex)
lines(x = data.frame(allchains2[[3]])$ATE, type = 'l', xlab = "Iterations, M", ylab = TeX("Treatment effect, $\\tau$"), col = 3, cex.lab=cex, cex.main=cex)
lines(x = data.frame(allchains2[[4]])$ATE, type = 'l', xlab = "Iterations, M", ylab = TeX("Treatment effect, $\\tau$"), col = 4, cex.lab=cex, cex.main=cex)


## Density of the estimated effect over chains
plot(density(x = data.frame(allchains2[[1]])$ATE, bw = 0.875), lwd = 2,
     col = 1, cex.lab=cex, cex.main=cex, main = "(b)")
dev.off()

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

