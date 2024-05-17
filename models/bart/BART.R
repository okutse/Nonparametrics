## Bayesian Additive Regression Trees (BART)

## clean the environment
rm(list = ls())
## load the data 
load("~/Desktop/PhD/Spring 2024/Bayesian/FinalProject/three.RData")

library(tidyverse)
library(kableExtra)
library(latex2exp)
library(foreach)
library(doParallel)

## Simulation Setting Three Observed/Complete Case Data Analysis

### Model 1: Adjustment for All Covariates
model_oneb <- function(df = NULL){
  df = dplyr::filter(df, R == 1)
  bart_two <- parsnip::bart(
    mode = "regression") %>%
    set_engine("dbarts") %>% 
    fit(formula = y ~ A + x1 + x2 + x3 + x4, data = df)
  ## set A = 0 and generate predictions for everyone
  df_A0 <- df
  df_A0$A <- 0
  pred_A0 <- predict(bart_two, df_A0)
  ## set A = 1 and generate predictions for everyone
  df_A1 <- df
  df_A1$A <- 1
  pred_A1 <- predict(bart_two, df_A1)
  ## compute the ATE
  ATE_adjusted = mean(pred_A1$.pred - pred_A0$.pred)
  ## compute the bias
  bias_adjusted = ATE_adjusted - 50
  ## return the results as a data frame
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
time 
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



## Model 2: Exclude X1 in Covariate Adjustment
##------------------------------------------------------------------------------
model_twob <- function(df = NULL){
  ## filter the data to have only individuals with R = 1
  df = dplyr::filter(df, R == 1)
  # fit random forest model for all individuals
  bart_two <- parsnip::bart(
    mode = "regression") %>%
    set_engine("dbarts") %>% 
    fit(formula = y ~ A + x2 + x3 + x4, data = df)
  ## set A = 0 and generate predictions for everyone
  df_A0 <- df
  df_A0$A <- 0
  pred_A0 <- predict(bart_two, df_A0)
  ## set A = 1 and generate predictions for everyone
  df_A1 <- df
  df_A1$A <- 1
  pred_A1 <- predict(bart_two, df_A1)
  ## compute the ATE
  ATE_adjusted = mean(pred_A1$.pred - pred_A0$.pred)
  ## compute the bias
  bias_adjusted = ATE_adjusted - 50
  ## return the results as a data frame
  rslt = data.frame("ATE_adjusted" = ATE_adjusted, "bias_adjusted" = bias_adjusted)
  return(rslt)
}

## apply the function under model one to each data set in the list
cores <- detectCores()-1
registerDoParallel(cores)
cl = makeCluster(cores)
pkgs = clusterEvalQ(cl, c(library(tidyverse), library(tidymodels), library(magrittr), set.seed(123)))
# get the results under each  simulated data set
twoaa = parLapply(cl, s3a, model_twob) 
twobb = parLapply(cl, s3b, model_twob) 
tt1 = Sys.time()
twocc = parLapply(cl, s3c, model_twob) 
tt2 = Sys.time()
tt2-tt1
tt2 = Sys.time()
twodd = parLapply(cl, s3d, model_twob) 
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

bart <- rbind(sim3_model_one, sim3_model_two)

## save the CMLR model [first 4 rows are model 1 then rest are 2]
write.csv(bart, file = "/Users/aokutse/Desktop/PhD/Spring 2024/Bayesian/FinalProject/models/bart/bart.csv", row.names = FALSE)
