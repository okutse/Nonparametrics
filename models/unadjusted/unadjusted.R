## Unadjusted results
##--------------------------------------------------------------

## clean the environment
rm(list = ls())
## load the data 
load("~/Desktop/PhD/Spring 2024/Bayesian/FinalProject/three.RData")
## settings correspond to observed data analysis under Model 1 where all observed covariates X are
## included in analyses and Model 2 where X1 is excluded. X1 informs missing data in this case.




## model one under full data analysis ~ linear regression function with correct specification of adjustment covariates
model_onea <- function(df = NULL){
  ## since this is based on the full data set, then use the full data set
  full_unadjusted = mean(df$y[df$A == 1]) - mean(df$y[df$A == 0])
  full_bias_unadjusted = full_unadjusted - 50
  
  ## subset the data to only subjects with R == 1
  df2 <- dplyr::filter(df, R == 1)
  observed_unadjusted = mean(df2$y[df2$A == 1]) - mean(df2$y[df2$A == 0])
  observed_bias = observed_unadjusted - 50
  return(data.frame(full_unadjusted, full_bias_unadjusted, observed_unadjusted, observed_bias))
}

## apply the function under model one to each data set in the list
cores <- detectCores()-1
registerDoParallel(cores)
cl = makeCluster(cores)
pkgs = clusterEvalQ(cl, c(library(tidyverse), library(tidymodels), library(magrittr)))
parallel::clusterSetRNGStream(cl, 123)
# get the results under each  simulated data set
onea = parLapply(cl, s3a, model_onea)
oneb = parLapply(cl, s3b, model_onea)
onec = parLapply(cl, s3c, model_onea)
oned = parLapply(cl, s3d, model_onea)
stopCluster(cl)

## convert the estimates to a data frame
onea <- onea %>% map_dfr(data.frame)
oneb <- oneb %>% map_dfr(data.frame)
onec <- onec %>% map_dfr(data.frame)
oned <- oned %>% map_dfr(data.frame)



#### Results

######################################
## Simulation Setting Three.         #
######################################

options(scipen = 999)
##------------------------------------------------------------------------------
## Case 1: n = 500, SD = 1
##------------------------------------------------------------------------------
## Full analysis
full = c(n = nrow(s_threea), ate = mean(onea$full_unadjusted), sd = sd(onea$full_unadjusted), bias = mean(onea$full_bias_unadjusted), sd_bias = sd(onea$full_bias_unadjusted))
full
## observed
obs = c(n = nrow(base::subset(s_threea, R == 1)), ate = mean(onea$observed_unadjusted), sd = sd(onea$observed_unadjusted), bias = mean(onea$observed_bias), sd_bias = sd(onea$observed_bias))
obs

##------------------------------------------------------------------------------
## Case 2: n = 500, SD = 45
##------------------------------------------------------------------------------
## Full analysis
full2 = c(n = nrow(s_threeb), ate = mean(oneb$full_unadjusted), sd = sd(oneb$full_unadjusted), bias = mean(oneb$full_bias_unadjusted), sd_bias = sd(oneb$full_bias_unadjusted))
full2
## observed
obs2 = c(n = nrow(base::subset(s_threeb, R == 1)), ate = mean(oneb$observed_unadjusted), sd = sd(oneb$observed_unadjusted), bias = mean(oneb$observed_bias), sd_bias = sd(oneb$observed_bias))
obs2

##------------------------------------------------------------------------------
## Case 3: n = 2000, SD = 1
##------------------------------------------------------------------------------

full3 = c(n = nrow(s_threec), ate = mean(onec$full_unadjusted), sd = sd(onec$full_unadjusted), bias = mean(onec$full_bias_unadjusted), sd_bias = sd(onec$full_bias_unadjusted))
full3
## observed
obs3 = c(n = nrow(base::subset(s_threec, R == 1)), ate = mean(onec$observed_unadjusted), sd = sd(onec$observed_unadjusted), bias = mean(onec$observed_bias), sd_bias = sd(onec$observed_bias))
obs3


##------------------------------------------------------------------------------
## Case 4: n = 2000, SD = 45
##------------------------------------------------------------------------------

full4 = c(n = nrow(s_threed), ate = mean(oned$full_unadjusted), sd = sd(oned$full_unadjusted), bias = mean(oned$full_bias_unadjusted), sd_bias = sd(oned$full_bias_unadjusted))
full4
## observed
obs4 = c(n = nrow(base::subset(s_threed, R == 1)), ate = mean(oned$observed_unadjusted), sd = sd(oned$observed_unadjusted), bias = mean(oned$observed_bias), sd_bias = sd(oned$observed_bias))
obs4

##------------------------------------------------------------------------------
## create final table of results
sim3_unadjusted = bind_rows(list("n = 500, SD = 1" = full, "n = 500, SD = 1" = obs, 
                                 "n = 500, SD = 45" = full2, "n = 500, SD = 45" = obs2, 
                                 "n = 2000, SD = 1" = full3, "n = 2000, SD = 1" = obs3, 
                                 "n = 2000, SD = 45" = full4, "n = 2000, SD = 45"=  obs4), 
                            .id = "Data generating values") 
kable(sim3_unadjusted, format = "latex", caption = "Unadjusted estimates of the average treatment effect using the correct outcome model across n = 1000 datasets under full and observed data analysis")


## the order of the rows starts with n = 500 
write.csv(sim3_unadjusted, file = "/Users/aokutse/Desktop/PhD/Spring 2024/Bayesian/FinalProject/models/unadjusted/unadjusted.csv", row.names = FALSE)
