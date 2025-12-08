### ### ### ### ### ### ### ### ### ### ###
# Season by sex iSSA models
### ### ### ### ### ### ### ### ### ### ###
#Wirsing et al. in review
#Code by Calum Cunningham; Updated by Beth Gardner
#November 20, 2025
#This script reads in the data and fits the sex by season iSSSA models
#The bottom of the script writes the results to 4 .csv files


#######install packages
pacman::p_load(tidyverse,sp, parallel, mgcv, glmmTMB)


### ### ### ### ### ### ### ### ### ### ###
# Load pre-processed data ----
# Scale continuous covariates
### ### ### ### ### ### ### ### ### ### ###

coug_data <- readRDS("Data/Wirsingetal_data.rds") %>% data.frame() %>%
  mutate(rugged_scale = scale(rugged), rd_dens_600_scale = scale(rd_dens_600), elev_scale = scale(elev), 
         wolfUD_scale = wolfUD - mean(wolfUD), forest_scale = forest_pc - mean(forest_pc), 
         human_scale = human_footprint - mean(human_footprint))

###Check correlations between main effects
cor(coug_data[,5:11])
 
##############################################
##### Create season by sex data sets
##############################################

coug_ssf_nonwinter_female <- coug_data %>%
  filter(season == "summer") %>%
  filter(sex == "F")

coug_ssf_winter_female <- coug_data %>%
  filter(season == "winter") %>%
  filter(sex == "F")

coug_ssf_nonwinter_male <- coug_data %>%
  filter(season == "summer") %>%
  filter(sex == "M")

coug_ssf_winter_male <- coug_data %>%
  filter(season == "winter") %>%
  filter(sex == "M")



####FIT THE MODELS ############################

######################################################################
###################   WINTER FEMALE   ################################

set.seed(121)
nt <- 10  #set number of cores
start.time=Sys.time()
# fit conditional logistic regression using glmmTMB following method of Muff et al
ssf_fit_winter_female <-  glmmTMB(case_ ~ -1 + # conditional logistic regression has no intercept
                       wolfUD_scale + (0 + wolfUD_scale | ID) +   
                       I(wolfUD_scale^2) + (0 + I(wolfUD_scale^2) | ID) + 

                       human_scale + (0 + human_scale | ID) + 
			     I(human_scale^2) + (0 + I(human_scale^2) | ID) + 
                       human_scale:wolfUD_scale + (0 + human_scale:wolfUD_scale | ID) +
                       
                       forest_scale + (0 + forest_scale | ID) + 
                       I(forest_scale^2) + (0 + I(forest_scale^2) | ID) + 
                       #forest_scale:wolfUD_scale + (0 + forest_scale:wolfUD_scale | ID) + 
                       forest_scale:human_scale + (0 + forest_scale:human_scale | ID) + 

			     elev_scale + (0 + elev_scale | ID) + 
                       I(elev_scale^2) + (0 + I(elev_scale^2) | ID) + 
                       elev_scale:wolfUD_scale + (0 + elev_scale:wolfUD_scale | ID) + 
                       elev_scale:human_scale + (0 + elev_scale:human_scale | ID) + 

			     rd_dens_600_scale + (0 + rd_dens_600_scale | ID) +
                       I(rd_dens_600_scale^2) + (0 + I(rd_dens_600_scale^2) | ID) +
                       rd_dens_600_scale:wolfUD_scale + (0 + rd_dens_600_scale:wolfUD_scale | ID) + 
                       rd_dens_600_scale:human_scale + (0 + rd_dens_600_scale:human_scale | ID) + 

			     rugged_scale + (0 + rugged_scale | ID) +
			     I(rugged_scale^2) + (0 + I(rugged_scale^2) | ID) +
 			     rugged_scale:wolfUD_scale + (0 + rugged_scale:wolfUD_scale | ID) + 
                       rugged_scale:human_scale + (0 + rugged_scale:human_scale | ID) + 

				# control for step length 
                       log_sl_ + (0 + log_sl_ | ID) +
			     log_sl_:wolfUD_scale + (0 + log_sl_:wolfUD_scale | ID) +
			     log_sl_:human_scale + (0 + log_sl_:human_scale | ID) +

                       # used/available strata
                       (1|step_id_),
                     family=poisson, doFit=T,
                     data = coug_ssf_winter_female, 
                     # Tell glmmTMB not to change the last standard deviation (for strata); all other values are freely estimated 
                     map = list(theta = factor(c(1:23, NA))),
                     # Set the value of the standard deviation of the strata (the last ranef) to large value 
                     start = list(theta = c(rep(0, times = 23),log(1000))),
                     control = glmmTMBControl(parallel = nt, optimizer=optim,
                                 optArgs=list(method="BFGS")))
end.time=Sys.time() 
summary(ssf_fit_winter_female)
end.time-start.time


##########################################################################
###################   NON-WINTER FEMALES  ################################

set.seed(121)
nt <- 10# min(parallel::detectCores()) - 2
start.time=Sys.time()
# fit conditional logistic regression using glmmTMB following method of Muff et al
ssf_fit_nonwinter_female <-  glmmTMB(case_ ~ -1 + # conditional logistic regression has no intercept
                       wolfUD_scale + (0 + wolfUD_scale | ID) +   
                       I(wolfUD_scale^2) + (0 + I(wolfUD_scale^2) | ID) + 

                       human_scale + (0 + human_scale | ID) + 
			     I(human_scale^2) + (0 + I(human_scale^2) | ID) + 
                       human_scale:wolfUD_scale + (0 + human_scale:wolfUD_scale | ID) +
                       
                       forest_scale + (0 + forest_scale | ID) + 
                       I(forest_scale^2) + (0 + I(forest_scale^2) | ID) + 
                       forest_scale:wolfUD_scale + (0 + forest_scale:wolfUD_scale | ID) + 
                       forest_scale:human_scale + (0 + forest_scale:human_scale | ID) + 

			     elev_scale + (0 + elev_scale | ID) + 
                       I(elev_scale^2) + (0 + I(elev_scale^2) | ID) + 
                       elev_scale:wolfUD_scale + (0 + elev_scale:wolfUD_scale | ID) + 
                       elev_scale:human_scale + (0 + elev_scale:human_scale | ID) + 

			     rd_dens_600_scale + (0 + rd_dens_600_scale | ID) +
                       I(rd_dens_600_scale^2) + (0 + I(rd_dens_600_scale^2) | ID) +
                       rd_dens_600_scale:wolfUD_scale + (0 + rd_dens_600_scale:wolfUD_scale | ID) + 
                       rd_dens_600_scale:human_scale + (0 + rd_dens_600_scale:human_scale | ID) + 

			     rugged_scale + (0 + rugged_scale | ID) +
			     I(rugged_scale^2) + (0 + I(rugged_scale^2) | ID) +
 			     rugged_scale:wolfUD_scale + (0 + rugged_scale:wolfUD_scale | ID) + 
                       rugged_scale:human_scale + (0 + rugged_scale:human_scale | ID) + 

				# control for step length 
                       log_sl_ + (0 + log_sl_ | ID) +
			     log_sl_:wolfUD_scale + (0 + log_sl_:wolfUD_scale | ID) +
			     log_sl_:human_scale + (0 + log_sl_:human_scale | ID) +

                       # used/available strata
                       (1|step_id_),
                     family=poisson, doFit=T,
                     data = coug_ssf_nonwinter_female, 
                     # Tell glmmTMB not to change the last standard deviation (for strata); all other values are freely estimated 
                     map = list(theta = factor(c(1:24, NA))),
                     # Set the value of the standard deviation of the strata (the last ranef) to large value 
                     start = list(theta = c(rep(0, times = 24),log(1000))),
                     control = glmmTMBControl(parallel = nt)) #, optimizer=optim,
                                 #optArgs=list(method="BFGS")))
end.time=Sys.time() 
summary(ssf_fit_nonwinter_female)
end.time-start.time




#######################################################################
###################   NONWINTER MALES  ################################

set.seed(121)
nt <- 10# min(parallel::detectCores()) - 2
start.time=Sys.time()
# fit conditional logistic regression using glmmTMB following method of Muff et al
ssf_fit_nonwinter_male <-  glmmTMB(case_ ~ -1 + # conditional logistic regression has no intercept
                       wolfUD_scale + (0 + wolfUD_scale | ID) +   
                       I(wolfUD_scale^2) + (0 + I(wolfUD_scale^2) | ID) + 

                       human_scale + (0 + human_scale | ID) + 
			     I(human_scale^2) + (0 + I(human_scale^2) | ID) + 
                       human_scale:wolfUD_scale + (0 + human_scale:wolfUD_scale | ID) +
                       
                       forest_scale + (0 + forest_scale | ID) + 
                       I(forest_scale^2) + (0 + I(forest_scale^2) | ID) + 
                       forest_scale:wolfUD_scale + (0 + forest_scale:wolfUD_scale | ID) + 
                       forest_scale:human_scale + (0 + forest_scale:human_scale | ID) + 

			     elev_scale + (0 + elev_scale | ID) + 
                       I(elev_scale^2) + (0 + I(elev_scale^2) | ID) + 
                       elev_scale:wolfUD_scale + (0 + elev_scale:wolfUD_scale | ID) + 
                       elev_scale:human_scale + (0 + elev_scale:human_scale | ID) + 

			     rd_dens_600_scale + (0 + rd_dens_600_scale | ID) +
                       I(rd_dens_600_scale^2) + (0 + I(rd_dens_600_scale^2) | ID) +
                       rd_dens_600_scale:wolfUD_scale + (0 + rd_dens_600_scale:wolfUD_scale | ID) + 
                       rd_dens_600_scale:human_scale + (0 + rd_dens_600_scale:human_scale | ID) + 

			     rugged_scale + (0 + rugged_scale | ID) +
			     I(rugged_scale^2) + (0 + I(rugged_scale^2) | ID) +
 			     rugged_scale:wolfUD_scale + (0 + rugged_scale:wolfUD_scale | ID) + 
                       rugged_scale:human_scale + (0 + rugged_scale:human_scale | ID) + 

				# control for step length 
                       log_sl_ + (0 + log_sl_ | ID) +
			     log_sl_:wolfUD_scale + (0 + log_sl_:wolfUD_scale | ID) +
			     log_sl_:human_scale + (0 + log_sl_:human_scale | ID) +

                       # used/available strata
                       (1|step_id_),
                     family=poisson, doFit=T,
                     data = coug_ssf_nonwinter_male, 
                     # Tell glmmTMB not to change the last standard deviation (for strata); all other values are freely estimated 
                     map = list(theta = factor(c(1:24, NA))),
                     # Set the value of the standard deviation of the strata (the last ranef) to large value 
                     start = list(theta = c(rep(0, times = 24),log(1000))),
                     control = glmmTMBControl(parallel = nt, optimizer=optim,
                                 optArgs=list(method="BFGS")))
end.time=Sys.time() 
summary(ssf_fit_nonwinter_male)
end.time-start.time


#######################################################################
###################   WINTER MALE   ###################################

set.seed(121)
start.time=Sys.time()
# fit conditional logistic regression using glmmTMB following method of Muff et al
ssf_fit_winter_male <-  glmmTMB(case_ ~ -1 + # conditional logistic regression has no intercept
                       wolfUD_scale + (0 + wolfUD_scale | ID) +   
                       I(wolfUD_scale^2) + (0 + I(wolfUD_scale^2) | ID) + 

                       human_scale + (0 + human_scale | ID) + 
			     I(human_scale^2) + (0 + I(human_scale^2) | ID) + 
                       human_scale:wolfUD_scale + (0 + human_scale:wolfUD_scale | ID) +
                       
                       forest_scale + (0 + forest_scale | ID) + 
                       I(forest_scale^2) + (0 + I(forest_scale^2) | ID) + 
                       forest_scale:wolfUD_scale + (0 + forest_scale:wolfUD_scale | ID) + 
                       forest_scale:human_scale + (0 + forest_scale:human_scale | ID) + 

			     elev_scale + (0 + elev_scale | ID) + 
                       I(elev_scale^2) + (0 + I(elev_scale^2) | ID) + 
                       elev_scale:wolfUD_scale + (0 + elev_scale:wolfUD_scale | ID) + 
                       elev_scale:human_scale + (0 + elev_scale:human_scale | ID) + 

			     rd_dens_600_scale + (0 + rd_dens_600_scale | ID) +
                       I(rd_dens_600_scale^2) + (0 + I(rd_dens_600_scale^2) | ID) +
                       rd_dens_600_scale:wolfUD_scale + (0 + rd_dens_600_scale:wolfUD_scale | ID) + 
                       rd_dens_600_scale:human_scale + (0 + rd_dens_600_scale:human_scale | ID) + 

			     rugged_scale + (0 + rugged_scale | ID) +
			     I(rugged_scale^2) + (0 + I(rugged_scale^2) | ID) +
 			     rugged_scale:wolfUD_scale + (0 + rugged_scale:wolfUD_scale | ID) + 
                       rugged_scale:human_scale + (0 + rugged_scale:human_scale | ID) + 

				# control for step length 
                       log_sl_ + (0 + log_sl_ | ID) +
			     log_sl_:wolfUD_scale + (0 + log_sl_:wolfUD_scale | ID) +
			     log_sl_:human_scale + (0 + log_sl_:human_scale | ID) +

                       # used/available strata
                       (1|step_id_),
                     family=poisson, doFit=T,
                     data = coug_ssf_winter_male, 
                     # Tell glmmTMB not to change the last standard deviation (for strata); all other values are freely estimated 
                     map = list(theta = factor(c(1:24, NA))),
                     # Set the value of the standard deviation of the strata (the last ranef) to large value 
                     start = list(theta = c(rep(0, times = 24),log(1000))),
                     control = glmmTMBControl(parallel = nt)) #, optimizer=optim, optArgs=list(method="BFGS")))

end.time=Sys.time() 
summary(ssf_fit_winter_male)
end.time-start.time

write.csv(coef(summary(ssf_fit_nonwinter_female))$cond, "results.table.nonwinterfemalewolf.csv")
write.csv(coef(summary(ssf_fit_nonwinter_male))$cond, "results.table.nonwintermale.wolf.csv")
write.csv(coef(summary(ssf_fit_winter_female))$cond, "results.table.winterfemalewolf.csv")
write.csv(coef(summary(ssf_fit_winter_male))$cond, "results.table.wintermalewolf.csv")

































