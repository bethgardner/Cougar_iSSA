### ### ### ### ### ### ### ### ### ### ###
# Finer scale distance to wolf iSSA model
### ### ### ### ### ### ### ### ### ### ###
#Wirsing et al. in review
#Code by Calum Cunningham; Updated by Beth Gardner
#November 20, 2025
#This script reads in the data and fits the distance to wolf iSSSA model
#The bottom of the script writes the results to a .csv file


#######install packages
pacman::p_load(tidyverse,sp, parallel, mgcv, glmmTMB)


### ### ### ### ### ### ### ### ### ### ###
# Load pre-processed data ----
### ### ### ### ### ### ### ### ### ### ###

coug_data <- readRDS("Data/Wirsingetal_data.rds") %>% data.frame() 

### ### ### ### ### ### ### ### ### ### ###
# Processed data for distance to wolf model
# remove individuals with less than 100 relocations
# Scale continuous covariates
### ### ### ### ### ### ### ### ### ### ###

# 24 hour
coug_data_subset <- coug_data %>% 
  filter(minDistStrata_24hrs_end < 10000) %>% 
  filter(!is.na(dist2wolf_24hrs_start)) %>%
  group_by(ID) %>%  
  mutate(n1 = n()/11) %>%# dividing by 11 because of 1:10 ratio
  ungroup() %>% 
  #only keeping animals with more than 100 used locations
  filter(n1 > 99)%>% 
  mutate(rugged_scale = scale(rugged), rd_dens_600_scale = scale(rd_dens_600), elev_scale = scale(elev), 
         dist2wolf_24hrs_end_scale = scale(dist2wolf_24hrs_end), 
         dist2wolf_24hrs_start_scale = scale(dist2wolf_24hrs_start),
	   wolfUD_scale = wolfUD - mean(wolfUD), forest_scale = forest_pc - mean(forest_pc), 
         human_scale = human_footprint - mean(human_footprint))
 

###Check correlations between main effects
cor(coug_data_subset[,5:12])

#### FIT THE MODEL ############################

set.seed(121)
nt <- 10  #set number of cores
start.time=Sys.time()

ssf_fit_d2wolf <-  glmmTMB(case_ ~ -1 + # conditional logistic regression has no intercept
                                   #main effects
                                   wolfUD_scale + (0 + wolfUD_scale | ID) +  I(wolfUD_scale^2) + (0 + I(wolfUD_scale^2) | ID) +
                                   human_scale + (0 + human_scale | ID) + I(human_scale^2) + (0 + I(human_scale^2) | ID) + 
                                   forest_scale + (0 + forest_scale | ID) + I(forest_scale^2) + (0 + I(forest_scale^2) | ID) + 
                                   elev_scale + (0 + elev_scale | ID) + I(elev_scale^2) + (0 + I(elev_scale^2) | ID) + 
                                   rugged_scale + (0 + rugged_scale | ID) + I(rugged_scale^2) + (0 + I(rugged_scale^2) | ID) + 
                                   rd_dens_600_scale + (0 + rd_dens_600_scale | ID) +I(rd_dens_600_scale^2) + (0 + I(rd_dens_600_scale^2) | ID) +
                                   
					                        # selection for dist2wolf & interactions
                                   dist2wolf_24hrs_end_scale + (0 + dist2wolf_24hrs_end_scale | ID) + I(dist2wolf_24hrs_end_scale^2) + (0 + I(dist2wolf_24hrs_end_scale^2) | ID) +
                           
                                  # control for step length 
                                    log_sl_ + (0 + log_sl_ | ID) + 
                                    log_sl_:dist2wolf_24hrs_end_scale + (0 + log_sl_:dist2wolf_24hrs_end_scale | ID) + 

                                  # used/available strata
                                    (1|step_id_),
                                  family=poisson, doFit=T,
                                  data = coug_data_subset,
                                  # Tell glmmTMB not to change the last standard deviation (for strata); all other values are freely estimated 
                                  map = list(theta = factor(c(1:16, NA))),
                                  # Set the value of the standard deviation of the strata (the last ranef) to large value 
                                  start = list(theta = c(rep(0, times = 16),log(1000))),
                                  control = glmmTMBControl(parallel = nt)) #, optimizer=optim, optArgs=list(method="BFGS"))) 
end.time=Sys.time() 
summary(ssf_fit_d2wolf)
end.time-start.time

write.csv(coef(summary(ssf_fit_d2wolf))$cond, "results.table.distancetowolf.csv")
