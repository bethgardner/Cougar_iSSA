# Cougar movements
# Calum Cunningham

### ### ### ### ### ### ### ### ### ### ###
# Reduced model incl wolf proximity ----
### ### ### ### ### ### ### ### ### ### ###

pacman::p_load(tidyverse, sp, leaflet, amt, readr, dplyr, tidyr, purrr, 
               atlastools, spdplyr, data.table, pals, raster, RStoolbox, sf, cowplot, ggmap,
               plotly, future.apply, geosphere, furrr, GGally, MuMIn, suncalc, survival, 
               raster, magick, rphylopic, parallel, mgcv, glmmTMB, gratia, scales)

### ### ### ### ### ### ### ### ### ### ###
# Load pre-processed data ----
### ### ### ### ### ### ### ### ### ### ###

# load 24-h data
# read file with extracted covariates
coug_ssf_dat <- readRDS("coug_ssf_dat_d2wolf3_20230913.rds") %>% data.frame() %>%
  # correctly specify strata separately for each animal ID
  mutate(step_id_ = paste(ID, step_id_, sep = "_")) %>%
  # overwrite log_sl_ because previously incorrect
  mutate(log_sl_ = log(sl_))%>%
  mutate(sex = substring(ID, 7,7)) %>%
  # find closest location within strata
  group_by(step_id_) %>% 
  mutate(minDistStrata_4hrs_end = min(dist2wolf_4hrs_end), minDistStrata_4hrs_start = min(dist2wolf_4hrs_start),
         minDistStrata_1.2hrs_end = min(dist2wolf_1.2hrs_end), minDistStrata_1.2hrs_start = min(dist2wolf_1.2hrs_start),
         minDistStrata_24hrs_end = min(dist2wolf_24hrs_end), minDistStrata_24hrs_start = min(dist2wolf_24hrs_start)) %>%
  ungroup()

# add sample size for each animal so that we can omit animals with few records and limit cougar relocations within WA
coug_ssf_dat1 <- coug_ssf_dat %>% group_by(ID) %>% 
  mutate(n = n()/11) %>%# dividing by 11 because of 1:10 ratio
  ungroup() %>% 
  filter (y1_ < 49)



### ### ### ### ### ### ### ### ### ### ###
# Processed data for distance to wolf model
### ### ### ### ### ### ### ### ### ### ###

# 24 hour
coug_ssf_dat1_subset_24hrs <- coug_ssf_dat1 %>% 
  filter(minDistStrata_24hrs_end < 10000) %>% 
  filter(!is.na(dist2wolf_24hrs_start)) %>%
  group_by(ID) %>% 
  mutate(n1 = n()/11) %>%# dividing by 11 because of 1:10 ratio
  ungroup() %>% 
  # only keeping animals with more than 100 used locations (change if desired)
  filter(n1 > 99)%>% 
  mutate(rugged_scale = scale(rugged), rd_dens_600_scale = scale(rd_dens_600), elev_scale = scale(elev), 
         dist2wolf_24hrs_end_scale = scale(dist2wolf_24hrs_end), 
         dist2wolf_24hrs_start_scale = scale(dist2wolf_24hrs_start),
	   wolfUD_c = wolfUD - mean(wolfUD), forest = forest_pc - mean(forest_pc), human_footprint = human_footprint - mean(human_footprint)) 


table(coug_ssf_dat1_subset_24hrs$ID)/11

cor(coug_ssf_dat1_subset_24hrs[,c(70:73, 76, 77, 22, 20)]


### ### ### ### ### ### ### ### ### ### ###
# Fit the model
### ### ### ### ### ### ### ### ### ### ###


nt <- 10

ssf_fit_d2wolf <-  glmmTMB(case_ ~ -1 + # conditional logistic regression has no intercept
                                   #main effects
                                   wolfUD_c + (0 + wolfUD_c | ID) +  I(wolfUD_c^2) + (0 + I(wolfUD_c^2) | ID) +
                                   human_footprint + (0 + human_footprint | ID) + I(human_footprint^2) + (0 + I(human_footprint^2) | ID) + 
                                   forest + (0 + forest | ID) + I(forest^2) + (0 + I(forest^2) | ID) + 
                                   elev_scale + (0 + elev_scale | ID) + I(elev_scale^2) + (0 + I(elev_scale^2) | ID) + 
                                   rugged_scale + (0 + rugged_scale | ID) + I(rugged_scale^2) + (0 + I(rugged_scale^2) | ID) + 
                                   rd_dens_600_scale + (0 + rd_dens_600_scale | ID) +I(rd_dens_600_scale^2) + (0 + I(rd_dens_600_scale^2) | ID) +
                                   
					    # selection for dist2wolf & interactions
                                   dist2wolf_24hrs_end_scale + (0 + dist2wolf_24hrs_end_scale | ID) + I(dist2wolf_24hrs_end_scale^2) + (0 + I(dist2wolf_24hrs_end_scale^2) | ID) +
                                   #forest_pc:dist2wolf_24hrs_end_scale + (0 + forest_pc:dist2wolf_24hrs_end_scale | ID) +
                                   #elev_scale:dist2wolf_24hrs_end_scale + (0 + elev_scale:dist2wolf_24hrs_end_scale | ID) +
                                   #rugged_scale:dist2wolf_24hrs_end_scale + (0 + rugged_scale:dist2wolf_24hrs_end_scale | ID) +

                                  # control for step length 
                                    log_sl_ + (0 + log_sl_ | ID) + 
                                    log_sl_:dist2wolf_24hrs_end_scale + (0 + log_sl_:dist2wolf_24hrs_end_scale | ID) + 

                                  # used/available strata
                                    (1|step_id_),
                                  family=poisson, doFit=T,
                                  data = coug_ssf_dat1_subset_24hrs,
                                  # Tell glmmTMB not to change the last standard deviation (for strata); all other values are freely estimated 
                                  map = list(theta = factor(c(1:16, NA))),
                                  # Set the value of the standard deviation of the strata (the last ranef) to large value 
                                  start = list(theta = c(rep(0, times = 16),log(1000))),
                                  control = glmmTMBControl(parallel = nt)) #, optimizer=optim, optArgs=list(method="BFGS"))) 
summary(ssf_fit_d2wolf)



write.csv(coef(summary(ssf_fit_d2wolf))$cond, "results.table.distancetowolf.csv")

#################################################################################################################################

### ### ### ### ### ### ### ### ### ### ###
# Create plots
### ### ### ### ### ### ### ### ### ### ###

source("Functions for cougar movement analysis V1.R")

# plot RSS in response to dist2wolf
d2wolf_24_dat <- find_observed_range_singleVar(fitted_data = coug_ssf_dat1_subset_24hrs, xvar = "dist2wolf_24hrs_end_scale", 
                                               xvar_point_of_comparison = mean, model = ssf_fit_d2wolf, predatorSp = "wolf") 

d2wolf_24_RSS <- calc_logRSS_glmmTMB(x1 = d2wolf_24_dat$x1, x2 = d2wolf_24_dat$x2, model = ssf_fit_d2wolf, ci_level = 0.95, model_type = "SSF")

# put d2wolf back on original scale
d2wolf_24_RSS$dist2wolf_24hrs_end <- d2wolf_24_RSS$dist2wolf_24hrs_end_scale * attributes(coug_ssf_dat1_subset_24hrs$dist2wolf_24hrs_end_scale)$`scaled:scale` + attributes(coug_ssf_dat1_subset_24hrs$dist2wolf_24hrs_end_scale)$`scaled:center`


max(d2wolf_24_RSS$logRSS)

# plot RSS
d2wolfplot<-ggplot(d2wolf_24_RSS, aes(dist2wolf_24hrs_end, exp(logRSS))) +
  theme_minimal() + ylim(.65, 1.25) + 
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr)), alpha = 0.3, colour = NA) +
  geom_hline(yintercept = 1, linetype = "dashed")  +
  labs(x = "Distance to wolf (m)", y = "Relative selection strength") +
  #ylim(c(0, exp(max(d2wolf_24_RSS$upr))))+
  #annotate("text", x = 2500, y = 1.2, label = "Males") +
  theme(plot.margin = margin(0.25,0.5,0.25,0.25, "cm"))



d2wolfplot





# plot RSS in response to forest
step_24_dat_near <- find_observed_range_singleVar(fitted_data = coug_ssf_dat1_subset_24hrs, xvar = "log_sl_", 
                                               xvar_point_of_comparison = mean, model = ssf_fit_d2wolf, predatorSp = "NA") 
step_24_dat_far <- step_24_dat_near


step_24_dat_near$x1$dist2wolf_24hrs_end_scale<-rep(quantile(coug_ssf_dat1_subset_24hrs$dist2wolf_24hrs_end_scale, .1),101)
step_24_dat_near$x2$dist2wolf_24hrs_end_scale<-quantile(coug_ssf_dat1_subset_24hrs$dist2wolf_24hrs_end_scale, .1)

step_24_dat_far$x1$dist2wolf_24hrs_end_scale<-rep(quantile(coug_ssf_dat1_subset_24hrs$dist2wolf_24hrs_end_scale, .9),101)
step_24_dat_far$x2$dist2wolf_24hrs_end_scale<-quantile(coug_ssf_dat1_subset_24hrs$dist2wolf_24hrs_end_scale, .9)

step_24_near <- calc_logRSS_glmmTMB(x1 = step_24_dat_near$x1, x2 =step_24_dat_near$x2, model = ssf_fit_d2wolf, ci_level = 0.95, model_type = "SSF")
step_24_far<- calc_logRSS_glmmTMB(x1 = step_24_dat_far$x1, x2 = step_24_dat_far$x2, model = ssf_fit_d2wolf, ci_level = 0.95, model_type = "SSF")



step_24_near$DistanceToWolf<-"Closer"
step_24_far$DistanceToWolf<-"Farther"

log_sl_RSS<-rbind(step_24_near, step_24_far)

wolf_colors <- c("Closer" = "#D55E00", "Farther" = "#0072B2")


# plot RSS
log_sl_plot_int<-ggplot(log_sl_RSS, aes(log_sl_, exp(logRSS), group=DistanceToWolf, colour=DistanceToWolf)) +
  theme_minimal() + ylim(0.65, 1.25) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr), fill=DistanceToWolf), alpha = 0.3,, colour = NA) +
  geom_hline(yintercept = 1, linetype = "dashed")  +
    labs(x = "log(step length)", y = "Relative selection strength")+
  scale_colour_manual(values = wolf_colors) +
  scale_fill_manual(values = wolf_colors)




library(patchwork)

d2wolfplot + log_sl_plot_int + plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(face = "bold", size = 14))




######################

# plot RSS in response to wolfUD
wolfUD_24_dat <- find_observed_range_singleVar(fitted_data = coug_ssf_dat1_subset_24hrs, xvar = "wolfUD", 
                                               xvar_point_of_comparison = mean, model = ssf_fit_d2wolf_24hrs, predatorSp = "wolf") 
wolfUD_24_RSS <- calc_logRSS_glmmTMB(x1 = wolfUD_24_dat$x1, x2 = wolfUD_24_dat$x2, model = ssf_fit_d2wolf_24hrs, ci_level = 0.95, model_type = "SSF")

# plot RSS
ggplot(wolfUD_24_RSS, aes(wolfUD, exp(logRSS))) +
  theme_minimal() +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr)), alpha = 0.3, colour = NA) +
  geom_hline(yintercept = 1, linetype = "dashed")  +
  labs(x = "wolf UD", y = "Relative selection strength") +
  ylim(c(0, exp(max(d2wolf_24_RSS$upr))))

x11()
# plot RSS in response to HFI
human_24_dat <- find_observed_range_singleVar(fitted_data = coug_ssf_dat1_subset_24hrs, xvar = "human_footprint", 
                                               xvar_point_of_comparison = min, model = ssf_fit_d2wolf_24hrs, predatorSp = "wolf") 
human_24_RSS <- calc_logRSS_glmmTMB(x1 = human_24_dat$x1, x2 = human_24_dat$x2, model = ssf_fit_d2wolf_24hrs, ci_level = 0.95, model_type = "SSF")

# plot RSS
ggplot(human_24_RSS, aes(human_footprint, exp(logRSS))) +
  theme_minimal() +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr)), alpha = 0.3, colour = NA) +
  geom_hline(yintercept = 1, linetype = "dashed")  +
  labs(x = "HFI", y = "Relative selection strength") +
  ylim(c(0, exp(max(d2wolf_24_RSS$upr))))


# plot RSS in response to forest
forest_24_dat <- find_observed_range_singleVar(fitted_data = coug_ssf_dat1_subset_24hrs, xvar = "forest_pc", 
                                               xvar_point_of_comparison = mean, model = ssf_fit_d2wolf_24hrs, predatorSp = "NA") 
forest_24_dat_far <- forest_24_dat
forest_24_dat_far$x1$dist2wolf_24hrs_end_scale <- max(coug_ssf_dat1_subset_24hrs$dist2wolf_24hrs_end_scale)
forest_24_dat_far$x2$dist2wolf_24hrs_end_scale <- max(coug_ssf_dat1_subset_24hrs$dist2wolf_24hrs_end_scale)
forest_24_dat_near <- forest_24_dat
forest_24_dat_near$x1$dist2wolf_24hrs_end_scale <- min(coug_ssf_dat1_subset_24hrs$dist2wolf_24hrs_end_scale)
forest_24_dat_near$x2$dist2wolf_24hrs_end_scale <- min(coug_ssf_dat1_subset_24hrs$dist2wolf_24hrs_end_scale)

forest_24_RSS_far <- calc_logRSS_glmmTMB(x1 = forest_24_dat_far$x1, x2 = forest_24_dat_far$x2, model = ssf_fit_d2wolf_24hrs, ci_level = 0.95, model_type = "SSF")
forest_24_RSS_near <- calc_logRSS_glmmTMB(x1 = forest_24_dat_near$x1, x2 = forest_24_dat_near$x2, model = ssf_fit_d2wolf_24hrs, ci_level = 0.95, model_type = "SSF")

# plot RSS
ggplot() +
  theme_minimal() +
  geom_line(data = forest_24_RSS_far, aes(forest_pc, exp(logRSS)), linewidth = 1, colour = "blue") +
  geom_line(data = forest_24_RSS_near, aes(forest_pc, exp(logRSS)), linewidth = 1, colour = "red") +
#  geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr)), alpha = 0.3, colour = NA) +
  geom_hline(yintercept = 1, linetype = "dashed")  


# ns
# plot RSS in response to human footprint
d2wolf_24_dat$x1$fit <- predict(ssf_fit_d2wolf_24hrs, d2wolf_24_dat$x1, se.fit = T, re.form = NA)[[1]]
d2wolf_24_dat$x1$se <- predict(ssf_fit_d2wolf_24hrs, d2wolf_24_dat$x1, se.fit = T, re.form = NA)[[2]]
# put d2wolf back on original scale
d2wolf_24_dat$x1$dist2wolf_24hrs_end <- d2wolf_24_dat$x1$dist2wolf_24hrs_end_scale * attributes(coug_ssf_dat1_subset_24hrs$dist2wolf_24hrs_end_scale)$`scaled:scale` + attributes(coug_ssf_dat1_subset_24hrs$dist2wolf_24hrs_end_scale)$`scaled:center`

# plot RSS
ggplot(d2wolf_24_dat$x1, aes(dist2wolf_24hrs_end, exp(fit))) +
  theme_minimal() +
  geom_line(linewidth = 1) +
#  geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr)), alpha = 0.3, colour = NA) +
  geom_hline(yintercept = 1, linetype = "dashed")  +
  labs(x = "Distance to wolf", y = "Relative selection strength")




##############################################################################################
##############################################################################################
##############################################################################################
##  BY SEX ANALYSIS


### ### ### ### ### ### ### ### ### ### ###
# Processed data for distance to wolf model
### ### ### ### ### ### ### ### ### ### ###

# 24 hour
coug_ssf_dat1_subset_24hrs <- coug_ssf_dat1 %>% 
  filter(minDistStrata_24hrs_end < 10000) %>% 
  # scale d2wolf
  mutate(dist2wolf_24hrs_end_scale = scale(dist2wolf_24hrs_end), 
         dist2wolf_24hrs_start_scale = scale(dist2wolf_24hrs_start)) %>%
  filter(!is.na(dist2wolf_24hrs_start)) 

cor(coug_ssf_dat1_subset_24hrs$wolfUD, coug_ssf_dat1_subset_24hrs$dist2wolf_24hrs_end)
hist(coug_ssf_dat1_subset_24hrs$dist2wolf_24hrs_end)


coug_ssf_24hr_male <- coug_ssf_dat1 %>%
    filter(sex == "M")%>%
filter(minDistStrata_24hrs_end < 10000) %>% 
  # scale d2wolf
  mutate(dist2wolf_24hrs_end_scale = scale(dist2wolf_24hrs_end), 
         dist2wolf_24hrs_start_scale = scale(dist2wolf_24hrs_start)) %>%
  filter(!is.na(dist2wolf_24hrs_start)) %>%
  mutate(rugged_scale = scale(rugged), rd_dens_600_scale = scale(rd_dens_600), elev_scale = scale(elev))


coug_ssf_24hr_female <- coug_ssf_dat1_subset_24hrs %>%
  filter(sex == "F")%>%
  filter(minDistStrata_24hrs_end < 10000) %>% 
  # scale d2wolf
  mutate(dist2wolf_24hrs_end_scale = scale(dist2wolf_24hrs_end), 
         dist2wolf_24hrs_start_scale = scale(dist2wolf_24hrs_start)) %>%
  filter(!is.na(dist2wolf_24hrs_start)) %>%
  mutate(rugged_scale = scale(rugged), rd_dens_600_scale = scale(rd_dens_600), elev_scale = scale(elev))




### ### ### ### ### ### ### ### ### ### ###
# Fit the model
### ### ### ### ### ### ### ### ### ### ###


nt <- 10

ssf_fit_d2wolf_f <-  glmmTMB(case_ ~ -1 + # conditional logistic regression has no intercept
                                   #main effects
                                   wolfUD + (0 + wolfUD | ID) +  I(wolfUD^2) + (0 + I(wolfUD^2) | ID) +
                                   human_footprint + (0 + human_footprint | ID) + I(human_footprint^2) + (0 + I(human_footprint^2) | ID) + 
                                   forest_pc + (0 + forest_pc | ID) + I(forest_pc^2) + (0 + I(forest_pc^2) | ID) + 
                                   elev_scale + (0 + elev_scale | ID) + I(elev_scale^2) + (0 + I(elev_scale^2) | ID) + 
                                   rugged_scale + (0 + rugged_scale | ID) + I(rugged_scale^2) + (0 + I(rugged_scale^2) | ID) + 
                                   rd_dens_600_scale + (0 + rd_dens_600_scale | ID) +
                                   
					    # selection for dist2wolf & interactions
                                   dist2wolf_24hrs_end_scale + (0 + dist2wolf_24hrs_end_scale | ID) + I(dist2wolf_24hrs_end_scale^2) + (0 + I(dist2wolf_24hrs_end_scale^2) | ID) +
                                   forest_pc:dist2wolf_24hrs_end_scale + (0 + forest_pc:dist2wolf_24hrs_end_scale | ID) +
                                  # human_footprint:dist2wolf_24hrs_end_scale + (0 + human_footprint:dist2wolf_24hrs_end_scale | ID) +
                                   rugged_scale:dist2wolf_24hrs_end_scale + (0 + rugged_scale:dist2wolf_24hrs_end_scale | ID) +

                                  # control for step length 
                                    log_sl_ + (0 + log_sl_ | ID) + 
                                    log_sl_:dist2wolf_24hrs_end_scale + (0 + log_sl_:dist2wolf_24hrs_end_scale | ID) + 

                                  # used/available strata
                                    (1|step_id_),
                                  family=poisson, doFit=T,
                                  data = coug_ssf_24hr_female,
                                  # Tell glmmTMB not to change the last standard deviation (for strata); all other values are freely estimated 
                                  map = list(theta = factor(c(1:17, NA))),
                                  # Set the value of the standard deviation of the strata (the last ranef) to large value 
                                  start = list(theta = c(rep(0, times = 17),log(1000))),
                                  control = glmmTMBControl(parallel = nt)) 
summary(ssf_fit_d2wolf_f)





ssf_fit_d2wolf_m <-  glmmTMB(case_ ~ -1 + # conditional logistic regression has no intercept
                                   #main effects
                                   wolfUD + (0 + wolfUD | ID) +  I(wolfUD^2) + (0 + I(wolfUD^2) | ID) +
                                   human_footprint + (0 + human_footprint | ID) + I(human_footprint^2) + (0 + I(human_footprint^2) | ID) + 
                                   forest_pc + (0 + forest_pc | ID) + I(forest_pc^2) + (0 + I(forest_pc^2) | ID) + 
                                   elev_scale + (0 + elev_scale | ID) + I(elev_scale^2) + (0 + I(elev_scale^2) | ID) + 
                                   rugged_scale + (0 + rugged_scale | ID) + I(rugged_scale^2) + (0 + I(rugged_scale^2) | ID) + 
                                   rd_dens_600_scale + (0 + rd_dens_600_scale | ID) +
                                   
					    # selection for dist2wolf & interactions
                                   dist2wolf_24hrs_end_scale + (0 + dist2wolf_24hrs_end_scale | ID) + I(dist2wolf_24hrs_end_scale^2) + (0 + I(dist2wolf_24hrs_end_scale^2) | ID) +
                                   forest_pc:dist2wolf_24hrs_end_scale + (0 + forest_pc:dist2wolf_24hrs_end_scale | ID) +
                                  # human_footprint:dist2wolf_24hrs_end_scale + (0 + human_footprint:dist2wolf_24hrs_end_scale | ID) +
                                   rugged_scale:dist2wolf_24hrs_end_scale + (0 + rugged_scale:dist2wolf_24hrs_end_scale | ID) +

                                  # control for step length 
                                    log_sl_ + (0 + log_sl_ | ID) + 
                                    log_sl_:dist2wolf_24hrs_end_scale + (0 + log_sl_:dist2wolf_24hrs_end_scale | ID) + 

                                  # used/available strata
                                    (1|step_id_),
                                  family=poisson, doFit=T,
                                  data = coug_ssf_24hr_male, 
                                  # Tell glmmTMB not to change the last standard deviation (for strata); all other values are freely estimated 
                                  map = list(theta = factor(c(1:17, NA))),
                                  # Set the value of the standard deviation of the strata (the last ranef) to large value 
                                  start = list(theta = c(rep(0, times = 17),log(1000))),
                                 control = glmmTMBControl(parallel = nt, optimizer=optim,
                                 optArgs=list(method="BFGS"))) 
summary(ssf_fit_d2wolf_m)




#################################################################################################################################

### ### ### ### ### ### ### ### ### ### ###
# Create plots
### ### ### ### ### ### ### ### ### ### ###

source("Functions for cougar movement analysis V1.R")

# plot RSS in response to dist2wolf
d2wolf_24_dat <- find_observed_range_singleVar(fitted_data = coug_ssf_24hr_male, xvar = "dist2wolf_24hrs_end_scale", 
                                               xvar_point_of_comparison = mean, model = ssf_fit_d2wolf_m, predatorSp = "wolf") 
d2wolf_24_RSS <- calc_logRSS_glmmTMB(x1 = d2wolf_24_dat$x1, x2 = d2wolf_24_dat$x2, model = ssf_fit_d2wolf_m, ci_level = 0.95, model_type = "SSF")
# put d2wolf back on original scale
d2wolf_24_RSS$dist2wolf_24hrs_end <- d2wolf_24_RSS$dist2wolf_24hrs_end_scale * attributes(coug_ssf_24hr_male$dist2wolf_24hrs_end_scale)$`scaled:scale` + attributes(coug_ssf_24hr_male$dist2wolf_24hrs_end_scale)$`scaled:center`


# plot RSS
male24<-ggplot(d2wolf_24_RSS, aes(dist2wolf_24hrs_end, exp(logRSS))) +
  theme_minimal() + ylim(.5, 1.35) + 
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr)), alpha = 0.3, colour = NA) +
  geom_hline(yintercept = 1, linetype = "dashed")  +
  labs(x = "Distance to wolf (m)", y = "Relative selection strength") +
  #ylim(c(0, exp(max(d2wolf_24_RSS$upr))))+
  annotate("text", x = 2500, y = 1.2, label = "Males") +
  theme(plot.margin = margin(0.25,0.5,0.25,0.25, "cm"))




# plot RSS in response to dist2wolf
d2wolf_24_dat <- find_observed_range_singleVar(fitted_data = coug_ssf_24hr_female, xvar = "dist2wolf_24hrs_end_scale", 
                                               xvar_point_of_comparison = mean, model = ssf_fit_d2wolf_f, predatorSp = "wolf") 
d2wolf_24_RSS <- calc_logRSS_glmmTMB(x1 = d2wolf_24_dat$x1, x2 = d2wolf_24_dat$x2, model = ssf_fit_d2wolf_f, ci_level = 0.95, model_type = "SSF")
# put d2wolf back on original scale
d2wolf_24_RSS$dist2wolf_24hrs_end <- d2wolf_24_RSS$dist2wolf_24hrs_end_scale * attributes(coug_ssf_24hr_female$dist2wolf_24hrs_end_scale)$`scaled:scale` + attributes(coug_ssf_24hr_female$dist2wolf_24hrs_end_scale)$`scaled:center`

# plot RSS
female24<-ggplot(d2wolf_24_RSS, aes(dist2wolf_24hrs_end, exp(logRSS))) +
  theme_minimal() + ylim(.5, 1.35) + 
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr)), alpha = 0.3, colour = NA) +
  geom_hline(yintercept = 1, linetype = "dashed")  +
  labs(x = "Distance to wolf (m)", y = "Relative selection strength") +
  #ylim(c(0, exp(max(d2wolf_24_RSS$upr))))+
  annotate("text", x = 2500, y = 1.2, label = "Females") +
  theme(plot.margin = margin(0.25,0.5,0.25,0.25, "cm"))


female24 + male24

# plot RSS in response to wolfUD
wolfUD_24_dat <- find_observed_range_singleVar(fitted_data = coug_ssf_dat1_subset_24hrs, xvar = "wolfUD", 
                                               xvar_point_of_comparison = mean, model = ssf_fit_d2wolf_24hrs, predatorSp = "wolf") 
wolfUD_24_RSS <- calc_logRSS_glmmTMB(x1 = wolfUD_24_dat$x1, x2 = wolfUD_24_dat$x2, model = ssf_fit_d2wolf_24hrs, ci_level = 0.95, model_type = "SSF")

# plot RSS
ggplot(wolfUD_24_RSS, aes(wolfUD, exp(logRSS))) +
  theme_minimal() +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr)), alpha = 0.3, colour = NA) +
  geom_hline(yintercept = 1, linetype = "dashed")  +
  labs(x = "wolf UD", y = "Relative selection strength") +
  ylim(c(0, exp(max(d2wolf_24_RSS$upr))))

x11()
# plot RSS in response to HFI
human_24_dat <- find_observed_range_singleVar(fitted_data = coug_ssf_dat1_subset_24hrs, xvar = "human_footprint", 
                                               xvar_point_of_comparison = min, model = ssf_fit_d2wolf_24hrs, predatorSp = "wolf") 
human_24_RSS <- calc_logRSS_glmmTMB(x1 = human_24_dat$x1, x2 = human_24_dat$x2, model = ssf_fit_d2wolf_24hrs, ci_level = 0.95, model_type = "SSF")

# plot RSS
ggplot(human_24_RSS, aes(human_footprint, exp(logRSS))) +
  theme_minimal() +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr)), alpha = 0.3, colour = NA) +
  geom_hline(yintercept = 1, linetype = "dashed")  +
  labs(x = "HFI", y = "Relative selection strength") +
  ylim(c(0, exp(max(d2wolf_24_RSS$upr))))


# plot RSS in response to forest
forest_24_dat <- find_observed_range_singleVar(fitted_data = coug_ssf_dat1_subset_24hrs, xvar = "forest_pc", 
                                               xvar_point_of_comparison = mean, model = ssf_fit_d2wolf_24hrs, predatorSp = "NA") 
forest_24_dat_far <- forest_24_dat
forest_24_dat_far$x1$dist2wolf_24hrs_end_scale <- max(coug_ssf_dat1_subset_24hrs$dist2wolf_24hrs_end_scale)
forest_24_dat_far$x2$dist2wolf_24hrs_end_scale <- max(coug_ssf_dat1_subset_24hrs$dist2wolf_24hrs_end_scale)
forest_24_dat_near <- forest_24_dat
forest_24_dat_near$x1$dist2wolf_24hrs_end_scale <- min(coug_ssf_dat1_subset_24hrs$dist2wolf_24hrs_end_scale)
forest_24_dat_near$x2$dist2wolf_24hrs_end_scale <- min(coug_ssf_dat1_subset_24hrs$dist2wolf_24hrs_end_scale)

forest_24_RSS_far <- calc_logRSS_glmmTMB(x1 = forest_24_dat_far$x1, x2 = forest_24_dat_far$x2, model = ssf_fit_d2wolf_24hrs, ci_level = 0.95, model_type = "SSF")
forest_24_RSS_near <- calc_logRSS_glmmTMB(x1 = forest_24_dat_near$x1, x2 = forest_24_dat_near$x2, model = ssf_fit_d2wolf_24hrs, ci_level = 0.95, model_type = "SSF")

# plot RSS
ggplot() +
  theme_minimal() +
  geom_line(data = forest_24_RSS_far, aes(forest_pc, exp(logRSS)), linewidth = 1, colour = "blue") +
  geom_line(data = forest_24_RSS_near, aes(forest_pc, exp(logRSS)), linewidth = 1, colour = "red") +
#  geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr)), alpha = 0.3, colour = NA) +
  geom_hline(yintercept = 1, linetype = "dashed")  


# ns
# plot RSS in response to human footprint
d2wolf_24_dat$x1$fit <- predict(ssf_fit_d2wolf_24hrs, d2wolf_24_dat$x1, se.fit = T, re.form = NA)[[1]]
d2wolf_24_dat$x1$se <- predict(ssf_fit_d2wolf_24hrs, d2wolf_24_dat$x1, se.fit = T, re.form = NA)[[2]]
# put d2wolf back on original scale
d2wolf_24_dat$x1$dist2wolf_24hrs_end <- d2wolf_24_dat$x1$dist2wolf_24hrs_end_scale * attributes(coug_ssf_dat1_subset_24hrs$dist2wolf_24hrs_end_scale)$`scaled:scale` + attributes(coug_ssf_dat1_subset_24hrs$dist2wolf_24hrs_end_scale)$`scaled:center`

# plot RSS
ggplot(d2wolf_24_dat$x1, aes(dist2wolf_24hrs_end, exp(fit))) +
  theme_minimal() +
  geom_line(linewidth = 1) +
#  geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr)), alpha = 0.3, colour = NA) +
  geom_hline(yintercept = 1, linetype = "dashed")  +
  labs(x = "Distance to wolf", y = "Relative selection strength")



