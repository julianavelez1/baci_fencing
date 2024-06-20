# 01a_glmm_animals.R: this script runs Bayesian generalized linear mixed models to estimate encounter probabilities of different animal species of the Colombian Orinoquía. Data were collected using a before-after control-impact (BACI) design. Separate models were fit for each species and season (dry and rainy). 

rm(list=ls()) 

library(tidyverse)
library(here)
library(rstan)
library(brms)
library(tidybayes)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

set.seed(537)

# Load data ---------------------------------------------------------------

# Covariates
covariates <- readRDS(here("data", "processed_data", "covariates.RDS"))

# Independent records from camera traps
before_after_recs <- readRDS(here("data", "processed_data", "before_after_recs.RDS"))


# Camera operational days for before and after periods
sum_ops <- readRDS(here("data", "processed_data", "sum_ops.RDS"))


# Format covariates ---------------------------------------------------------------

covariates <- covariates[order(covariates$placename),]
covariates <- replace(covariates, is.na(covariates), 0)

covariates[,3:ncol(covariates)] <- scale(covariates[,3:ncol(covariates)], center = T, scale = T) # center and scale covariates


# Summarise number of detections and operational days of each sampling station --------
  
# days with >= 1 detection

days_detected <- before_after_recs %>%
  mutate(yearjulian = format(date_time, '%Y-%j')) %>%
  group_by(placename, species, before_after, 
           ci, season, patch, yearjulian)  %>%
  summarise(detections = n()) %>% # num of detections per day
  group_by(placename, species, before_after, 
           ci, season, patch) %>% 
  summarise(days_with_detections = n()) %>% # num of days detected
  mutate(patch_temp = paste0(patch, "_", before_after))


# Create non-detection days  ----------------------------------------------

days_det_nodet <- days_detected %>% 
  left_join(sum_ops, by = c("placename", "before_after")) %>% 
  arrange(species, patch, before_after, season) %>% 
  ungroup() %>% 
  mutate(nondetect_days = op_days - days_with_detections)

# Join covariates to sampling stations ------------------------------------

dat <- days_det_nodet %>% 
  left_join(covariates, by = "placename") %>% 
  filter(season %in% "Dry") %>% # Specify "Rainy" or "Dry" season
  mutate(before_after = factor(before_after, levels = c('before', 'after')))

num_recs <- dat %>% 
  group_by(species) %>% 
  summarise(n = n())

op_days_season <- dat %>% 
  group_by(season) %>% 
  summarise(mean(op_days))


# Fit models for each species---------------------------------------------------------------

bform <- bf(days_with_detections | trials(op_days) ~ 0 + 
              Intercept +
              before_after*ci +
              dist_edge_500 + dist_water_500 +
              (1|patch/placename), 
            family = binomial(link = "logit"))

get_prior(bform, data = dat)

prior_vec <- c(set_prior("normal(0,5)", class = "b"))


iter = 150000
warmup = 90000

# Each model takes ~ 15 minutes to run on a MacBook Pro (M1) with 10 cores and 16 Gb of memory.

bos <- dat %>% filter(species %in% "Bos Species")
mod_bos <- brm(bform, data = bos, cores = 4, prior = prior_vec, control = list(adapt_delta = 0.999, max_treedepth = 20), iter = iter, warmup = warmup, thin = 50)

agouti <- dat %>% filter(species %in% "Black Agouti")
mod_agouti <- brm(bform, data = agouti, cores = 4, prior = prior_vec, control = list(adapt_delta = 0.9999, max_treedepth = 20), iter = iter, warmup = warmup, thin = 50) 

collared_pec <- dat %>% filter(species %in% "Collared Peccary")
mod_collared <- brm(bform, data = collared_pec, cores = 4, prior = prior_vec, control = list(adapt_delta = 0.99, max_treedepth = 20), iter = iter, warmup = warmup, thin = 50)

tapir <- dat %>% filter(species %in% "Lowland Tapir")
mod_tapir <- brm(bform, data = tapir, cores = 4, prior = prior_vec, control = list(adapt_delta = 0.99, max_treedepth = 20), iter = iter, warmup = warmup, thin = 50)

coati <- dat %>% filter(species %in% "South American Coati")
mod_coati <- brm(bform, data = coati, cores = 4, prior = prior_vec, control = list(adapt_delta = 0.9999, max_treedepth = 20), iter = iter, warmup = warmup, thin = 50) 

paca <- dat %>% filter(species %in% "Spotted Paca")
mod_paca <- brm(bform, data = paca, cores = 4, prior = prior_vec, control = list(adapt_delta = 0.99, max_treedepth = 20), iter = iter, warmup = warmup, thin = 50)

white_pec <- dat %>% filter(species %in% "White-lipped Peccary")
mod_white_pec <- brm(bform, data = white_pec, cores = 4, prior = prior_vec, control = list(adapt_delta = 0.9999, max_treedepth = 20), iter = iter, warmup = warmup, thin = 50)

deer <- dat %>% filter(species %in% "White-tailed Deer")
mod_deer <- brm(bform, data = deer, cores = 4, prior = prior_vec, control = list(adapt_delta = 0.99, max_treedepth = 20), iter = iter, warmup = warmup, thin = 50) 


mod_sp <- list(mod_agouti,
                  mod_bos,
                  mod_collared,
                  mod_tapir,
                  mod_coati,
                  mod_paca,
                  mod_deer,
                  mod_white_pec)

names(mod_sp) <- c('Black agouti',
                  'Bos sp.',
                  'Collared peccary',
                  'Lowland tapir',
                  'South American coati',
                  'Spotted paca',
                  'White-tailed deer',
                  'White-lipped peccary')


# Save models -------------------------------------------------------------

#saveRDS(dat, here("data", "processed_data", "dat_mods_rainy.RDS"))
#saveRDS(mod_sp, file = "data/output_data/mod_sp_rainy.RDS")

#saveRDS(dat, here("data", "processed_data", "dat_mods_dry.RDS"))
#saveRDS(mod_sp, file = "data/output_data/mod_sp_dry.RDS")


