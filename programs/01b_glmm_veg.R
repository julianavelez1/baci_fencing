# 01b_glmm_veg.R: this script fits a beta regression model to measured areal ground cover proportion in 4-m2 plots located within 5 m of each camera-trap station. Data collection was performed using a before-after control-impact (BACI) design. This script also plots posterior median Time x Treatment coefficient and posterior distributions of predicted areal cover proportion.

rm(list=ls()) 

library(tidyverse)
library(here)
library(rstan)
library(brms)
library(tidybayes)
library(scales) # labels percent
library(ggthemes) # themes_clean()


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

set.seed(537)

# Load data ---------------------------------------------------------------

# Covariates
covariates <- readRDS(here("data", "processed_data", "covariates.RDS")) 

# Get control-impact sites

ci <- readRDS(here("data", "processed_data", "before_after_recs.RDS"))

ci <- ci %>% ungroup() %>% dplyr::select(placename, ci, patch) %>% 
  distinct(placename, .keep_all = TRUE)

# Format covariates ---------------------------------------------------------------
covariates <- covariates[order(covariates$placename),]
covariates <- covariates %>% relocate(geometry, .after = placename)
covariates <- replace(covariates, is.na(covariates), 0)

covariates[,3:10] <- scale(covariates[,3:10], center = T, scale = T) # center and scale covariates

filt_covs <- covariates %>% 
  left_join(ci, by = "placename") %>% 
  # get two sampling periods in same column
  pivot_longer(cols = c(all_shrubs_t1, all_shrubs_t2), 
               names_to = "before_after", values_to = "all_shrubs") %>% 
  mutate_if(is.character, as.factor) %>% 
  # to proportion
  mutate(all_shrubs = round(all_shrubs, 0),
         all_shrubs = ifelse(all_shrubs == 0, 0.001, all_shrubs),
         all_shrubs = all_shrubs/100)


# Explore data ------------------------------------------------------------

shrubs_densities <- ggplot(filt_covs, aes(x = all_shrubs, fill = before_after)) +
  geom_density(alpha = 0.6) +
  scale_x_continuous(labels = label_percent(accuracy = 1, scale = 100)) +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  labs(x = "Shrub percentage", y = "Density", fill = "Period") +
  theme_clean() +
  theme(legend.position = "bottom")


# Fit models---------------------------------------------------------------

bform <- bf(
  # mu (mean part) - logit scale
  all_shrubs ~ before_after*ci +
              (1|patch/placename), 
  # phi (precision part) - log scale
  phi ~ 1
)

get_pri <- get_prior(bform, 
                     data = filt_covs, 
                     family = Beta()
)

priors <- c(set_prior("student_t(3, 0, 2.5)", class = "Intercept"),
            set_prior("normal(0, 1)", class = "b"))

iter = 150000
warmup = 90000

veg_model <- brm(bform,
                 data = filt_covs,
                 family = Beta(),
                 cores = 4,
                 seed = 7337,
                 prior = priors,
                 init = 0,
                 control = list(adapt_delta = 0.999999,
                                max_treedepth = 20),
                 iter = iter, warmup = warmup,
                 thin = 50, chains = 4)

# Save model ---------------------------------------------------------------
#saveRDS(veg_model, file = here("data", "output_data", "mod_beta_mx.RDS"))

# Read model ---------------------------------------------

veg_model <- readRDS(file = here("data", "output_data", "mod_beta_mx.RDS"))

# Shrub proportion - predictions --------------------------------------------------------
new_data <- filt_covs %>% 
  group_by(before_after, ci) %>% 
  summarise()

# Predictions

beta_linpred <- veg_model %>% 
  linpred_draws(newdata = new_data, 
                re_formula = NA,
                transform = TRUE) # mu in the model in probability or proportion values

p_beta_linpred <- ggplot(beta_linpred, 
                         aes(x = .linpred, fill = ci)) +
  stat_halfeye(.width = c(0.5, 0.89), point_interval = "median_hdi",
               slab_alpha = 0.75) +
  scale_fill_viridis_d(option = "viridis", end = 0.6, name = NULL, labels = c("Control", "Fence"), direction = -1) +
  labs(x = "Proportion of ground cover", 
       y = "Density") +
  theme_clean() +
  theme(legend.position = "bottom") + 
  facet_wrap(~before_after, labeller = labeller(before_after = c("all_shrubs_t1" = "Before", "all_shrubs_t2" = "After"))) + 
  theme(plot.background = element_blank())



# Beta coefficients -------------------------------------------------------

model_fixed <- broom.mixed::tidy(veg_model, effects = "fixed", conf.level = 0.89)

coef_draws <- veg_model %>% 
             as_draws_df() %>% 
  rename(baci = `b_before_afterall_shrubs_t2:ciimpact`) %>% 
  median_qi(baci, .width = 0.89)


# Transform the veg_model object to extract and organize posterior betas
posterior_beta <- veg_model %>%
  # Extract draws for variables matching 'b_.*'
  gather_draws(`b_.*`, regex = TRUE) %>%
  
  # Add a new column 'component' to categorize variables as 'Precision' or 'Mean'
  mutate(component = case_when(
    str_detect(.variable, "phi_") ~ "Precision", # If variable name contains 'phi_', it's 'Precision'
    TRUE ~ "Mean"
  )) %>%
  
  # Add a new column 'intercept'
  mutate(intercept = str_detect(.variable, "Intercept")) %>%
  
  # Filter to keep only rows where 'component' is 'Mean'
  filter(component == "Mean") %>%
  
  # Rename the '.variable' values for better readability
  mutate(.variable = case_when(
    .variable == "b_Intercept" ~ "Intercept",
    .variable == "b_dist_water_500" ~ "Distance to water",
    .variable == "b_dist_edge_500" ~ "Distance to forest edge",
    .variable == "b_ciimpact" ~ "Treatment",
    .variable == "b_before_afterall_shrubs_t2:ciimpact" ~ "Time x Treatment",
    .variable == "b_before_afterall_shrubs_t2" ~ "Time",
    TRUE ~ .variable # Keep original name if no match
  ))



p_coeff_veg <- ggplot(posterior_beta, 
                         aes(x = .value, y = .variable, fill = component)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  stat_halfeye(aes(slab_alpha = intercept), .width = c(0.5, 0.89), point_interval = "median_hdi") +
  scale_fill_viridis_d(option = "viridis", end = 0.6) +
  scale_slab_alpha_discrete(range = c(0.5, 0.8)) +
  guides(fill = "none", slab_alpha = "none") +
  labs(x = "Coefficient", y = NULL) +
  theme_clean() + 
  theme(plot.background = element_blank())


# Save plots ---------------------------------------------------------------
saveRDS(p_beta_linpred, file = here("figures", "p_beta_linpred_veg.RDS"))
saveRDS(p_coeff_veg, file = here("figures", "p_coeff_veg.RDS"))

        