# 03_baci_preds.R: this script plots the posterior distributions of predicted mean encounter probability and 89% credible intervals for the animal species studied.

rm(list=ls()) 

library(tidyverse)
library(tidybayes) # linpred_draws
library(here)
library(brms)
library(ggthemes) # theme_clean()



# Load model output -------------------------------------------------------

mod_sp_dry <- readRDS(here("data", "output_data", "mod_sp_dry.RDS"))
names_dry <- paste0(names(mod_sp_dry), "_dry")
names(mod_sp_dry) <- names_dry

mod_sp_rainy <- readRDS(here("data", "output_data", "mod_sp_rainy.RDS"))
names_rainy <- paste0(names(mod_sp_rainy), "_rainy")
names(mod_sp_rainy) <- names_rainy

all_mods <- append(mod_sp_dry, mod_sp_rainy)

# Visualize draws ---------------------------------------------------------

coef_draws <- all_mods %>% 
  map_df(~.x %>% as_draws_df(), 
         .id = 'species_period') 


linpreds <- all_mods %>% 
  map_df(~.x %>% linpred_draws(newdata = .$data, 
                                   transform = TRUE), 
         .id = 'species') 

# Posteriors of encounter probability --------

# Posteriors of encounter probability across groups of sites, discarding site-to-site variability within groups

linpreds_mean <- linpreds %>% 
  group_by(species, before_after, ci, .draw) %>% 
  summarise(`Encounter probability` = mean(.linpred)) %>%
  mutate(Period = case_when(before_after %in% "before" ~ "Before", TRUE ~ "After")) %>% 
  mutate(Period = factor(Period, levels = c('Before', 'After')),
         ci = case_when(ci %in% "control" ~ "Control", 
                        TRUE ~ "Fence"),
         species = str_replace(species, "Bos sp.", "Cattle")) %>% 
  separate(species, into = c("species", "period"), sep = "_")


# Plots -------------------------------------------------------------------

linpreds_dry <- ggplot(linpreds_mean %>% 
         filter(period %in% "dry"), 
       aes(Period, `Encounter probability`, 
           group = interaction(before_after, ci), 
           colour = ci)) + 
  stat_pointinterval(position = position_dodge(),
                     point_size = 1.1,
                     .width = 0.89,
                     interval_size_range = c(0.5, 1)) + 
  facet_wrap(~species, ncol = 4) + 
  scale_color_viridis_d(option = "viridis", 
                        end = 0.6, name = NULL, 
                        direction = -1) +
  theme_clean() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank()) + 
  theme(plot.background = element_blank())

linpreds_rainy <- ggplot(linpreds_mean %>% 
                         filter(period %in% "rainy"), 
                       aes(Period, `Encounter probability`, 
                           group = interaction(before_after, ci), 
                           colour = ci)) + 
  stat_pointinterval(position = position_dodge(),
                     point_size = 1.1,
                     .width = 0.89,
                     interval_size_range = c(0.5, 1)) + 
  facet_wrap(~species, ncol = 4) + 
  scale_color_viridis_d(option = "viridis", 
                        end = 0.6, name = NULL, 
                        direction = -1) +
  theme_clean() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank()) + 
  theme(plot.background = element_blank())


# Save plots ------------------------------------------------------------


save(linpreds_dry, linpreds_rainy, file = here("figures", "p_linpreds_animals.rdata"))

           
