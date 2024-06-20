# 02_baci_plots.R: this script plots posterior median Time x Treatment coefficient with 50% and 89% credible intervals for the animal species studied.

rm(list=ls()) 

library(tidyverse)
library(here)
library(brms)
library(tidybayes) # statpoint_interval


# Load model output -------------------------------------------------------

mod_sp_dry <- readRDS("data/output_data/mod_sp_dry.RDS")
mod_sp_rainy <- readRDS("data/output_data/mod_sp_rainy.RDS")
mod_sp <- list(mod_sp_dry, mod_sp_rainy)

coef_draws_dry <- mod_sp_dry %>% 
  map_df(~.x %>% 
           as_draws_df(), .id = 'species') %>% 
  rename(baci = `b_before_afterafter:ciimpact`) %>% 
  mutate(period = "dry")

coef_draws_rainy <- mod_sp_rainy %>% 
  map_df(~.x %>% 
           as_draws_df(), .id = 'species') %>% 
  rename(baci = `b_before_afterafter:ciimpact`) %>%
  mutate(period = "rainy")

coef_bind <- bind_rows(coef_draws_dry, coef_draws_rainy) %>% 
  mutate(species = str_replace(species, "Bos sp.", "Cattle"))

sp_order <- coef_bind %>% group_by(species) %>% 
  summarise(baci = median(baci)) %>% 
  arrange(baci) %>% pull(species) %>% sort()

# Plot --------------------------------------------------------------------


mod_plot <- coef_bind %>% 
  mutate(species = factor(species, levels = sp_order)) %>%  # Factorize species names
  ggplot(aes(x = species, y = baci, color = period)) + 
  stat_pointinterval(position = "dodge", 
                     point_interval = "median_qi", 
                     .width = c(0.5, 0.89), 
                     point_size = 2) +  # Add point intervals
  coord_flip() +  # Flip coordinates
  geom_hline(yintercept = 0, linetype = "dashed") +  # Add vertical line in zero
  labs(y = "Time x Treatment interaction", x = NULL) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.title = element_blank()) + 
  scale_color_manual(values = c("#952EA0", "#F98477"),
                     name = "Period", 
                     labels = c("Dry", "Rainy")) +  
  geom_vline(xintercept = seq(1.5, 8 - 0.5, 1), colour = 'grey', linetype = 3) +
  guides(alpha = 'none')



# Save plot ---------------------------------------------------------------

#saveRDS(mod_plot, file = here("figures/p_baci_animals.RDS"))







