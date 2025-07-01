## baci_fencing

This repository contains data, code and results associated with:

Vélez, J., McShea, W., Pukazhenthi, B., Rodríguez, J.D., Suárez, M.F., Torres, J.M., Barrera, C. and J. Fieberg. 2024. Cattle exclusion increases encounters of wild herbivores in Neotropical forests. Journal of Applied Ecology, 61, 2444-2454. https://doi.org/10.1111/1365-2664.14751

This study implements a BACI experimental sampling design to quantify the effect of cattle exclusion on encounter probability of the native community of browsers and fruit consumers, and percent ground cover in multifunctional landscapes of the Colombian Orinoquía.

## Data

### Processed
- before_after_recs.RDS: data frame of species records.
- covariates.RDS: data frame of site covariates.
- dat_mods_dry:  data frame of species records in the dry season, annotated with covariates and summarized by site, time period, treatment, and patch.
- dat_mods_rainy:  data frame of species records in the rainy season annotated with covariates and summarized by site, time period, treatment, and patch.
- sum_ops.RDS: operational days of each camera-trap station.

### Output
- mod_beta_mx.RDS: beta regression model fit to areal percent ground cover. 
- mod_sp_dry.RDS: generalized linear mixed models for each animal species (dry season).
- mod_sp_rainy.RDS: generalized linear mixed models for each animal species (rainy season).


## Programs

- 01a_glmm_animals.R: runs Bayesian generalized linear mixed models to estimate encounter probabilities of different animal species.
- 02_baci_plots.R: plots posterior median Time x Treatment coefficients.
- 03_baci_preds.R: plots posterior distributions of predicted mean encounter probability.
- 04_diagnostics_ppc.Rmd: assesses model fit and convergence.

## Figures

- p_baci_animals.RDS: figure 2 in manuscript.
- p_coeff_veg.RDS: figure 3 in manuscript.
- p_beta_linpred_veg.RDS: figure 4 in manuscript.
- p_linpreds_animals.rdata: figure 3 in supplementary information.
