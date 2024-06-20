## Baci_fences

This repository contains data, code and results associated with:

Vélez, J., McShea, W., Pukazhenthi, B., Rodríguez, J.D., Suárez, M.F., Torres, J.M., Barrera, C. and J. Fieberg. 2024. Cattle exclusion increases encounters of wild herbivores in Neotropical forests. Journal of Applied Ecology.

This study implements a BACI experimental sampling design to quantify the effect of cattle exclusion on encounter probability of the native community of browsers and fruit consumers, and percent ground cover in multifunctional landscapes of the Colombian Orinoquía. Wildlife-permeable fences were built along forest edges in four forest patches (i.e., blocks) containing control and fenced (treatment) sites. We installed 33 camera traps to obtain information about wildlife and cattle encounter probabilities, before and after the fences were constructed. We fit Bayesian generalized linear mixed effects models to quantify the effect of fences via the interaction between the time period (before and after the fences were built) and treatment (control or fenced sites).


## Data

- processed_data/before_after_recs.RDS: data frame of species records.
- processed_data/covariates.RDS: data frame of site covariates.
- processed_data/dat_mods_dry:  data frame of species records in the dry season, annotated with covariates and summarized by site, time period, treatment, and patch.
- processed_data/dat_mods_rainy:  data frame of species records in the rainy season annotated with covariates and summarized by site, time period, treatment, and patch.
- processed_data/sum_ops.RDS: operational days of each camera-trap station.
- output_data/mod_beta_mx.RDS: beta regression model fit to areal percent ground cover. 
- output_data/mod_sp_dry.RDS: generalized linear mixed models for each animal species (dry season).
- output_data/mod_sp_rainy.RDS: generalized linear mixed models for each animal species (rainy season).


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

## Session information

R version 4.3.2 (2023-10-31)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS Monterey 12.2.1

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/Bogota
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] forcats_1.0.0       stringr_1.5.1       dplyr_1.1.4        
 [4] purrr_1.0.2         readr_2.1.4         tidyr_1.3.0        
 [7] tibble_3.2.1        ggplot2_3.5.1       tidyverse_1.3.2    
[10] tidybayes_3.0.6     brms_2.20.4         Rcpp_1.0.12        
[13] rstan_2.32.3        StanHeaders_2.26.28 here_1.0.1         

loaded via a namespace (and not attached):
  [1] tensorA_0.36.2       rstudioapi_0.15.0    jsonlite_1.8.8      
  [4] magrittr_2.0.3       estimability_1.4.1   farver_2.1.1        
  [7] rmarkdown_2.25       fs_1.6.3             vctrs_0.6.5         
 [10] base64enc_0.1-3      htmltools_0.5.7      distributional_0.3.2
 [13] curl_5.2.1           haven_2.5.4          broom_1.0.5         
 [16] cellranger_1.1.0     googlesheets4_1.1.1  parallelly_1.36.0   
 [19] htmlwidgets_1.6.4    plyr_1.8.9           emmeans_1.10.0      
 [22] zoo_1.8-12           lubridate_1.9.3      igraph_2.0.2        
 [25] mime_0.12            lifecycle_1.0.4      pkgconfig_2.0.3     
 [28] colourpicker_1.3.0   Matrix_1.6-1.1       R6_2.5.1            
 [31] fastmap_1.1.1        future_1.33.0        shiny_1.8.0         
 [34] digest_0.6.33        colorspace_2.1-0     furrr_0.3.1         
 [37] rprojroot_2.0.4      crosstalk_1.2.1      fansi_1.0.6         
 [40] timechange_0.3.0     httr_1.4.7           abind_1.4-5         
 [43] compiler_4.3.2       gargle_1.5.2         withr_3.0.0         
 [46] backports_1.4.1      inline_0.3.19        shinystan_2.6.0     
 [49] DBI_1.2.2            QuickJSR_1.0.8       pkgbuild_1.4.3      
 [52] broom.mixed_0.2.9.4  gtools_3.9.5         loo_2.6.0           
 [55] tools_4.3.2          googledrive_2.1.1    httpuv_1.6.13       
 [58] threejs_0.3.3        glue_1.7.0           nlme_3.1-163        
 [61] promises_1.2.1       grid_4.3.2           checkmate_2.3.1     
 [64] reshape2_1.4.4       generics_0.1.3       gtable_0.3.5        
 [67] tzdb_0.4.0           hms_1.1.3            xml2_1.3.6          
 [70] utf8_1.2.4           pillar_1.9.0         ggdist_3.3.2        
 [73] markdown_1.12        posterior_1.5.0      later_1.3.2         
 [76] splines_4.3.2        lattice_0.21-9       tidyselect_1.2.1    
 [79] miniUI_0.1.1.1       knitr_1.45           arrayhelpers_1.1-0  
 [82] gridExtra_2.3        V8_4.4.1             stats4_4.3.2        
 [85] xfun_0.43            bridgesampling_1.1-2 matrixStats_1.2.0   
 [88] DT_0.31              stringi_1.8.3        yaml_2.3.8          
 [91] evaluate_0.23        codetools_0.2-19     cli_3.6.2           
 [94] RcppParallel_5.1.7   shinythemes_1.2.0    xtable_1.8-4        
 [97] munsell_0.5.1        modelr_0.1.11        readxl_1.4.3        
[100] globals_0.16.2       dbplyr_2.4.0         coda_0.19-4         
[103] svUnit_1.0.6         parallel_4.3.2       rstantools_2.3.1.1  
[106] ellipsis_0.3.2       dygraphs_1.1.1.6     bayesplot_1.10.0    
[109] reprex_2.0.2         Brobdingnag_1.2-9    listenv_0.9.0       
[112] mvtnorm_1.2-4        scales_1.3.0         xts_0.13.1          
[115] crayon_1.5.2         rlang_1.1.3          rvest_1.0.3         
[118] shinyjs_2.1.0       