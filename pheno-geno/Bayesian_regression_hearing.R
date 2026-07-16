# Bayesian regression model to estimate phenotypic effects of different genotypes on hearing loss; 
# ESKD and ocular leisons were used the same method; 
# in Chinese AS cohort
# Author: Zhen Y

library(rstan)
library(brms)
library(bayesrules)
library(ggsci)
library(cowplot)
library(ggdist)
library(tidybayes)
library(tidyverse)
# hearing loss

# generate prior
moment_estimation <- function(params) {
  alpha <- params[1]
  beta <- params[2]
  
  # Theoretical moments
  theoretical_mean <- alpha / (alpha + beta)
  theoretical_var <- (alpha * beta) / ((alpha + beta)^2 * (alpha + beta + 1))
  
  # Sample moments
  sample_mean <- m / N
  sample_var <- (m * (N - m)) / (N^2 * (N - 1))
  
  # Objective function: sum of squared differences between theoretical and sample moments
  obj <- (theoretical_mean - sample_mean)^2 + (theoretical_var - sample_var)^2
  
  return(obj)
}


# load hearing loss phenotype & genotype data from test_data folder
hear_df = read.table("hearing_loss_data.xls", header = T, sep = "\t")
#
# convert hearing loss to integer: 0 for normal and 1 for hearing loss
hear_df$hearing = as.integer(factor(hear_df$hearing, levels = c("Normal", "Hearing loss")))

# remove genotypes with too few cases
hear_df %>% group_by(genotype) %>% mutate(freq = n()) %>% filter(freq > 10) %>% as.data.frame -> hear_df


# setting prior pars; N for total cases number and m for hearing loss number 
N <- 274
m <- 104

# make age a linear fixed effect
formula <- bf(hearing | trials(273) ~ Age + (1 | genotype))

#get_prior(formula = formula,
#          data = hear_df,
#          family = binomial(link = "logit"))

# non-informative prior with stan implementation
model <- brm(
    formula,
    family = binomial(link = "logit"),
    data = hear_df,
    prior = prior(beta(1, 1), class = b, lb = 0, ub = 1),
    chains = 4, warmup = 1000, iter = 2000, seed = 123, # user can increase the iter and chains for convergence
  refresh = 0
)


# plotting
plot1 <- model |> tidy_draws() |> pivot_longer(cols = starts_with("r_genotype")) |> select(name, value) |> 
mutate(name = str_remove(str_remove(name, "r_genotype\\["), ",Intercept]")) |> ggplot(aes(x = value, y = name, fill=name, color=name)) + 
stat_halfeye(point_interval = median_qi, .width = .95, color="black") + 
labs(title = "Posterior Distributions of Mutation Effects on Hearing Phenotype",
     x = "Effect on Hearing Phenotype",
     y = "Density") + 
theme_bw(base_size = 16, base_rect_size = 1, base_line_size = 1) + scale_fill_flatui() + scale_color_flatui() + theme(legend.position ="none") + 
theme(plot.margin = unit(c(1,1,1,1), "cm"))

# getting posterior distributions
draws_fit <- as_draws_array(model)
posterior_df <- as.data.frame(draws_fit)

# visualizing posterior bar plots
plot2 <- barplot(posterior_summary(model)[,1][5:12], col = pal_flatui()(8), xaxt="n")
axis(side = 1, at = 1:8, labels = rownames(posterior_summary(model))[5:12], las=2)




# sessionInfo:
'''
R version 4.3.1 (2023-06-16 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 11 x64 (build 26200)

Matrix products: default


locale:
[1] LC_COLLATE=English_Canada.utf8  LC_CTYPE=English_Canada.utf8    LC_MONETARY=English_Canada.utf8
[4] LC_NUMERIC=C                    LC_TIME=English_Canada.utf8    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] lubridate_1.9.3    forcats_1.0.0      stringr_1.5.1      dplyr_1.1.4        purrr_1.0.2       
 [6] readr_2.1.5        tidyr_1.3.0        tibble_3.2.1       ggplot2_3.5.1      tidyverse_2.0.0   
[11] tidybayes_3.0.6    ggdist_3.3.2       cowplot_1.1.3      ggsci_3.2.0        bayesrules_0.0.2  
[16] brms_2.21.0        Rcpp_1.0.11        rstan_2.32.6       StanHeaders_2.32.6

loaded via a namespace (and not attached):
  [1] gridExtra_2.3        inline_0.3.19        rlang_1.1.2          magrittr_2.0.3      
  [5] snakecase_0.11.1     matrixStats_1.2.0    e1071_1.7-14         compiler_4.3.1      
  [9] loo_2.8.0            vctrs_0.6.5          reshape2_1.4.4       arrayhelpers_1.1-0  
 [13] pkgconfig_2.0.3      fastmap_1.2.0        backports_1.4.1      labeling_0.4.3      
 [17] utf8_1.2.4           threejs_0.3.3        promises_1.2.1       rmarkdown_2.27      
 [21] tzdb_0.4.0           markdown_1.13        nloptr_2.0.3         xfun_0.46           
 [25] jsonlite_1.8.8       later_1.3.2          broom_1.0.6          parallel_4.3.1      
 [29] R6_2.5.1             dygraphs_1.1.1.6     stringi_1.8.3        car_3.1-2           
 [33] boot_1.3-30          knitr_1.48           zoo_1.8-12           base64enc_0.1-3     
 [37] bayesplot_1.11.1     httpuv_1.6.13        Matrix_1.6-4         splines_4.3.1       
 [41] igraph_1.6.0         timechange_0.2.0     tidyselect_1.2.1     yaml_2.3.10         
 [45] rstudioapi_0.16.0    abind_1.4-5          codetools_0.2-20     miniUI_0.1.1.1      
 [49] pkgbuild_1.4.4       lattice_0.22-5       plyr_1.8.9           withr_3.0.1         
 [53] shiny_1.9.1          bridgesampling_1.1-2 groupdata2_2.0.3     posterior_1.6.0     
 [57] coda_0.19-4.1        evaluate_0.24.0      survival_3.5-7       proxy_0.4-27        
 [61] RcppParallel_5.1.7   xts_0.13.2           pillar_1.9.0         carData_3.0-5       
 [65] tensorA_0.36.2.1     checkmate_2.3.1      DT_0.33              stats4_4.3.1        
 [69] shinyjs_2.1.0        distributional_0.4.0 generics_0.1.3       hms_1.1.3           
 [73] rstantools_2.4.0     munsell_0.5.1        scales_1.3.0         minqa_1.2.6         
 [77] gtools_3.9.5         xtable_1.8-4         class_7.3-22         glue_1.6.2          
 [81] janitor_2.2.0        tools_4.3.1          shinystan_2.6.0      lme4_1.1-35.1       
 [85] colourpicker_1.3.0   mvtnorm_1.2-4        grid_4.3.1           QuickJSR_1.1.3      
 [89] crosstalk_1.2.1      colorspace_2.1-0     patchwork_1.2.0      nlme_3.1-164        
 [93] cli_3.6.2            svUnit_1.0.6         fansi_1.0.6          Brobdingnag_1.2-9   
 [97] gtable_0.3.5         rstatix_0.7.2        digest_0.6.33        farver_2.1.2        
[101] htmlwidgets_1.6.4    htmltools_0.5.9      lifecycle_1.0.4      mime_0.12           
[105] rstanarm_2.32.1      shinythemes_1.2.0    MASS_7.3-60         
'''
