# Bayesian regression model to estimate phenotypic effects of different genotypes on hearing loss;
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



# eyesight is applied as the same method
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
'''
