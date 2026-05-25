# ESKD survival analysis
# ESKD data depoists at folder test_data
# Author: Zhen Y

library(ggfortify)
library(survminer)
library(ggpubr)
library(ggplot2)
library(survival)
library(xlsx)
library(dplyr)
library(cowplot)

# load data 
eskd <- read.table("./test_data/ESKD_data.xls", header = T, sep = "\t")[,1:3]
# columns: time status genotype
# status: 1 > normal; 2 > ESKD

# remove genotype with too few cases
eskd %>% group_by(genotype) %>% mutate(freq = n()) %>% 
                       filter(freq > 10) %>% as.data.frame -> eskd



# Cox proportional hazards regression models
cox1 <- coxph(Surv(time, status) ~ 1 + genotype, data = eskd)
cox1_sum <- summary(cox1)
print(cox1_sum$sctest[3]) # log-rank p
# P = 1.599728e-08
print(cox1_sum$coefficients[,1]) # HR
"""
genotypeCOL4A3_collagen_hom  genotypeCOL4A4_collagen_het           genotypeCOL4A4_hom 
                  0.93121161                   0.33043814                   0.03444387 
genotypeCOL4A5_collagen_hemi  genotypeCOL4A5_collagen_het      genotypeCOL4A5_NC1_hemi 
                  1.92766118                  -0.23594279                   2.06459783 
      genotypeCOL4A5_NC1_het 
                  1.02955536 
"""
cox1_newdata <- data.frame(genotype = cox1$xlevels, time = rep(mean(eskd$time), 8))
rownames(cox1_newdata) = cox1_newdata$genotype
cox1_fit <- survfit(cox1, newdata = cox1_newdata)

# strata is following the rownames of cox4_newdata

# plotting
ggsurvplot(
  cox1_fit,             
  data = eskd,         
  #risk.table = TRUE,     
  pval = TRUE,          
  conf.int = TRUE,     
  palette = "npg",
  xlim = c(0,65),         
  xlab = "Ages",  
  break.time.by = 10,    
  ggtheme = theme_cowplot(),
  #risk.table.y.text.col = T,
  #risk.table.height = 0.25, 
  #risk.table.y.text = FALSE,
  #ncensor.plot = TRUE,     
  #ncensor.plot.height = 0.25,
  conf.int.style = "step", 
  surv.median.line = "hv"
)




# sessionInfo:
'''
R version 4.3.1 (2023-06-16 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 11 x64 (build 26200)

Matrix products: default


locale:
[1] LC_COLLATE=English_Canada.utf8  LC_CTYPE=English_Canada.utf8    LC_MONETARY=English_Canada.utf8
[4] LC_NUMERIC=C                    LC_TIME=English_Canada.utf8    

time zone: Asia/Shanghai
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] xlsx_0.6.5         survival_3.5-7     survminer_0.4.9    ggpubr_0.6.0       ggfortify_0.4.17  
 [6] lubridate_1.9.3    forcats_1.0.0      stringr_1.5.1      dplyr_1.1.4        purrr_1.0.2       
[11] readr_2.1.5        tidyr_1.3.0        tibble_3.2.1       ggplot2_3.5.1      tidyverse_2.0.0   
[16] tidybayes_3.0.6    ggdist_3.3.2       cowplot_1.1.3      ggsci_3.2.0        bayesrules_0.0.2  
[21] brms_2.21.0        Rcpp_1.0.11        rstan_2.32.6       StanHeaders_2.32.6

loaded via a namespace (and not attached):
  [1] tensorA_0.36.2.1     rstudioapi_0.16.0    jsonlite_1.8.8       magrittr_2.0.3      
  [5] farver_2.1.2         nloptr_2.0.3         rmarkdown_2.27       vctrs_0.6.5         
  [9] minqa_1.2.6          base64enc_0.1-3      rstatix_0.7.2        janitor_2.2.0       
 [13] htmltools_0.5.9      distributional_0.4.0 broom_1.0.6          htmlwidgets_1.6.4   
 [17] plyr_1.8.9           zoo_1.8-12           igraph_1.6.0         mime_0.12           
 [21] lifecycle_1.0.4      pkgconfig_2.0.3      colourpicker_1.3.0   Matrix_1.6-4        
 [25] R6_2.5.1             fastmap_1.2.0        shiny_1.9.1          snakecase_0.11.1    
 [29] digest_0.6.33        colorspace_2.1-0     patchwork_1.2.0      crosstalk_1.2.1     
 [33] labeling_0.4.3       km.ci_0.5-6          fansi_1.0.6          timechange_0.2.0    
 [37] abind_1.4-5          compiler_4.3.1       proxy_0.4-27         withr_3.0.1         
 [41] backports_1.4.1      inline_0.3.19        shinystan_2.6.0      carData_3.0-5       
 [45] QuickJSR_1.1.3       pkgbuild_1.4.4       ggsignif_0.6.4       MASS_7.3-60         
 [49] gtools_3.9.5         loo_2.8.0            tools_4.3.1          httpuv_1.6.13       
 [53] threejs_0.3.3        glue_1.6.2           nlme_3.1-164         promises_1.2.1      
 [57] grid_4.3.1           checkmate_2.3.1      reshape2_1.4.4       generics_0.1.3      
 [61] gtable_0.3.5         KMsurv_0.1-5         tzdb_0.4.0           class_7.3-22        
 [65] data.table_1.15.4    hms_1.1.3            car_3.1-2            utf8_1.2.4          
 [69] pillar_1.9.0         markdown_1.13        posterior_1.6.0      later_1.3.2         
 [73] rJava_1.0-11         splines_4.3.1        lattice_0.22-5       tidyselect_1.2.1    
 [77] miniUI_0.1.1.1       knitr_1.48           arrayhelpers_1.1-0   gridExtra_2.3       
 [81] groupdata2_2.0.3     stats4_4.3.1         xfun_0.46            rstanarm_2.32.1     
 [85] bridgesampling_1.1-2 matrixStats_1.2.0    DT_0.33              stringi_1.8.3       
 [89] yaml_2.3.10          boot_1.3-30          xlsxjars_0.6.1       evaluate_0.24.0     
 [93] codetools_0.2-20     cli_3.6.2            RcppParallel_5.1.7   shinythemes_1.2.0   
 [97] xtable_1.8-4         munsell_0.5.1        survMisc_0.5.6       coda_0.19-4.1       
[101] svUnit_1.0.6         parallel_4.3.1       rstantools_2.4.0     dygraphs_1.1.1.6    
[105] bayesplot_1.11.1     Brobdingnag_1.2-9    lme4_1.1-35.1        mvtnorm_1.2-4       
[109] scales_1.3.0         xts_0.13.2           e1071_1.7-14         crayon_1.5.3        
[113] rlang_1.1.2          shinyjs_2.1.0      
'''




