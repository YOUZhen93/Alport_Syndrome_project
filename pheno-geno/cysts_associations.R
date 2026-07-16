# Glm model for genotype association analysis with renal cysts phenotype (Fig. 6c)
# Chinese AS cohort
# Author: Zhen Y

library(MASS)
library(vcdExtra)
library(tidyverse)
library(car)
library(RColorBrewer)
library(mosaic)
library(gridExtra)
library(ggpubr)
library(grid)
library(gtable)
library(ComplexHeatmap)
library(boot)
library(glmnet)
# load cysts genotype and phenotype data from test_data folder
cysts = read.table("./test_data/cyst_data.xls", sep = "\t", header = T)
cysts <- cysts[rep(seq_len(nrow(cysts)), cysts$patient.counts), , drop = FALSE]

# remove GT with few cases that lower the statitics power
cysts %>% group_by(GT) %>% mutate(freq = n()) %>% filter(freq > 10) %>% as.data.frame -> cysts

# reform phenotypic data; 1 for normal; 2 for 1-2 cysts; 3 for greater than 3 cysts
cysts$cysts = as.numeric(factor(cysts$cysts, levels = c("Normal", "1-2 cysts", "gt 3 cysts")))

# poisson linked glm model
model <- glm(cysts ~ GT - 1, data = cysts, family = "poisson")
x = model.matrix(~ GT - 1, data = cysts)
y = cysts$cysts
model2 = glmnet(x, y, family = "poisson", alpha = 1, lambda = 0)
summ2 = coef(model2)
"""
coef of genotypes: both glm and glmnet are consistent
COL4A5_NC1_het COL4A5_collagen_het COL4A3_collagen_NC1_hom COL4A4_NC1_het 
-0.12220210    -0.06637367         -0.06637366             -0.02211853 
COL4A3_collagen_hom COL4A4_hom COL4A4_collagen_het COL4A5_NC1_hemi COL4A5_collagen_hemi
-0.01028422         0.01343244 0.03483631          0.05425432       0.07042717
COL4A3_collagen_het (Intercept)      
0.25427542           0.23342774              
"""
summ = summary(model)

# table drawing
g <- tableGrob(signif(summ$coefficients[order(summ$coefficients[,4]),c(1,4)],2), theme = ttheme_minimal())
separators <- replicate(nrow(g),
                        segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"),
                                     y0 = unit(0, "npc"), y1 = unit(0, "npc"),
                                     gp=gpar(lwd=2)),
                        simplify=FALSE)

ggg = signif(summ$coefficients[order(summ$coefficients[,4]),c(1,4)],2)
rownames(ggg) = gsub("^GT", "", rownames(ggg))
g <- tableGrob(ggg, theme = ttheme_minimal())
separators <- replicate(nrow(g),
                        segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"),
                                     y0 = unit(0, "npc"), y1 = unit(0, "npc"),
                                     gp=gpar(lwd=2)),
                        simplify=FALSE)

g <- gtable::gtable_add_grob(g, grobs = separators, t = 1, r = ncol(g), b = seq_len(nrow(g)), l = 1)
grid.newpage()
grid.draw(g)

# plotting
plot1 = cysts %>% mutate(cysts=as.factor(cysts)) %>% group_by(cysts,GT) %>%
  tally %>% group_by(GT) %>% mutate(freq = n/sum(n)) %>% 
  ggplot(aes(x = GT, y=freq)) + geom_bar(aes(fill=cysts), position="stack", stat="identity", width = 0.5) +
  theme_bw() + geom_text(aes(label=GT), y=0, angle=90, hjust=0) +
  theme(legend.title = element_blank(), axis.text.x = element_blank()) +
  theme(axis.text.y = element_text(size = 12, color="black")) + xlab("") + ylab("") +
  theme(plot.margin = margin(0,1,1,1, "cm")) + theme(axis.ticks.x = element_blank()) +
  scale_fill_brewer(palette = "Blues")

# AIC and statistics for GLM model
a1 = anova(model, test = "LRT")
a1
stepAIC(model)

# coefficients of each genotype
coeff = coefficients(model); names(coeff) = gsub("GT", "", names(coeff))
cysts$coef = coeff[match(cysts$GT, names(coeff))]

plot2 = cysts %>% ggplot(aes(x = 1, y = GT, fill = coef)) + geom_tile(color="black") + geom_text(aes(label=round(coef,2)), size=4) +
  scale_fill_gradient2(low="white", mid = "#dbc5ca", high = "#dc143c", midpoint = 0.2) +
  theme_void() + theme(axis.title = element_blank(), axis.text = element_blank(), 
                       axis.ticks = element_blank(), plot.margin = margin(0.5,2.6,0,2.4, "cm")) + theme(legend.position = "none") +
  ggtitle(paste("P =", formatC(a1[5][2,1], format = "e", digits = 2))) + theme(plot.title = element_text(size = 12, face = "bold.italic")) + coord_flip()

ppp1 = summ$coefficients[order(summ$coefficients[,4]),c(1,4)][,2]
names(ppp1) = gsub("GT", "", names(ppp1))
cysts$pval = ppp1[match(cysts$GT, names(ppp1))]
plot3 = cysts %>% ggplot(aes(x = 1, y = GT, fill = -log10(pval))) + geom_tile(color="black") + geom_text(aes(label=signif(pval,2)), size=2) +
  scale_fill_fermenter(palette = "PRGn") +
  theme_void() + theme(axis.title = element_blank(), axis.text = element_blank(), 
                       axis.ticks = element_blank(), plot.margin = margin(0.5,2.6,0,2.4, "cm")) + theme(legend.position = "none") +
  ggtitle(paste("P =", formatC(a1[5][2,1], format = "e", digits = 2))) + theme(plot.title = element_text(size = 12, face = "bold.italic")) + coord_flip()

ggarrange(plot3, plot2, plot1, ncol=1, heights = c(1,1,6))


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
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] glmnet_4.1-8          boot_1.3-30           ComplexHeatmap_2.16.0 gtable_0.3.5         
 [5] ggpubr_0.6.0          gridExtra_2.3         mosaic_1.9.1          mosaicData_0.20.4    
 [9] ggformula_0.12.0      Matrix_1.6-4          lattice_0.22-5        RColorBrewer_1.1-3   
[13] car_3.1-2             carData_3.0-5         lubridate_1.9.3       forcats_1.0.0        
[17] stringr_1.5.1         dplyr_1.1.4           purrr_1.0.2           readr_2.1.5          
[21] tidyr_1.3.0           tibble_3.2.1          ggplot2_3.5.1         tidyverse_2.0.0      
[25] vcdExtra_0.8-5        gnm_1.1-5             vcd_1.4-12            MASS_7.3-60          

loaded via a namespace (and not attached):
  [1] shinythemes_1.2.0    splines_4.3.1        later_1.3.2          cellranger_1.1.0    
  [5] brms_2.21.0          janitor_2.2.0        xts_0.13.2           lifecycle_1.0.4     
  [9] rstatix_0.7.2        doParallel_1.0.17    rprojroot_2.0.4      StanHeaders_2.32.6  
 [13] crosstalk_1.2.1      backports_1.4.1      magrittr_2.0.3       rmarkdown_2.27      
 [17] yaml_2.3.10          httpuv_1.6.13        pkgbuild_1.4.4       minqa_1.2.6         
 [21] abind_1.4-5          BiocGenerics_0.48.1  nnet_7.3-19          tensorA_0.36.2.1    
 [25] circlize_0.4.16      labelled_2.13.0      S4Vectors_0.40.1     IRanges_2.36.0      
 [29] inline_0.3.19        bridgesampling_1.1-2 codetools_0.2-20     DT_0.33             
 [33] shape_1.4.6.1        tidyselect_1.2.1     rstanarm_2.32.1      bayesplot_1.11.1    
 [37] lme4_1.1-35.1        matrixStats_1.2.0    stats4_4.3.1         base64enc_0.1-3     
 [41] jsonlite_1.8.8       GetoptLong_1.0.5     e1071_1.7-14         ggridges_0.5.6      
 [45] survival_3.5-7       iterators_1.0.14     foreach_1.5.2        tools_4.3.1         
 [49] groupdata2_2.0.3     Rcpp_1.0.11          glue_1.6.2           xfun_0.46           
 [53] here_1.0.1           distributional_0.4.0 ca_0.71.1            loo_2.8.0           
 [57] withr_3.0.1          fastmap_1.2.0        fansi_1.0.6          shinyjs_2.1.0       
 [61] digest_0.6.33        timechange_0.2.0     R6_2.5.1             mime_0.12           
 [65] colorspace_2.1-0     bayesrules_0.0.2     gtools_3.9.5         markdown_1.13       
 [69] threejs_0.3.3        utf8_1.2.4           generics_0.1.3       class_7.3-22        
 [73] htmlwidgets_1.6.4    pkgconfig_2.0.3      dygraphs_1.1.1.6     lmtest_0.9-40       
 [77] htmltools_0.5.9      clue_0.3-65          scales_1.3.0         png_0.1-8           
 [81] posterior_1.6.0      snakecase_0.11.1     knitr_1.48           rstudioapi_0.16.0   
 [85] tzdb_0.4.0           reshape2_1.4.4       rjson_0.2.21         coda_0.19-4.1       
 [89] checkmate_2.3.1      nlme_3.1-164         nloptr_2.0.3         proxy_0.4-27        
 [93] zoo_1.8-12           GlobalOptions_0.1.2  relimp_1.0-5         parallel_4.3.1      
 [97] miniUI_0.1.1.1       pillar_1.9.0         vctrs_0.6.5          shinystan_2.6.0     
[101] promises_1.2.1       cluster_2.1.6        xtable_1.8-4         evaluate_0.24.0     
[105] mvtnorm_1.2-4        cli_3.6.2            compiler_4.3.1       rlang_1.1.2         
[109] crayon_1.5.3         rstantools_2.4.0     ggsignif_0.6.4       plyr_1.8.9          
[113] stringi_1.8.3        rstan_2.32.6         QuickJSR_1.1.3       munsell_0.5.1       
[117] colourpicker_1.3.0   Brobdingnag_1.2-9    mosaicCore_0.9.4.0   qvcalc_1.0.3        
[121] hms_1.1.3            patchwork_1.2.0      shiny_1.9.1          haven_2.5.4         
[125] igraph_1.6.0         broom_1.0.6          RcppParallel_5.1.7   readxl_1.4.3  
'''

