# Glm model for genotype association analysis with renal cysts phenotype
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

# load cysts genotype and phenotype data from test_data folder
cysts = read.table("./test_data/cyst_data.xls", sep = "\t", header = T)
cysts <- cysts[rep(seq_len(nrow(cysts)), cysts$patient.counts), , drop = FALSE]

# remove GT with few cases that lower the statitics power
cysts %>% group_by(GT) %>% mutate(freq = n()) %>% filter(freq > 10) %>% as.data.frame -> cysts

# reform phenotypic data; 1 for normal; 2 for 1-2 cysts; 3 for greater than 3 cysts
cysts$cysts = as.numeric(factor(cysts$cysts, levels = c("Normal", "1-2 cysts", "gt 3 cysts")))

# poisson linked glm model
model <- glm(cysts ~ GT - 1, data = cysts, family = "poisson")

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
 [1] boot_1.3-30           ComplexHeatmap_2.16.0 gtable_0.3.5          gridExtra_2.3        
 [5] mosaic_1.9.1          mosaicData_0.20.4     ggformula_0.12.0      Matrix_1.6-4         
 [9] lattice_0.22-5        RColorBrewer_1.1-3    car_3.1-2             carData_3.0-5        
[13] vcdExtra_0.8-5        gnm_1.1-5             vcd_1.4-12            MASS_7.3-60          
[17] survival_3.5-7        survminer_0.4.9       ggpubr_0.6.0          ggfortify_0.4.17     
[21] lubridate_1.9.3       forcats_1.0.0         stringr_1.5.1         dplyr_1.1.4          
[25] purrr_1.0.2           readr_2.1.5           tidyr_1.3.0           tibble_3.2.1         
[29] ggplot2_3.5.1         tidyverse_2.0.0       tidybayes_3.0.6       ggdist_3.3.2         
[33] cowplot_1.1.3         ggsci_3.2.0           bayesrules_0.0.2      brms_2.21.0          
[37] Rcpp_1.0.11           rstan_2.32.6          StanHeaders_2.32.6   

loaded via a namespace (and not attached):
  [1] svUnit_1.0.6         shinythemes_1.2.0    splines_4.3.1        later_1.3.2         
  [5] cellranger_1.1.0     janitor_2.2.0        xts_0.13.2           lifecycle_1.0.4     
  [9] rstatix_0.7.2        doParallel_1.0.17    rprojroot_2.0.4      crosstalk_1.2.1     
 [13] backports_1.4.1      magrittr_2.0.3       rmarkdown_2.27       yaml_2.3.10         
 [17] httpuv_1.6.13        pkgbuild_1.4.4       minqa_1.2.6          abind_1.4-5         
 [21] BiocGenerics_0.48.1  nnet_7.3-19          tensorA_0.36.2.1     xlsxjars_0.6.1      
 [25] circlize_0.4.16      labelled_2.13.0      S4Vectors_0.40.1     IRanges_2.36.0      
 [29] KMsurv_0.1-5         inline_0.3.19        bridgesampling_1.1-2 codetools_0.2-20    
 [33] DT_0.33              shape_1.4.6.1        tidyselect_1.2.1     rstanarm_2.32.1     
 [37] bayesplot_1.11.1     farver_2.1.2         lme4_1.1-35.1        matrixStats_1.2.0   
 [41] stats4_4.3.1         base64enc_0.1-3      jsonlite_1.8.8       GetoptLong_1.0.5    
 [45] e1071_1.7-14         iterators_1.0.14     ggridges_0.5.6       foreach_1.5.2       
 [49] tools_4.3.1          groupdata2_2.0.3     glue_1.6.2           xfun_0.46           
 [53] here_1.0.1           distributional_0.4.0 ca_0.71.1            loo_2.8.0           
 [57] withr_3.0.1          fastmap_1.2.0        fansi_1.0.6          shinyjs_2.1.0       
 [61] digest_0.6.33        timechange_0.2.0     R6_2.5.1             mime_0.12           
 [65] colorspace_2.1-0     gtools_3.9.5         markdown_1.13        threejs_0.3.3       
 [69] utf8_1.2.4           generics_0.1.3       data.table_1.15.4    class_7.3-22        
 [73] htmlwidgets_1.6.4    pkgconfig_2.0.3      dygraphs_1.1.1.6     rJava_1.0-11        
 [77] lmtest_0.9-40        survMisc_0.5.6       htmltools_0.5.9      clue_0.3-65         
 [81] scales_1.3.0         png_0.1-8            posterior_1.6.0      snakecase_0.11.1    
 [85] knitr_1.48           km.ci_0.5-6          rstudioapi_0.16.0    rjson_0.2.21        
 [89] tzdb_0.4.0           reshape2_1.4.4       coda_0.19-4.1        checkmate_2.3.1     
 [93] nlme_3.1-164         nloptr_2.0.3         GlobalOptions_0.1.2  proxy_0.4-27        
 [97] zoo_1.8-12           relimp_1.0-5         parallel_4.3.1       miniUI_0.1.1.1      
[101] pillar_1.9.0         vctrs_0.6.5          shinystan_2.6.0      promises_1.2.1      
[105] arrayhelpers_1.1-0   cluster_2.1.6        xtable_1.8-4         evaluate_0.24.0     
[109] mvtnorm_1.2-4        cli_3.6.2            compiler_4.3.1       rlang_1.1.2         
[113] crayon_1.5.3         xlsx_0.6.5           rstantools_2.4.0     ggsignif_0.6.4      
[117] labeling_0.4.3       plyr_1.8.9           stringi_1.8.3        QuickJSR_1.1.3      
[121] munsell_0.5.1        colourpicker_1.3.0   mosaicCore_0.9.4.0   Brobdingnag_1.2-9   
[125] qvcalc_1.0.3         hms_1.1.3            patchwork_1.2.0      shiny_1.9.1         
[129] haven_2.5.4          igraph_1.6.0         broom_1.0.6          RcppParallel_5.1.7  
[133] readxl_1.4.3        
'''

