Analysis for consensusDE paper - airway RNA-seq
================
Ashley J. Waardenberg
30/06/2019

-   [Notes](#notes)
-   [functions used for analysis](#functions-used-for-analysis)
-   [Analysis of RNA-seq data with consensusDE](#analysis-of-rna-seq-data-with-consensusde)
-   [LogFC goodness of fit - without RUV](#logfc-goodness-of-fit---without-ruv)
-   [LogFC goodness of fit - with RUV](#logfc-goodness-of-fit---with-ruv)
-   [Jaccard Coefficient and Set Sizes](#jaccard-coefficient-and-set-sizes)
-   [Mean LogFC standard deviation of sets](#mean-logfc-standard-deviation-of-sets)
-   [sessionInfo](#sessioninfo)
-   [functions used](#functions-used)

Notes
=====

30-06-2019 + updated for compatability with consensusDE version 1.3.3 + setting of norm\_method = "all\_defaults" for multi\_de\_pairs

``` r
library(consensusDE)
library(airway)
data(airway)
```

functions used for analysis
===========================

``` r
source("./functions/return_airway_stats.r")
```

Analysis of RNA-seq data with consensusDE
=========================================

``` r
# establishment of data for pair-wise analysis
colData(airway)$group <- colData(airway)$dex
colData(airway)$file <- rownames(colData(airway))

# read into a summarized experiment, with low counts filtered
airway_filter <- buildSummarized(summarized = airway,
                                 filter = TRUE)

# run DE analysis with RUV
airway_ruv <- multi_de_pairs(summarized = airway_filter,
                             norm_method = "all_defaults",
                             ruv_correct = TRUE)

# run DE analysis without RUV
airway_NOruv <- multi_de_pairs(summarized = airway_filter,
                               norm_method = "all_defaults",
                               ruv_correct = FALSE)
```

LogFC goodness of fit - without RUV
===================================

``` r
# extract logFC values for edgeR, voom and DESeq2
all_FC <- get_merged_logfc(airway_NOruv)
# plot results and obtain goodness of fit

# 1. Voom vs. DESeq2
plot(all_FC$logFC_VM, 
     all_FC$logFC_DE, 
     xlab="voom - logFC", 
     ylab="DESeq2 - logFC")
```

![](1_airway_analysis_files/figure-markdown_github/unnamed-chunk-3-1.png)

``` r
summary(lm(all_FC$logFC_VM~all_FC$logFC_DE))
```

    ## 
    ## Call:
    ## lm(formula = all_FC$logFC_VM ~ all_FC$logFC_DE)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.30193 -0.01724 -0.00804  0.00671  1.60902 
    ## 
    ## Coefficients:
    ##                  Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)     0.0132990  0.0006389   20.82   <2e-16 ***
    ## all_FC$logFC_DE 0.9822373  0.0011410  860.87   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.08062 on 15924 degrees of freedom
    ## Multiple R-squared:  0.979,  Adjusted R-squared:  0.979 
    ## F-statistic: 7.411e+05 on 1 and 15924 DF,  p-value: < 2.2e-16

``` r
# 2. Voom vs. edgeR
plot(all_FC$logFC_VM, 
     all_FC$logFC_ER,
     xlab="voom - logFC", 
     ylab="edgeR - logFC")
```

![](1_airway_analysis_files/figure-markdown_github/unnamed-chunk-3-2.png)

``` r
summary(lm(all_FC$logFC_VM~all_FC$logFC_ER))
```

    ## 
    ## Call:
    ## lm(formula = all_FC$logFC_VM ~ all_FC$logFC_ER)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.30755 -0.01756 -0.00821  0.00671  1.91259 
    ## 
    ## Coefficients:
    ##                  Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)     0.0099360  0.0006436   15.44   <2e-16 ***
    ## all_FC$logFC_ER 0.9843375  0.0011520  854.47   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.08122 on 15924 degrees of freedom
    ## Multiple R-squared:  0.9787, Adjusted R-squared:  0.9787 
    ## F-statistic: 7.301e+05 on 1 and 15924 DF,  p-value: < 2.2e-16

``` r
# 3. edgeR vs. DESeq2
plot(all_FC$logFC_ER, 
     all_FC$logFC_DE, 
     xlab="edgeR - logFC", 
     ylab="DESeq2 - logFC")
```

![](1_airway_analysis_files/figure-markdown_github/unnamed-chunk-3-3.png)

``` r
summary(lm(all_FC$logFC_ER~all_FC$logFC_DE))
```

    ## 
    ## Call:
    ## lm(formula = all_FC$logFC_ER ~ all_FC$logFC_DE)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.52743 -0.00084  0.00013  0.00106  0.14335 
    ## 
    ## Coefficients:
    ##                  Estimate Std. Error  t value Pr(>|t|)    
    ## (Intercept)     3.417e-03  4.651e-05    73.46   <2e-16 ***
    ## all_FC$logFC_DE 9.977e-01  8.306e-05 12011.37   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.005869 on 15924 degrees of freedom
    ## Multiple R-squared:  0.9999, Adjusted R-squared:  0.9999 
    ## F-statistic: 1.443e+08 on 1 and 15924 DF,  p-value: < 2.2e-16

LogFC goodness of fit - with RUV
================================

``` r
# extract logFC values for edgeR, voom and DESeq2
all_FC_ruv <- get_merged_logfc(airway_ruv)
# plot results and obtain goodness of fit

# 1. Voom vs. DESeq2
plot(all_FC_ruv$logFC_VM, 
     all_FC_ruv$logFC_DE, 
     xlab="voom - logFC", 
     ylab="DESeq2 - logFC")
```

![](1_airway_analysis_files/figure-markdown_github/unnamed-chunk-4-1.png)

``` r
summary(lm(all_FC_ruv$logFC_VM~all_FC_ruv$logFC_DE))
```

    ## 
    ## Call:
    ## lm(formula = all_FC_ruv$logFC_VM ~ all_FC_ruv$logFC_DE)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.07048 -0.00970 -0.00203  0.00656  1.40422 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)         0.0058481  0.0004775   12.25   <2e-16 ***
    ## all_FC_ruv$logFC_DE 0.9914890  0.0008559 1158.48   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.06025 on 15924 degrees of freedom
    ## Multiple R-squared:  0.9883, Adjusted R-squared:  0.9883 
    ## F-statistic: 1.342e+06 on 1 and 15924 DF,  p-value: < 2.2e-16

``` r
# 2. Voom vs. edgeR
plot(all_FC_ruv$logFC_VM, 
     all_FC_ruv$logFC_ER,
     xlab="voom - logFC", 
     ylab="edgeR - logFC")
```

![](1_airway_analysis_files/figure-markdown_github/unnamed-chunk-4-2.png)

``` r
summary(lm(all_FC_ruv$logFC_VM~all_FC_ruv$logFC_ER))
```

    ## 
    ## Call:
    ## lm(formula = all_FC_ruv$logFC_VM ~ all_FC_ruv$logFC_ER)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.07644 -0.00964 -0.00219  0.00627  1.94077 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error  t value Pr(>|t|)    
    ## (Intercept)         0.0024585  0.0004789    5.134 2.87e-07 ***
    ## all_FC_ruv$logFC_ER 0.9947762  0.0008611 1155.291  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.06042 on 15924 degrees of freedom
    ## Multiple R-squared:  0.9882, Adjusted R-squared:  0.9882 
    ## F-statistic: 1.335e+06 on 1 and 15924 DF,  p-value: < 2.2e-16

``` r
# 3. edgeR vs. DESeq2
plot(all_FC_ruv$logFC_ER, 
     all_FC_ruv$logFC_DE, 
     xlab="edgeR - logFC", 
     ylab="DESeq2 - logFC")
```

![](1_airway_analysis_files/figure-markdown_github/unnamed-chunk-4-3.png)

``` r
summary(lm(all_FC_ruv$logFC_ER~all_FC_ruv$logFC_DE))
```

    ## 
    ## Call:
    ## lm(formula = all_FC_ruv$logFC_ER ~ all_FC_ruv$logFC_DE)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.54069 -0.00075  0.00019  0.00106  0.19056 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)         0.0034086  0.0000644   52.93   <2e-16 ***
    ## all_FC_ruv$logFC_DE 0.9965569  0.0001154 8633.86   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.008126 on 15924 degrees of freedom
    ## Multiple R-squared:  0.9998, Adjusted R-squared:  0.9998 
    ## F-statistic: 7.454e+07 on 1 and 15924 DF,  p-value: < 2.2e-16

Jaccard Coefficient and Set Sizes
=================================

``` r
non_ruv <- return_cde_stats(airway_NOruv)
ruv <- return_cde_stats(airway_ruv)

# Jaccard coefficients
JC_stats <- rbind(unlist(non_ruv$props),
                  unlist(ruv$props))
rownames(JC_stats) <- c("non_RUV", "RUV")

#print results
JC_stats
```

    ##                ER        DE        VM
    ## non_RUV 0.8174078 0.6290499 0.9201278
    ## RUV     0.8565588 0.7524693 0.9229075

``` r
# set sizes
non_ruv_set_size <- c(length(non_ruv$venn_input$VM),
                      length(non_ruv$venn_input$DE),
                      length(non_ruv$venn_input$ER),
                      length(non_ruv$intersect))

ruv_set_size <- c(length(ruv$venn_input$VM),
                  length(ruv$venn_input$DE),
                  length(ruv$venn_input$ER),
                  length(ruv$intersect))

# percetage increases (nonRUV to RUV)
percent_diff <- c((length(ruv$venn_input$VM)-length(non_ruv$venn_input$VM))/length(ruv$venn_input$VM),
                  (length(ruv$venn_input$DE)-length(non_ruv$venn_input$DE))/length(ruv$venn_input$DE),
                  (length(ruv$venn_input$ER)-length(non_ruv$venn_input$ER))/length(ruv$venn_input$ER),
                  (length(ruv$intersect)-length(non_ruv$intersect))/length(ruv$intersect))

set_sizes <- rbind(non_ruv_set_size,
                   ruv_set_size,
                   percent_diff)

colnames(set_sizes) <- c("voom", "DESeq2", "edgeR", "intersect")
rownames(set_sizes) <- c("non_RUV", "RUV", "percent_diff")

# print results
set_sizes
```

    ##                      voom       DESeq2        edgeR    intersect
    ## non_RUV      1878.0000000 2747.0000000 2114.0000000 1728.0000000
    ## RUV          2724.0000000 3341.0000000 2935.0000000 2514.0000000
    ## percent_diff    0.3105727    0.1777911    0.2797274    0.3126492

Mean LogFC standard deviation of sets
=====================================

``` r
mean_sd <- rbind(c(mean(abs(non_ruv$ER$LogFC_sd)),
                   mean(abs(non_ruv$VM$LogFC_sd)),
                   mean(abs(non_ruv$DE$LogFC_sd))),
                 c(mean(abs(ruv$ER$LogFC_sd)),
                   mean(abs(ruv$VM$LogFC_sd)),
                   mean(abs(ruv$DE$LogFC_sd))))
colnames(mean_sd) <- c("edgeR", "voom", "DESeq2")
rownames(mean_sd) <- c("non_RUV", "RUV")

# print results
mean_sd
```

    ##              edgeR       voom     DESeq2
    ## non_RUV 0.02133190 0.01693218 0.01795723
    ## RUV     0.01562589 0.01318950 0.01321216

``` r
# of the non-intersecting sets for each algorithm
mean_sd_not_intersect <- rbind(c(mean(abs(non_ruv$ER[!non_ruv$ER$ID %in% non_ruv$intersect,]$LogFC_sd)),
                                 mean(abs(non_ruv$VM[!non_ruv$VM$ID %in% non_ruv$intersect,]$LogFC_sd)),
                                 mean(abs(non_ruv$DE[!non_ruv$DE$ID %in% non_ruv$intersect,]$LogFC_sd))),
                               c(mean(abs(ruv$ER[!ruv$ER$ID %in% ruv$intersect,]$LogFC_sd)),
                                 mean(abs(ruv$VM[!ruv$VM$ID %in% ruv$intersect,]$LogFC_sd)),
                                 mean(abs(ruv$DE[!ruv$DE$ID %in% ruv$intersect,]$LogFC_sd))))
colnames(mean_sd_not_intersect) <- c("edgeR", "voom", "DESeq2")
rownames(mean_sd_not_intersect) <- c("non_RUV", "RUV")

# print results
mean_sd_not_intersect
```

    ##              edgeR       voom     DESeq2
    ## non_RUV 0.03929638 0.01247605 0.01903953
    ## RUV     0.02747661 0.00778035 0.01190751

``` r
# overall intersect
mean_sd_intersect <- c(mean(abs(non_ruv$DE[non_ruv$DE$ID %in% non_ruv$intersect,]$LogFC_sd)),
                       mean(abs(ruv$DE[ruv$DE$ID %in% ruv$intersect,]$LogFC_sd)))
names(mean_sd_intersect) <- c("non_RUV", "RUV")

# print results
mean_sd_intersect
```

    ##    non_RUV        RUV 
    ## 0.01731900 0.01364134

``` r
# overall union
mean_sd_union <- c(mean(abs(unique(rbind(non_ruv$ER,non_ruv$VM,non_ruv$DE)$LogFC_sd))),
                   mean(abs(unique(rbind(ruv$ER,ruv$VM,ruv$DE)$LogFC_sd))))
names(mean_sd_union) <- c("non_RUV", "RUV")

# print results
mean_sd_union
```

    ##    non_RUV        RUV 
    ## 0.01869146 0.01413187

sessionInfo
===========

``` r
sessionInfo()
```

    ## R version 3.5.1 (2018-07-02)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS High Sierra 10.13.6
    ## 
    ## Matrix products: default
    ## BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_AU.UTF-8/en_AU.UTF-8/en_AU.UTF-8/C/en_AU.UTF-8/en_AU.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] airway_0.114.0              SummarizedExperiment_1.10.1
    ##  [3] DelayedArray_0.6.6          BiocParallel_1.16.6        
    ##  [5] matrixStats_0.54.0          Biobase_2.40.0             
    ##  [7] GenomicRanges_1.32.7        GenomeInfoDb_1.16.0        
    ##  [9] IRanges_2.14.12             S4Vectors_0.18.3           
    ## [11] consensusDE_1.3.3           BiocGenerics_0.26.0        
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] colorspace_1.3-2                         
    ##   [2] hwriter_1.3.2                            
    ##   [3] htmlTable_1.12                           
    ##   [4] XVector_0.20.0                           
    ##   [5] base64enc_0.1-3                          
    ##   [6] rstudioapi_0.7                           
    ##   [7] bit64_0.9-7                              
    ##   [8] AnnotationDbi_1.44.0                     
    ##   [9] splines_3.5.1                            
    ##  [10] R.methodsS3_1.7.1                        
    ##  [11] DESeq_1.32.0                             
    ##  [12] geneplotter_1.58.0                       
    ##  [13] knitr_1.23                               
    ##  [14] Formula_1.2-3                            
    ##  [15] Rsamtools_1.34.1                         
    ##  [16] annotate_1.58.0                          
    ##  [17] cluster_2.0.7-1                          
    ##  [18] R.oo_1.22.0                              
    ##  [19] compiler_3.5.1                           
    ##  [20] httr_1.3.1                               
    ##  [21] backports_1.1.2                          
    ##  [22] assertthat_0.2.1                         
    ##  [23] Matrix_1.2-14                            
    ##  [24] lazyeval_0.2.2                           
    ##  [25] limma_3.36.5                             
    ##  [26] acepack_1.4.1                            
    ##  [27] htmltools_0.3.6                          
    ##  [28] prettyunits_1.0.2                        
    ##  [29] tools_3.5.1                              
    ##  [30] bindrcpp_0.2.2                           
    ##  [31] gtable_0.3.0                             
    ##  [32] glue_1.3.1                               
    ##  [33] GenomeInfoDbData_1.1.0                   
    ##  [34] dplyr_0.7.8                              
    ##  [35] ShortRead_1.38.0                         
    ##  [36] Rcpp_0.12.19                             
    ##  [37] TxDb.Dmelanogaster.UCSC.dm3.ensGene_3.2.2
    ##  [38] Biostrings_2.48.0                        
    ##  [39] rtracklayer_1.40.6                       
    ##  [40] xfun_0.7                                 
    ##  [41] stringr_1.4.0                            
    ##  [42] ensembldb_2.6.7                          
    ##  [43] XML_3.98-1.16                            
    ##  [44] dendextend_1.12.0                        
    ##  [45] edgeR_3.22.5                             
    ##  [46] zlibbioc_1.26.0                          
    ##  [47] MASS_7.3-50                              
    ##  [48] scales_1.0.0                             
    ##  [49] aroma.light_3.10.0                       
    ##  [50] pcaMethods_1.72.0                        
    ##  [51] hms_0.4.2                                
    ##  [52] ProtGenerics_1.14.0                      
    ##  [53] AnnotationFilter_1.6.0                   
    ##  [54] RColorBrewer_1.1-2                       
    ##  [55] yaml_2.2.0                               
    ##  [56] curl_3.2                                 
    ##  [57] memoise_1.1.0                            
    ##  [58] RUVSeq_1.16.1                            
    ##  [59] gridExtra_2.3                            
    ##  [60] ggplot2_3.1.1                            
    ##  [61] biomaRt_2.36.1                           
    ##  [62] rpart_4.1-13                             
    ##  [63] latticeExtra_0.6-28                      
    ##  [64] stringi_1.4.3                            
    ##  [65] RSQLite_2.1.1                            
    ##  [66] genefilter_1.62.0                        
    ##  [67] checkmate_1.8.5                          
    ##  [68] GenomicFeatures_1.32.3                   
    ##  [69] rlang_0.3.4                              
    ##  [70] pkgconfig_2.0.2                          
    ##  [71] bitops_1.0-6                             
    ##  [72] evaluate_0.14                            
    ##  [73] lattice_0.20-35                          
    ##  [74] purrr_0.2.5                              
    ##  [75] bindr_0.1.1                              
    ##  [76] GenomicAlignments_1.16.0                 
    ##  [77] htmlwidgets_1.2                          
    ##  [78] bit_1.1-14                               
    ##  [79] tidyselect_0.2.5                         
    ##  [80] plyr_1.8.4                               
    ##  [81] magrittr_1.5                             
    ##  [82] DESeq2_1.20.0                            
    ##  [83] R6_2.2.2                                 
    ##  [84] Hmisc_4.1-1                              
    ##  [85] DBI_1.0.0                                
    ##  [86] pillar_1.3.1                             
    ##  [87] foreign_0.8-71                           
    ##  [88] survival_2.42-6                          
    ##  [89] RCurl_1.95-4.11                          
    ##  [90] nnet_7.3-12                              
    ##  [91] tibble_2.1.1                             
    ##  [92] EDASeq_2.14.1                            
    ##  [93] crayon_1.3.4                             
    ##  [94] rmarkdown_1.13                           
    ##  [95] viridis_0.5.1                            
    ##  [96] progress_1.2.0                           
    ##  [97] locfit_1.5-9.1                           
    ##  [98] grid_3.5.1                               
    ##  [99] data.table_1.12.2                        
    ## [100] blob_1.1.1                               
    ## [101] digest_0.6.18                            
    ## [102] xtable_1.8-3                             
    ## [103] R.utils_2.7.0                            
    ## [104] munsell_0.5.0                            
    ## [105] viridisLite_0.3.0

functions used
==============

printing for static reference

``` r
return_cde_stats
```

    ## function (input_pairs, p_cut = 0.05) 
    ## {
    ##     merged_data <- input_pairs$merged[[1]]
    ##     ER <- merged_data[merged_data$edger_adj_p <= p_cut, ]
    ##     DE <- merged_data[merged_data$deseq_adj_p <= p_cut, ]
    ##     VM <- merged_data[merged_data$voom_adj_p <= p_cut, ]
    ##     ER_ID <- as.character(ER$ID)
    ##     DE_ID <- as.character(DE$ID)
    ##     VM_ID <- as.character(VM$ID)
    ##     ER.DE.VM_ID <- (intersect(intersect(ER_ID, DE_ID), VM_ID))
    ##     ER_prop <- length(intersect(ER_ID, ER.DE.VM_ID))/length(ER_ID)
    ##     DE_prop <- length(intersect(DE_ID, ER.DE.VM_ID))/length(DE_ID)
    ##     VM_prop <- length(intersect(VM_ID, ER.DE.VM_ID))/length(VM_ID)
    ##     props <- list(ER = ER_prop, DE = DE_prop, VM = VM_prop)
    ##     venn.inputs <- list(ER_ID, DE_ID, VM_ID)
    ##     names(venn.inputs) <- c("ER", "DE", "VM")
    ##     return(list(venn_input = venn.inputs, intersect = ER.DE.VM_ID, 
    ##         props = props, ER = ER, DE = DE, VM = VM))
    ## }
    ## <bytecode: 0x7f8f98f827b0>

``` r
get_merged_logfc
```

    ## function (cde_object) 
    ## {
    ##     DE.short <- cde_object$deseq$short_results[[1]]
    ##     DE.short <- data.frame(logFC_DE = DE.short$logFC, ID = rownames(DE.short))
    ##     VM.short <- cde_object$voom$short_results[[1]]
    ##     VM.short <- data.frame(logFC_VM = VM.short$logFC, ID = rownames(VM.short))
    ##     ER.short <- cde_object$edger$short_results[[1]]
    ##     ER.short <- data.frame(logFC_ER = ER.short$logFC, ID = rownames(ER.short))
    ##     all_FC <- merge(merge(DE.short, VM.short, by = "ID"), ER.short, 
    ##         by = "ID")
    ##     return(all_FC)
    ## }
    ## <bytecode: 0x7f8f9e040430>
