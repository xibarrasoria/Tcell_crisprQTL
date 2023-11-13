---
title: "<span style='font-size: 16px; color: grey'> crisprQTL proof-of-concept experiment in T cells </span> <br> <span style='font-size: 36px'> Visualisation of perturbation effects </span>"
author: "Ximena Ibarra-Soria"
date: '12 December, 2022'
output:
  html_document:
    theme: paper
    keep_md: true
    fig_width: 5
    fig_height: 5
    fig_caption: yes
    df_print: paged
    code_folding: hide
    toc: true
    toc_depth: 4
    toc_float: 
      collapsed: false
---



We have determined the set of genes that change expression upon perturbation of each target locus. Here, we plot the expression changes for all significant results.


```r
## data
sce <- readRDS(paste0(dir, "results/03_sce_POC_Tcells.goodQual.calls.NORM.Rds"))

## gRNA annotation
gRNA_ann <- read.table(paste0(dir, "data/gRNA_library.tsv"), sep="\t", header = TRUE)

## DEGs
degs <- read.table(paste0(dir, "results/06_DEgenes.tsv"), header = TRUE)
# all DE testing results
res <- read.table(paste0(dir, "results/06_DEresults_targetLevel.tsv"), header = TRUE)
res$tier <- degs[match(res$pair_tgt, degs$pair_tgt),'tier']
res$distance <- degs[match(res$pair_tgt, degs$pair_tgt),'distance']

# define NT cells for plotting
nt_guides <- gRNA_ann[gRNA_ann$class=="NT",]$ID
cells <- colnames(sce[,sce$call == "unique"])
cells.neg <- apply(counts(altExp(sce, 'gRNA_calls'))[nt_guides, cells], 1, function(x) names(which(x)))
cells.neg <- unlist(cells.neg)
# use 5K
set.seed(749)
cells.neg <- sample(cells.neg, size=5e3)
```


### True positive effects {.tabset}

First, we plot the expression of the genes that are expected to be downregulated, from positive control perturbations.

#### TSS


```r
selected <- setdiff(unique(gRNA_ann[gRNA_ann$class=="TSS",]$target), c("IL12A_TSS", "TBX1_TSS")) # IL12A and TBX1 are not expressed

## plot every target
plots <- list()
for(target in selected){
  # expected gene
  gene <- unique(gRNA_ann[gRNA_ann$target==target,]$expected_DE_gene)

  # DE results
  hit <- res[res$target==target & res$gene_id==gene,]

  # plot
  plots[[target]] <- plot_gene_sc(sce = sce, cells.neg = cells.neg, hit = hit, thr=0.05)
}
ggarrange(plotlist = plots, ncol=2, nrow=17, align="hv")
```

![](07_visualise_perturbation_effects_files/figure-html/TP_tss-1.png)<!-- -->


#### LCR


```r
selected <- unique(gRNA_ann[gRNA_ann$class=="LCR",]$target)

## plot every target
plots <- list()
for(target in selected){
  # expected gene
  gene <- unique(gRNA_ann[gRNA_ann$target==target,]$expected_DE_gene)
  
  # DE results
  hit <- res[res$target==target & res$gene_id==gene,]
    
  # plot
  plots[[target]] <- plot_gene_sc(sce = sce, cells.neg = cells.neg, hit = hit, thr=0.05)
}
ggarrange(plotlist = plots, ncol=2, nrow=2, align="hv")
```

![](07_visualise_perturbation_effects_files/figure-html/TP_LCR-1.png)<!-- -->


#### ENH


```r
selected <- setdiff(unique(gRNA_ann[gRNA_ann$class=="ENH",]$target), c("USP6NL_ENH", "TUBB2A_ENH")) # USP6NL and TUBB2A are not expressed

## plot every target
plots <- list()
for(target in selected){
  # expected gene
  gene <- unique(gRNA_ann[gRNA_ann$target==target,]$expected_DE_gene)

  # DE results
  hit <- res[res$target==target & res$gene_id==gene,]
    
  # plot
  plots[[target]] <- plot_gene_sc(sce = sce, cells.neg = cells.neg, hit = hit, thr=0.05)
}
ggarrange(plotlist = plots, ncol=2, nrow=13, align="hv")
```

![](07_visualise_perturbation_effects_files/figure-html/TP_enh-1.png)<!-- -->



### Additional hits {.tabset}

Below are any other **high- and medium-confidence hits**, that are not the expected gene for intergenic perturbations.

#### ENH


```r
## plot every target
selected <- degs[degs$class=="ENH" & !degs$expected & degs$tier!="low",]

plots <- list()
for(i in 1:nrow(selected)){
  hit <- selected[i,]

  plots[[i]] <- plot_gene_sc(sce = sce, cells.neg = cells.neg, hit = hit, thr=0.05)
}
ggarrange(plotlist = plots, ncol=2, nrow=21, align="hv")
```

![](07_visualise_perturbation_effects_files/figure-html/other_enh-1.png)<!-- -->


#### INTRON


```r
## plot every target
selected <- degs[degs$class=="INTRON" & degs$tier!="low",]

plots <- list()
for(i in 1:nrow(selected)){
  hit <- selected[i,]

  plots[[i]] <- plot_gene_sc(sce = sce, cells.neg = cells.neg, hit = hit, thr=0.05)
}
ggarrange(plotlist = plots, ncol=2, nrow=12, align="hv")
```

![](07_visualise_perturbation_effects_files/figure-html/other_intron-1.png)<!-- -->

#### INTERGENIC


```r
## plot every target
selected <- degs[degs$class=="INTERGENIC" & degs$tier!="low",]

plots <- list()
for(i in 1:nrow(selected)){
  hit <- selected[i,]

  plots[[i]] <- plot_gene_sc(sce = sce, cells.neg = cells.neg, hit = hit, thr=0.05)
}
ggarrange(plotlist = plots, ncol=2, nrow=6, align="hv")
```

![](07_visualise_perturbation_effects_files/figure-html/other_intergenic-1.png)<!-- -->

### Overlapping gene for intronic perturbations

For intronic elements, we check directly the expression of the overlapping gene. *LCA5L* is not detected in at least 5% of cells and was not tested.


```r
## plot every target
df <- data.frame(target = unique(gRNA_ann[gRNA_ann$class=="INTRON",]$target),
                 gene = unlist(lapply(strsplit(unique(gRNA_ann[gRNA_ann$class=="INTRON",]$target), "_"),'[[',1)))
df <- df[-which(df$target == "LCA5L_INTRON"),]

plots <- list()
for(i in 1:nrow(df)){
  target <- df[i,1]
  gene <- df[i,2]
  hit <- res[res$target == target & res$gene_id == gene,]
  
  plots[[i]] <- plot_gene_sc(sce = sce, cells.neg = cells.neg, hit = hit, thr=0.05)
}
ggarrange(plotlist = plots, ncol=2, nrow=5, align="hv")
```

![](07_visualise_perturbation_effects_files/figure-html/intron_gene-1.png)<!-- -->


### Nearby genes for intergenic perturbations

We also plot the expression of the four closest genes (detected in at least 5% of the cells) for `INTERGENIC` targets.


```r
## neighbourhood
gene_distance_target <- read.table(paste0(dir, "results/06_gene_target_distance.tsv"), header=TRUE)

## plot every target
df <- gene_distance_target[grep("INTERGENIC", gene_distance_target$target),]

plots <- list()
for(target in unique(gRNA_ann[gRNA_ann$class=="INTERGENIC",]$target)){
  selected <- df[which(df$target==target & df$expressed),][1:4,]
  for(i in 1:nrow(selected)){
    hit <- res[res$target == target & res$gene_id == selected[i,'external_gene_name'],]
    plots[[paste0(target,i)]] <- plot_gene_sc(sce = sce, cells.neg = cells.neg, hit = hit, thr=0.05)
  }
}
ggarrange(plotlist = plots, ncol=2, nrow=6, align="hv")
```

![](07_visualise_perturbation_effects_files/figure-html/intergenic_nearby-1.png)<!-- -->



```r
sessionInfo()
```

```
## R version 4.0.2 (2020-06-22)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: CentOS Linux 7 (Core)
## 
## Matrix products: default
## BLAS:   /usr/lib64/libblas.so.3.4.2
## LAPACK: /hpc/apps/2018/R/v4.0.2.app/lib/R/lib/libRlapack.so
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] ggpubr_0.4.0                ggplot2_3.3.5              
##  [3] scales_1.1.1                scran_1.16.0               
##  [5] SingleCellExperiment_1.10.1 SummarizedExperiment_1.18.2
##  [7] DelayedArray_0.14.1         matrixStats_0.60.1         
##  [9] Biobase_2.48.0              GenomicRanges_1.40.0       
## [11] GenomeInfoDb_1.24.2         IRanges_2.22.2             
## [13] S4Vectors_0.26.1            BiocGenerics_0.34.0        
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-7              tools_4.0.2              
##  [3] backports_1.2.1           bslib_0.3.0              
##  [5] utf8_1.2.2                R6_2.5.1                 
##  [7] irlba_2.3.3               vipor_0.4.5              
##  [9] DBI_1.1.1                 colorspace_2.0-2         
## [11] withr_2.4.2               tidyselect_1.1.1         
## [13] gridExtra_2.3             curl_4.3.2               
## [15] compiler_4.0.2            BiocNeighbors_1.6.0      
## [17] labeling_0.4.2            sass_0.4.0               
## [19] stringr_1.4.0             digest_0.6.27            
## [21] foreign_0.8-81            rmarkdown_2.7            
## [23] rio_0.5.26                XVector_0.28.0           
## [25] scater_1.16.2             pkgconfig_2.0.3          
## [27] htmltools_0.5.2           highr_0.8                
## [29] fastmap_1.1.0             limma_3.44.3             
## [31] readxl_1.3.1              rlang_0.4.11             
## [33] DelayedMatrixStats_1.10.1 farver_2.1.0             
## [35] jquerylib_0.1.4           generics_0.1.0           
## [37] jsonlite_1.7.2            BiocParallel_1.22.0      
## [39] zip_2.1.1                 dplyr_1.0.7              
## [41] car_3.0-10                RCurl_1.98-1.3           
## [43] magrittr_2.0.1            BiocSingular_1.4.0       
## [45] GenomeInfoDbData_1.2.3    Matrix_1.3-2             
## [47] Rcpp_1.0.7                ggbeeswarm_0.6.0         
## [49] munsell_0.5.0             fansi_0.5.0              
## [51] abind_1.4-5               viridis_0.6.1            
## [53] lifecycle_1.0.0           stringi_1.5.3            
## [55] yaml_2.2.1                edgeR_3.30.3             
## [57] carData_3.0-4             zlibbioc_1.34.0          
## [59] grid_4.0.2                dqrng_0.3.0              
## [61] forcats_0.5.1             crayon_1.4.1             
## [63] lattice_0.20-41           cowplot_1.1.1            
## [65] haven_2.3.1               hms_1.1.0                
## [67] locfit_1.5-9.4            knitr_1.31               
## [69] pillar_1.6.2              igraph_1.2.6             
## [71] ggsignif_0.6.1            glue_1.4.2               
## [73] evaluate_0.14             data.table_1.14.0        
## [75] vctrs_0.3.8               cellranger_1.1.0         
## [77] gtable_0.3.0              purrr_0.3.4              
## [79] tidyr_1.1.3               assertthat_0.2.1         
## [81] openxlsx_4.2.3            xfun_0.22                
## [83] rsvd_1.0.5                broom_0.7.6              
## [85] rstatix_0.7.0             viridisLite_0.4.0        
## [87] tibble_3.1.4              beeswarm_0.4.0           
## [89] statmod_1.4.36            ellipsis_0.3.2
```

