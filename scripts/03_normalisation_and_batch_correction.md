---
title: "<span style='font-size: 16px; color: grey'> crisprQTL proof-of-concept experiment in T cells </span> <br> <span style='font-size: 36px'> Normalisation and batch correction </span>"
author: "Ximena Ibarra-Soria"
date: '22 July, 2022'
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
    toc_depth: 3
    toc_float: 
      collapsed: false
---




We start from the normalised data for ~150K good-quality cells with at least one gRNA detected.


```r
sce <- readRDS(paste0(dir, "results/02_sce_POC_Tcells.goodQual.calls.Rds"))
ncol(sce)
```

```
## [1] 152403
```

### Normalisation

We use `scran` to normalise differences in sequencing depth between cells and other biases.


```r
# pre-cluster the data to protect the size factor estimation from too many DEGs
clusters  <- quickCluster(sce, min.size = 100, method = "igraph", BPPARAM = bpp)
# estimate size factors
sce <- computeSumFactors(sce, clusters = clusters, min.mean = 0.1, BPPARAM = bpp)

par(mfrow=c(1,2))
plot(sce$sizeFactor, 10^sce$libSize/1e6, 
     pch=16, 
     xlab="size factors", ylab="library size (millions)", 
     bty="l")
abline(lm((10^sce$libSize/1e6)~sce$sizeFactor))

## normalise
sce <- logNormCounts(sce)
```

![](03_normalisation_and_batch_correction_files/figure-html/normalisation-1.png)<!-- -->


### Visualisation

To assess the success of normalisation and QC we perform dimensionality reduction using UMAP. We observe some effects from library size, but these are likely correlated with differences in cell cycle phase, instead of failures in the normalisation. 


```r
## HVGs
var_decomp <- modelGeneVar(sce, BPPARAM = bpp)
hvgs <- getTopHVGs(var_decomp, n=2000)
rowData(sce)$hvg <- row.names(sce) %in% hvgs

## PCA
set.seed(497)
sce <- runPCA(sce, subset_row = hvgs, BSPARAM=IrlbaParam())

## UMAP
sce <- runUMAP(sce, dimred = "PCA", BPPARAM = bpp)

## plot
umap <- as.data.frame(reducedDim(sce, 'UMAP'))
colnames(umap) <- c("x", "y")
umap <- cbind(umap, as.data.frame(colData(sce)))

o <- sample(1:nrow(umap))
plots <- list()
plots[[1]] <- ggplot(umap[o,], aes(x, y, colour=libSize)) +
  geom_point(size=0.2) +
  scale_color_gradientn(colours = rev(brewer.pal(n=9, "RdYlBu"))) +
  xlab("UMAP1") + ylab("UMAP2") +
  ggtitle("library size") +
  th + theme(axis.text.x = element_blank(),
             axis.text.y = element_blank(),
             axis.ticks.x = element_blank(),
             axis.ticks.y = element_blank(),
             legend.position = "bottom")
plots[[2]] <- ggplot(umap[o,], aes(x, y, colour=mt)) +
  geom_point(size=0.2) +
  scale_color_gradientn(colours = rev(brewer.pal(n=9, "RdYlBu"))) +
  xlab("UMAP1") + ylab("UMAP2") +
  ggtitle("mt") +
  th + theme(axis.text.x = element_blank(),
             axis.text.y = element_blank(),
             axis.ticks.x = element_blank(),
             axis.ticks.y = element_blank(),
             legend.position = "bottom")

ggarrange(plotlist = plots, ncol=2, nrow=1, align="hv", legend = "bottom")
```

![](03_normalisation_and_batch_correction_files/figure-html/umap-1.png)<!-- -->
 
To confirm this, we use `tricycle` to infer the cell cycle phase of each cell. Indeed, the effects from library size observed on the UMAP correspond to differences of the cells along the cell cycle, with smaller libraries mapping to cells in G1.


```r
## project data to reference
sce <- project_cycle_space(sce, gname.type = 'SYMBOL', species = 'human')

## estimate cell cycle score
sce <- estimate_cycle_position(sce)

## plot
p <- plot_emb_circle_scale(sce, 
                           dimred = 'UMAP',
                           point.size = 0.5, point.alpha = 0.9)
legend <- circle_scale_legend(text.size = 5, alpha = 0.9)
cowplot::plot_grid(p, legend, ncol = 2, rel_widths = c(1, 0.4))
```

![](03_normalisation_and_batch_correction_files/figure-html/tricycle-1.png)<!-- -->

We do, however, observe strong batch effects from the two different experiments.


```r
## plot by batch+sample
plots <- list()
# batch1
tmp <- umap[umap$batch == "batch1",]
o <- sample(1:nrow(tmp))
plots[[1]] <- ggplot(tmp[o,], aes(x, y, colour=Sample)) +
  geom_point(size=0.2) +
  xlab("UMAP1") + ylab("UMAP2") +
  ggtitle("batch 1") +
  th + theme(axis.text.x = element_blank(),
             axis.text.y = element_blank(),
             axis.ticks.x = element_blank(),
             axis.ticks.y = element_blank(),
             legend.position = "none")
# batch2
tmp <- umap[umap$batch == "batch2",]
o <- sample(1:nrow(tmp))
plots[[2]] <- ggplot(tmp[o,], aes(x, y, colour=Sample)) +
  geom_point(size=0.2) +
  xlab("UMAP1") + ylab("UMAP2") +
  ggtitle("batch 2") +
  th + theme(axis.text.x = element_blank(),
             axis.text.y = element_blank(),
             axis.ticks.x = element_blank(),
             axis.ticks.y = element_blank(),
             legend.position = "none")

ggarrange(plotlist = plots, ncol=2, nrow=1, align="hv")
```

![](03_normalisation_and_batch_correction_files/figure-html/batch-1.png)<!-- -->

### Batch correction

We perform batch correction using `harmony` to account for the differences between the two batches. We recompute the UMAP and observe successful batch correction. 


```r
## remove batch effect from the two experiments
harmony_embeddings <- HarmonyMatrix(data_mat = reducedDim(sce, 'PCA'),
                                    meta_data = sce$batch,
                                    do_pca = FALSE,
                                    plot_convergence = FALSE,
                                    verbose = FALSE) # 11.6 min
## add to sce and recompute UMAP
reducedDim(sce, 'harmony') <- harmony_embeddings
sce <- runUMAP(sce, dimred = "harmony", name = 'UMAP', BPPARAM = bpp) # 13min

## get data
umap <- as.data.frame(reducedDim(sce, 'UMAP'))
colnames(umap) <- c("x", "y")
umap <- cbind(umap, as.data.frame(colData(sce)))

## plot
o <- sample(1:nrow(umap))
ggplot(umap[o,], aes(x, y, colour=batch)) +
  geom_point(size=0.2) +
  xlab("UMAP1") + ylab("UMAP2") +
  th + theme(axis.text.x = element_blank(),
             axis.text.y = element_blank(),
             axis.ticks.x = element_blank(),
             axis.ticks.y = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=1)))
```

![](03_normalisation_and_batch_correction_files/figure-html/harmony-1.png)<!-- -->

We now have normalised data and a batch-corrected embedding for visualisation in downstream analyses. We will have to account for the differences in batches in all the perturbation analyses.


```r
## save
saveRDS(sce, paste0(dir, "results/03_sce_POC_Tcells.goodQual.calls.NORM.Rds"))
```





```r
sessionInfo()
```

```
## R version 4.1.1 (2021-08-10)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: CentOS Linux 7 (Core)
## 
## Matrix products: default
## BLAS:   /usr/lib64/libblas.so.3.4.2
## LAPACK: /hpc/apps/2018/R/v4.1.app/lib/R/lib/libRlapack.so
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
##  [1] grid      parallel  stats4    stats     graphics  grDevices utils    
##  [8] datasets  methods   base     
## 
## other attached packages:
##  [1] ComplexHeatmap_2.8.0        inlmisc_0.5.2              
##  [3] RColorBrewer_1.1-2          ggpubr_0.4.0               
##  [5] harmony_0.1.0               Rcpp_1.0.7                 
##  [7] BiocSingular_1.8.1          BiocParallel_1.26.1        
##  [9] tricycle_1.0.0              scran_1.20.1               
## [11] scater_1.20.1               ggplot2_3.3.5              
## [13] scuttle_1.2.1               SingleCellExperiment_1.14.1
## [15] SummarizedExperiment_1.22.0 Biobase_2.52.0             
## [17] GenomicRanges_1.44.0        GenomeInfoDb_1.28.1        
## [19] IRanges_2.26.0              S4Vectors_0.30.0           
## [21] BiocGenerics_0.38.0         MatrixGenerics_1.4.2       
## [23] matrixStats_0.60.0         
## 
## loaded via a namespace (and not attached):
##   [1] circlize_0.4.13           readxl_1.3.1             
##   [3] backports_1.2.1           igraph_1.2.6             
##   [5] sp_1.4-5                  scattermore_0.7          
##   [7] digest_0.6.27             foreach_1.5.1            
##   [9] htmltools_0.5.1.1         viridis_0.6.1            
##  [11] fansi_0.5.0               magrittr_2.0.1           
##  [13] memoise_2.0.0             ScaledMatrix_1.0.0       
##  [15] cluster_2.1.2             doParallel_1.0.16        
##  [17] openxlsx_4.2.4            limma_3.48.3             
##  [19] Biostrings_2.60.2         colorspace_2.0-2         
##  [21] blob_1.2.2                haven_2.4.3              
##  [23] xfun_0.25                 dplyr_1.0.7              
##  [25] rgdal_1.5-23              crayon_1.4.1             
##  [27] RCurl_1.98-1.4            jsonlite_1.7.2           
##  [29] iterators_1.0.13          glue_1.4.2               
##  [31] gtable_0.3.0              zlibbioc_1.38.0          
##  [33] XVector_0.32.0            GetoptLong_1.0.5         
##  [35] DelayedArray_0.18.0       car_3.0-11               
##  [37] shape_1.4.6               abind_1.4-5              
##  [39] scales_1.1.1              mvtnorm_1.1-2            
##  [41] DBI_1.1.1                 edgeR_3.34.0             
##  [43] rstatix_0.7.0             viridisLite_0.4.0        
##  [45] clue_0.3-59               dqrng_0.3.0              
##  [47] foreign_0.8-81            bit_4.0.4                
##  [49] rsvd_1.0.5                metapod_1.0.0            
##  [51] httr_1.4.2                ellipsis_0.3.2           
##  [53] farver_2.1.0              pkgconfig_2.0.3          
##  [55] uwot_0.1.10               sass_0.4.0               
##  [57] locfit_1.5-9.4            utf8_1.2.2               
##  [59] labeling_0.4.2            tidyselect_1.1.1         
##  [61] rlang_0.4.11              AnnotationDbi_1.54.1     
##  [63] munsell_0.5.0             cellranger_1.1.0         
##  [65] tools_4.1.1               cachem_1.0.5             
##  [67] generics_0.1.0            RSQLite_2.2.7            
##  [69] broom_0.7.9               evaluate_0.14            
##  [71] stringr_1.4.0             fastmap_1.1.0            
##  [73] yaml_2.2.1                knitr_1.33               
##  [75] bit64_4.0.5               zip_2.2.0                
##  [77] purrr_0.3.4               KEGGREST_1.32.0          
##  [79] sparseMatrixStats_1.4.2   compiler_4.1.1           
##  [81] beeswarm_0.4.0            curl_4.3.2               
##  [83] png_0.1-7                 ggsignif_0.6.2           
##  [85] tibble_3.1.3              statmod_1.4.36           
##  [87] bslib_0.2.5.1             stringi_1.7.3            
##  [89] highr_0.9                 RSpectra_0.16-0          
##  [91] forcats_0.5.1             lattice_0.20-44          
##  [93] circular_0.4-93           bluster_1.2.1            
##  [95] Matrix_1.3-4              vctrs_0.3.8              
##  [97] pillar_1.6.2              lifecycle_1.0.0          
##  [99] GlobalOptions_0.1.2       jquerylib_0.1.4          
## [101] RcppAnnoy_0.0.19          BiocNeighbors_1.10.0     
## [103] data.table_1.14.0         cowplot_1.1.1            
## [105] bitops_1.0-7              irlba_2.3.3              
## [107] raster_3.4-13             R6_2.5.1                 
## [109] gridExtra_2.3             rio_0.5.27               
## [111] vipor_0.4.5               codetools_0.2-18         
## [113] boot_1.3-28               assertthat_0.2.1         
## [115] rjson_0.2.20              withr_2.4.2              
## [117] GenomeInfoDbData_1.2.6    hms_1.1.0                
## [119] beachmat_2.8.1            tidyr_1.1.3              
## [121] rmarkdown_2.10            DelayedMatrixStats_1.14.2
## [123] carData_3.0-4             Cairo_1.5-12.2           
## [125] ggbeeswarm_0.6.0
```

