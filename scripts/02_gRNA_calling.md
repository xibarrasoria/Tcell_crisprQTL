---
title: "<span style='font-size: 16px; color: grey'> crisprQTL proof-of-concept experiment in T cells </span> <br> <span style='font-size: 36px'> gRNA calling </span>"
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



We start from the set of ~250K good-quality cells retained after quality control. 


```r
sce <- readRDS(paste0(dir, "results/01_sce_POC_Tcells.goodQual.RAW.Rds"))

c(n_cells=ncol(sce))
```

```
## n_cells 
##  250195
```

As the next step, we need to determine which gRNAs are present in each cell.

### gRNA detection in cDNA library

Before performing gRNA calling, we check the concordance in gRNA transcript detection from the cDNA and gRNA libraries, to ensure the PCR enrichment process doesn't introduce significant biases. The cDNA data was mapped against a modified reference genome containing the gRNA sequences from the plasmid library, allowing us to recover reads from gRNA transcripts directly.

We compare the readouts from both the cDNA and gRNA fractions. As expected, the vast majority of cells have no detection in the cDNA fraction. However, the majority of the cells with detection in both libraries identify the same gRNA as the most abundant one. The discordant calls are an upper bound, since some of these will be concordant after accounting for multiple gRNAs per cell (currently only looking at whether the highest detected is the same in both libraries).


```r
## matched counts for gRNAs in the `gRNA` and `cDNA` fractions in the two altExp slots
## find which is the most abundant guide in each
## if there are zero counts, 1 is returned by which.max
df <- data.frame(top_gRNA = apply(counts(altExp(sce, 'CRISPR')), 2, which.max),
                 top_cDNA = apply(counts(altExp(sce, 'gRNA_in_cDNA')), 2, which.max))
df$total_gRNA <- colSums(counts(altExp(sce, 'CRISPR')))
df$total_cDNA <- colSums(counts(altExp(sce, 'gRNA_in_cDNA')))
df$class <- ifelse(df$top_cDNA==1 & df$total_cDNA==0, "not_detected_cDNA",
                   ifelse(df$top_gRNA==1 & df$top_gRNA==0, "not_detected_gRNA",
                          ifelse(df$top_cDNA == df$top_gRNA, "concordant", "discordant")))
df[df$class=="not_detected_cDNA",]$top_cDNA <- -10

plots <- list()
plots[[1]] <- ggplot(df, aes(top_gRNA, top_cDNA, colour=class)) +
  geom_point(size=1, alpha = 0.1) +
  scale_color_manual(values = as.character(GetColors(3, scheme = "light"))) +
  xlab("gRNAs 1 --> 355 - gRNA library") +
  ylab("gRNAs 1 --> 355 - cDNA library") +
  guides(colour = guide_legend(override.aes = list(alpha=1))) +
  th + theme(axis.text.x = element_blank(),
             axis.text.y = element_blank(),
             axis.ticks = element_blank())

plots[[2]] <- ggplot(df, aes(1, fill=class)) +
  geom_bar() +
  scale_fill_manual(values = as.character(GetColors(3, scheme = "light"))) +
  xlab("") + ylab("number of cells") +
  annotate("text", x=1, y=c(0.69e5, 1.44e5, 2e5), 
           label=c("not_detected_cDNA", 
                   "discordant", "concordant"), size=4) +
  th + theme(axis.text.x = element_blank(),
             axis.ticks.x = element_blank())

ggarrange(plotlist = plots, ncol=2, widths = c(0.65, 0.35), 
          legend="none", align = "h")
```

![](02_gRNA_calling_files/figure-html/compare_grna_counts-1.png)<!-- -->


```r
round(prop.table(table(df$class))*100, 2)
```

```
## 
##        concordant        discordant not_detected_cDNA 
##             40.14              4.56             55.30
```

This indicates that the PCR amplification step doesn't introduce spurious signal, but it greatly enhances gRNA transcript detection.

### gRNA calling 

We use a binomial test to determine which gRNAs are likely present in each cell. 


```r
call_gRNAs_binom <- function(sce=NULL, pval_threshold=0.001, 
                             max_to_test=5, total_guides=NULL,
                             prob=NULL,
                             min_count=3, min_fraction=0.05, 
                             plot=TRUE, res=200,
                             ann=NULL, columns=NULL, add_stats=TRUE){
  if(is.null(total_guides)){
    stop("Please provide the total number of gRNAs to use in the binomial test. 
         This should be the total number of gRNAs present in the library.")
  }
  
  ## create data table with detected gRNA-cell combinations
  gRNA_calls <- reshape2::melt(as.matrix(counts(altExp(sce, 'CRISPR'))))
  ## remove guide-cell pairs with no counts
  gRNA_calls <- gRNA_calls[gRNA_calls$value>0,]
  colnames(gRNA_calls) <- c("gRNA", "cell", "UMI_count")
  
  ## remove guides with too few UMIs
  gRNA_calls <- gRNA_calls[gRNA_calls$UMI_count > min_count,]
  
  ## annotate the targets associated with each gRNA
  if(!is.null(ann)){
    if(length(setdiff(columns, colnames(ann)))>0){
      stop(paste("Column(s)", setdiff(columns, colnames(ann)), "are not present in annotation data.frame"))
    }
    for(i in seq_along(columns)){
      gRNA_calls[,columns[i]] <- ann[match(gRNA_calls$gRNA, ann$ID), columns[i]]
    }
  }
  
  ## add total UMIs per cell, to restrict analysis to the largest counts per cell
  totals <- colSums(counts(altExp(sce, 'CRISPR')))
  gRNA_calls$cell_total_UMIs <- totals[gRNA_calls$cell]
  # compute % of total UMIs taken by each gRNA
  gRNA_calls$fraction_total_UMIs <- gRNA_calls$UMI_count/gRNA_calls$cell_total_UMIs
  
  ## convert to data.table
  gRNA_calls <- as.data.table(gRNA_calls)
  
  ## order decreasingly by UMI_count
  setorder(gRNA_calls, -UMI_count)
  # add column with order per cell based on UMI abundance
  gRNA_calls[, order := 1:.N, by = c("cell")] 
  
  # test top gRNAs per cell as determined by max_to_test
  # if more than max_to_test integrations are true, the rest won't be called
  # higher number in max_to_test results in higher computation time
  
  ## use a binomial test to determine if observed counts are greater than background
  # if no probabilitis provided, use equal probabilities
  if(is.null(prob)){
    gRNA_calls$prob <- 1/total_guides
  }else{
    # add proportions of guides in guide library to gRNA_calls
    # provided as a named vector, with names being gRNA ids
    prob <- as.data.table(prob, keep.rownames = TRUE)
    gRNA_calls <- merge(gRNA_calls, prob, by.x="gRNA", by.y="rn")
  }
  gRNA_calls[order <= max_to_test,
             binom := binom.test(UMI_count, cell_total_UMIs, prob, "greater")$p.value, 
             by = seq_len(nrow(gRNA_calls[order <= max_to_test]))]
  # bonferroni correction, per cell
  gRNA_calls[, p_adj := binom * .N,  by = "cell"] 
  gRNA_calls[p_adj > 1, p_adj := 1] # reset pvals that are larger than 1 after correction
  
  # count significant results
  gRNA_calls[, n_detected_gRNAs := length(which(p_adj <= pval_threshold)), by = "cell"]
  
  ## check if there are cells with at least max_to_test called hits
  # if so, retest a larger number
  while(max(gRNA_calls$n_detected_gRNAs) >= max_to_test){
    # get cells that need higher threshold
    done <- gRNA_calls[gRNA_calls$n_detected_gRNAs < max_to_test,]
    missing <- gRNA_calls[gRNA_calls$n_detected_gRNAs >= max_to_test,]
    
    # compute pvals for those cells with twice as many cells tested
    max_to_test <- max_to_test*2
    missing[order <= max_to_test,
            binom := binom.test(UMI_count, cell_total_UMIs, (1 / total_guides), "greater")$p.value, 
            by = seq_len(nrow(missing[order <= max_to_test]))]
    missing[, p_adj := binom * .N,  by = "cell"] 
    missing[p_adj > 1, p_adj := 1]
    
    # join with previous calls
    gRNA_calls <- rbind(done, missing)
    gRNA_calls[, n_detected_gRNAs := length(which(p_adj <= pval_threshold)), by = "cell"]
  }
  
  # for unique calls, also remove those where the fraction from total UMIs is too low
  gRNA_calls[gRNA_calls$n_detected_gRNAs == 1 &
               gRNA_calls$fraction_total_UMIs <= min_fraction,]$p_adj <- 1
  # count significant results again after adjustments above
  gRNA_calls[, n_detected_gRNAs := length(which(p_adj <= pval_threshold)), by = "cell"]
  
  # assign call based on number of significant gRNAs
  gRNA_calls$call <- ifelse(gRNA_calls$n_detected_gRNAs == 0, "unassigned",
                            ifelse(gRNA_calls$n_detected_gRNAs == 1, "unique", "multiple"))
  
  # annotate significant results
  gRNA_calls[, detected_gRNAs := .(paste0(gRNA[which(p_adj <= pval_threshold)], collapse="|")), by = "cell"]
  gRNA_calls[, detected_UMIs := .(paste0(UMI_count[which(p_adj <= pval_threshold)], collapse="|")), by = "cell"]
  gRNA_calls[, detected_fraction := .(paste0(fraction_total_UMIs[which(p_adj <= pval_threshold)], collapse="|")), by = "cell"]
  gRNA_calls[, detected_padj := .(paste0(p_adj[which(p_adj <= pval_threshold)], collapse="|")), by = "cell"]
  
  ## add annotation if required
  if(!is.null(ann)){
    for(i in seq_along(columns)){
      gRNA_calls[, tmp := .(paste0(get(columns[i])[which(p_adj <= pval_threshold)], collapse="|")), by = "cell"]
      colnames(gRNA_calls)[ncol(gRNA_calls)] <- paste0('detected_', columns[i])
    }
  }
  
  ## QC plot
  if(plot){
    ## remove overlapping points
    idx <- subsetPointsByGrid(X = gRNA_calls$fraction_total_UMIs,
                              Y = log10(gRNA_calls$UMI_count),
                              grouping = gRNA_calls$col,
                              resolution = res)
    df <- gRNA_calls[idx,]
    df$call <- factor(df$call, levels=c("unassigned", "multiple", "unique"))
    df <- df[order(df$call),]
    
    ## assign colors to unique/multiple/unassigned
    cols <- c(multiple="indianred3", unassigned="grey", unique="grey10")
    df$col <- cols[match(df$call, names(cols))]
    
    ## plot 
    p <- ggplot(df, aes(fraction_total_UMIs, log10(UMI_count), colour=call)) +
      scale_color_manual(values = cols) +
      geom_point(size=0.25) +
      xlab("% of total UMIs in cell") +
      ylab(expression('log'[10]*' UMIs in gRNA')) +
      geom_hline(yintercept = log10(min_count), lty=2, col="grey10") +
      geom_vline(xintercept = min_fraction, lty=2, col="grey10") +
      th
  }
  
  ## produce sparse matrix of cell-gRNA pairs that are significant
  # subset to significant
  gRNA_calls_sig <- gRNA_calls[gRNA_calls$p_adj < pval_threshold,]
  # gRNA indices
  rows <- match(gRNA_calls_sig$gRNA, rownames(altExp(sce, 'CRISPR')))
  # cell indices
  cols <- match(gRNA_calls_sig$cell, colnames(altExp(sce, 'CRISPR')))
  # create sparse matrix
  gRNA_cell_matrix <- sparseMatrix(i = rows, j = cols, dims = dim(altExp(sce, 'CRISPR')))
  row.names(gRNA_cell_matrix) <- row.names(altExp(sce, 'CRISPR'))
  colnames(gRNA_cell_matrix) <- colnames(altExp(sce, 'CRISPR'))
  # add to sce obect, copying over col/rowData from CRISPR counts
  altExp(sce, 'gRNA_calls') <- altExp(sce, 'CRISPR')
  assay(altExp(sce, 'gRNA_calls')) <- gRNA_cell_matrix
  
  ## produce matrix of assignments
  gRNA_calls <- gRNA_calls[order == 1, ] # keep only one entry per cell
  gRNA_calls <- gRNA_calls[, c("gRNA", "UMI_count", "cell_total_UMIs", "fraction_total_UMIs", "order", "binom", "p_adj"):=NULL]
  if(!is.null(columns)){
    gRNA_calls <- gRNA_calls[, (columns):=NULL]
  }
  ## add to sce colData
  if(add_stats == FALSE){
    # only add number of gRNAs called, and call annotation
    sel <- c("n_detected_gRNAs", "call", "detected_gRNAs")
    colData(sce) <- cbind(colData(sce), gRNA_calls[match(colnames(sce), gRNA_calls$cell), ..sel])
  }else{
    # add everything
    colData(sce) <- cbind(colData(sce), gRNA_calls[match(colnames(sce), gRNA_calls$cell),-1])
  }
  # remove NAs produced by cells with _no_ gRNA counts at all (missing from analysis above)
  if(sum(is.na(sce$n_detected_gRNAs))>0){
    colData(sce[,is.na(sce$n_detected_gRNAs)])[,'call'] <- "unassigned"
    colData(sce[,is.na(sce$n_detected_gRNAs)])[,'n_detected_gRNAs'] <- 0
  }
  
  ## return results
  if(plot == TRUE){
    return(list(sce = sce,
                calls = gRNA_calls,
                plot = p))
  }else{
    return(list(sce = sce,
                calls = gRNA_calls))
  }
}
```

We perform the test using the probability of selecting each gRNA based on their frequencies in the initial plasmid library. The vast majority of cells have a single gRNA detected, consistent with this being a low MOI experiment.


```r
## gRNA anntation - same library as previous pilots
gRNA_ann <- read.table(paste0(dir, "data/gRNA_library.tsv"), sep="\t", header = TRUE)

## gRNA representation in initial plasmid library
gRNA_prob <- read.table(paste0(dir, "data/gRNA_count_plasmid_library.tsv"), header = TRUE)
probs <- gRNA_prob$Count/sum(gRNA_prob$Count)
names(probs) <- gRNA_prob$gRNA_ID

## call gRNAs for each sample
sce.list <- sapply(unique(sce$Sample), function(x) sce[,sce$Sample == x])
gRNA_calls <- bplapply(sce.list, 
                       function(x) call_gRNAs_binom(sce = x, 
                                                    pval_threshold = 0.001,
                                                    prob = probs,
                                                    max_to_test = 5,
                                                    total_guides = nrow(altExp(sce)),
                                                    plot = TRUE, 
                                                    ann = gRNA_ann, 
                                                    columns = c("target", "expected_DE_gene", "class"),
                                                    min_count = 3, min_fraction = 0.05,
                                                    add_stats = FALSE),
                       BPPARAM = bpp)

## merge all samples
sce <- do.call(cbind, lapply(gRNA_calls, '[[', 1))
plots <- lapply(gRNA_calls, '[[', 3)
gRNA_calls <- lapply(gRNA_calls, '[[', 2)
gRNA_calls <- do.call(rbind, gRNA_calls)

## plot
ggarrange(plotlist = plots, ncol=5, nrow=9,
          labels = unique(sce$Sample), font.label = list(size=8),
          common.legend = TRUE, legend = "bottom")
```

![](02_gRNA_calling_files/figure-html/call_gRNAs_cells-1.png)<!-- -->

Overall, around 63% of cells have a gRNA detected.


```r
# proportion of cells with each call type
df <- round(prop.table(table(sce$Sample, sce$call),1)*100, 2)

par(mar=c(6,4,2,2))
barplot(t(df[,c(3,1,2)]), col=c("grey10", "indianred3", "grey"), 
        las=2, cex.names = 0.7,
        width = 0.8, space = 0.2, xlim=c(0,40))
legend("topright", legend = c("unassigned", "multiple", "unique"),
       col = c("grey", "indianred3", "grey10"), pch = 15, cex=0.7)
```

![](02_gRNA_calling_files/figure-html/prop_unique-1.png)<!-- -->

Moving forward we only use cells with a gRNA assignment. This leaves us with this many cells:


```r
sce <- sce[,sce$call != "unassigned"]
c(n_cells=ncol(sce))
```

```
## n_cells 
##  152403
```



```r
saveRDS(sce, paste0(dir, "results/02_sce_POC_Tcells.goodQual.calls.Rds"))
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
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] iSEE_2.4.0                  inlmisc_0.5.2              
##  [3] viridis_0.6.1               viridisLite_0.4.0          
##  [5] BuenColors_0.5.6            MASS_7.3-54                
##  [7] RColorBrewer_1.1-2          ggpubr_0.4.0               
##  [9] BiocSingular_1.8.1          BiocParallel_1.26.1        
## [11] Matrix_1.3-4                data.table_1.14.0          
## [13] scran_1.20.1                scater_1.20.1              
## [15] ggplot2_3.3.5               scuttle_1.2.1              
## [17] SingleCellExperiment_1.14.1 SummarizedExperiment_1.22.0
## [19] Biobase_2.52.0              GenomicRanges_1.44.0       
## [21] GenomeInfoDb_1.28.1         IRanges_2.26.0             
## [23] S4Vectors_0.30.0            BiocGenerics_0.38.0        
## [25] MatrixGenerics_1.4.2        matrixStats_0.60.0         
## 
## loaded via a namespace (and not attached):
##   [1] readxl_1.3.1              backports_1.2.1          
##   [3] circlize_0.4.13           igraph_1.2.6             
##   [5] shinydashboard_0.7.1      splines_4.1.1            
##   [7] sp_1.4-5                  digest_0.6.27            
##   [9] foreach_1.5.1             htmltools_0.5.1.1        
##  [11] fansi_0.5.0               checkmate_2.0.0          
##  [13] magrittr_2.0.1            ScaledMatrix_1.0.0       
##  [15] cluster_2.1.2             doParallel_1.0.16        
##  [17] openxlsx_4.2.4            limma_3.48.3             
##  [19] ComplexHeatmap_2.8.0      colorspace_2.0-2         
##  [21] ggrepel_0.9.1             haven_2.4.3              
##  [23] xfun_0.25                 dplyr_1.0.7              
##  [25] rgdal_1.5-23              crayon_1.4.1             
##  [27] RCurl_1.98-1.4            jsonlite_1.7.2           
##  [29] iterators_1.0.13          glue_1.4.2               
##  [31] gtable_0.3.0              zlibbioc_1.38.0          
##  [33] XVector_0.32.0            GetoptLong_1.0.5         
##  [35] DelayedArray_0.18.0       car_3.0-11               
##  [37] shape_1.4.6               abind_1.4-5              
##  [39] scales_1.1.1              DBI_1.1.1                
##  [41] edgeR_3.34.0              rstatix_0.7.0            
##  [43] miniUI_0.1.1.1            Rcpp_1.0.7               
##  [45] xtable_1.8-4              clue_0.3-59              
##  [47] dqrng_0.3.0               foreign_0.8-81           
##  [49] rsvd_1.0.5                DT_0.18                  
##  [51] metapod_1.0.0             htmlwidgets_1.5.3        
##  [53] shinyAce_0.4.1            ellipsis_0.3.2           
##  [55] farver_2.1.0              pkgconfig_2.0.3          
##  [57] sass_0.4.0                locfit_1.5-9.4           
##  [59] utf8_1.2.2                labeling_0.4.2           
##  [61] tidyselect_1.1.1          rlang_0.4.11             
##  [63] later_1.3.0               munsell_0.5.0            
##  [65] cellranger_1.1.0          tools_4.1.1              
##  [67] generics_0.1.0            rintrojs_0.3.0           
##  [69] broom_0.7.9               fastmap_1.1.0            
##  [71] evaluate_0.14             stringr_1.4.0            
##  [73] yaml_2.2.1                knitr_1.33               
##  [75] zip_2.2.0                 purrr_0.3.4              
##  [77] nlme_3.1-152              sparseMatrixStats_1.4.2  
##  [79] mime_0.11                 compiler_4.1.1           
##  [81] beeswarm_0.4.0            curl_4.3.2               
##  [83] png_0.1-7                 ggsignif_0.6.2           
##  [85] tibble_3.1.3              statmod_1.4.36           
##  [87] bslib_0.2.5.1             stringi_1.7.3            
##  [89] highr_0.9                 forcats_0.5.1            
##  [91] lattice_0.20-44           bluster_1.2.1            
##  [93] shinyjs_2.0.0             vctrs_0.3.8              
##  [95] pillar_1.6.2              lifecycle_1.0.0          
##  [97] jquerylib_0.1.4           GlobalOptions_0.1.2      
##  [99] BiocNeighbors_1.10.0      cowplot_1.1.1            
## [101] bitops_1.0-7              irlba_2.3.3              
## [103] raster_3.4-13             httpuv_1.6.2             
## [105] R6_2.5.1                  promises_1.2.0.1         
## [107] gridExtra_2.3             rio_0.5.27               
## [109] vipor_0.4.5               codetools_0.2-18         
## [111] colourpicker_1.1.0        assertthat_0.2.1         
## [113] rjson_0.2.20              shinyWidgets_0.6.0       
## [115] withr_2.4.2               GenomeInfoDbData_1.2.6   
## [117] mgcv_1.8-36               hms_1.1.0                
## [119] grid_4.1.1                beachmat_2.8.1           
## [121] tidyr_1.1.3               rmarkdown_2.10           
## [123] DelayedMatrixStats_1.14.2 carData_3.0-4            
## [125] Cairo_1.5-12.2            shiny_1.6.0              
## [127] ggbeeswarm_0.6.0
```

