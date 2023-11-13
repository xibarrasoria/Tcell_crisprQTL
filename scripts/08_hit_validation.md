---
title: "<span style='font-size: 16px; color: grey'> crisprQTL proof-of-concept experiment in T cells </span> <br> <span style='font-size: 36px'> Hit validation </span>"
author: "Ximena Ibarra-Soria"
date: '10 November, 2023'
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



We have identified deferentially expressed genes after silencing intergenic elements in CD4 T cells. To validate the results, new experiments are performed for four different targets:


```r
## bulk RNA-seq validation of selected hits
files <- list.files(paste0(dir, "data/validation_experiments/Gene"), pattern = "short.txt")

data <- lapply(files, function(x) { read.table(paste0(dir, "data/validation_experiments/Gene/", x),header=1,row.names=1) })
data <- do.call(cbind, data)
colnames(data) <- sub(".trimmomatic.genome.deduplicated.gene_id.exon.ct.short.txt", "", files)

# metadata
meta <- read.table(paste0(dir, "data/validation_experiments/validation_metadata.tsv"), header = TRUE)
meta <- meta[meta$IRMS_Sample_ID %in% colnames(data),]
meta <- meta[match(colnames(data), meta$IRMS_Sample_ID),]

meta$target <- as.factor(meta$target)
meta$target <- relevel(meta$target, ref = "NT")
levels(meta$target)[-1]
```

```
## [1] "CTSC_ENH1"  "HPCAL1_ENH" "KCNN4_ENH"  "PKM_ENH"
```

```r
## annotate gene names
# take from the annotation used by STAR
# zcat Homo_sapiens.GRCh38.96.gtf.gz | grep -w gene | cut -f 9 | cut -d " " -f2,6 > Homo_sapiens.GRCh38.96.ann
gene_ann <- read.table(paste0(dir, "data/validation_experiments/Homo_sapiens.GRCh38.96.ann"))[,c(1,3)]
colnames(gene_ann) <- c("id", "name")
```

Instead of silencing the elements, this time we use CRISPRn with two gRNAs flanking the putative enhancer to induce a deletion. The perturbed cells are then sequenced in bulk to assess transcriptional changes. Experiments are done using cells from two different donors. Two non-targeting controls are also included to serve as controls.

### Differential expression analysis

We use `edgeR` to test for differences in expression levels between the edited cells and those with NT control gRNAs.


```r
## edgeR object
y <- DGEList(data, samples=meta)
# remove lowly expressed genes
y <- y[filterByExpr(y, group=y$samples$target),]
# estimate size factors
y <- calcNormFactors(y)
# normalised data
data.cpm <- cpm(y, log = TRUE)

# design
design <- model.matrix(~target, y$samples)

# fit model
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust = TRUE)

## test each perturbation against NT controls
# also consider only genes within 1Mb as done in the screen
gene_gRNA_pairs <- readRDS(paste0(dir, "results/04_gene_gRNA_pairs.Rds"))
gRNA_ann <- read.table(paste0(dir, "data/gRNA_library.tsv"), sep="\t", header = TRUE)
gene_gRNA_pairs$target <- gRNA_ann[match(gene_gRNA_pairs$gRNA_id, gRNA_ann$ID),'target']
gene_gRNA_pairs <- unique(gene_gRNA_pairs[,c(1,3)])
gene_gRNA_pairs$id <- gene_ann[match(gene_gRNA_pairs$gene_id, gene_ann$name), 'id']

## edgeR object
y <- DGEList(data, samples=meta)
# remove lowly expressed genes
y <- y[filterByExpr(y, group=y$samples$target),]
# estimate size factors
y <- calcNormFactors(y)
# normalised data
data.cpm <- cpm(y, log = TRUE)


edger <- list()
plots <- list()
for(target in levels(y$samples$target)[-1]){
  # subset 
  y.sub <- y[,y$samples$target %in% c(target, "NT")]
  y.sub$samples$target <- droplevels(y.sub$samples$target)
  y.sub$samples$target <- relevel(y.sub$samples$target, ref="NT")

  # design
  design <- model.matrix(~target, y.sub$samples)

  # fit model
  y.sub <- estimateDisp(y.sub, design)
  fit <- glmQLFit(y.sub, design, robust = TRUE)

  edger[[target]] <- as.data.frame(topTags(glmQLFTest(fit), n=Inf))
  edger[[target]]$target <- target

  # add gene names
  edger[[target]]$gene <- gene_ann[match(row.names(edger[[target]]), gene_ann$id),'name']

  # flag genes in the +-1Mb neighbourhood
  genes <- intersect(gene_gRNA_pairs[gene_gRNA_pairs$target == target,]$id, row.names(y$counts))
  edger[[target]]$cis <- row.names(edger[[target]]) %in% genes
  
  # reorder columns and sort
  edger[[target]] <- edger[[target]][,c(6:8, 1:5)]
  edger[[target]] <- edger[[target]][order(-edger[[target]]$cis, edger[[target]]$FDR),]

  # volcano plot
  plots[[target]] <- ggplot(edger[[target]], aes(logFC, -log10(FDR), colour=FDR < 0.05)) +
    geom_point(size=0.5, alpha=0.5) +
    scale_color_manual(values = c("grey", "indianred")) +
    ggtitle(target) +
    xlab("log fold-change") +
    ylab("-log10 FDR") +
    th + labs(colour="DEG") +
    guides(colour = guide_legend(override.aes = list(size=0.75)))
}
ggarrange(plotlist = plots, ncol=4, nrow=1, legend = "none")
```

![](08_hit_validation_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

### Validation of expected genes

We check whether the gene we identified as DE form the single-cell experiments is recapitulated. Three out of the four perturbations show significant dowregulation of the expected gene.


```r
expected <- rbind(edger[['CTSC_ENH1']][edger[['CTSC_ENH1']]$gene == "CTSC",],
                  edger[['HPCAL1_ENH']][edger[['HPCAL1_ENH']]$gene == "HPCAL1",],
                  edger[['KCNN4_ENH']][edger[['KCNN4_ENH']]$gene == "KCNN4",],
                  edger[['PKM_ENH']][edger[['PKM_ENH']]$gene == "PKM",])
expected[order(expected$FDR),]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["target"],"name":[1],"type":["chr"],"align":["left"]},{"label":["gene"],"name":[2],"type":["chr"],"align":["left"]},{"label":["cis"],"name":[3],"type":["lgl"],"align":["right"]},{"label":["logFC"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["logCPM"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["F"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["PValue"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["FDR"],"name":[8],"type":["dbl"],"align":["right"]}],"data":[{"1":"CTSC_ENH1","2":"CTSC","3":"TRUE","4":"-0.3971774","5":"7.225263","6":"32.440575","7":"0.001643771","8":"0.01888646","_rn_":"ENSG00000109861"},{"1":"PKM_ENH","2":"PKM","3":"TRUE","4":"-0.4589847","5":"11.430490","6":"20.197600","7":"0.004172798","8":"0.03698556","_rn_":"ENSG00000067225"},{"1":"HPCAL1_ENH","2":"HPCAL1","3":"TRUE","4":"-0.4720598","5":"6.337765","6":"20.026932","7":"0.004410062","8":"0.04017257","_rn_":"ENSG00000115756"},{"1":"KCNN4_ENH","2":"KCNN4","3":"TRUE","4":"-0.1376707","5":"4.524061","6":"1.172519","7":"0.317744841","8":"0.49608763","_rn_":"ENSG00000104783"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

Although `KCNN4_ENH` doesn't reach statistical significance, we do observe the expected downregulation but the variability between donors is larger than the perturbation effect. It is likely that with a larger number of donors this effect would be validated.


```r
plots <- list()
for(target in levels(meta$target)[-1]){
  id <- gene_ann[gene_ann$name == meta[meta$target==target,'expected_gene'],'id']
  df <- meta[meta$target %in% c(target, "NT"),]
  df$expr <- as.numeric(data.cpm[id, df$IRMS_Sample_ID])

  plots[[target]] <- ggplot(df, aes(target, expr, colour=donor)) +
    geom_point() +
    scale_color_manual(values = as.character(GetColors(n=2, scheme="bright"))) +
    xlab("") +
    ylab("log2 CPM") +
    ggtitle(target) +
    th
}

ggarrange(plotlist = plots, ncol=4, nrow=1, common.legend = TRUE, legend = "bottom")
```

![](08_hit_validation_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

### Conclusions

Based on these results, we can conclude that the effects discovered in the single-cell crisprQTL screen are indeed the result of inactivating the putative enhancers. The changes are reproducible in two independent donors.


```r
## save
write.table(data.cpm, paste0(dir, "results/08_validation_experiment_CPMexpr.tsv"), quote = FALSE, sep="\t")

saveRDS(edger, paste0(dir, "results/08_validation_experiment_edgeR.Rds"))
for(target in names(edger)){
  write.table(edger[[target]], paste0(dir, "results/08_validation_experiment_edgeR_", target, ".tsv"), quote = FALSE, sep="\t")
}
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
##  [1] inlmisc_0.5.2               RColorBrewer_1.1-2         
##  [3] ggpubr_0.4.0                scales_1.1.1               
##  [5] edgeR_3.34.0                limma_3.48.3               
##  [7] scran_1.20.1                scater_1.20.1              
##  [9] ggplot2_3.3.5               scuttle_1.2.1              
## [11] SingleCellExperiment_1.14.1 SummarizedExperiment_1.22.0
## [13] Biobase_2.52.0              GenomicRanges_1.44.0       
## [15] GenomeInfoDb_1.28.1         IRanges_2.26.0             
## [17] S4Vectors_0.30.0            BiocGenerics_0.38.0        
## [19] MatrixGenerics_1.4.2        matrixStats_0.60.0         
## 
## loaded via a namespace (and not attached):
##  [1] ggbeeswarm_0.6.0          colorspace_2.0-2         
##  [3] ggsignif_0.6.2            ellipsis_0.3.2           
##  [5] rio_0.5.27                rgdal_1.5-23             
##  [7] bluster_1.2.1             XVector_0.32.0           
##  [9] BiocNeighbors_1.10.0      farver_2.1.0             
## [11] fansi_0.5.0               splines_4.1.1            
## [13] codetools_0.2-18          sparseMatrixStats_1.4.2  
## [15] knitr_1.33                jsonlite_1.7.2           
## [17] broom_0.7.9               cluster_2.1.2            
## [19] compiler_4.1.1            dqrng_0.3.0              
## [21] backports_1.2.1           assertthat_0.2.1         
## [23] Matrix_1.3-4              BiocSingular_1.8.1       
## [25] htmltools_0.5.1.1         tools_4.1.1              
## [27] rsvd_1.0.5                igraph_1.2.6             
## [29] gtable_0.3.0              glue_1.4.2               
## [31] GenomeInfoDbData_1.2.6    dplyr_1.0.7              
## [33] Rcpp_1.0.7                carData_3.0-4            
## [35] raster_3.4-13             cellranger_1.1.0         
## [37] jquerylib_0.1.4           vctrs_0.3.8              
## [39] DelayedMatrixStats_1.14.2 xfun_0.25                
## [41] stringr_1.4.0             openxlsx_4.2.4           
## [43] beachmat_2.8.1            lifecycle_1.0.0          
## [45] irlba_2.3.3               statmod_1.4.36           
## [47] rstatix_0.7.0             zlibbioc_1.38.0          
## [49] hms_1.1.0                 yaml_2.2.1               
## [51] curl_4.3.2                gridExtra_2.3            
## [53] sass_0.4.0                stringi_1.7.3            
## [55] highr_0.9                 ScaledMatrix_1.0.0       
## [57] checkmate_2.0.0           zip_2.2.0                
## [59] BiocParallel_1.26.1       rlang_0.4.11             
## [61] pkgconfig_2.0.3           bitops_1.0-7             
## [63] evaluate_0.14             lattice_0.20-44          
## [65] purrr_0.3.4               labeling_0.4.2           
## [67] cowplot_1.1.1             tidyselect_1.1.1         
## [69] magrittr_2.0.1            R6_2.5.1                 
## [71] generics_0.1.0            metapod_1.0.0            
## [73] DelayedArray_0.18.0       DBI_1.1.1                
## [75] pillar_1.6.2              haven_2.4.3              
## [77] foreign_0.8-81            withr_2.4.2              
## [79] sp_1.4-5                  abind_1.4-5              
## [81] RCurl_1.98-1.4            tibble_3.1.3             
## [83] crayon_1.4.1              car_3.0-11               
## [85] utf8_1.2.2                rmarkdown_2.10           
## [87] viridis_0.6.1             locfit_1.5-9.4           
## [89] grid_4.1.1                readxl_1.3.1             
## [91] data.table_1.14.0         forcats_0.5.1            
## [93] digest_0.6.27             tidyr_1.1.3              
## [95] munsell_0.5.0             beeswarm_0.4.0           
## [97] viridisLite_0.4.0         vipor_0.4.5              
## [99] bslib_0.2.5.1
```

