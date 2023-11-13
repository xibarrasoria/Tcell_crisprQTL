---
title: "<span style='font-size: 16px; color: grey'> crisprQTL proof-of-concept experiment in T cells </span> <br> <span style='font-size: 36px'> Quality control </span>"
author: "Ximena Ibarra-Soria"
date: '21 July, 2022'
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



This is a proof-of-concept experiment to enable CRISPRi of intergenic elements in primary CD4+ T cells, with a single-cell RNA-seq readout (CROP-seq). To assess the performance of the assay we use a library of positive control gRNAs including:

- 35 gene **transcription start sites** (TSS). For these perturbations, we expect silencing of the targeted genes
- 3 **locus control regions** (LCR) of CD2. Again, silencing these loci should lead to downregulation of *CD2* expression.
- 28 **intergenic loci** that have been **linked** to target genes by [Gasperini et al.](https://doi.org/10.1016/j.cell.2018.11.029) We expect decreased expression, but the caveat is that Gasperini et al. performed their screen in K562 cells.
- 14 additional loci that overlap candidate cis regulatory elements, often located within introns. For these, we do not have a clear expectation of the effects, but are annotated based on the nearest gene.

Each gene/locus is targeted by **four different sgRNAs**. 

We also have **35 non-targeting control sgRNAs** as negative controls.

### Experimental design

Primary CD4+ T cells were run through the experimental protocol to induce CRISPRi with the library detailed above. The perturbed cells were then loaded onto a 10X chromium chip for cell capture and cDNA conversion, using all **8 channels as technical replicates**. The remaining cells were frozen down.

A second 10X run was performed using the frozen cells, running four independent chips to process **32 additional technical replicates** to increase cell numbers. The gRNA library from one sample failed sequencing, and will not be considered further.

### Data pre-processing

Sequencing data was processed using `cellranger v4.0`, with a modified reference genome containing an artificial chromosome composed of all the gRNA sequences used, to capture any gRNA transcripts present in the cDNA library.

All samples showed good mapping statistics and adequate barcode-rank plots. There was one sample from the first experiment with very low cell recovery, likely due to a clog in the 10X channel. This sample will not be included in the data analysis.

We start downstream analyses from `cellranger` filtered count matrices.

### Dataset

Here, we will combine and analyse the data from both experiments; we refer to these two experiments as two different batches. We have ~7K +/- 1.5K cells per sample and 30 samples (technical replicates) in total.


```r
## read in raw data
# ignore sample 5 which had a blockage
# sample 9 failed sequencing and is missing
dirs <- c(paste0(data.dir, "pilot_CROPseq_CRISPRi_Tcells_batch1/",
                 "Tcell_CRISPRi_", c(1:4,6:8), "_filtered_feature_bc_matrix"),
          paste0(data.dir, "pilot_CROPseq_CRISPRi_Tcells_batch2/",
                 "Tcell_CRISPRi_", 10:40, ".filtered_feature_bc_matrix"))

sce <- read10xCounts(samples = dirs, 
                     sample.names = paste0("Tcell_CRISPRi_", c(1:4,6:8,10:40)),
                     col.names=TRUE, 
                     BPPARAM = bpp)
  
## separate cDNA counts from gRNA counts
sce <- splitAltExps(sce, f = rowData(sce)$Type, ref="Gene Expression")
# rename the alternative experiment
altExpNames(sce)[altExpNames(sce) == "CRISPR Guide Capture"] <- "CRISPR"
  
## also separate gRNA counts captured from the cDNA library
sce <- splitAltExps(sce, f = grepl("guide", rowData(sce)[,'ID']))
altExpNames(sce)[2] <- "gRNA_in_cDNA"

## use gene names as row.names instead of Ensembl IDs
rownames(sce) <- uniquifyFeatureNames(rowData(sce)[,'ID'], rowData(sce)[,'Symbol'])

## plot
# relevel samples to correctly order them
sce$Sample = factor(sce$Sample, levels = paste0("Tcell_CRISPRi_", c(1:4, 6:8, 10:40)))

# add batch info
sce$batch <- ifelse(sce$Sample %in% paste0("Tcell_CRISPRi_", c(1:4,6:8)), "batch1", "batch2")
# add chip information (groups of 8 samples)
chip <- paste0("chip", 1:5)
chip <- c(rep(chip[1:2], each=7), rep(chip[3:5], each=8))
names(chip) <- unique(sce$Sample)
sce$chip <- chip[sce$Sample]

# n cells per sample
df <- as.data.frame(table(sce$Sample))
df$batch <- ifelse(df$Var1 %in% paste0("Tcell_CRISPRi_", c(1:4,6:8)), "batch1", "batch2")
df$chip <- chip

ggplot(df, aes(Var1, Freq, fill=chip)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = as.character(GetColors(n=5, scheme = "light"))) +
  facet_grid(~batch, scales = "free_x", space = "free_x") +
  xlab("technical replicates") +
  ylab("number of cells recovered") +
  th + theme(axis.text.x = element_blank(),
             axis.ticks.x = element_blank())
```

![](01_quality_control_files/figure-html/readData-1.png)<!-- -->


### Quality control

To remove barcodes associated with poor quality data we check the library size, number of detected genes and fraction of reads from mitochondrial genes in each sample. All three metrics show unimodal distributions and are comparable across batches and chips. 


```r
## compute QC metrics and find outliers
# define mitochondrial genes
mito <- grepl("^MT-", rowData(sce)$Symbol)
# compute metrics
qc_metrics <- perCellQCMetrics(sce, subsets = list(Mt = mito), BPPARAM = bpp)

# add to sce
sce$libSize <- log10(qc_metrics$sum)
sce$n_genes <- log10(qc_metrics$detected)
sce$mt <- qc_metrics$subsets_Mt_percent
# also for gRNA features
sce$libSize_gRNA <- log10(qc_metrics$altexps_CRISPR_sum+1)
sce$n_gRNAs <- qc_metrics$altexps_CRISPR_detected
sce$pct_gRNA_total <- qc_metrics$altexps_CRISPR_percent
    
## define outliers
outliers.libSize <- isOutlier(qc_metrics$sum, 
                              type = "lower", 
                              log = TRUE, 
                              nmads = 3,
                              batch = sce$Sample)
outliers.ngenes <- isOutlier(qc_metrics$detected, 
                             type = "lower", 
                             log = TRUE, 
                             nmads = 3,
                             batch = sce$Sample)
outliers.mt <- isOutlier(qc_metrics$subsets_Mt_percent, 
                         type = "higher",
                         log = FALSE,
                         nmads = 3,
                         batch = sce$Sample)
outliers <- outliers.libSize | outliers.ngenes | outliers.mt
qc_metrics$passQC <- !outliers
sce$passQC <- !outliers

# and thresholds
thresholds <- data.frame(libSize = attr(outliers.libSize, "thresholds")[1,],
                         nGenes = attr(outliers.ngenes, "thresholds")[1,],
                         mt = attr(outliers.mt, "thresholds")[2,])

## plot
plots <- list()

# add sample to the qc dataframe
qc_metrics$sample <- sce$Sample
qc_metrics <- as.data.frame(qc_metrics)

# add coordinates to plot the thresholds
thresholds$start <- 1:nrow(thresholds)-0.4
thresholds$end <- 1:nrow(thresholds)+0.4
    
# colour samples by chip
qc_metrics$colour <- sce$chip

# plot library size
plots[['libSize']] <- ggplot(qc_metrics, aes(sample, sum)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(side="l") +
  geom_violin(aes(colour = colour)) + 
  geom_boxplot(aes(colour = colour), width=0.05) + 
  geom_segment(data = thresholds, 
               aes(x=start, xend=end, y=libSize, yend=libSize),
               colour="grey20", lwd=0.5, lty=2) +
  ggtitle("library size") +
  xlab("") + 
  ylab(expression('log'[10]*' total UMIs')) + 
  th +
  theme(legend.position = "none",
        axis.text.x = element_blank())
# plot number of genes
plots[['nGenes']] <- ggplot(qc_metrics, aes(sample, detected)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(side="l") +
  geom_violin(aes(colour = colour)) + 
  geom_boxplot(aes(colour = colour), width=0.05) + 
  geom_segment(data = thresholds, 
               aes(x=start, xend=end, y=nGenes, yend=nGenes),
               colour="grey20", lwd=0.5, lty=2) +
  ggtitle("total genes per cell") +
  xlab("") + 
  ylab(expression('log'[10]*' genes detected')) +
  th + 
  theme(legend.position = "none",
        axis.text.x = element_blank())
# plot mitochondrial proportion
plots[['mt']] <- ggplot(qc_metrics, aes(sample, subsets_Mt_percent)) + 
  geom_violin(aes(colour = colour)) + 
  geom_boxplot(aes(colour = colour), width=0.05) + 
  geom_segment(data = thresholds, 
               aes(x=start, xend=end, y=mt, yend=mt),
               colour="grey20", lwd=0.5, lty=2) +
  ggtitle("% reads in mitochondrial genes") +
  xlab("") + 
  ylab(expression('% reads')) +
  th + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggarrange(plotlist = plots, ncol=1, nrow=length(plots), align = "v", heights = c(0.275,0.275,0.4))
```

![](01_quality_control_files/figure-html/qc_metrics-1.png)<!-- -->

Based on these distributions, we define thresholds to determine outliers as any barcodes that deviate by more than 3 median absolute deviations from the median of each sample, as indicated below.


```r
thresholds$start <- NULL
thresholds$end <- NULL
thresholds
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["libSize"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["nGenes"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["mt"],"name":[3],"type":["dbl"],"align":["right"]}],"data":[{"1":"2022.297","2":"1094.9442","3":"11.342294","_rn_":"Tcell_CRISPRi_1"},{"1":"2290.710","2":"1308.0186","3":"8.433789","_rn_":"Tcell_CRISPRi_10"},{"1":"2638.967","2":"1383.0450","3":"8.545091","_rn_":"Tcell_CRISPRi_11"},{"1":"2681.244","2":"1396.7176","3":"8.687898","_rn_":"Tcell_CRISPRi_12"},{"1":"1988.392","2":"1121.1051","3":"8.648437","_rn_":"Tcell_CRISPRi_13"},{"1":"2444.164","2":"1276.8105","3":"8.962162","_rn_":"Tcell_CRISPRi_14"},{"1":"2782.909","2":"1398.8462","3":"8.887379","_rn_":"Tcell_CRISPRi_15"},{"1":"2324.670","2":"1269.7217","3":"8.672490","_rn_":"Tcell_CRISPRi_16"},{"1":"1791.367","2":"1041.1036","3":"7.577636","_rn_":"Tcell_CRISPRi_17"},{"1":"1786.801","2":"1034.0001","3":"7.563547","_rn_":"Tcell_CRISPRi_18"},{"1":"1673.168","2":"982.1276","3":"7.624784","_rn_":"Tcell_CRISPRi_19"},{"1":"2173.663","2":"1124.6084","3":"11.731268","_rn_":"Tcell_CRISPRi_2"},{"1":"1865.144","2":"1067.0251","3":"7.881745","_rn_":"Tcell_CRISPRi_20"},{"1":"1688.746","2":"963.7544","3":"7.405660","_rn_":"Tcell_CRISPRi_21"},{"1":"2284.644","2":"1262.9033","3":"7.950886","_rn_":"Tcell_CRISPRi_22"},{"1":"2415.452","2":"1321.9865","3":"8.105254","_rn_":"Tcell_CRISPRi_23"},{"1":"2294.659","2":"1311.8791","3":"8.237508","_rn_":"Tcell_CRISPRi_24"},{"1":"2770.913","2":"1493.7837","3":"8.394258","_rn_":"Tcell_CRISPRi_25"},{"1":"2783.955","2":"1499.9020","3":"8.284758","_rn_":"Tcell_CRISPRi_26"},{"1":"2322.060","2":"1281.8092","3":"8.421382","_rn_":"Tcell_CRISPRi_27"},{"1":"2493.758","2":"1334.2304","3":"8.365254","_rn_":"Tcell_CRISPRi_28"},{"1":"2565.891","2":"1385.9129","3":"8.521719","_rn_":"Tcell_CRISPRi_29"},{"1":"2050.579","2":"1104.0481","3":"11.325982","_rn_":"Tcell_CRISPRi_3"},{"1":"2512.466","2":"1342.8854","3":"8.487700","_rn_":"Tcell_CRISPRi_30"},{"1":"2236.550","2":"1234.6779","3":"8.546640","_rn_":"Tcell_CRISPRi_31"},{"1":"2517.157","2":"1376.4516","3":"8.543261","_rn_":"Tcell_CRISPRi_32"},{"1":"2601.540","2":"1302.7560","3":"8.496854","_rn_":"Tcell_CRISPRi_33"},{"1":"2390.586","2":"1191.8632","3":"8.388476","_rn_":"Tcell_CRISPRi_34"},{"1":"2099.322","2":"1099.3597","3":"8.507051","_rn_":"Tcell_CRISPRi_35"},{"1":"1756.146","2":"921.7005","3":"8.448994","_rn_":"Tcell_CRISPRi_36"},{"1":"1924.587","2":"1010.1203","3":"8.513033","_rn_":"Tcell_CRISPRi_37"},{"1":"2139.746","2":"1091.9896","3":"8.548180","_rn_":"Tcell_CRISPRi_38"},{"1":"2249.594","2":"1114.3858","3":"8.640327","_rn_":"Tcell_CRISPRi_39"},{"1":"2231.225","2":"1135.3083","3":"11.599120","_rn_":"Tcell_CRISPRi_4"},{"1":"2588.720","2":"1271.3284","3":"8.222638","_rn_":"Tcell_CRISPRi_40"},{"1":"2174.841","2":"1125.6007","3":"11.615065","_rn_":"Tcell_CRISPRi_6"},{"1":"2162.024","2":"1125.0741","3":"11.798252","_rn_":"Tcell_CRISPRi_7"},{"1":"2477.015","2":"1257.4834","3":"11.766294","_rn_":"Tcell_CRISPRi_8"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

Barcodes with library sizes or number of detected genes **lower** than these thresholds, or with **higher** mitochondrial fraction are flagged as low-quality. This corresponds to ~4% of the data.


```r
round(prop.table(table(pass_QC = sce$passQC))*100, 2)
```

```
## pass_QC
## FALSE  TRUE 
##  4.07 95.93
```

Finally, there are a few barcodes that show almost no mitochondrial gene expression. These usually represent stripped nuclei instead of cells.


```r
plot(density(sce$mt), 
     xlab = "% reads in mitochondrial genes",
     main = "", bty="l")
```

![](01_quality_control_files/figure-html/mt-1.png)<!-- -->

To identify these potential nuclei, we define outliers with very **low** mitochondrial content. Consistent with these representing nuclei, we observe much lower number of genes detected. 


```r
outliers.nuclei <- isOutlier(qc_metrics$subsets_Mt_percent, 
                               nmads = 3, type = "lower", log = FALSE, 
                               batch = qc_metrics$sample)
qc_metrics$nuclei <- outliers.nuclei
sce$nuclei <- outliers.nuclei

ggplot(qc_metrics, aes(nuclei, detected, colour=nuclei)) +
  geom_violin() +
  geom_boxplot(width=0.1) +
  xlab("is nuclei") +
  ylab("number of genes detected") +
  th + theme(legend.position = "none")
```

![](01_quality_control_files/figure-html/nuclei-1.png)<!-- -->

### Doublet detection

To identify potential doublets we use the `scDblFinder` package, under default settings.


```r
set.seed(9378)
sce <- scDblFinder(sce, samples="Sample", BPPARAM=bpp)
prop.table(table(sce$Sample, sce$scDblFinder.class), 1)*100
```

```
##                   
##                      singlet   doublet
##   Tcell_CRISPRi_1  89.521023 10.478977
##   Tcell_CRISPRi_2  90.241113  9.758887
##   Tcell_CRISPRi_3  89.798432 10.201568
##   Tcell_CRISPRi_4  90.706320  9.293680
##   Tcell_CRISPRi_6  90.246603  9.753397
##   Tcell_CRISPRi_7  89.289816 10.710184
##   Tcell_CRISPRi_8  89.942016 10.057984
##   Tcell_CRISPRi_10 91.085103  8.914897
##   Tcell_CRISPRi_11 92.306417  7.693583
##   Tcell_CRISPRi_12 91.801334  8.198666
##   Tcell_CRISPRi_13 91.110511  8.889489
##   Tcell_CRISPRi_14 91.015070  8.984930
##   Tcell_CRISPRi_15 92.581953  7.418047
##   Tcell_CRISPRi_16 91.978941  8.021059
##   Tcell_CRISPRi_17 92.766104  7.233896
##   Tcell_CRISPRi_18 92.750267  7.249733
##   Tcell_CRISPRi_19 91.784339  8.215661
##   Tcell_CRISPRi_20 90.512312  9.487688
##   Tcell_CRISPRi_21 92.731433  7.268567
##   Tcell_CRISPRi_22 92.097313  7.902687
##   Tcell_CRISPRi_23 92.666552  7.333448
##   Tcell_CRISPRi_24 92.839973  7.160027
##   Tcell_CRISPRi_25 92.189291  7.810709
##   Tcell_CRISPRi_26 91.562952  8.437048
##   Tcell_CRISPRi_27 89.509269 10.490731
##   Tcell_CRISPRi_28 91.540785  8.459215
##   Tcell_CRISPRi_29 89.996362 10.003638
##   Tcell_CRISPRi_30 89.938986 10.061014
##   Tcell_CRISPRi_31 91.315713  8.684287
##   Tcell_CRISPRi_32 91.261014  8.738986
##   Tcell_CRISPRi_33 91.249165  8.750835
##   Tcell_CRISPRi_34 90.249066  9.750934
##   Tcell_CRISPRi_35 90.961343  9.038657
##   Tcell_CRISPRi_36 91.608168  8.391832
##   Tcell_CRISPRi_37 90.184295  9.815705
##   Tcell_CRISPRi_38 90.430562  9.569438
##   Tcell_CRISPRi_39 91.423695  8.576305
##   Tcell_CRISPRi_40 91.666667  8.333333
```

As expected, barcodes labelled as doublets have higher library sizes and many more detected genes.


```r
## add to QC stats
qc_metrics$scDblFinder.class <- sce$scDblFinder.class

## plot libSize and n_genes
plots <- list()
plots[[1]] <- ggplot(qc_metrics[qc_metrics$passQC,], aes(scDblFinder.class, log10(sum), colour=scDblFinder.class)) +
  geom_violin() +
  geom_boxplot(width=0.1) +
  annotation_logticks(side="l") +
  xlab("") +
  ylab("log10 total UMI counts") +
  ggtitle("library size") +
  th + theme(legend.position = "none")
plots[[2]] <- ggplot(qc_metrics[qc_metrics$passQC,], aes(scDblFinder.class, log10(detected), colour=scDblFinder.class)) +
  geom_violin() +
  geom_boxplot(width=0.1) +
  annotation_logticks(side="l") +
  xlab("") +
  ylab("log10 number of genes") +
  ggtitle("genes detected") +
  th + theme(legend.position = "none")

ggarrange(plotlist = plots, ncol=2)
```

![](01_quality_control_files/figure-html/doublet_stats-1.png)<!-- -->


### Good-quality data

In summary, we remove barcodes that:

- Have too low library sizes or number of detected genes.
- Have too high fraction of reads in mitochondrial genes.
- Have very low fraction of reads in mitochondrial genes, since this is indicative of nuclei.
- Are likely to be doublets.

We retain around 250K cells for downstream analyses.


```r
## clean up metadata
colData(sce) <- colData(sce)[c(1:12,15:16)]

## only keep good-quality cells
sce <- sce[,sce$passQC]
## remove nuclei
sce <- sce[,!sce$nuclei]
## remove doublets
sce <- sce[,sce$scDblFinder.class == "singlet"]

## remove non-expressed genes
sce <- sce[rowMeans(counts(sce))>0,]

ncol(sce)
```

```
## [1] 250195
```

We save the good-quality data for downstream analyses.


```r
# save QC metrics
write.table(qc_metrics, paste0(dir, "results/01_QCmetrics_POC_Tcells.tsv"),
            quote = FALSE, sep="\t")

# and save raw data after removing low-quality samples
saveRDS(sce, paste0(dir, "results/01_sce_POC_Tcells.goodQual.RAW.Rds"))
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
##  [5] BiocSingular_1.8.1          BiocParallel_1.26.1        
##  [7] scDblFinder_1.6.0           scran_1.20.1               
##  [9] scater_1.20.1               ggplot2_3.3.5              
## [11] scuttle_1.2.1               DropletUtils_1.12.2        
## [13] SingleCellExperiment_1.14.1 SummarizedExperiment_1.22.0
## [15] Biobase_2.52.0              GenomicRanges_1.44.0       
## [17] GenomeInfoDb_1.28.1         IRanges_2.26.0             
## [19] S4Vectors_0.30.0            BiocGenerics_0.38.0        
## [21] MatrixGenerics_1.4.2        matrixStats_0.60.0         
## 
## loaded via a namespace (and not attached):
##   [1] ggbeeswarm_0.6.0          colorspace_2.0-2         
##   [3] ggsignif_0.6.2            ellipsis_0.3.2           
##   [5] rio_0.5.27                rgdal_1.5-23             
##   [7] bluster_1.2.1             XVector_0.32.0           
##   [9] BiocNeighbors_1.10.0      farver_2.1.0             
##  [11] fansi_0.5.0               codetools_0.2-18         
##  [13] R.methodsS3_1.8.1         sparseMatrixStats_1.4.2  
##  [15] knitr_1.33                jsonlite_1.7.2           
##  [17] broom_0.7.9               cluster_2.1.2            
##  [19] R.oo_1.24.0               HDF5Array_1.20.0         
##  [21] compiler_4.1.1            dqrng_0.3.0              
##  [23] backports_1.2.1           assertthat_0.2.1         
##  [25] Matrix_1.3-4              limma_3.48.3             
##  [27] htmltools_0.5.1.1         tools_4.1.1              
##  [29] rsvd_1.0.5                igraph_1.2.6             
##  [31] gtable_0.3.0              glue_1.4.2               
##  [33] GenomeInfoDbData_1.2.6    dplyr_1.0.7              
##  [35] Rcpp_1.0.7                carData_3.0-4            
##  [37] raster_3.4-13             cellranger_1.1.0         
##  [39] jquerylib_0.1.4           vctrs_0.3.8              
##  [41] rhdf5filters_1.4.0        DelayedMatrixStats_1.14.2
##  [43] xfun_0.25                 stringr_1.4.0            
##  [45] openxlsx_4.2.4            beachmat_2.8.1           
##  [47] lifecycle_1.0.0           irlba_2.3.3              
##  [49] statmod_1.4.36            rstatix_0.7.0            
##  [51] edgeR_3.34.0              zlibbioc_1.38.0          
##  [53] hms_1.1.0                 rhdf5_2.36.0             
##  [55] yaml_2.2.1                curl_4.3.2               
##  [57] gridExtra_2.3             sass_0.4.0               
##  [59] stringi_1.7.3             highr_0.9                
##  [61] ScaledMatrix_1.0.0        checkmate_2.0.0          
##  [63] zip_2.2.0                 rlang_0.4.11             
##  [65] pkgconfig_2.0.3           bitops_1.0-7             
##  [67] evaluate_0.14             lattice_0.20-44          
##  [69] purrr_0.3.4               Rhdf5lib_1.14.2          
##  [71] labeling_0.4.2            cowplot_1.1.1            
##  [73] tidyselect_1.1.1          magrittr_2.0.1           
##  [75] R6_2.5.1                  generics_0.1.0           
##  [77] metapod_1.0.0             DelayedArray_0.18.0      
##  [79] DBI_1.1.1                 pillar_1.6.2             
##  [81] haven_2.4.3               foreign_0.8-81           
##  [83] withr_2.4.2               sp_1.4-5                 
##  [85] abind_1.4-5               RCurl_1.98-1.4           
##  [87] tibble_3.1.3              crayon_1.4.1             
##  [89] car_3.0-11                xgboost_1.4.1.1          
##  [91] utf8_1.2.2                rmarkdown_2.10           
##  [93] viridis_0.6.1             locfit_1.5-9.4           
##  [95] grid_4.1.1                readxl_1.3.1             
##  [97] data.table_1.14.0         forcats_0.5.1            
##  [99] digest_0.6.27             tidyr_1.1.3              
## [101] R.utils_2.10.1            munsell_0.5.0            
## [103] beeswarm_0.4.0            viridisLite_0.4.0        
## [105] vipor_0.4.5               bslib_0.2.5.1
```

