---
title: "<span style='font-size: 16px; color: grey'> crisprQTL proof-of-concept experiment in T cells </span> <br> <span style='font-size: 36px'> Perturbation analysis </span>"
author: "Ximena Ibarra-Soria"
date: '01 August, 2022'
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



### Analysis strategy

Using the normalised data for good-quality cells with at least one gRNA assigned, we test the effects of each perturbation using several differential expression approaches:

- [`limma-voom`](https://doi.org/10.1186/gb-2014-15-2-r29): well-established methodology and extended to work with single-cell data.
- [`MAST`](https://doi.org/10.1186/s13059-015-0844-5): specifically developed for single-cell data.
- [`SCEPTRE`](https://doi.org/10.1186/s13059-015-0844-5): specifically developed for *high MOI* crisprQTL data.
- `wilcoxon rank sum test`: non-parametric test.


```r
sce <- readRDS(paste0(dir, "results/03_sce_POC_Tcells.goodQual.calls.NORM.Rds"))
c(n_cells = ncol(sce))
```

```
## n_cells 
##  152403
```

```r
# genes detected in at least 5% of the cells
genes_to_keep <- rowSums(counts(sce)>0) > ncol(sce)*0.05
genes_detected <- rownames(sce)[genes_to_keep]

## gRNA annotation
gRNA_ann <- read.table(paste0(dir, "data/gRNA_library.tsv"), sep="\t", header = TRUE)
```

The majority of the perturbations included in the library correspond to positive controls, where we have an expected effect:

- 35 gene **transcription start sites** (TSS). For these perturbations, we expect silencing of the targeted genes. Denoted `TSS`.
- 3 **locus control regions** (LCR) of CD2. Again, silencing these loci should lead to downregulation of *CD2* expression. Denoted `LCR`.
- 28 **intergenic loci** that have been **linked** to target genes by Gasperini et al. We expect decreased expression, but the caveat is that Gasperini et al. performed their screen in K562 cells, so we do not expect to recover all effects. Denoted `ENH`.

We detect most of the genes with an expected effect, except for:


```r
## check if all expected genes were detected
not_expr <- setdiff(unique(gRNA_ann$expected_DE_gene), c(genes_detected, NA))
not_expr <- unique(gRNA_ann[gRNA_ann$expected_DE_gene %in% not_expr,c(3:5)])
not_expr[order(not_expr$class, decreasing = TRUE),]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["class"],"name":[1],"type":["chr"],"align":["left"]},{"label":["target"],"name":[2],"type":["chr"],"align":["left"]},{"label":["expected_DE_gene"],"name":[3],"type":["chr"],"align":["left"]}],"data":[{"1":"TSS","2":"TBX1_TSS","3":"TBX1","_rn_":"217"},{"1":"TSS","2":"IL12A_TSS","3":"IL12A","_rn_":"277"},{"1":"ENH","2":"USP6NL_ENH","3":"USP6NL","_rn_":"81"},{"1":"ENH","2":"TUBB2A_ENH","3":"TUBB2A","_rn_":"85"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

We additionally have 23 **intergenic loci** that overlap ENCODE cCREs. For these there is no hypothesis.

The library contains 355 gRNAs, each detected in a median of 282 cells. 


```r
ncells_per_gRNA <- rowSums(counts(altExp(sce, 'gRNA_calls')))

ggplot(as.data.frame(ncells_per_gRNA), aes(1, ncells_per_gRNA)) +
  scale_y_log10() +
  geom_violin() +
  geom_boxplot(width=0.1) +
  annotation_logticks(side="l") +
  xlab("") +
  ylab("number of cells per gRNA") +
  th + theme(axis.text.x = element_blank(),
             axis.ticks.x = element_blank())
```

![](05_perturbation_analysis_files/figure-html/cells_per_guide-1.png)<!-- -->

```r
summary(ncells_per_gRNA)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##     7.0   152.5   282.0   450.3   522.5  5245.0
```

If we aggregate all cells with any of the 4 gRNAs for each target, we have a median of 1,435 cells per target.


```r
df <- data.frame(gRNA = names(ncells_per_gRNA),
                 n_cells = ncells_per_gRNA,
                 target = gRNA_ann[match(names(ncells_per_gRNA), gRNA_ann$ID),]$target)
df <- df %>% group_by(target) %>%
  summarise(ncells_per_target = sum(n_cells))

ggplot(df, aes(1, ncells_per_target)) +
  scale_y_log10() +
  geom_violin() +
  geom_boxplot(width=0.1) +
  annotation_logticks(side="l") +
  xlab("") +
  ylab("number of cells per target") +
  th + theme(axis.text.x = element_blank(),
             axis.ticks.x = element_blank())
```

![](05_perturbation_analysis_files/figure-html/cells_per_target-1.png)<!-- -->

```r
summary(df$ncells_per_target)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##     293    1007    1435    1974    2182   10280
```

For the implementation of the algorithms:

- We focus the analysis on genes within 1Mb (up/downstream) of the target and that are detected in at least 5% of the cells. 
- Cells with a given gRNA are tested against a random sample of five thousand cells containing only non-targeting gRNAs, which represent the *wild-type* transcriptional state.
- For `limma` and `MAST` the model fit and dispersion estimates are performed using all genes (detected in at least 5% of cells), and all cells containing *any* of the four gRNAs from each target. The effect of each gRNA is then tested with the appropriate contrast.
- We include the number of genes detected per cell as a covariate.
- We block on the batch (the two independent experiments performed; see `03_normalisation_and_batch_correction.Rmd` for details).
- To assess p-value calibration, we also test sets of 4 non-targeting gRNAs against any gene within 1Mb of a target (i.e. 4 NT gRNAs create a mock new target for which we expect no differentially expressed genes). 


### Differential expression testing

For details on the implementation of each method, please see:

- `sceptre`: using the original release of the package (`v0.0.1`).
  + `04.1_gene_gRNA_pairs.R`: defining the genes in the vicinity of each target to be tested.
  + `04.2.1_sceptre_data_prep.R`: reformats data in the object types required.
  + `04.2.2_param_file_sceptre.R`: parameters used to run `sceptre`.
- `MAST`: 
  + `04.3_run_MAST.R`: run in parallel for subsets of targets.
- `limma-voom`:
  + `04.4_run_limma.R`
- `wilcoxon rank sum test`: non-parametric test.
  + `04.5_run_wilcoxong.R`: using the implementation from `scran`.

### Results

For all algorithms, we read in the results for all targets, and adjust the p-values for multiple testing (Benjamini & Hochberg) if not done already.


```r
sceptre <- readRDS(paste0(dir, "results/04_sceptre/results/all_results.rds"))
## adjust for multiple testing (per perturbation)
sceptre$FDR <- NA
for(guide in unique(sceptre$gRNA_id)){
  sceptre[sceptre$gRNA_id==guide,]$FDR <- p.adjust(sceptre[sceptre$gRNA_id==guide,]$p_value, 'fdr')
}
  
## add true positives
sceptre$target <- gRNA_ann[match(sceptre$gRNA_id, gRNA_ann$ID),]$target
sceptre$expected_gene <- gRNA_ann[match(sceptre$target, gRNA_ann$target),]$expected_DE_gene
sceptre$class <- gRNA_ann[match(sceptre$target, gRNA_ann$target),]$class
sceptre$sig <- sceptre$FDR < 0.05

# separate NT data
sceptre.nt <- sceptre[sceptre$target == "NT",]
sceptre <- sceptre[sceptre$target != "NT",]

# reorder columns
sceptre <- sceptre[,c(1,2,12:14,3,11,15, 4:10)]
sceptre.nt <- sceptre.nt[,c(1,2,12:14,3,11,15, 4:10)]
```


```r
files <- list.files(paste0(dir, "results/04_MAST/"))
mast <- list()
for(file in grep("neighbourhood_seed204", files, value=TRUE)){
  mast[[file]] <- readRDS(paste0(dir, "results/04_MAST/", file))
  mast[[file]] <- do.call(rbind, mast[[file]])
}
mast <- as.data.frame(do.call(rbind, mast))

# rename columns
colnames(mast)[1] <- "gene_id"
colnames(mast)[2] <- "p_value"
colnames(mast)[3] <- "logFC"
colnames(mast)[4] <- "logFC.hi"
colnames(mast)[5] <- "logFC.lo"
colnames(mast)[8] <- "gRNA_id"

## add true positives
mast$expected_gene <- gRNA_ann[match(mast$target, gRNA_ann$target),]$expected_DE_gene
mast$class <- gRNA_ann[match(mast$target, gRNA_ann$target),]$class
mast$sig <- mast$FDR <0.05

# separate NT data
mast.nt <- mast[grep("^NT_", mast$target),]
mast <- mast[grep("^NT_", mast$target, invert=TRUE),]

# reorder columns
mast <- mast[,c(1,8,7,9:10,2,6,11, 3:5)]
mast.nt <- mast.nt[,c(1,8,7,9:10,2,6,11, 3:5)]
```


```r
limma <- readRDS(paste0(dir, "results/04.4_results_limma_neighbourhood.Rds"))
limma <- do.call(rbind, limma)

# rename columns
colnames(limma)[2] <- "gene_id"
colnames(limma)[7] <- "p_value"
colnames(limma)[8] <- "FDR"
colnames(limma)[11] <- "gRNA_id"

## add true positives
limma$expected_gene <- gRNA_ann[match(limma$target, gRNA_ann$target),]$expected_DE_gene
limma$class <- gRNA_ann[match(limma$target, gRNA_ann$target),]$class
limma$sig <- limma$FDR < 0.05

# separate NT data
limma.nt <- limma[grep("^NT_", limma$target),]
limma <- limma[grep("^NT_", limma$target, invert=TRUE),]

# reorder columns
limma <- limma[,c(2,11,10,12:13,7:8,14, 1,4:6,9)]
limma.nt <- limma.nt[,c(2,11,10,12:13,7:8,14, 1,4:6,9)]
```


```r
wilcox <- readRDS(paste0(dir, "results/04.5_results_wilcoxon_neighbourhood.Rds"))
for(i in 1:length(wilcox)){
  wilcox[[i]]$gene_id <- row.names(wilcox[[i]])
}
wilcox <- as.data.frame(do.call(rbind, wilcox))

# rename columns
colnames(wilcox)[2] <- "p_value"
colnames(wilcox)[4] <- "gRNA_id"

## add true positives
wilcox$target <- gRNA_ann[match(wilcox$gRNA, gRNA_ann$ID),]$target
wilcox$expected_gene <- gRNA_ann[match(wilcox$target, gRNA_ann$target),]$expected_DE_gene
wilcox$class <- gRNA_ann[match(wilcox$target, gRNA_ann$target),]$class
wilcox$sig <- wilcox$FDR < 0.05

# separate NT data
wilcox.nt <- wilcox[wilcox$target == "NT",]
wilcox <- wilcox[wilcox$target != "NT",]

# reorder columns
wilcox <- wilcox[,c(5,4,6:8,2:3,9, 1)]
wilcox.nt <- wilcox.nt[,c(5,4,6:8,2:3,9, 1)]
```


#### Detection of true positives {.tabset}

First, we look at whether algorithms are able to identify the expected downregulation of the genes associated with `TSS`, `LCR` and `ENH` perturbations.

##### sceptre

`sceptre` identifies significant differential downregulation for all `LCR` and all but two `TSS` targets, with at least one gRNA. For `ENH` perturbations, 16 (out of 26) also lead to a significant result for at least one gRNA per target.


```r
## assess TP recovery
sceptre.exp <- sceptre[!is.na(sceptre$expected_gene),] # separate targets with expected effect
sceptre.tp <- sceptre.exp[sceptre.exp$gene_id == sceptre.exp$expected_gene,]

## visualise as heatmap
pvals <- matrix(-log10(sceptre.tp$FDR), ncol=4, nrow=length(unique(sceptre.tp$target)), byrow = TRUE)
colnames(pvals) <- paste0("gRNA", 1:4)
row.names(pvals) <- unique(sceptre.tp$target)

# annotate by perturbation type
cols <- as.character(GetColors(n=5, scheme = "muted")[c(1,3,5)])
names(cols) <- unique(sceptre.tp$class)
ha  <- rowAnnotation(class = sceptre.tp[seq(1,nrow(sceptre.tp),4),]$class,
                     n_DE = anno_barplot(sapply(seq(1,nrow(sceptre.tp),4), 
                                                   function(x) sum(sceptre.tp[x:(x+3),]$sig & sceptre.tp[x:(x+3),]$z_value < 0))),
                     col = list(class = cols))
Heatmap(pvals,
        col = circlize::colorRamp2(breaks = c(0,-log10(0.05),3,5,10,15),
                                   colors = c("white", brewer.pal(n=9, "BuPu")[3:7])),
        cluster_columns = FALSE, 
        right_annotation = ha,
        split = sceptre.tp[seq(1,nrow(sceptre.tp),4),]$class,
        row_names_gp = gpar(fontsize = 8),
        heatmap_legend_param = list(title = "-log10(FDR)"),
        column_title = "sceptre")
```

![](05_perturbation_analysis_files/figure-html/sceptre_tp-1.png)<!-- -->

The vast majority of perturbation effects are supported by at least 2 independent gRNAs, and over 60% are significant with three or four gRNAs.


```r
prop.table(table(n_gRNAs=table(sceptre.tp[sceptre.tp$sig & sceptre.tp$z_value < 0, c('target')])))*100
```

```
## n_gRNAs
##  1  2  3  4 
## 20 18 28 34
```

As expected, almost all significant results correspond to gene downregulation, compared to *wild-type* cells. However, a small number of gRNAs induce significant upregulation instead, almost all from `ENH` targets. Generally, a single gRNA shows this effect, but we observe consistent upregulation of *TMSB4X* (with three independent enhancers and 5 different gRNAs) suggesting a possible difference in regulation between K562 and T cells.


```r
sceptre.tp[sceptre.tp$sig & sceptre.tp$z_value > 0,]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["gene_id"],"name":[1],"type":["chr"],"align":["left"]},{"label":["gRNA_id"],"name":[2],"type":["chr"],"align":["left"]},{"label":["target"],"name":[3],"type":["chr"],"align":["left"]},{"label":["expected_gene"],"name":[4],"type":["chr"],"align":["left"]},{"label":["class"],"name":[5],"type":["chr"],"align":["left"]},{"label":["p_value"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["FDR"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["sig"],"name":[8],"type":["lgl"],"align":["right"]},{"label":["skew_t_fit_success"],"name":[9],"type":["lgl"],"align":["right"]},{"label":["xi"],"name":[10],"type":["dbl"],"align":["right"]},{"label":["omega"],"name":[11],"type":["dbl"],"align":["right"]},{"label":["alpha"],"name":[12],"type":["dbl"],"align":["right"]},{"label":["nu"],"name":[13],"type":["dbl"],"align":["right"]},{"label":["z_value"],"name":[14],"type":["dbl"],"align":["right"]},{"label":["n_successful_resamples"],"name":[15],"type":["int"],"align":["right"]}],"data":[{"1":"YPEL3","2":"gRNA_2","3":"YPEL3_ENH","4":"YPEL3","5":"ENH","6":"7.277490e-04","7":"5.943283e-03","8":"TRUE","9":"TRUE","10":"-0.5614570","11":"1.2829109","12":"0.62526052","13":"3.882327e+02","14":"3.978056","15":"500","_rn_":"78"},{"1":"BTG2","2":"gRNA_6","3":"BTG2_ENH","4":"BTG2","5":"ENH","6":"1.439288e-04","7":"1.727146e-03","8":"TRUE","9":"TRUE","10":"-1.3552872","11":"1.8999448","12":"1.71853796","13":"5.196821e+06","14":"5.867355","15":"500","_rn_":"217"},{"1":"BTG2","2":"gRNA_7","3":"BTG2_ENH","4":"BTG2","5":"ENH","6":"1.590611e-08","7":"1.908733e-07","8":"TRUE","9":"TRUE","10":"-1.0630217","11":"1.8006192","12":"1.29647608","13":"8.340444e+01","14":"10.209088","15":"500","_rn_":"229"},{"1":"TMSB4X","2":"gRNA_17","3":"TMSB4X_ENH1","4":"TMSB4X","5":"ENH","6":"9.093171e-12","7":"5.455902e-11","8":"TRUE","9":"TRUE","10":"0.6111629","11":"1.1418317","12":"-0.92380583","13":"7.993131e+06","14":"7.176340","15":"500","_rn_":"326"},{"1":"TMSB4X","2":"gRNA_19","3":"TMSB4X_ENH1","4":"TMSB4X","5":"ENH","6":"8.357975e-04","7":"4.428816e-03","8":"TRUE","9":"TRUE","10":"0.2615502","11":"0.9247799","12":"-0.45298928","13":"1.136306e+01","14":"3.993284","15":"500","_rn_":"338"},{"1":"TMSB4X","2":"gRNA_20","3":"TMSB4X_ENH1","4":"TMSB4X","5":"ENH","6":"7.304288e-04","7":"4.382573e-03","8":"TRUE","9":"TRUE","10":"0.2419068","11":"1.0196237","12":"-0.37488784","13":"9.567716e+04","14":"3.314377","15":"500","_rn_":"344"},{"1":"CTSC","2":"gRNA_120","3":"CTSC_ENH1","4":"CTSC","5":"ENH","6":"2.220446e-16","7":"2.220446e-16","8":"TRUE","9":"TRUE","10":"5.2904418","11":"1.2059163","12":"0.05968157","13":"3.133879e+01","14":"24.499790","15":"500","_rn_":"21504"},{"1":"IER3","2":"gRNA_126","3":"IER3_ENH","4":"IER3","5":"ENH","6":"8.426961e-04","7":"4.297750e-02","8":"TRUE","9":"TRUE","10":"-1.8351381","11":"2.4729590","12":"2.31363237","13":"1.184801e+01","14":"9.127130","15":"500","_rn_":"21580"},{"1":"PHF19","2":"gRNA_134","3":"PHF19_ENH","4":"PHF19","5":"ENH","6":"1.223878e-05","7":"1.101490e-04","8":"TRUE","9":"TRUE","10":"-0.6847198","11":"1.4874317","12":"0.59084660","13":"4.144043e+05","14":"5.978550","15":"500","_rn_":"21770"},{"1":"BTG1","2":"gRNA_159","3":"BTG1_ENH","4":"BTG1","5":"ENH","6":"6.082714e-04","7":"2.433086e-03","8":"TRUE","9":"TRUE","10":"-0.8308141","11":"1.4093155","12":"0.76016487","13":"9.233569e+05","14":"4.256226","15":"500","_rn_":"31529"},{"1":"TMSB4X","2":"gRNA_165","3":"TMSB4X_ENH2","4":"TMSB4X","5":"ENH","6":"1.481661e-02","7":"2.963322e-02","8":"TRUE","9":"TRUE","10":"0.2719050","11":"1.0065552","12":"-0.43275100","13":"1.112062e+02","14":"2.369224","15":"500","_rn_":"31566"},{"1":"TMSB4X","2":"gRNA_171","3":"TMSB4X_ENH3","4":"TMSB4X","5":"ENH","6":"7.564077e-04","7":"4.538446e-03","8":"TRUE","9":"TRUE","10":"0.6545680","11":"1.0161264","12":"-0.67392033","13":"2.148159e+01","14":"3.584905","15":"500","_rn_":"31602"},{"1":"CD2","2":"gRNA_209","3":"CD2_LCR2","4":"CD2","5":"LCR","6":"1.048339e-02","7":"3.145017e-02","8":"TRUE","9":"TRUE","10":"-0.4039381","11":"1.1791033","12":"0.41753097","13":"3.479273e+01","14":"3.038853","15":"500","_rn_":"32135"},{"1":"GATA3","2":"gRNA_222","3":"GATA3_TSS","4":"GATA3","5":"TSS","6":"2.809802e-03","7":"8.429407e-03","8":"TRUE","9":"TRUE","10":"-1.2357384","11":"2.2243819","12":"1.18273863","13":"3.091898e+03","14":"5.488201","15":"500","_rn_":"32275"},{"1":"ETS1","2":"gRNA_248","3":"ETS1_TSS","4":"ETS1","5":"TSS","6":"3.738004e-05","7":"3.738004e-05","8":"TRUE","9":"TRUE","10":"0.7058531","11":"1.0917334","12":"-0.74205092","13":"2.769787e+06","14":"3.995251","15":"500","_rn_":"41243"},{"1":"PTS","2":"gRNA_271","3":"PTS_TSS","4":"PTS","5":"TSS","6":"4.094406e-07","7":"1.364802e-06","8":"TRUE","9":"TRUE","10":"2.0216966","11":"0.9874595","12":"0.53879933","13":"1.847931e+01","14":"9.534380","15":"500","_rn_":"41682"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

##### MAST

`MAST` also identifies significant downregulation of all `LCR` targets, and all but one `TSS` perturbations. 17 out of 26 `ENH` targets have at least one significant gRNA result.


```r
## assess TP recovery
mast.exp <- mast[!is.na(mast$expected_gene),] # separate targets with expected effect
mast.tp <- mast.exp[mast.exp$gene_id == mast.exp$expected_gene,]

## visualise as heatmap
pvals <- matrix(-log10(mast.tp$FDR), ncol=4, nrow=length(unique(mast.tp$target)), byrow = TRUE)
colnames(pvals) <- paste0("gRNA", 1:4)
row.names(pvals) <- unique(mast.tp$target)
# some p-values are 0; assign smallest non-zero value to avoid Inf
pvals[pvals==Inf] <- -log10(min(mast.tp[mast.tp$FDR > 0,]$FDR))

# annotate by perturbation type
cols <- as.character(GetColors(n=5, scheme = "muted")[c(1,3,5)])
names(cols) <- unique(mast.tp$class)
# and count number of significnat gRNAs
# consider that if all the perturbed cells have 0 count a fold-change is not estimated, but it should be considered as a negative (downregulation) effect
ha  <- rowAnnotation(class = mast.tp[seq(1,nrow(mast.tp),4),]$class,
                     n_DE = anno_barplot(sapply(seq(1,nrow(mast.tp),4), 
                                                   function(x) sum(mast.tp[x:(x+3),]$sig & 
                                                                     (mast.tp[x:(x+3),]$logFC < 0 | is.na(mast.tp[x:(x+3),]$logFC))))),
                     col = list(class = cols))
Heatmap(pvals,
        col = circlize::colorRamp2(breaks = c(0,-log10(0.05),30,50,100,150),
                                   colors = c("white", brewer.pal(n=9, "BuPu")[3:7])),
        cluster_columns = FALSE, 
        right_annotation = ha,
        split = mast.tp[seq(1,nrow(mast.tp),4),]$class,
        row_names_gp = gpar(fontsize = 8),
        heatmap_legend_param = list(title = "-log10(FDR)"),
        column_title = "MAST")
```

![](05_perturbation_analysis_files/figure-html/mast_tp-1.png)<!-- -->

Again, the vast majority of perturbation effects are supported by at least 2 independent gRNAs, and in this case over 65% are significant with three or four gRNAs. A larger proportion of targets are supported by higher number of gRNAs compared to `sceptre`, suggesting higher sensitivity.


```r
round(prop.table(table(n_gRNAs=table(mast.tp[mast.tp$sig & mast.tp$logFC < 0, c('target')])))*100, 2)
```

```
## n_gRNAs
##     1     2     3     4 
## 19.23 15.38 25.00 40.38
```

`MAST` also reports a few cases of gene upregulation, and all were also identified by `sceptre`. 


```r
tmp <- mast.tp[mast.tp$sig,]
tmp[which(tmp$logFC > 0),]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["gene_id"],"name":[1],"type":["chr"],"align":["left"]},{"label":["gRNA_id"],"name":[2],"type":["chr"],"align":["left"]},{"label":["target"],"name":[3],"type":["chr"],"align":["left"]},{"label":["expected_gene"],"name":[4],"type":["chr"],"align":["left"]},{"label":["class"],"name":[5],"type":["chr"],"align":["left"]},{"label":["p_value"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["FDR"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["sig"],"name":[8],"type":["lgl"],"align":["right"]},{"label":["logFC"],"name":[9],"type":["dbl"],"align":["right"]},{"label":["logFC.hi"],"name":[10],"type":["dbl"],"align":["right"]},{"label":["logFC.lo"],"name":[11],"type":["dbl"],"align":["right"]}],"data":[{"1":"BTG2","2":"gRNA_6","3":"BTG2_ENH","4":"BTG2","5":"ENH","6":"5.707114e-04","7":"6.848537e-03","8":"TRUE","9":"0.21550031","10":"0.3470031","11":"0.08399756","_rn_":"211"},{"1":"BTG2","2":"gRNA_7","3":"BTG2_ENH","4":"BTG2","5":"ENH","6":"2.636554e-06","7":"3.163865e-05","8":"TRUE","9":"0.25528053","10":"0.3710722","11":"0.13948888","_rn_":"223"},{"1":"TMSB4X","2":"gRNA_17","3":"TMSB4X_ENH1","4":"TMSB4X","5":"ENH","6":"8.860866e-09","7":"5.316520e-08","8":"TRUE","9":"0.16954518","10":"0.2280087","11":"0.11108161","_rn_":"329"},{"1":"TMSB4X","2":"gRNA_19","3":"TMSB4X_ENH1","4":"TMSB4X","5":"ENH","6":"2.468397e-05","7":"1.481038e-04","8":"TRUE","9":"0.12898358","10":"0.1870768","11":"0.07089035","_rn_":"341"},{"1":"CTSC","2":"gRNA_120","3":"CTSC_ENH1","4":"CTSC","5":"ENH","6":"1.222312e-147","7":"1.222312e-147","8":"TRUE","9":"0.51223203","10":"0.5651701","11":"0.45929398","_rn_":"2504"},{"1":"IER3","2":"gRNA_126","3":"IER3_ENH","4":"IER3","5":"ENH","6":"9.208517e-06","7":"4.696343e-04","8":"TRUE","9":"0.28049582","10":"0.4384122","11":"0.12257939","_rn_":"2585"},{"1":"PHF19","2":"gRNA_134","3":"PHF19_ENH","4":"PHF19","5":"ENH","6":"2.392646e-05","7":"1.076691e-04","8":"TRUE","9":"0.19350272","10":"0.2966677","11":"0.09033773","_rn_":"2770"},{"1":"BTG1","2":"gRNA_159","3":"BTG1_ENH","4":"BTG1","5":"ENH","6":"4.313190e-03","7":"1.725276e-02","8":"TRUE","9":"0.35416864","10":"0.5661038","11":"0.14223343","_rn_":"3030"},{"1":"ETS1","2":"gRNA_248","3":"ETS1_TSS","4":"ETS1","5":"TSS","6":"3.056358e-03","7":"3.056358e-03","8":"TRUE","9":"0.18648830","10":"0.2943088","11":"0.07866777","_rn_":"4199"},{"1":"PTS","2":"gRNA_271","3":"PTS_TSS","4":"PTS","5":"TSS","6":"3.771283e-11","7":"9.428208e-11","8":"TRUE","9":"0.08434343","10":"0.1120872","11":"0.05659965","_rn_":"4635"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>


##### limma

`limma-voom` also identifies significant differential downregulation for all `LCR` and all but two `TSS` targets, with at least one gRNA. For `ENH` perturbations, 16 (out of 26) lead to a significant result for at least one gRNA per target.


```r
## assess TP recovery
limma.exp <- limma[!is.na(limma$expected_gene),] # separate targets with expected effect
limma.tp <- limma.exp[limma.exp$gene_id == limma.exp$expected_gene,]


## visualise as heatmap
pvals <- matrix(-log10(limma.tp$FDR), ncol=4, nrow=length(unique(limma.tp$target)), byrow = TRUE)
colnames(pvals) <- paste0("gRNA", 1:4)
row.names(pvals) <- unique(limma.tp$target)
# some p-values are 0; assign smallest non-zero value to avoid Inf
pvals[pvals==Inf] <- -log10(min(limma.tp[limma.tp$FDR > 0,]$FDR))

# annotate by perturbation type
cols <- as.character(GetColors(n=5, scheme = "muted")[c(1,3,5)])
names(cols) <- unique(limma.tp$class)
ha  <- rowAnnotation(class = limma.tp[seq(1,nrow(limma.tp),4),]$class,
                     n_DE = anno_barplot(sapply(seq(1,nrow(limma.tp),4), 
                                                   function(x) sum(limma.tp[x:(x+3),]$sig & limma.tp[x:(x+3),]$logFC < 0))),
                     col = list(class = cols))
Heatmap(pvals,
        col = circlize::colorRamp2(breaks = c(0,-log10(0.05),30,50,100,150),
                                   colors = c("white", brewer.pal(n=9, "BuPu")[3:7])),
        cluster_columns = FALSE, 
        right_annotation = ha,
        split = limma.tp[seq(1,nrow(limma.tp),4),]$class,
        row_names_gp = gpar(fontsize = 8),
        heatmap_legend_param = list(title = "-log10(FDR)"),
        column_title = "limma-voom")
```

![](05_perturbation_analysis_files/figure-html/limma_tp-1.png)<!-- -->

The vast majority of perturbation effects are supported by at least 2 independent gRNAs, and over 60% are significant with three or four gRNAs.


```r
round(prop.table(table(n_gRNAs=table(limma.tp[limma.tp$sig & limma.tp$logFC < 0, c('target')])))*100, 2)
```

```
## n_gRNAs
##     1     2     3     4 
## 16.36 20.00 21.82 41.82
```

Again, almost all significant results correspond to gene downregulation, but a few genes instead show higher expression. All of these *hits* were also identified by `sceptre` and `MAST`. 


```r
limma.tp[limma.tp$sig & limma.tp$logFC > 0,]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["gene_id"],"name":[1],"type":["chr"],"align":["left"]},{"label":["gRNA_id"],"name":[2],"type":["chr"],"align":["left"]},{"label":["target"],"name":[3],"type":["chr"],"align":["left"]},{"label":["expected_gene"],"name":[4],"type":["chr"],"align":["left"]},{"label":["class"],"name":[5],"type":["chr"],"align":["left"]},{"label":["p_value"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["FDR"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["sig"],"name":[8],"type":["lgl"],"align":["right"]},{"label":["ID"],"name":[9],"type":["chr"],"align":["left"]},{"label":["logFC"],"name":[10],"type":["dbl"],"align":["right"]},{"label":["AveExpr"],"name":[11],"type":["dbl"],"align":["right"]},{"label":["t"],"name":[12],"type":["dbl"],"align":["right"]},{"label":["B"],"name":[13],"type":["dbl"],"align":["right"]}],"data":[{"1":"BTG2","2":"gRNA_6","3":"BTG2_ENH","4":"BTG2","5":"ENH","6":"3.309843e-03","7":"3.971811e-02","8":"TRUE","9":"ENSG00000159388","10":"0.23339195","11":"6.827293","12":"2.938512","13":"-2.2047022","_rn_":"BTG2_ENH.gRNA_6.BTG2"},{"1":"BTG2","2":"gRNA_7","3":"BTG2_ENH","4":"BTG2","5":"ENH","6":"2.701992e-05","7":"3.242391e-04","8":"TRUE","9":"ENSG00000159388","10":"0.29387236","11":"6.827293","12":"4.200323","13":"2.2315742","_rn_":"BTG2_ENH.gRNA_7.BTG2"},{"1":"TMSB4X","2":"gRNA_17","3":"TMSB4X_ENH1","4":"TMSB4X","5":"ENH","6":"2.334388e-08","7":"1.400633e-07","8":"TRUE","9":"ENSG00000205542","10":"0.15897651","11":"13.718622","12":"5.592044","13":"8.6010247","_rn_":"TMSB4X_ENH1.gRNA_17.TMSB4X"},{"1":"TMSB4X","2":"gRNA_19","3":"TMSB4X_ENH1","4":"TMSB4X","5":"ENH","6":"2.326662e-05","7":"1.395997e-04","8":"TRUE","9":"ENSG00000205542","10":"0.12679750","11":"13.718622","12":"4.234034","13":"2.1270609","_rn_":"TMSB4X_ENH1.gRNA_19.TMSB4X"},{"1":"TMSB4X","2":"gRNA_20","3":"TMSB4X_ENH1","4":"TMSB4X","5":"ENH","6":"7.533209e-03","7":"4.519925e-02","8":"TRUE","9":"ENSG00000205542","10":"0.05694787","11":"13.718622","12":"2.673133","13":"-3.5422361","_rn_":"TMSB4X_ENH1.gRNA_20.TMSB4X"},{"1":"CTSC","2":"gRNA_120","3":"CTSC_ENH1","4":"CTSC","5":"ENH","6":"2.587120e-158","7":"2.587120e-158","8":"TRUE","9":"ENSG00000109861","10":"0.58736241","11":"9.233058","12":"27.371519","13":"350.5787926","_rn_":"CTSC_ENH1.gRNA_120"},{"1":"IER3","2":"gRNA_126","3":"IER3_ENH","4":"IER3","5":"ENH","6":"2.882939e-04","7":"1.470299e-02","8":"TRUE","9":"ENSG00000137331","10":"0.37872028","11":"5.842344","12":"3.627864","13":"0.1876787","_rn_":"IER3_ENH.gRNA_126.IER3"},{"1":"PHF19","2":"gRNA_134","3":"PHF19_ENH","4":"PHF19","5":"ENH","6":"2.555124e-03","7":"5.749030e-03","8":"TRUE","9":"ENSG00000119403","10":"0.16628009","11":"7.001105","12":"3.017997","13":"-2.1273973","_rn_":"PHF19_ENH.gRNA_134.PHF19"},{"1":"BTG1","2":"gRNA_159","3":"BTG1_ENH","4":"BTG1","5":"ENH","6":"1.691491e-03","7":"6.765962e-03","8":"TRUE","9":"ENSG00000133639","10":"0.44085781","11":"8.498388","12":"3.141236","13":"-1.3199596","_rn_":"BTG1_ENH.gRNA_159.BTG1"},{"1":"ETS1","2":"gRNA_248","3":"ETS1_TSS","4":"ETS1","5":"TSS","6":"2.359702e-03","7":"2.359702e-03","8":"TRUE","9":"ENSG00000134954","10":"0.21925051","11":"8.402600","12":"3.041860","13":"-1.9002481","_rn_":"ETS1_TSS.gRNA_248.ETS1"},{"1":"PTS","2":"gRNA_271","3":"PTS_TSS","4":"PTS","5":"TSS","6":"9.856758e-56","7":"9.856758e-55","8":"TRUE","9":"ENSG00000150787","10":"0.29791674","11":"5.721802","12":"15.810330","13":"114.9449145","_rn_":"PTS_TSS.gRNA_271.PTS"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

##### wilcoxon

Results are also very similar when using a `wilcoxon rank sum` test instead. Significant downregulation for at least one gRNA is detected for all `LCR` and all but four `TSS` targets. In this case, 19 out of 26 `ENH` targets have a significant result.


```r
## assess TP recovery
wilcox.exp <- wilcox[!is.na(wilcox$expected_gene),] # separate targets with expected effect
wilcox.tp <- wilcox.exp[wilcox.exp$gene_id == wilcox.exp$expected_gene,]

## visualise as heatmap
pvals <- matrix(-log10(wilcox.tp$FDR), ncol=4, nrow=length(unique(wilcox.tp$target)), byrow = TRUE)
colnames(pvals) <- paste0("gRNA", 1:4)
row.names(pvals) <- unique(wilcox.tp$target)

# annotate by perturbation type
cols <- as.character(GetColors(n=5, scheme = "muted")[c(1,3,5)])
names(cols) <- unique(wilcox.tp$class)
ha  <- rowAnnotation(class = wilcox.tp[seq(1,nrow(wilcox.tp),4),]$class,
                     n_DE = anno_barplot(sapply(seq(1,nrow(wilcox.tp),4), 
                                                   function(x) sum(wilcox.tp[x:(x+3),]$sig & wilcox.tp[x:(x+3),]$AUC < 0.5))),
                     col = list(class = cols))
Heatmap(pvals,
        col = circlize::colorRamp2(breaks = c(0,-log10(0.05),30,50,100,150),
                                   colors = c("white", brewer.pal(n=9, "BuPu")[3:7])),
        cluster_columns = FALSE, 
        right_annotation = ha,
        split = wilcox.tp[seq(1,nrow(wilcox.tp),4),]$class,
        row_names_gp = gpar(fontsize = 8),
        heatmap_legend_param = list(title = "-log10(FDR)"),
        column_title = "wilcoxon rank sum test")
```

![](05_perturbation_analysis_files/figure-html/wilcox_tp-1.png)<!-- -->

Most targets are supported by several gRNAs, although the cases with three and four gRNAs is a little lower. 


```r
round(prop.table(table(n_gRNAs=table(wilcox.tp[wilcox.tp$sig & wilcox.tp$AUC < 0.5, c('target')])))*100, 2)
```

```
## n_gRNAs
##     1     2     3     4 
## 21.57 23.53 23.53 31.37
```

Again, for the few cases where significant upregulation is detected, most are also reported by the other methods, but there are a few additional *hits*. 


```r
wilcox.tp[wilcox.tp$sig & wilcox.tp$AUC > 0.5,]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["gene_id"],"name":[1],"type":["chr"],"align":["left"]},{"label":["gRNA_id"],"name":[2],"type":["chr"],"align":["left"]},{"label":["target"],"name":[3],"type":["chr"],"align":["left"]},{"label":["expected_gene"],"name":[4],"type":["chr"],"align":["left"]},{"label":["class"],"name":[5],"type":["chr"],"align":["left"]},{"label":["p_value"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["FDR"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["sig"],"name":[8],"type":["lgl"],"align":["right"]},{"label":["AUC"],"name":[9],"type":["dbl"],"align":["right"]}],"data":[{"1":"YPEL3","2":"gRNA_3","3":"YPEL3_ENH","4":"YPEL3","5":"ENH","6":"2.119084e-04","7":"2.595878e-03","8":"TRUE","9":"0.5487925","_rn_":"YPEL3.2"},{"1":"BTG2","2":"gRNA_6","3":"BTG2_ENH","4":"BTG2","5":"ENH","6":"1.739364e-03","7":"1.237603e-02","8":"TRUE","9":"0.5666159","_rn_":"BTG2.1"},{"1":"BTG2","2":"gRNA_7","3":"BTG2_ENH","4":"BTG2","5":"ENH","6":"2.157472e-07","7":"2.588966e-06","8":"TRUE","9":"0.6032796","_rn_":"BTG2.2"},{"1":"PTPRC","2":"gRNA_16","3":"PTPRC_ENH","4":"PTPRC","5":"ENH","6":"8.316206e-03","7":"2.079051e-02","8":"TRUE","9":"0.5431736","_rn_":"PTPRC.3"},{"1":"TMSB4X","2":"gRNA_17","3":"TMSB4X_ENH1","4":"TMSB4X","5":"ENH","6":"2.824343e-09","7":"1.694606e-08","8":"TRUE","9":"0.6147604","_rn_":"TMSB4X"},{"1":"TMSB4X","2":"gRNA_19","3":"TMSB4X_ENH1","4":"TMSB4X","5":"ENH","6":"6.759201e-04","7":"4.055521e-03","8":"TRUE","9":"0.5710100","_rn_":"TMSB4X.2"},{"1":"CTSC","2":"gRNA_120","3":"CTSC_ENH1","4":"CTSC","5":"ENH","6":"3.367782e-163","7":"3.367782e-163","8":"TRUE","9":"0.6918790","_rn_":"CTSC.3"},{"1":"IER3","2":"gRNA_126","3":"IER3_ENH","4":"IER3","5":"ENH","6":"4.186548e-04","7":"2.135140e-02","8":"TRUE","9":"0.6126283","_rn_":"IER3.5"},{"1":"BTG1","2":"gRNA_159","3":"BTG1_ENH","4":"BTG1","5":"ENH","6":"2.616682e-05","7":"1.046673e-04","8":"TRUE","9":"0.6811255","_rn_":"BTG1.2"},{"1":"ETS1","2":"gRNA_245","3":"ETS1_TSS","4":"ETS1","5":"TSS","6":"9.464197e-05","7":"1.892839e-04","8":"TRUE","9":"0.6112453","_rn_":"ETS1"},{"1":"ETS1","2":"gRNA_248","3":"ETS1_TSS","4":"ETS1","5":"TSS","6":"5.454655e-05","7":"6.738503e-05","8":"TRUE","9":"0.6040628","_rn_":"ETS1.3"},{"1":"PTS","2":"gRNA_271","3":"PTS_TSS","4":"PTS","5":"TSS","6":"8.296396e-79","7":"8.296396e-78","8":"TRUE","9":"0.6127314","_rn_":"PTS.10"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

####

---

Overall, all methods are able to recover true positives well, and report largely consistent results. The `wilcoxon rank sum test` recovers slightly fewer significant expected results (147), and `limma-voom` recovers the most (170).

#### Additional DE genes around the target{.tabset}

Next, we check how many other genes are found significantly changing expression in the 2Mb interval around the target.

##### sceptre

- 40% of all gRNAs result in a single DE gene, and it corresponds to the expected gene in most cases (75%).
- Another 40% of all gRNAs result in 2 to 4 DE genes and, again, most of these include the expected gene (79%).
- A small number of gRNAs induce expression changes in a larger (>10) number of genes.


```r
## total DEGs per gRNA (in positive control perturbations)
n_DE <- as.data.frame(table(sceptre.exp[sceptre.exp$sig,]$gRNA))
colnames(n_DE) <- c("gRNA", "n_DEGs_all")

# exclude TP
tmp <- sceptre.exp[which(sceptre.exp$gene_id != sceptre.exp$expected_gene),]
tmp <- tmp[tmp$sig,]
tmp <- as.data.frame(table(tmp$gRNA_id))
n_DE$other <- tmp[match(as.character(n_DE$gRNA), as.character(tmp$Var1)),2]
n_DE[is.na(n_DE$other),]$other <- 0
n_DE$expected <- n_DE$n_DEGs_all - n_DE$other

df <- reshape2::melt(n_DE[,-2])
df$gRNA <- factor(df$gRNA, levels=unique(df[order(df$value, decreasing = TRUE),]$gRNA))
df$class <- factor(gRNA_ann[match(df$gRNA, gRNA_ann$ID),]$class, levels=c("TSS", "LCR", "ENH"))
colnames(df)[2] <- "DEgene"

ggplot(df, aes(gRNA, value, fill=DEgene)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("grey60", "indianred")) +
  facet_grid(~class, scale="free_x", space="free_x") +
  xlab("all gRNAs with at least one DE gene") +
  ylab("number of DE genes") +
  ggtitle("sceptre") +
  th + theme(axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             panel.border = element_rect(fill=NA, size=0.5),
             panel.grid.major.y = element_line(),
             legend.position = "bottom")
```

![](05_perturbation_analysis_files/figure-html/sceptre_total_degs-1.png)<!-- -->

The gRNAs that result in 10 or more DE genes correspond to `TSS` perturbations and one `ENH` target. Given the `TSS` controls involve proteins important for T cell biology, it is expected that these might have more widespread effects in the transcriptome.


```r
many <- n_DE[n_DE$n_DEGs_all >= 10,-3]
many <- cbind(many, gRNA_ann[match(many$gRNA, gRNA_ann$ID),2:4])
many
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["gRNA"],"name":[1],"type":["fct"],"align":["left"]},{"label":["n_DEGs_all"],"name":[2],"type":["int"],"align":["right"]},{"label":["expected"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Symbol"],"name":[4],"type":["chr"],"align":["left"]},{"label":["class"],"name":[5],"type":["chr"],"align":["left"]},{"label":["target"],"name":[6],"type":["chr"],"align":["left"]}],"data":[{"1":"gRNA_1","2":"19","3":"1","4":"ENH_YPEL3_1","5":"ENH","6":"YPEL3_ENH","_rn_":"1"},{"1":"gRNA_2","2":"10","3":"1","4":"ENH_YPEL3_2","5":"ENH","6":"YPEL3_ENH","_rn_":"59"},{"1":"gRNA_230","2":"10","3":"1","4":"TSS_JUNB_2","5":"TSS","6":"JUNB_TSS","_rn_":"83"},{"1":"gRNA_241","2":"10","3":"1","4":"TSS_IRF1_1","5":"TSS","6":"IRF1_TSS","_rn_":"93"},{"1":"gRNA_242","2":"10","3":"1","4":"TSS_IRF1_2","5":"TSS","6":"IRF1_TSS","_rn_":"94"},{"1":"gRNA_243","2":"10","3":"1","4":"TSS_IRF1_3","5":"TSS","6":"IRF1_TSS","_rn_":"95"},{"1":"gRNA_267","2":"21","3":"1","4":"TSS_SIRT7_3","5":"TSS","6":"SIRT7_TSS","_rn_":"115"},{"1":"gRNA_294","2":"10","3":"1","4":"TSS_TMSB10_2","5":"TSS","6":"TMSB10_TSS","_rn_":"137"},{"1":"gRNA_301","2":"14","3":"1","4":"TSS_FTL_1","5":"TSS","6":"FTL_TSS","_rn_":"145"},{"1":"gRNA_302","2":"15","3":"1","4":"TSS_FTL_2","5":"TSS","6":"FTL_TSS","_rn_":"146"},{"1":"gRNA_311","2":"11","3":"1","4":"TSS_STAT3_3","5":"TSS","6":"STAT3_TSS","_rn_":"155"},{"1":"gRNA_314","2":"10","3":"1","4":"TSS_RUNX3_2","5":"TSS","6":"RUNX3_TSS","_rn_":"157"},{"1":"gRNA_76","2":"12","3":"1","4":"TSS_TRIM26_4","5":"TSS","6":"TRIM26_TSS","_rn_":"178"},{"1":"gRNA_77","2":"15","3":"1","4":"TSS_SHARPIN_1","5":"TSS","6":"SHARPIN_TSS","_rn_":"179"},{"1":"gRNA_80","2":"15","3":"1","4":"TSS_SHARPIN_4","5":"TSS","6":"SHARPIN_TSS","_rn_":"182"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

Importantly, the set of *additional* DE genes are most often only identified with one out of the four gRNAs. In contrast, 80% of the expected DE genes are supported by at least 2 independent gRNAs.


```r
degs <- table(paste0(sceptre.exp[sceptre.exp$sig,]$gene_id, ".", sceptre.exp[sceptre.exp$sig,]$target))
degs.exp <- table(paste0(sceptre.tp[sceptre.tp$sig,]$gene_id, ".", sceptre.tp[sceptre.tp$sig,]$target))

df <- cbind(1:4,
            round(prop.table(table(degs))*100, 2),
            round(prop.table(table(degs.exp))*100, 2))
colnames(df) <- c("n_sig_gRNAs", "all", "expected")
df
```

```
##   n_sig_gRNAs   all expected
## 1           1 60.16    19.64
## 2           2 21.37    17.86
## 3           3 12.14    30.36
## 4           4  6.33    32.14
```

This suggests that at least a fraction of the additional DE genes are false positives and/or off-target effects.

##### MAST

Similar results observed with `MAST`.

- 40% of all gRNAs result in a single DE gene, and it corresponds to the expected gene in most cases (72%).
- Another 40% of all gRNAs result in 2 to 4 DE genes and, again, most of these include the expected gene (82%).

However, `MAST` identifies fewer additional *hits* compared to `sceptre`, especially for larger number of DE genes (except for one outlier `ENH` gRNA that returns a very large number of DEGs, the number of gRNAs with 10 or more DEGs is only 7 compared to 15 for `sceptre`).


```r
## total DEGs per gRNA (in positive control perturbations)
n_DE <- as.data.frame(table(mast.exp[mast.exp$sig,]$gRNA))
colnames(n_DE) <- c("gRNA", "n_DEGs_all")

# exclude TP
tmp <- mast.exp[which(mast.exp$gene_id != mast.exp$expected_gene),]
tmp <- tmp[tmp$sig,]
tmp <- as.data.frame(table(tmp$gRNA_id))
n_DE$other <- tmp[match(as.character(n_DE$gRNA), as.character(tmp$Var1)),2]
n_DE[is.na(n_DE$other),]$other <- 0
n_DE$expected <- n_DE$n_DEGs_all - n_DE$other

df <- reshape2::melt(n_DE[,-2])
df$gRNA <- factor(df$gRNA, levels=unique(df[order(df$value, decreasing = TRUE),]$gRNA))
df$class <- factor(gRNA_ann[match(df$gRNA, gRNA_ann$ID),]$class, levels=c("TSS", "LCR", "ENH"))
colnames(df)[2] <- "DEgene"

ggplot(df, aes(gRNA, value, fill=DEgene)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("grey60", "indianred")) +
  facet_grid(~class, scale="free_x", space="free_x") +
  xlab("all gRNAs with at least one DE gene") +
  ylab("number of DE genes") +
  ggtitle("MAST") +
  th + theme(axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             panel.border = element_rect(fill=NA, size=0.5),
             panel.grid.major.y = element_line(),
             legend.position = "bottom")
```

![](05_perturbation_analysis_files/figure-html/mast_total_degs-1.png)<!-- -->

Again, the gRNAs that result in 10 or more DE genes correspond to `TSS` perturbations and one `ENH` target.


```r
many <- n_DE[n_DE$n_DEGs_all >= 10,-3]
many <- cbind(many, gRNA_ann[match(many$gRNA, gRNA_ann$ID),2:4])
many
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["gRNA"],"name":[1],"type":["fct"],"align":["left"]},{"label":["n_DEGs_all"],"name":[2],"type":["int"],"align":["right"]},{"label":["expected"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Symbol"],"name":[4],"type":["chr"],"align":["left"]},{"label":["class"],"name":[5],"type":["chr"],"align":["left"]},{"label":["target"],"name":[6],"type":["chr"],"align":["left"]}],"data":[{"1":"gRNA_1","2":"36","3":"1","4":"ENH_YPEL3_1","5":"ENH","6":"YPEL3_ENH","_rn_":"1"},{"1":"gRNA_242","2":"10","3":"1","4":"TSS_IRF1_2","5":"TSS","6":"IRF1_TSS","_rn_":"102"},{"1":"gRNA_243","2":"10","3":"1","4":"TSS_IRF1_3","5":"TSS","6":"IRF1_TSS","_rn_":"103"},{"1":"gRNA_267","2":"18","3":"1","4":"TSS_SIRT7_3","5":"TSS","6":"SIRT7_TSS","_rn_":"124"},{"1":"gRNA_76","2":"14","3":"1","4":"TSS_TRIM26_4","5":"TSS","6":"TRIM26_TSS","_rn_":"187"},{"1":"gRNA_77","2":"15","3":"1","4":"TSS_SHARPIN_1","5":"TSS","6":"SHARPIN_TSS","_rn_":"188"},{"1":"gRNA_80","2":"11","3":"1","4":"TSS_SHARPIN_4","5":"TSS","6":"SHARPIN_TSS","_rn_":"191"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

Again, many of the additional hits are only supported by a single gRNA, whereas expected genes are usually significant with independent gRNAs.


```r
degs <- table(paste0(mast.exp[mast.exp$sig,]$gene_id, ".", mast.exp[mast.exp$sig,]$target))
degs.exp <- table(paste0(mast.tp[mast.tp$sig,]$gene_id, ".", mast.tp[mast.tp$sig,]$target))

df <- cbind(1:4,
            round(prop.table(table(degs))*100, 2),
            round(prop.table(table(degs.exp))*100, 2))
colnames(df) <- c("n_sig_gRNAs", "all", "expected")
df
```

```
##   n_sig_gRNAs   all expected
## 1           1 60.45    17.86
## 2           2 22.56    17.86
## 3           3  8.36    23.21
## 4           4  8.64    41.07
```

##### limma

Results with `limma-voom` are very similar to those from `MAST`.

- 39% of all gRNAs result in a single DE gene, and it corresponds to the expected gene in most cases (77%).
- Another 43% of all gRNAs result in 2 to 4 DE genes and, again, most of these include the expected gene (82%).
- Similar distribution of number of additional genes as seen with `MAST`.


```r
## total DEGs per gRNA (in positive control perturbations)
n_DE <- as.data.frame(table(limma.exp[limma.exp$sig,]$gRNA))
colnames(n_DE) <- c("gRNA", "n_DEGs_all")

# exclude TP
tmp <- limma.exp[which(limma.exp$gene_id != limma.exp$expected_gene),]
tmp <- tmp[tmp$sig,]
tmp <- as.data.frame(table(tmp$gRNA_id))
n_DE$other <- tmp[match(as.character(n_DE$gRNA), as.character(tmp$Var1)),2]
n_DE[is.na(n_DE$other),]$other <- 0
n_DE$expected <- n_DE$n_DEGs_all - n_DE$other

df <- reshape2::melt(n_DE[,-2])
df$gRNA <- factor(df$gRNA, levels=unique(df[order(df$value, decreasing = TRUE),]$gRNA))
df$class <- factor(gRNA_ann[match(df$gRNA, gRNA_ann$ID),]$class, levels=c("TSS", "LCR", "ENH"))
colnames(df)[2] <- "DEgene"

ggplot(df, aes(gRNA, value, fill=DEgene)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("grey60", "indianred")) +
  facet_grid(~class, scale="free_x", space="free_x") +
  xlab("all gRNAs with at least one DE gene") +
  ylab("number of DE genes") +
  ggtitle("limma-voom") +
  th + theme(axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             panel.border = element_rect(fill=NA, size=0.5),
             panel.grid.major.y = element_line(),
             legend.position = "bottom")
```

![](05_perturbation_analysis_files/figure-html/limma_total_degs-1.png)<!-- -->

Again, the gRNAs that result in 10 or more DE genes correspond to `TSS` perturbations and one `ENH` target.


```r
many <- n_DE[n_DE$n_DEGs_all >= 10,-3]
many <- cbind(many, gRNA_ann[match(many$gRNA, gRNA_ann$ID),2:4])
many
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["gRNA"],"name":[1],"type":["fct"],"align":["left"]},{"label":["n_DEGs_all"],"name":[2],"type":["int"],"align":["right"]},{"label":["expected"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Symbol"],"name":[4],"type":["chr"],"align":["left"]},{"label":["class"],"name":[5],"type":["chr"],"align":["left"]},{"label":["target"],"name":[6],"type":["chr"],"align":["left"]}],"data":[{"1":"gRNA_1","2":"32","3":"1","4":"ENH_YPEL3_1","5":"ENH","6":"YPEL3_ENH","_rn_":"1"},{"1":"gRNA_230","2":"10","3":"1","4":"TSS_JUNB_2","5":"TSS","6":"JUNB_TSS","_rn_":"96"},{"1":"gRNA_242","2":"11","3":"1","4":"TSS_IRF1_2","5":"TSS","6":"IRF1_TSS","_rn_":"108"},{"1":"gRNA_243","2":"10","3":"1","4":"TSS_IRF1_3","5":"TSS","6":"IRF1_TSS","_rn_":"109"},{"1":"gRNA_267","2":"17","3":"1","4":"TSS_SIRT7_3","5":"TSS","6":"SIRT7_TSS","_rn_":"129"},{"1":"gRNA_302","2":"11","3":"1","4":"TSS_FTL_2","5":"TSS","6":"FTL_TSS","_rn_":"159"},{"1":"gRNA_80","2":"11","3":"1","4":"TSS_SHARPIN_4","5":"TSS","6":"SHARPIN_TSS","_rn_":"197"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

Again, many of the additional hits are only supported by a single gRNA, whereas expected genes are usually significant with independent gRNAs.


```r
degs <- table(paste0(limma.exp[limma.exp$sig,]$gene_id, ".", limma.exp[limma.exp$sig,]$target))
degs.exp <- table(paste0(limma.tp[limma.tp$sig,]$gene_id, ".", limma.tp[limma.tp$sig,]$target))

df <- cbind(1:4,
            round(prop.table(table(degs))*100, 2),
            round(prop.table(table(degs.exp))*100, 2))
colnames(df) <- c("n_sig_gRNAs", "all", "expected")
df
```

```
##   n_sig_gRNAs   all expected
## 1           1 62.10    13.79
## 2           2 20.70    22.41
## 3           3  8.33    20.69
## 4           4  8.87    43.10
```

##### wilcoxon

Similar results with the `wilcoxon rank sum` test but, again, lower number of additional hits.

- 43% of all gRNAs result in a single DE gene, and it corresponds to the expected gene in most cases (70%).
- Another 42% of all gRNAs result in 2 to 4 DE genes and, again, most of these include the expected gene (81%).


```r
## total DEGs per gRNA (in positive control perturbations)
n_DE <- as.data.frame(table(wilcox.exp[wilcox.exp$sig,]$gRNA))
colnames(n_DE) <- c("gRNA", "n_DEGs_all")

# exclude TP
tmp <- wilcox.exp[which(wilcox.exp$gene_id != wilcox.exp$expected_gene),]
tmp <- tmp[tmp$sig,]
tmp <- as.data.frame(table(tmp$gRNA_id))
n_DE$other <- tmp[match(as.character(n_DE$gRNA), as.character(tmp$Var1)),2]
n_DE[is.na(n_DE$other),]$other <- 0
n_DE$expected <- n_DE$n_DEGs_all - n_DE$other

df <- reshape2::melt(n_DE[,-2])
df$gRNA <- factor(df$gRNA, levels=unique(df[order(df$value, decreasing = TRUE),]$gRNA))
df$class <- factor(gRNA_ann[match(df$gRNA, gRNA_ann$ID),]$class, levels=c("TSS", "LCR", "ENH"))
colnames(df)[2] <- "DEgene"

ggplot(df, aes(gRNA, value, fill=DEgene)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("grey60", "indianred")) +
  facet_grid(~class, scale="free_x", space="free_x") +
  xlab("all gRNAs with at least one DE gene") +
  ylab("number of DE genes") +
  ggtitle("wilcoxon rank sum") +
  th + theme(axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             panel.border = element_rect(fill=NA, size=0.5),
             panel.grid.major.y = element_line(),
             legend.position = "bottom")
```

![](05_perturbation_analysis_files/figure-html/wilcox_total_degs-1.png)<!-- -->



```r
many <- n_DE[n_DE$n_DEGs_all >= 10,-3]
many <- cbind(many, gRNA_ann[match(many$gRNA, gRNA_ann$ID),2:4])
many
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["gRNA"],"name":[1],"type":["fct"],"align":["left"]},{"label":["n_DEGs_all"],"name":[2],"type":["int"],"align":["right"]},{"label":["expected"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Symbol"],"name":[4],"type":["chr"],"align":["left"]},{"label":["class"],"name":[5],"type":["chr"],"align":["left"]},{"label":["target"],"name":[6],"type":["chr"],"align":["left"]}],"data":[{"1":"gRNA_1","2":"40","3":"1","4":"ENH_YPEL3_1","5":"ENH","6":"YPEL3_ENH","_rn_":"1"},{"1":"gRNA_243","2":"10","3":"1","4":"TSS_IRF1_3","5":"TSS","6":"IRF1_TSS","_rn_":"92"},{"1":"gRNA_267","2":"11","3":"1","4":"TSS_SIRT7_3","5":"TSS","6":"SIRT7_TSS","_rn_":"112"},{"1":"gRNA_271","2":"10","3":"1","4":"TSS_PTS_3","5":"TSS","6":"PTS_TSS","_rn_":"116"},{"1":"gRNA_304","2":"13","3":"1","4":"TSS_FTL_4","5":"TSS","6":"FTL_TSS","_rn_":"145"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

As mentioned before, the support from 3/4 gRNAs of expected genes is a little lower compared to other methods, and we still observe that many of the additional hits come from a single gRNA.


```r
degs <- table(paste0(wilcox.exp[wilcox.exp$sig,]$gene_id, ".", wilcox.exp[wilcox.exp$sig,]$target))
degs.exp <- table(paste0(wilcox.tp[wilcox.tp$sig,]$gene_id, ".", wilcox.tp[wilcox.tp$sig,]$target))

df <- cbind(1:4,
            round(prop.table(table(degs))*100, 2),
            round(prop.table(table(degs.exp))*100, 2))
colnames(df) <- c("n_sig_gRNAs", "all", "expected")
df
```

```
##   n_sig_gRNAs   all expected
## 1           1 59.64    18.18
## 2           2 23.80    27.27
## 3           3  9.34    23.64
## 4           4  7.23    30.91
```

#### P-value calibration {.tabset}

To get an idea on how likely it is that some of the additional hits are false positives, we check how often do NT gRNAs result in significant results. We test the 35 NT gRNAs against all genes within 1MB of a targeted locus (40250 tests).

##### sceptre

`sceptre` reports 3415 gene-NTgRNA pairs with a significant p-value, under a FDR of 5%. This indicates `sceptre` does not achieve adequate FDR control. Indeed, the p-values for NT gRNAs are heavily inflated, indicating that `sceptre` is not appropriate for this dataset, where the MOI is 1.


```r
## code adapted from https://github.com/Katsevich-Lab/sceptre-manuscript/blob/master/sceptre_paper/
## check p-value callibration
ci <- 0.95
p_thresh <- 1e-8

pvals_nt <- sceptre.nt[,c('gene_id', 'gRNA_id', 'p_value')]

## compute expected p-values for QQplots
df <- pvals_nt %>%
    mutate(r = rank(p_value), expected = ppoints(n())[r],
           clower = qbeta(p=(1-ci)/2, shape1 = r, shape2 = n()+1-r),
           cupper = qbeta(p=(1+ci)/2, shape1 = r, shape2 = n()+1-r)) %>%
  ungroup() %>% mutate()

df %>% 
  mutate(p_value = ifelse(p_value < p_thresh, p_thresh, p_value)) %>%
  ggplot(aes(x = expected, y = p_value, ymin = clower, ymax = cupper)) +
    geom_point(size = 1, alpha = 0.5) +
    geom_ribbon(alpha = 0.2) +
    geom_abline(intercept = 0, slope = 1) +
    scale_x_continuous(trans = revlog_trans(base = 10)) + 
    scale_y_continuous(trans = revlog_trans(base = 10)) +
    xlab(expression(paste("Expected null p-value"))) +
    ylab(expression(paste("Observed p-value"))) +
    ggtitle("sceptre") +
    th
```

![](05_perturbation_analysis_files/figure-html/sceptre_nt-1.png)<!-- -->

##### MAST

`MAST` reports 2811 gene-NTgRNA pairs with a significant p-value, under a FDR of 5%. A bit better than `sceptre` but still showing lack of adequate FDR control, and heavily inflated p-values.


```r
## check p-value callibration
ci <- 0.95
p_thresh <- 1e-8

pvals_nt <- mast.nt[,c('gene_id', 'gRNA_id', 'p_value')]

## compute expected p-values for QQplots
df <- pvals_nt %>%
    mutate(r = rank(p_value), expected = ppoints(n())[r],
           clower = qbeta(p=(1-ci)/2, shape1 = r, shape2 = n()+1-r),
           cupper = qbeta(p=(1+ci)/2, shape1 = r, shape2 = n()+1-r)) %>%
  ungroup() %>% mutate()

df %>% 
  mutate(p_value = ifelse(p_value < p_thresh, p_thresh, p_value)) %>%
  ggplot(aes(x = expected, y = p_value, ymin = clower, ymax = cupper)) +
    geom_point(size = 1, alpha = 0.5) +
    geom_ribbon(alpha = 0.2) +
    geom_abline(intercept = 0, slope = 1) +
    scale_x_continuous(trans = revlog_trans(base = 10)) + 
    scale_y_continuous(trans = revlog_trans(base = 10)) +
    xlab(expression(paste("Expected null p-value"))) +
    ylab(expression(paste("Observed p-value"))) +
    ggtitle("MAST") +
    th
```

![](05_perturbation_analysis_files/figure-html/mast_nt-1.png)<!-- -->

##### limma

`limma-voom` reports 3109 gene-NTgRNA pairs with a significant p-value, and also shows heavily inflated p-values.


```r
## check p-value callibration
ci <- 0.95
p_thresh <- 1e-8

pvals_nt <- limma.nt[,c('gene_id', 'gRNA_id', 'p_value')]

## compute expected p-values for QQplots
df <- pvals_nt %>%
    mutate(r = rank(p_value), expected = ppoints(n())[r],
           clower = qbeta(p=(1-ci)/2, shape1 = r, shape2 = n()+1-r),
           cupper = qbeta(p=(1+ci)/2, shape1 = r, shape2 = n()+1-r)) %>%
  ungroup() %>% mutate()

df %>% 
  mutate(p_value = ifelse(p_value < p_thresh, p_thresh, p_value)) %>%
  ggplot(aes(x = expected, y = p_value, ymin = clower, ymax = cupper)) +
    geom_point(size = 1, alpha = 0.5) +
    geom_ribbon(alpha = 0.2) +
    geom_abline(intercept = 0, slope = 1) +
    scale_x_continuous(trans = revlog_trans(base = 10)) + 
    scale_y_continuous(trans = revlog_trans(base = 10)) +
    xlab(expression(paste("Expected null p-value"))) +
    ylab(expression(paste("Observed p-value"))) +
    ggtitle("limma-voom") +
    th
```

![](05_perturbation_analysis_files/figure-html/limma_nt-1.png)<!-- -->

##### wilcoxon

The `wilcoxon rank sum` test reports 2887 gene-NTgRNA pairs with a significant p-value, the second lowest total. 


```r
## check p-value callibration
ci <- 0.95
p_thresh <- 1e-8

pvals_nt <- wilcox.nt[,c('gene_id', 'gRNA_id', 'p_value')]

## compute expected p-values for QQplots
df <- pvals_nt %>%
    mutate(r = rank(p_value), expected = ppoints(n())[r],
           clower = qbeta(p=(1-ci)/2, shape1 = r, shape2 = n()+1-r),
           cupper = qbeta(p=(1+ci)/2, shape1 = r, shape2 = n()+1-r)) %>%
  ungroup() %>% mutate()

df %>% 
  mutate(p_value = ifelse(p_value < p_thresh, p_thresh, p_value)) %>%
  ggplot(aes(x = expected, y = p_value, ymin = clower, ymax = cupper)) +
    geom_point(size = 1, alpha = 0.5) +
    geom_ribbon(alpha = 0.2) +
    geom_abline(intercept = 0, slope = 1) +
    scale_x_continuous(trans = revlog_trans(base = 10)) + 
    scale_y_continuous(trans = revlog_trans(base = 10)) +
    xlab(expression(paste("Expected null p-value"))) +
    ylab(expression(paste("Observed p-value"))) +
    ggtitle("wilcoxon rank sum") +
    th
```

![](05_perturbation_analysis_files/figure-html/wilcox_nt-1.png)<!-- -->


### Method comparison

#### True positive recovery


```r
compare <- gRNA_ann[,c(1,4,3,5)]
compare$sceptre <- sceptre[match(paste0(compare$ID, ".", compare$expected_DE_gene),
                                 paste0(sceptre$gRNA_id, ".", sceptre$gene_id)),]$FDR
compare$sceptre_sig <- ifelse(compare$sceptre < 0.05, 1, 0)

compare$mast <- mast[match(paste0(compare$ID, ".", compare$expected_DE_gene),
                           paste0(mast$gRNA_id, ".", mast$gene_id)),]$FDR
compare$mast_sig <- ifelse(compare$mast < 0.05, 1, 0)

compare$limma <- limma[match(paste0(compare$ID, ".", compare$expected_DE_gene),
                             paste0(limma$gRNA_id, ".", limma$gene_id)),]$FDR
compare$limma_sig <- ifelse(compare$limma < 0.05, 1, 0)

compare$wilcoxon <- wilcox[match(paste0(compare$ID, ".", compare$expected_DE_gene),
                                 paste0(wilcox$gRNA_id, ".", wilcox$gene_id)),]$FDR
compare$wilcoxon_sig <- ifelse(compare$wilcoxon < 0.05, 1, 0)

compare <- compare[!is.na(compare$sceptre),]
compare$sum <- rowSums(compare[,seq(6,12,2)])
```

All methods are able to recover true positives and report similar results. From **248 different gRNA-expected_gene pairs**:

- 187 are deemed significant by at least one method.
  + From the 61 pairs without any significant results, the majority correspond to `ENH` perturbations, but there are also a number of `TSS` pairs.
  

```r
table(compare[compare$sum==0,]$class)
```

```
## 
## ENH LCR TSS 
##  41   2  18
```
  
- 
  + The `ENH` gRNAs without detected effects correspond to 21 different targets. From these, two thirds have significant results for other gRNAs targeting the same locus, suggesting that the lack of results stem from non-functional gRNAs and/or small effects that are harder to detect. Instead, 6 targets have at most a single significant gRNA and might correspond to non-functional enhancers in T cells, or very small effects.
  

```r
tmp <- table(n_gRNAs_without_DEGs=compare[compare$sum==0 & compare$class == "ENH",]$target)
tmp[order(tmp)]
```

```
## n_gRNAs_without_DEGs
##           BTG2_ENH          CTSC_ENH1          KLF6_ENH2          PHF19_ENH 
##                  1                  1                  1                  1 
##            PKM_ENH            TKT_ENH        TMSB4X_ENH3          YPEL3_ENH 
##                  1                  1                  1                  1 
##          EPB41_ENH         GIGYF1_ENH           IER3_ENH          KLF6_ENH1 
##                  2                  2                  2                  2 
##          PTPRC_ENH         SLC2A3_ENH        TMSB4X_ENH2           BBC3_ENH 
##                  2                  2                  2                  3 
##           BTG1_ENH           CD69_ENH          CERS2_ENH          SOCS3_ENH 
##                  3                  3                  3                  3 
## CYSTM1_UBE2D2D_ENH 
##                  4
```

- 
  + The `TSS` gRNAs without detected DE results correspond to 14 different targets. In this case, all but one of the targets have significant results for other gRNAs, suggesting the lack of hits derive from non-functional gRNAs.


```r
tmp <- table(n_gRNAs_without_DEGs=compare[compare$sum==0 & compare$class == "TSS",]$target)
tmp[order(tmp)]
```

```
## n_gRNAs_without_DEGs
##        ETS1_TSS        LCP1_TSS        MTX3_TSS     PPP2R1B_TSS        RORA_TSS 
##               1               1               1               1               1 
##       RTKN2_TSS     SHARPIN_TSS        SIK2_TSS      TRIM26_TSS WDR78_MIER1_TSS 
##               1               1               1               1               1 
##        XBP1_TSS     PABPC1L_TSS       RUNX3_TSS       GATA3_TSS 
##               1               2               2               3
```

- Each method reports ~150 significant results.
- 127 gRNA-gene pairs are detected by all four methods.
  + A further 14 are identified by all except the `wilcoxon rank sum test`. 
- `sceptre` misses 18 hits recovered by other methods.
- All methods except `MAST` show a small number of hits that are not reported by any other method.


```r
upset(compare[, c(seq(6,12,2),3)], order.by = "freq",
      queries = list(
        list(query = elements, 
             params = list("class", c("TSS","ENH","LCR")), query.name = "LCR",
             color = GetColors(n=5, scheme = "muted")[5], active = TRUE),
        list(query = elements, 
             params = list("class", c("TSS","ENH")), query.name = "ENH",
             color = GetColors(n=5, scheme = "muted")[1], active = TRUE),
        list(query = elements, 
             params = list("class", "TSS"), query.name = "TSS",
             color = GetColors(n=5, scheme = "muted")[3], active = TRUE)),
      query.legend="bottom")
```

![](05_perturbation_analysis_files/figure-html/compare_methods_upset-1.png)<!-- -->

When looking at the groups of genes that are detected by only a subset of methods, compared to those detected by all (or none), we observe that:

- As expected, pairs that are not detected by any method tend to involve genes expressed at lower levels and gRNAs detected in lower number of cells.
- Pairs missed by the `wilcoxon rank sum` test involve genes with low expression and variability.
- Pairs that are not reported by `sceptre` involve genes with higher technical variability.

- `MAST` and `limma-voom` are able to report significant findings for genes expressed at very low levels, also observed for pairs reported by `limma-voom` only.
- Pairs reported only by `limma-voom` involve gRNAs detected in smaller numbers of cells. It is likely that with better gRNA recovery the other methods can pick them up too.
- Pairs reported only by `sceptre` involve genes expressed at higher levels, but with lower variability.
- Pairs reported only by the `wilcoxon rank sum` test also show high variability.


```r
## annotate the hits specific to some methods
compare$group <- NA
compare[compare$sum==0,]$group <- "none"
compare[compare$sum==4,]$group <- "all"
compare[compare$sum==1 & compare$sceptre_sig,]$group <- "sceptre_only"
compare[compare$sum==1 & compare$limma_sig,]$group <- "limma_only"
compare[compare$sum==1 & compare$wilcoxon_sig,]$group <- "wilcox_only"
compare[compare$sum==3 & !compare$sceptre_sig,]$group <- "sceptre_miss"
compare[compare$sum==3 & !compare$wilcoxon_sig,]$group <- "wilcox_miss"
compare[compare$sum==2 & compare$mast_sig & compare$limma_sig,]$group <- "mast_limma"

## correlate to features that could explain (lack of) detection
var_decomp <- modelGeneVar(sce, BPPARAM = bpp)

df <- compare[,c(1:4,ncol(compare))]
df <- cbind(df, var_decomp[match(df$expected_DE_gene, row.names(var_decomp)), 1:4])
df <- df[!is.na(df$group),] # for simplicity, remove the groups that have very low counts

df$group <- factor(df$group, levels=c("none", "all", "wilcox_miss", "sceptre_miss", "mast_limma", 
                                      paste0(c("limma", "sceptre", "wilcox"), "_only")))
df$n_cells <- ncells_per_gRNA[df$ID]

plots <- list()
plots[[1]] <- ggplot(df, aes(group, n_cells)) +
  scale_y_log10() +
  geom_violin(scale="width") +
  geom_boxplot(width=0.1) +
  facet_wrap(~group, scale="free_x", ncol=8) +
  xlab("") +
  ylab("cells per gRNA") +
  th + theme(axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             panel.grid.major.y = element_line())
plots[[2]] <- ggplot(df, aes(group, mean)) +
  geom_violin(scale="width") +
  geom_boxplot(width=0.1) +
  facet_wrap(~group, scale="free_x", ncol=8) +
  xlab("") +
  ylab(expression('log'[2]*' mean expression')) +
  th + theme(axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             strip.background = element_blank(),
             strip.text.x = element_blank(),
             panel.grid.major.y = element_line())
tmp <- reshape2::melt(df[,c(5,7:9)])
plots[[3]] <- ggplot(tmp, aes(group, value, colour=variable)) +
  geom_violin(scale="width") +
  geom_boxplot(width=0.1, position = position_dodge(0.9)) +
  scale_color_manual(values = c(total="grey30", tech="grey80", bio="indianred2")) +
  facet_wrap(~group, scale="free_x", ncol=8) +
  xlab("") +
  ylab("variation") +
  labs(colour="") +
  th + theme(axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             panel.grid.major.y = element_line(),
             strip.background = element_blank(),
             strip.text.x = element_blank(), 
             legend.position = "bottom")
ggarrange(plotlist = plots, ncol=1, align="v", heights = c(0.35,0.3,0.4))
```

![](05_perturbation_analysis_files/figure-html/compare_subsets-1.png)<!-- -->

#### Other DE genes


```r
compare <- data.frame(pair = unique(c(paste(sceptre.exp[sceptre.exp$sig,]$gRNA_id, sceptre.exp[sceptre.exp$sig,]$gene_id, sep="|"),
                                      paste(mast.exp[mast.exp$sig,]$gRNA_id, mast.exp[mast.exp$sig,]$gene_id, sep="|"),
                                      paste(limma.exp[limma.exp$sig,]$gRNA_id, limma.exp[limma.exp$sig,]$gene_id, sep="|"),
                                      paste(wilcox.exp[wilcox.exp$sig,]$gRNA_id, wilcox.exp[wilcox.exp$sig,]$gene_id, sep="|"))))
compare$gRNA_id <- unlist(lapply(strsplit(compare$pair, "|", fixed=TRUE), '[[', 1))
compare$gene_id <- unlist(lapply(strsplit(compare$pair, "|", fixed=TRUE), '[[', 2))
compare$target <- gRNA_ann[match(compare$gRNA_id, gRNA_ann$ID),]$target
compare$class <- gRNA_ann[match(compare$gRNA_id, gRNA_ann$ID),]$class
compare$expected_DE_gene <- gRNA_ann[match(compare$gRNA_id, gRNA_ann$ID),]$expected_DE_gene
compare$expected <- compare$gene_id == compare$expected_DE_gene

compare$sceptre <- sceptre[match(compare$pair, paste0(sceptre$gRNA_id, "|", sceptre$gene_id)),]$FDR
compare$sceptre_sig <- ifelse(compare$sceptre < 0.05, 1, 0)

compare$mast <- mast[match(compare$pair, paste0(mast$gRNA_id, "|", mast$gene_id)),]$FDR
compare$mast_sig <- ifelse(compare$mast < 0.05, 1, 0)

compare$limma <- limma[match(compare$pair, paste0(limma$gRNA_id, "|", limma$gene_id)),]$FDR
compare$limma_sig <- ifelse(compare$limma < 0.05, 1, 0)

compare$wilcoxon <- wilcox[match(compare$pair, paste0(wilcox$gRNA_id, "|", wilcox$gene_id)),]$FDR
compare$wilcoxon_sig <- ifelse(compare$wilcoxon < 0.05, 1, 0)

compare$sum <- rowSums(compare[,seq(9,15,2)])
```

Overall, there are 994 gRNA-gene pairs that are deemed significant by at least one method. From these, 187 involve the expected gene, leaving a large number of *additional* DE genes in the vicinity of the targets. This set of additional DEGs contains a subset of 170 pairs that are reported by all methods. However, in contrast to the expected effects, much higher numbers of pairs are reported by a single method, indicative of likely false positives. `sceptre` and the `wilcoxon rank sum` test report substantially higher number of these *private* hits, compared to `limma-voom` and `MAST`.


```r
upset(compare[!compare$expected, c(seq(9,15,2),5)], order.by = "freq",
      queries = list(
        list(query = elements, 
             params = list("class", c("TSS","ENH","LCR")), query.name = "LCR",
             color = GetColors(n=5, scheme = "muted")[5], active = TRUE),
        list(query = elements, 
             params = list("class", c("TSS","ENH")), query.name = "ENH",
             color = GetColors(n=5, scheme = "muted")[1], active = TRUE),
        list(query = elements, 
             params = list("class", "TSS"), query.name = "TSS",
             color = GetColors(n=5, scheme = "muted")[3], active = TRUE)),
      query.legend="bottom")
```

![](05_perturbation_analysis_files/figure-html/compare_methods_other_upset-1.png)<!-- -->

Consistently we can see that

- DE genes called by several methods are mostly supported by at least 2 independent gRNAs when grouped by target.
- In contrast, DE genes that are identified by a single method have a large fraction of hits identified with a single gRNA, thus more likely to be false positive or off-target effects. 
- DE genes reported by only two methods still show a large fraction of cases supported by independent gRNAs, but have lower proportion of 3 and 4-gRNA support, indicating these might include smaller effects that are harder to identify.


```r
## add number of gRNAs supporting each hit
# group by target
compare$parir_tgt <- paste(compare$target, compare$gene_id, sep='|')

# use the maximum achieved by any method
tmp <- rbind(sceptre = as.data.frame(table(compare[compare$sceptre_sig==1,]$parir_tgt)),
            mast = as.data.frame(table(compare[compare$mast_sig==1,]$parir_tgt)),
            limma = as.data.frame(table(compare[compare$limma_sig==1,]$parir_tgt)),
            wilcox = as.data.frame(table(compare[compare$wilcoxon_sig==1,]$parir_tgt)))

tmp <- tmp %>% group_by(Var1) %>% summarise(max = max(Freq))
compare$n_gRNAs <- tmp[match(compare$parir_tgt, tmp$Var1),]$max

# replot
upset(compare[!compare$expected, c(seq(9,15,2),18)], order.by = "freq",
      queries = list(
        list(query = elements, 
             params = list("n_gRNAs", c('1','2','3','4')), query.name = '4',
             color = GetColors(n=11, scheme = "sunset")[11], active = TRUE),
        list(query = elements, 
             params = list("n_gRNAs", c('1','2','3')), query.name = '3',
             color = GetColors(n=11, scheme = "sunset")[10], active = TRUE),
        list(query = elements, 
             params = list("n_gRNAs", c('1','2')), query.name = '2',
             color = GetColors(n=11, scheme = "sunset")[8], active = TRUE),
        list(query = elements, 
             params = list("n_gRNAs", c('1')), query.name = '1',
             color = GetColors(n=11, scheme = "sunset")[6], active = TRUE)),
      query.legend="bottom")
```

![](05_perturbation_analysis_files/figure-html/compare_methods_other_ngrnas-1.png)<!-- -->

### Conclusions

Overall, we can conclude that

- all methods are good at detecting true positives, and do so consistently.
- all methods also report many additional hits. A fraction of these are likely to be false positives since p-values are miscalibrated.



```r
write.table(sceptre, paste0(dir, "results/05_DE_results.sceptre.tsv"), quote = FALSE, sep="\t", row.names = FALSE)
write.table(mast, paste0(dir, "results/05_DE_results.mast.tsv"), quote = FALSE, sep="\t", row.names = FALSE)
write.table(limma, paste0(dir, "results/05_DE_results.limma.tsv"), quote = FALSE, sep="\t", row.names = FALSE)
write.table(wilcox, paste0(dir, "results/05_DE_results.wilcox.tsv"), quote = FALSE, sep="\t", row.names = FALSE)

write.table(sceptre.nt, paste0(dir, "results/05_DE_results_NT.sceptre.tsv"), quote = FALSE, sep="\t", row.names = FALSE)
write.table(mast.nt, paste0(dir, "results/05_DE_results_NT.mast.tsv"), quote = FALSE, sep="\t", row.names = FALSE)
write.table(limma.nt, paste0(dir, "results/05_DE_results_NT.limma.tsv"), quote = FALSE, sep="\t", row.names = FALSE)
write.table(wilcox.nt, paste0(dir, "results/05_DE_results_NT.wilcox.tsv"), quote = FALSE, sep="\t", row.names = FALSE)

write.table(compare, paste0(dir, "results/05_DEGs_compareMethods.tsv"), quote = FALSE, sep="\t", row.names = FALSE)
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
##  [1] UpSetR_1.4.0                inlmisc_0.5.2              
##  [3] RColorBrewer_1.1-2          ggpubr_0.4.0               
##  [5] ComplexHeatmap_2.8.0        BiocSingular_1.8.1         
##  [7] BiocParallel_1.26.1         scales_1.1.1               
##  [9] dplyr_1.0.7                 scran_1.20.1               
## [11] scater_1.20.1               ggplot2_3.3.5              
## [13] scuttle_1.2.1               SingleCellExperiment_1.14.1
## [15] SummarizedExperiment_1.22.0 Biobase_2.52.0             
## [17] GenomicRanges_1.44.0        GenomeInfoDb_1.28.1        
## [19] IRanges_2.26.0              S4Vectors_0.30.0           
## [21] BiocGenerics_0.38.0         MatrixGenerics_1.4.2       
## [23] matrixStats_0.60.0         
## 
## loaded via a namespace (and not attached):
##   [1] readxl_1.3.1              backports_1.2.1          
##   [3] circlize_0.4.13           plyr_1.8.6               
##   [5] igraph_1.2.6              sp_1.4-5                 
##   [7] digest_0.6.27             foreach_1.5.1            
##   [9] htmltools_0.5.1.1         viridis_0.6.1            
##  [11] fansi_0.5.0               magrittr_2.0.1           
##  [13] checkmate_2.0.0           ScaledMatrix_1.0.0       
##  [15] cluster_2.1.2             doParallel_1.0.16        
##  [17] openxlsx_4.2.4            limma_3.48.3             
##  [19] colorspace_2.0-2          haven_2.4.3              
##  [21] xfun_0.25                 rgdal_1.5-23             
##  [23] crayon_1.4.1              RCurl_1.98-1.4           
##  [25] jsonlite_1.7.2            iterators_1.0.13         
##  [27] glue_1.4.2                gtable_0.3.0             
##  [29] zlibbioc_1.38.0           XVector_0.32.0           
##  [31] GetoptLong_1.0.5          DelayedArray_0.18.0      
##  [33] car_3.0-11                shape_1.4.6              
##  [35] abind_1.4-5               DBI_1.1.1                
##  [37] edgeR_3.34.0              rstatix_0.7.0            
##  [39] Rcpp_1.0.7                viridisLite_0.4.0        
##  [41] clue_0.3-59               dqrng_0.3.0              
##  [43] foreign_0.8-81            rsvd_1.0.5               
##  [45] metapod_1.0.0             ellipsis_0.3.2           
##  [47] pkgconfig_2.0.3           farver_2.1.0             
##  [49] sass_0.4.0                locfit_1.5-9.4           
##  [51] utf8_1.2.2                tidyselect_1.1.1         
##  [53] labeling_0.4.2            rlang_0.4.11             
##  [55] reshape2_1.4.4            munsell_0.5.0            
##  [57] cellranger_1.1.0          tools_4.1.1              
##  [59] generics_0.1.0            broom_0.7.9              
##  [61] evaluate_0.14             stringr_1.4.0            
##  [63] yaml_2.2.1                knitr_1.33               
##  [65] zip_2.2.0                 purrr_0.3.4              
##  [67] sparseMatrixStats_1.4.2   compiler_4.1.1           
##  [69] beeswarm_0.4.0            curl_4.3.2               
##  [71] png_0.1-7                 ggsignif_0.6.2           
##  [73] tibble_3.1.3              statmod_1.4.36           
##  [75] bslib_0.2.5.1             stringi_1.7.3            
##  [77] highr_0.9                 forcats_0.5.1            
##  [79] lattice_0.20-44           bluster_1.2.1            
##  [81] Matrix_1.3-4              vctrs_0.3.8              
##  [83] pillar_1.6.2              lifecycle_1.0.0          
##  [85] jquerylib_0.1.4           GlobalOptions_0.1.2      
##  [87] BiocNeighbors_1.10.0      cowplot_1.1.1            
##  [89] data.table_1.14.0         bitops_1.0-7             
##  [91] irlba_2.3.3               raster_3.4-13            
##  [93] R6_2.5.1                  gridExtra_2.3            
##  [95] rio_0.5.27                vipor_0.4.5              
##  [97] codetools_0.2-18          assertthat_0.2.1         
##  [99] rjson_0.2.20              withr_2.4.2              
## [101] GenomeInfoDbData_1.2.6    hms_1.1.0                
## [103] beachmat_2.8.1            tidyr_1.1.3              
## [105] rmarkdown_2.10            DelayedMatrixStats_1.14.2
## [107] carData_3.0-4             Cairo_1.5-12.2           
## [109] ggbeeswarm_0.6.0
```

