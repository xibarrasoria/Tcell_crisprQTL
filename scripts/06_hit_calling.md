---
title: "<span style='font-size: 16px; color: grey'> crisprQTL proof-of-concept experiment in T cells </span> <br> <span style='font-size: 36px'> Hit calling </span>"
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



We have compared four different approaches to detect differentially expressed genes in perturbed cells. All methods were good at detecting expected changes in expression for positive-control perturbations. However, all showed miscalibrated p-values, which is expected since we do not have biological replicates in these data. **`MAST` is our method of choice for this dataset since it showed the best trade-off between true and false positives.**


```r
# data
sce <- readRDS(paste0(dir, "results/03_sce_POC_Tcells.goodQual.calls.NORM.Rds"))
# genes detected in at least 5% of the cells
genes_to_keep <- rowSums(counts(sce)>0) > ncol(sce)*0.05
genes_detected <- rownames(sce)[genes_to_keep]

# gRNA annotation
gRNA_ann <- read.table(paste0(dir, "data/gRNA_library.tsv"), sep="\t", header = TRUE)
targets <- setdiff(gRNA_ann$target, "NT")

# define NT cells for plotting
nt_guides <- gRNA_ann[gRNA_ann$class=="NT",]$ID
cells <- colnames(sce[,sce$call == "unique"])
cells_neg <- apply(counts(altExp(sce, 'gRNA_calls'))[nt_guides, cells], 1, function(x) names(which(x)))
cells_neg <- unlist(cells_neg)
```

### Target-level results

In the previous script we observed p-values are inflated and thus we have a *large* fraction of false positives. We also observed that DE genes that are reproduced by several DE algorithms tend to be supported by independent gRNAs. Thus, the combined evidence from separate perturbations for the same target is a good approach to enrich for robust hits.

We perform hit calling at the target level, taking into account the combined information from all 4 gRNAs from each target. The raw p-values for each gRNA are integrated into a single p-value using Fisher's method. We correct for multiple testing at the screen level, considering all tests performed for genes in the vicinity of their targets. 

Considering all targets together, we recover  a little under 400 significant element-to-gene (E2G) pairs.


```r
## MAST results
mast <- read.table(paste0(dir, "results/05_DE_results.mast.tsv"), header=TRUE)
mast$pair_tgt <- paste(mast$target, mast$gene_id, sep="|")

## integrate results at target level for reproducible hits
mast.tgt <- mast %>% 
  group_by(pair_tgt) %>%
  mutate(gRNAs = paste(gRNA_id, collapse=","),
         p_values = paste(p_value, collapse=","),
         logFCs = paste(logFC, collapse=","),
         n_gRNAs = sum(p_value < 0.05),
         n_down = sum(ifelse(is.na(logFC[p_value < 0.05] < 0), TRUE, logFC[p_value < 0.05] < 0)), # if all cells are 0 FC is NA (count as down)
         n_up = sum(logFC[p_value < 0.05] > 0, na.rm=TRUE))
mast.tgt$direction <- ifelse(mast.tgt$n_down > 0 & mast.tgt$n_up > 0, "discordant",
                            ifelse(mast.tgt$n_down > 0, "down", ifelse(mast.tgt$n_up > 0, "up", NA)))

## for a summary logFC, we take the mean across all gRNAs that have a significant effect
mast.tgt <- mast.tgt %>% 
  group_by(pair_tgt) %>% 
  mutate(logFC_summary = mean(logFC[p_value < 0.05], na.rm=TRUE))

## for a summary p-value, we integrate results across gRNAs using Fisher's method
pvals.tgt <- combineGroupedPValues(p.values = mast$p_value, grouping = mast$pair_tgt, method="fisher")$p.value
mast.tgt$pval_fisher <- pvals.tgt[match(mast.tgt$pair_tgt, names(pvals.tgt))]

# keep relevant columns and uniquify
mast.tgt <- as.data.frame(mast.tgt[,c("target","gene_id","pair_tgt", "class","expected_gene", "pval_fisher","logFC_summary", 
                                      "n_gRNAs","n_down","n_up","direction", "gRNAs","p_values","logFCs")])
mast.tgt <- unique(mast.tgt)
# FDR correction
mast.tgt$FDR <- p.adjust(mast.tgt$pval_fisher, method = 'fdr')
mast.tgt <- mast.tgt[,c(1:6,15,7:14)]

mast.sig <- mast.tgt[mast.tgt$FDR < 0.05,]
c(total_DEGs = nrow(mast.sig))
```

```
## total_DEGs 
##        378
```

The magnitude of the per-target adjusted p-value nicely reflects the amount of evidence provided by the independent gRNAs. 


```r
ggplot(mast.sig, aes(as.factor(n_gRNAs), (FDR))) +
  geom_violin(scale="width") +
  geom_boxplot(width=0.1) +
  xlab("number of gRNAs < 0.05 (raw gRNA-level p-values)") +
  ylab("target-level FDR") +
  th
```

![](06_hit_calling_files/figure-html/fisher_pval_vs_n_grnas-1.png)<!-- -->

### Differentially expressed genes

The library contains 80 targets:


```r
table(unique(gRNA_ann[,c('target', 'class')])$class)[c(6,4,1,3,2)]
```

```
## 
##        TSS        LCR        ENH     INTRON INTERGENIC 
##         35          3         28         11          3
```

From these, only four have no DE genes within 1Mb. 


```r
## targets with no hits
no_degs <- setdiff(targets, mast.sig$target)
no_degs
```

```
## [1] "IL2RA_INTRON" "BTG1_ENH"     "TMSB4X_ENH2"  "IL12A_TSS"
```

The gRNAs for these four targets are generally detected in fewer cells compared to the rest of the library. Thus, it is possible that some of these perturbations do induce expression changes but we are underpowered to detect them.


```r
n_cells_per_gRNA <- as.data.frame(rowSums(counts(altExp(sce, 'gRNA_calls'))))
colnames(n_cells_per_gRNA) <- "n_cells"
n_cells_per_gRNA$target <- gRNA_ann[match(row.names(n_cells_per_gRNA), gRNA_ann$ID),'target']
n_cells_per_gRNA$no_DEGs <- n_cells_per_gRNA$target %in% no_degs

ggplot(n_cells_per_gRNA, aes(no_DEGs, n_cells)) +
  scale_y_log10() +
  annotation_logticks(side="l") +
  geom_violin() +
  geom_boxplot(width=0.1) +
  ylab("number of cells per gRNA") +
  xlab("zero DE genes detected within 1 Mb") +
  th
```

![](06_hit_calling_files/figure-html/no_hits_ncells-1.png)<!-- -->

Most perturbations result in a small number of DE genes within 1 Mb of the perturbation site, with a median of 4 DE genes.


```r
## number of hits per target
n_DE_per_tgt <- as.data.frame(table(mast.sig$target))
colnames(n_DE_per_tgt) <- c("target", "n_DE")
n_DE_per_tgt$class <- gRNA_ann[match(n_DE_per_tgt$target, gRNA_ann$target),'class']
n_DE_per_tgt$class <- factor(n_DE_per_tgt$class, levels=c("TSS", "LCR", "ENH", "INTRON", "INTERGENIC"))

summary(n_DE_per_tgt$n_DE)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   1.000   2.000   4.000   4.974   6.000  31.000
```

A large number of targets result in 1-2 DE genes only, and perturbing the `TSS` generally induces more changes than intergenic perturbations.


```r
ggplot(n_DE_per_tgt, aes(n_DE, fill=class)) +
  geom_bar() +
  scale_fill_manual(values = colors) +
  xlab("number of DE genes within 1Mb") +
  ylab("number of targets") +
  th
```

![](06_hit_calling_files/figure-html/n_DE_plot-1.png)<!-- -->

#### Tier classification

We classify DE genes into confidence tiers based on the amount of support from independent gRNAs.

- `low` tier E2G pairs: only **one** gRNA-level raw p-value < 0.05.
- `medium` tier E2G pairs: **two** gRNA-level raw p-value < 0.05, **both in the same direction**.
- `high` tier E2G pairs: **three or four** gRNA-level raw p-value < 0.05, **all in the same direction**.


```r
## confidence tiers
mast.sig$tier <- factor(ifelse(mast.sig$n_gRNAs == 1, "low",
                               ifelse(mast.sig$n_up == 2 | mast.sig$n_down == 2, "medium", 
                                       ifelse(mast.sig$n_up > 2 | mast.sig$n_down > 2, "high", "low"))), 
                        levels=c("high","medium","low"))
table(mast.sig$tier)
```

```
## 
##   high medium    low 
##     87    140    151
```

As highlighted above, 76 targets have at least one DE gene within 1Mb of the target site. Almost all of these (92%) have at least one `medium` hit, and the majority (72%) have at least one `high` hit. This indicates that the vast majority of perturbations result in a change in gene expression in *cis* that is detectable in the single-cell data.


```r
c(any_tier = length(unique(mast.sig$target)),
  at_least_medium = length(unique(mast.sig[mast.sig$tier!="low",]$target)),
  at_least_1_high = length(unique(mast.sig[mast.sig$tier=="high",]$target)))
```

```
##        any_tier at_least_medium at_least_1_high 
##              76              70              55
```

The hits in the `low` tier that are only identified with one gRNA are more likely to be off-target or false-positive effects. Many of these have small raw gRNA-level p-values, but also include gRNAs that are detected in much larger number of cells. Together, these data suggest that many are likely false positives that stem from increased statistical power due to large cell numbers. Thus, we discard these hits from downstream analyses.


```r
mast$tier <- mast.sig[match(mast$pair_tgt, mast.sig$pair_tgt),'tier']

df <- mast[mast$p_value < 0.05 & !(is.na(mast$tier)),]
df$n_cells <- rowSums(counts(altExp(sce, 'gRNA_calls')))[df$gRNA_id]

plots <- list()
plots[[1]] <- ggplot(df, aes(tier, p_value)) +
  geom_violin() +
  geom_boxplot(width=0.1) +
  xlab("") +
  ylab("raw gRNA-level p-value") +
  th
plots[[2]] <- ggplot(df, aes(tier, n_cells)) +
  geom_violin() +
  geom_boxplot(width=0.1) +
  xlab("") +
  ylab("number of cells per gRNA") +
  th
ggarrange(plotlist = plots)
```

![](06_hit_calling_files/figure-html/low_tier-1.png)<!-- -->

Ignoring `low` tier hits, 88% of all targets have 5 or fewer DE genes in the vicinity of the perturbation site, and nearly half (33 out of 70, 47%) result in only one or two DE genes.

- The proportion of targets with only 1/2 hits is twice as high for intergenic (59%) compared to `TSS` perturbations (36%). This might be related to the smaller changes in expression induced by intergenic perturbations, which in turn lead to smaller downstream effects, compared to a `TSS` target that more closely resembles a KO.


```r
## number of hits per target
n_DE_per_tgt <- as.data.frame(table(mast.sig[mast.sig$tier!="low",]$target))
colnames(n_DE_per_tgt) <- c("target", "n_DE")
n_DE_per_tgt$class <- gRNA_ann[match(n_DE_per_tgt$target, gRNA_ann$target),'class']
n_DE_per_tgt$class <- factor(n_DE_per_tgt$class, levels=c("TSS", "LCR", "ENH", "INTRON", "INTERGENIC"))

summary(n_DE_per_tgt$n_DE)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   1.000   1.000   3.000   3.243   4.000  15.000
```

```r
ggplot(n_DE_per_tgt, aes(n_DE, fill=class)) +
  geom_bar() +
  scale_fill_manual(values = colors) +
  xlab("number of DE genes within 1Mb") +
  ylab("number of targets") +
  th
```

![](06_hit_calling_files/figure-html/n_DE_plot_no_low-1.png)<!-- -->


```r
plots <- list()

n_DE_per_tgt <- as.data.frame(table(mast.sig[mast.sig$tier=="high",]$target))
colnames(n_DE_per_tgt) <- c("target", "n_DE")
n_DE_per_tgt$class <- gRNA_ann[match(n_DE_per_tgt$target, gRNA_ann$target),'class']
n_DE_per_tgt$class <- factor(n_DE_per_tgt$class, levels=c("TSS", "LCR", "ENH", "INTRON", "INTERGENIC"))
plots[[1]] <- ggplot(n_DE_per_tgt, aes(n_DE, fill=class)) +
  geom_bar() +
  scale_fill_manual(values = colors) +
  xlim(c(0,13)) +
  ggtitle("high-confidence hits") +
  xlab("number of DE genes within 1Mb") +
  ylab("number of targets") +
  th

n_DE_per_tgt <- as.data.frame(table(mast.sig[mast.sig$tier=="medium",]$target))
colnames(n_DE_per_tgt) <- c("target", "n_DE")
n_DE_per_tgt$class <- gRNA_ann[match(n_DE_per_tgt$target, gRNA_ann$target),'class']
n_DE_per_tgt$class <- factor(n_DE_per_tgt$class, levels=c("TSS", "LCR", "ENH", "INTRON", "INTERGENIC"))
plots[[2]] <- ggplot(n_DE_per_tgt, aes(n_DE, fill=class)) +
  geom_bar() +
  scale_fill_manual(values = colors) +
  xlim(c(0,13)) +
  ggtitle("medium-confidence hits") +
  xlab("number of DE genes within 1Mb") +
  ylab("number of targets") +
  th
ggarrange(plotlist = plots, common.legend = TRUE, legend = "right")
```

![](06_hit_calling_files/figure-html/n_DE_plot_tier-1.png)<!-- -->

#### Expected genes for positive control perturbations

As discussed previously and consistent with results at the gRNA level, positive control perturbations result in differential expression of the expected gene most of the time:

- 94% for `TSS` targets.
- 100% for `LCR` targets.
- 73% for `ENH` targets.


```r
pos.ctrls <- mast.sig[!is.na(mast.sig$expected_gene),]

df <- as.data.frame(table(unique(pos.ctrls[,c('target','class')])$class))
colnames(df) <- c("class", "n_targets")
df$n_expected_DE <- as.data.frame(table(pos.ctrls[pos.ctrls$expected_gene==pos.ctrls$gene_id,]$class))$Freq
df
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["class"],"name":[1],"type":["fct"],"align":["left"]},{"label":["n_targets"],"name":[2],"type":["int"],"align":["right"]},{"label":["n_expected_DE"],"name":[3],"type":["int"],"align":["right"]}],"data":[{"1":"ENH","2":"26","3":"19"},{"1":"LCR","2":"3","3":"3"},{"1":"TSS","2":"34","3":"32"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

Only 4 of the targets have `low` confidence effects.


```r
table(pos.ctrls[pos.ctrls$expected_gene==pos.ctrls$gene_id,]$class,
      pos.ctrls[pos.ctrls$expected_gene==pos.ctrls$gene_id,]$tier)
```

```
##      
##       high medium low
##   ENH    9      8   2
##   LCR    2      1   0
##   TSS   28      2   2
```

And almost all the changes correspond to downregulation effects, as expected.


```r
table(pos.ctrls[pos.ctrls$expected_gene==pos.ctrls$gene_id,]$class,
      pos.ctrls[pos.ctrls$expected_gene==pos.ctrls$gene_id,]$direction)
```

```
##      
##       discordant down up
##   ENH          3   14  2
##   LCR          0    3  0
##   TSS          2   30  0
```


### Limitations on perturbation effect detection

Gene expression levels are one of the key limitations of single-cell technologies, and this greatly impacts our ability to detect perturbation effects. 

- We only tested genes that are detected in at least 5% of the cells in the dataset (~10K genes), corresponding to the top third from all genes detected in the experiment. 
- `TSS` controls were selected to span different expression levels. When tested, all but one result in significant differential expression suggesting strong effects can be detected effectively as long as the gene can be profiled by single-cell techniques.
- Instead, perturbation of intergenic `regulatory elements` result in significant differential expression in a more limited range of higher expression levels (top 6K genes). Enhancers regulating genes expressed at lower levels do not produce detectable expression changes, likely because the effects are much smaller. Increased power for these events will require larger numbers of cells and/or targeted transcriptomic approaches.
  + A few expected genes expressed at high levels are not detected as DE genes. Three out of four correspond to `ENH` elements and thus likely correspond to enhancers that are not functional in T cells. 


```r
## mean expr
means <- rowMeans(logcounts(sce))
mast.sig$mean <- means[mast.sig$gene_id]

## expression ranks
ranks <- means[order(means, decreasing = TRUE)]
ranks <- data.frame(mean=ranks,
                    rank=1:length(ranks))
ranks$detected <- row.names(ranks) %in% genes_detected

## annotate perturbation status
ranks$perturbed <- row.names(ranks) %in% gRNA_ann$expected_DE_gene
# flag which are detected
ranks$sig <- row.names(ranks) %in% mast.sig$gene_id
ranks$class <- gRNA_ann[match(row.names(ranks), gRNA_ann$expected_DE_gene),]$class

## plot all genes in a line
ranks$dummy <- 1
## separate significant genes
ranks[ranks$sig,]$dummy <- 2
## separate genes targeted
ranks[ranks$perturbed,]$dummy <- 6
## separate significant
ranks[which(ranks$sig & ranks$class != "TSS"),]$dummy <- 4
ranks[which(ranks$sig & ranks$class == "TSS"),]$dummy <- 5

## set colors
ranks$status <- ifelse(ranks$sig & ranks$perturbed, "targeted + DE",
                    ifelse(ranks$sig, "non-targeted + DE",
                           ifelse(ranks$perturbed, "targeted + non-DE", 
                                  ifelse(ranks$detected, "tested", "not-tested"))))
ranks$group <- "not-targeted"
ranks[which(ranks$class == "TSS"),]$group <- "TSS"
ranks[which(ranks$class != "TSS"),]$group <- "INTERGENIC"
ranks$group <- factor(ranks$group, levels=c('not-targeted', 'TSS', 'INTERGENIC'))

ggplot(ranks, aes(rank, dummy, shape=group, colour=status)) +
  geom_point(size=2) +
  scale_color_manual(values = c('tested'='grey50', 'not-tested'='grey80', 'non-targeted + DE'='indianred2', 
                                'targeted + DE'='indianred3', 'targeted + non-DE'='pink2')) +
  ylim(0,7) +
  ylab("") +
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         shape=guide_legend(nrow=3, byrow=TRUE)) +
  labs(colour="", shape = "") +
  th + theme(axis.text.x = element_blank(),
             axis.text.y = element_blank(),
             axis.ticks.y = element_blank(),
             legend.position = "bottom")
```

![](06_hit_calling_files/figure-html/mean_rank_expr-1.png)<!-- -->

- As expected, perturbations supported by higher number of independent gRNAs correspond to genes expressed at higher levels, which are easier to detect. Plot below shows only tested genes, split by the number of gRNAs with a raw gRNA p-value < 0.05.


```r
max_support <- mast.sig %>% 
  group_by(gene_id) %>% 
  summarise(n_grna = max(n_gRNAs))
ranks$n_grna <- max_support[match(row.names(ranks), max_support$gene_id),]$n_grna
ranks[is.na(ranks$n_grna),]$n_grna <- 0

df <- ranks[ranks$detected,]
df <- df[order(df$n_grna),]

df$dummy <- 2
df[which(df$class != "TSS"),]$dummy <- 4
df[which(df$class == "TSS"),]$dummy <- 5

plots <- list()
plots[[1]] <- ggplot(df, aes(rank, dummy, shape=group, colour=as.factor(n_grna))) +
  geom_point(size=2) +
  scale_color_manual(values = c('grey60',GetColors(n=11, scheme = "sunset")[c(6,8,10,11)]),
                     guide="none") +
  xlab("") +
  ylab("") +
  ylim(1,6) +
  facet_wrap(~n_grna, ncol=5) +
  labs(shape = "") +
  th + theme(axis.text.x = element_blank(),
             axis.text.y = element_blank(),
             axis.ticks.y = element_blank(),
             strip.background = element_rect(fill=NA, size=1),
             strip.text = element_text(size=10, face="bold"),
             panel.border = element_rect(fill=NA),
             legend.position = "top")
plots[[2]] <- ggplot(df, aes(rank, colour=as.factor(n_grna))) +
  geom_density() +
  scale_color_manual(values = c('grey60',GetColors(n=11, scheme = "sunset")[c(6,8,10,11)])) +
  facet_wrap(~n_grna, ncol=5) +
  ylab("") +
  th + theme(axis.text.x = element_blank(),
             axis.text.y = element_blank(),
             axis.ticks.y = element_blank(),
             strip.background = element_blank(),
             strip.text = element_blank(),
             panel.border = element_rect(fill=NA),
             legend.position = "none")
ggarrange(plotlist = plots, ncol=1, heights = c(0.6,0.4))
```

![](06_hit_calling_files/figure-html/mean_rank_expr_detected-1.png)<!-- -->

Similarly, the magnitude of the perturbation effect has a strong influence on our detection power. Although fold-change estimates are unreliable when the amount of zeroes is large, we clearly observe that perturbations with larger fold-changes are more easily detected and have support by more gRNAs.


```r
ggplot(mast.sig, aes(as.factor(n_gRNAs), abs(logFC_summary))) +
  geom_violin() +
  geom_boxplot(width=0.1) +
  xlab("number of independent gRNAs with raw p-value < 0.05") +
  ylab("absolute fold-change") +
  th
```

![](06_hit_calling_files/figure-html/fold_change_abs-1.png)<!-- -->

The majority of the expected genes from positive control perturbations are downregulated (negative fold-changes); and the effect sizes are much larger when targeting the `TSS`, compared to `intergenic` regulatory elements. 


```r
mast.sig$expected <- mast.sig$gene_id == mast.sig$expected_gene
ggplot(mast.sig[order(mast.sig$expected),], aes(as.factor(n_gRNAs), logFC_summary)) +
  geom_violin() +
  geom_boxplot(width=0.1) +
  geom_jitter(stat="identity", aes(colour=expected),
              position=position_jitter(0.1), size=1, alpha=0.5) +
  scale_color_manual(values = c("grey80", "indianred2")) +
  geom_hline(yintercept = 0, lwd=1, col="grey40", lty=2) +
  facet_wrap(~ifelse(class=="TSS", "TSS", "intergenic")) +
  xlab("number of gRNAs with raw p-value < 0.05") +
  ylab("log2 fold-change") +
  th + theme(strip.background = element_rect(fill=NA, size=1),
             strip.text = element_text(size=10, face="bold"),
             panel.border = element_rect(fill=NA))
```

![](06_hit_calling_files/figure-html/fold_change-1.png)<!-- -->

### Gene-target distance

We compute the distance for all significant E2G links. High-confidence pairs are generally close to the perturbation site (55% within 50kb; 46% within 25 kb). In contrast, low-confidence hits tend to span much longer distances, supporting some of these might be false positives.


```r
## annotate gene neighbourhood
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version=98)

## for each target, retrieve the genes in the neighbourhood
## rank them by distance
## annotate whether they are expressed
distance <- 1e6
gene_distance_target <- sapply(unique(mast.sig$target), function(t){
  ## get coordinates of target
  coords <- unique(gRNA_ann[gRNA_ann$target==t, c('chr','start','end')])
  # remmove 'chr' prefix to match ensembl chr names
  coords[1] <- substr(coords[1], 4,5)
  
  ## retreive nearby genes
  genes <- getBM(attributes = c('external_gene_name','chromosome_name','start_position','end_position','strand'),
                 filters = c('chromosome_name','start','end'),
                 values = list(coords$chr, 
                               coords$start-distance, 
                               coords$end+distance),
                 mart = ensembl)
  genes$target <- t
  
  ## distance from target midpoint to gene _start_ (as proxy for TSS)
  midpoint <- ceiling( coords$start+((coords$end-coords$start)/2))
  genes$distance <- ifelse(genes$strand==1, abs(genes$start_position-midpoint), abs(genes$end_position-midpoint))
                          
  genes <- genes[order(genes$distance),]
  genes$dist_rank <- 1:nrow(genes)
  genes$expressed <- genes_to_keep[genes$external_gene_name]  ## FALSE for those in less than 5%; NA for non-detected
  return(genes)
}, simplify = FALSE)
gene_distance_target <- do.call(rbind, gene_distance_target)
gene_distance_target$pair_tgt <- paste(gene_distance_target$target, gene_distance_target$external_gene_name, sep="|")
gene_distance_target$DE <- gene_distance_target$pair_tgt %in% mast.sig$pair_tgt
gene_distance_target$tier <- mast.sig[match(gene_distance_target$pair_tgt, mast.sig$pair_tgt),]$tier

## add distance to results
mast.sig$distance <- gene_distance_target[match(mast.sig$pair_tgt, gene_distance_target$pair_tgt),]$distance

ggplot(mast.sig, aes(distance/1e3, colour=tier)) +
  geom_density(size=1) +
  scale_color_manual(values = GetColors(n=11, scheme = "sunset")[c(10,8,6)]) +
  facet_wrap(~ifelse(class=="TSS", "TSS", "intergenic"), scales="free_y") +
  xlab("gene-target distance (kb)") +
  th
```

![](06_hit_calling_files/figure-html/distance-1.png)<!-- -->

For `intergenic` perturbations, we check how often the DE genes include the nearest **detected** gene (at least 5% of cells). 

- Over two thirds (69%) of all intergenic perturbations include the nearest detected gene.


```r
## annotate number of _expressed_ genes that are nearer to the target compared to DE gene
df <- gene_distance_target[gene_distance_target$DE,]
df$nearer_expr <- sapply(1:nrow(df), function(i) 
  sum(gene_distance_target[gene_distance_target$target == df$target[i] &
                             gene_distance_target$dist_rank < df$dist_rank[i],]$expressed, na.rm = TRUE))
df$nearer_expr_l <- df$nearer_expr == 0
df$class <- gRNA_ann[match(df$target, gRNA_ann$target),]$class

# add to results df
mast.sig$nearest_expr <- df[match(mast.sig$pair_tgt, df$pair_tgt),'nearer_expr_l']

tmp <- t(table(df[df$class != "TSS",]$target, df[df$class != "TSS",]$nearer_expr_l))
tmp <- tmp[,order(colSums(tmp))]

# plot
tmp <- as.data.frame(tmp)
tmp$class <- gRNA_ann[match(tmp$Var2, gRNA_ann$target),]$class
ggplot(tmp, aes(Var2, Freq, fill=Var1)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("grey", "indianred2")) +
  facet_grid(~class, scales = "free_x", space = "free_x") +
  labs(fill = "is nearest detected") +
  xlab("") +
  ylab("# DE genes") +
  th + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
```

![](06_hit_calling_files/figure-html/nearest_gene-1.png)<!-- -->

- Almost all DE genes that are the nearest detected are `high` or `medium` confidence hits.


```r
round(prop.table(table(df[df$nearer_expr_l,]$tier))*100, 2)
```

```
## 
##   high medium    low 
##  59.38  28.12  12.50
```

- From `high` confidence hits, 40% are the nearest gene, compared to only 12% of `medium` and 6% of `low` hits.


```r
round(prop.table(table(df$tier, is_nearest=df$nearer_expr_l), 1)*100, 2)
```

```
##         is_nearest
##          FALSE  TRUE
##   high   56.32 43.68
##   medium 87.14 12.86
##   low    94.70  5.30
```

Thus, robust detection of expression changes is more effective if the affected genes are close to the putative regulatory element, perhaps because shorter distances result in stronger effects that are easier to detect.


```r
ggplot(mast.sig[mast.sig$class != "TSS",], aes(distance/1e3, logFC_summary, colour=tier)) +
  geom_point() +
  scale_x_log10() +
  scale_color_manual(values = GetColors(n=11, scheme = "sunset")[c(10,8,6)]) +
  geom_smooth(method="lm") +
  facet_wrap(~tier) +
  xlab("DE gene - target distance (kb)") +
  ylab("log2 fold-change") +
  th + theme(legend.position = "none")
```

![](06_hit_calling_files/figure-html/dsitance_tier-1.png)<!-- -->

### Comparison to Gasperini et al.

The `ENH` perturbations profiled are hits from the Gasperini et al. (2009) study conducted in K562 cells. We use the gene linked by Gasperini et al. as the *expected gene*.

From 28 `ENH` targets, 26 result in at least one DE gene within 1 Mb. From these, we detect the expected gene in 19 cases (73.08%), with most being confident E2Gs. Thus, many of these regulatory interactions are active across diverse cell types.


```r
enh <- mast.sig[mast.sig$class=="ENH",]
table(enh[enh$expected,]$tier)
```

```
## 
##   high medium    low 
##      9      8      2
```

For the seven cases where the expected gene is not detected as DE, 


```r
## 'missed' expected genes
missing <- mast.tgt[mast.tgt$class == "ENH" &
                      mast.tgt$gene_id == mast.tgt$expected_gene &
                      mast.tgt$FDR > 0.05,]
plots <- list()
for(i in 1:nrow(missing)){
  target <- missing$target[i]
  gene <- missing$gene_id[i]
  pair <- paste(target, gene, sep="|")
  pvals <- mast[mast$pair_tgt==pair,]$p_value
  names(pvals) <- mast[mast$pair_tgt==pair,]$gRNA_id

  plots[[i]] <- plot_target_expr(sce = sce,
                 gene = gene,
                 guides = gRNA_ann[gRNA_ann$target==target,]$ID, 
                 target = target,
                 nt_cells = cells_neg,
                 per_guide = TRUE, 
                 fdr_thr = 0.05,
                 ann = pvals)
}
```

- The enhancer associated with regulation of *CD69* expression is a false negative. We observe downregulation with all four gRNAs, but given the small number of cells with each gRNA, the difference is not statistically significant. This hints to detrimental effects from silencing *CD69* expression.


```r
plots[[which(missing$gene_id=="CD69")]]
```

![](06_hit_calling_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

- Similarly, the perturbations involving *CERS2* and *CYSTM1* result in lower fraction of perturbed cells with detected expression: *CERS2* expression is detected in 44% of wild-type cells, versus ~30% of perturbed cells; *CYSTM1* expression is detected in 15% of wild-type cells, versus ~8% of perturbed cells. However, these genes are expressed at low levels and thus contain a large number of dropouts, making it challenging to confidently detect small downregulation effects.


```r
ggarrange(plotlist = plots[which(missing$gene_id %in% c("CERS2", "CYSTM1"))])
```

![](06_hit_calling_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

- In contrast, perturbations affecting *BTG1* and *TMSB4X*, which are both expressed at high levels, indicate slight gene upregulation in perturbed cells. However, the effect is too small to be detected confidently. These might represent differences in the regulatory wiring in T cells compared to K562 cells.


```r
## remove outliers from TMSB4X plot to be able to see changes more clearly
p <- plots[[which(missing$gene_id=="TMSB4X")]] + ylim(c(6,10))
ggarrange(plotlist = list(plots[[which(missing$gene_id=="BTG1")]], p))
```

![](06_hit_calling_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

- Finally, the perturbations against *PTPRC* and *SOCS3* do not show any robust effects. These might represent enhancers that are not active in T cells or include non-functional gRNAs.


```r
ggarrange(plotlist = plots[which(missing$gene_id %in% c("PTPRC", "SOCS3"))])
```

![](06_hit_calling_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

In sum, most perturbations profiled in Gasperini et al. result in the same changes in expression when profiled in T cells.

----

Several of the `ENH` perturbations result in additional confident changes in expression on top of the gene reported in Gasperini et al.


```r
additional <- enh[!enh$expected & enh$tier != "low",]
additional$tier <- droplevels(additional$tier)
table(additional$target, additional$tier)
```

```
##                     
##                      high medium
##   BBC3_ENH              1      2
##   CYSTM1_UBE2D2D_ENH    0      1
##   GIGYF1_ENH            1      4
##   HPCAL1_ENH            0      2
##   IER3_ENH              1      3
##   KCNN4_ENH             1      1
##   PHF19_ENH             0      2
##   RNF213_ENH            3      2
##   SOCS3_ENH             1      3
##   TKT_ENH               2      1
##   TUBB2A_ENH            1      0
##   YPEL3_ENH             3      7
```

The majority of these additional hits are located hundreds of kilobases away from the perturbation site, and might not be *direct* effects.


```r
summary(additional$distance)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   20208  126186  457640  449284  742707  997523
```

```r
ggplot(additional, aes(distance/1e3, colour=tier)) +
  geom_density(size=1) +
  scale_color_manual(values = GetColors(n=11, scheme = "sunset")[c(10,8)]) +
  xlab("gene-target distance (kb)") +
  th
```

![](06_hit_calling_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

For most cases, the additional hits are farther away than the *expected* gene reported by Gasperini et al. 


```r
matched <- enh[enh$expected,]
additional$dist_expected <- matched[match(additional$target, matched$target),'distance']

ggplot(additional, aes(dist_expected, distance, colour=tier)) +
  geom_point() +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks() +
  geom_abline(slope = 1, intercept = 0, lty=2) +
  scale_color_manual(values = GetColors(n=11, scheme = "sunset")[c(10,8)]) +
  ggtitle("gene-target distance (bp)") +
  xlab("expected gene") +
  ylab("additional hits") +
  th
```

![](06_hit_calling_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

And additional hits generally have smaller effects than *expected* hits. Also, many more additional DE genes show upregulation effects, again suggesting these might not be direct effects.


```r
ggplot(enh[enh$tier != "low",], aes(tier, logFC_summary, colour=expected)) +
  geom_violin() +
  geom_boxplot(width=0.1, position = position_dodge(0.9)) +
  geom_hline(yintercept = 0, lty=2) +
  scale_color_manual(values = c("grey", "indianred")) +
  xlab("") + 
  ylab("log2 fold-change") +
  th
```

![](06_hit_calling_files/figure-html/unnamed-chunk-10-1.png)<!-- -->


### cCRE perturbations

Beside positive control perturbations, the library also includes 14 targets that overlap ENCODE cCREs: 11 are within gene introns, and 3 are intergenic. 

All but two of the perturbations targeting gene introns result in at least one confident DE gene within 1 Mb, although the majority are `medium` rather than `high` tier.


```r
intron <- mast.sig[mast.sig$class == "INTRON" & mast.sig$tier != "low",]
intron$tier <- droplevels(intron$tier)

table(intron$target, intron$tier)
```

```
##               
##                high medium
##   ACAP1_INTRON    1      2
##   DVL1_INTRON1    1      3
##   DVL1_INTRON2    1      4
##   GLB1_INTRON     0      1
##   LCA5L_INTRON    0      1
##   LSP1_INTRON1    2      2
##   LSP1_INTRON2    0      1
##   LSP1_INTRON3    1      0
##   SMAD3_INTRON    0      3
```

From these 9, 5 include the overlapping gene.


```r
intron[intron$gene_id %in% c("ACAP1", "DVL1", "GLB1", "LCAS5L", "LSP1", "SMAD3"), c(1:2,6:9,12,16:17,19)]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["target"],"name":[1],"type":["chr"],"align":["left"]},{"label":["gene_id"],"name":[2],"type":["chr"],"align":["left"]},{"label":["pval_fisher"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["FDR"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["logFC_summary"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["n_gRNAs"],"name":[6],"type":["int"],"align":["right"]},{"label":["direction"],"name":[7],"type":["chr"],"align":["left"]},{"label":["tier"],"name":[8],"type":["fct"],"align":["left"]},{"label":["mean"],"name":[9],"type":["dbl"],"align":["right"]},{"label":["distance"],"name":[10],"type":["dbl"],"align":["right"]}],"data":[{"1":"DVL1_INTRON1","2":"DVL1","3":"0.0053235322","4":"0.021726693","5":"-0.08022955","6":"2","7":"down","8":"medium","9":"0.08658387","10":"11589","_rn_":"360","_row":"DVL1_INTRON1|DVL1"},{"1":"LSP1_INTRON1","2":"LSP1","3":"0.0036377277","4":"0.015830633","5":"-0.09335539","6":"2","7":"down","8":"medium","9":"3.00699605","10":"15091","_rn_":"719","_row":"LSP1_INTRON1|LSP1"},{"1":"SMAD3_INTRON","2":"SMAD3","3":"0.0015156154","4":"0.007555456","5":"-0.06803584","6":"2","7":"down","8":"medium","9":"0.55259002","10":"86517","_rn_":"1076","_row":"SMAD3_INTRON|SMAD3"},{"1":"ACAP1_INTRON","2":"ACAP1","3":"0.0000000000","4":"0.000000000","5":"-1.51868658","6":"4","7":"down","8":"high","9":"2.07964646","10":"553","_rn_":"1115","_row":"ACAP1_INTRON|ACAP1"},{"1":"LSP1_INTRON2","2":"LSP1","3":"0.0003886146","4":"0.002287176","5":"0.12497002","6":"2","7":"up","8":"medium","9":"3.00699605","10":"16741","_rn_":"3187","_row":"LSP1_INTRON2|LSP1"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

The most robust change is observed for *ACAP1*, but the target site is quite close to the gene TSS, and thus we are likely silencing the gene's promoter (making this a `TSS` perturbation rather than intergenic).


```r
target <- "ACAP1_INTRON"
gene <- "ACAP1"
pair <- paste(target, gene, sep="|")
pvals <- mast[mast$pair_tgt==pair,]$p_value
names(pvals) <- mast[mast$pair_tgt==pair,]$gRNA_id

plot_target_expr(sce = sce,
                 gene = gene,
                 guides = gRNA_ann[gRNA_ann$target==target,]$ID, 
                 target = target,
                 nt_cells = cells_neg,
                 per_guide = TRUE, 
                 fdr_thr = 0.05,
                 ann = pvals)
```

![](06_hit_calling_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

For the introns of *DVL1* and *LSP1* we have several targets, all for the same intron. For these, the `high` confidence hits are generally recapitulated across targets.


```r
tmp <- intron[grep("LSP", intron$target),c(1:2,6:9,12,16:17)]
tmp[order(tmp$tier,tmp$FDR),]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["target"],"name":[1],"type":["chr"],"align":["left"]},{"label":["gene_id"],"name":[2],"type":["chr"],"align":["left"]},{"label":["pval_fisher"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["FDR"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["logFC_summary"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["n_gRNAs"],"name":[6],"type":["int"],"align":["right"]},{"label":["direction"],"name":[7],"type":["chr"],"align":["left"]},{"label":["tier"],"name":[8],"type":["fct"],"align":["left"]},{"label":["mean"],"name":[9],"type":["dbl"],"align":["right"]}],"data":[{"1":"LSP1_INTRON1","2":"CTSD","3":"5.418625e-14","4":"9.281133e-13","5":"0.21618384","6":"3","7":"up","8":"high","9":"1.4031346","_rn_":"717","_row":"LSP1_INTRON1|CTSD"},{"1":"LSP1_INTRON1","2":"TSPAN32","3":"1.183271e-13","4":"1.927866e-12","5":"-0.08621026","6":"4","7":"discordant","8":"high","9":"0.4574127","_rn_":"723","_row":"LSP1_INTRON1|TSPAN32"},{"1":"LSP1_INTRON3","2":"TSPAN32","3":"3.571061e-09","4":"4.458820e-08","5":"-0.15763340","6":"3","7":"down","8":"high","9":"0.4574127","_rn_":"3143","_row":"LSP1_INTRON3|TSPAN32"},{"1":"LSP1_INTRON1","2":"C11orf21","3":"1.397114e-07","4":"1.414049e-06","5":"-0.06080552","6":"3","7":"discordant","8":"medium","9":"0.2794156","_rn_":"714","_row":"LSP1_INTRON1|C11orf21"},{"1":"LSP1_INTRON2","2":"LSP1","3":"3.886146e-04","4":"2.287176e-03","5":"0.12497002","6":"2","7":"up","8":"medium","9":"3.0069960","_rn_":"3187","_row":"LSP1_INTRON2|LSP1"},{"1":"LSP1_INTRON1","2":"LSP1","3":"3.637728e-03","4":"1.583063e-02","5":"-0.09335539","6":"2","7":"down","8":"medium","9":"3.0069960","_rn_":"719","_row":"LSP1_INTRON1|LSP1"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>


```r
tmp <- intron[grep("DVL", intron$target),c(1:2,6:9,12,16:17)]
tmp[order(tmp$tier,tmp$FDR),]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["target"],"name":[1],"type":["chr"],"align":["left"]},{"label":["gene_id"],"name":[2],"type":["chr"],"align":["left"]},{"label":["pval_fisher"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["FDR"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["logFC_summary"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["n_gRNAs"],"name":[6],"type":["int"],"align":["right"]},{"label":["direction"],"name":[7],"type":["chr"],"align":["left"]},{"label":["tier"],"name":[8],"type":["fct"],"align":["left"]},{"label":["mean"],"name":[9],"type":["dbl"],"align":["right"]}],"data":[{"1":"DVL1_INTRON2","2":"ISG15","3":"1.446252e-172","4":"1.288128e-170","5":"-0.27946142","6":"4","7":"down","8":"high","9":"1.22095171","_rn_":"501","_row":"DVL1_INTRON2|ISG15"},{"1":"DVL1_INTRON1","2":"ISG15","3":"2.048274e-12","4":"3.219405e-11","5":"-0.25968922","6":"3","7":"down","8":"high","9":"1.22095171","_rn_":"365","_row":"DVL1_INTRON1|ISG15"},{"1":"DVL1_INTRON2","2":"ACAP3","3":"1.746729e-29","4":"4.667260e-28","5":"-0.06909625","6":"2","7":"down","8":"medium","9":"0.15086657","_rn_":"485","_row":"DVL1_INTRON2|ACAP3"},{"1":"DVL1_INTRON2","2":"PUSL1","3":"2.833398e-08","4":"3.102803e-07","5":"-0.02533839","6":"2","7":"down","8":"medium","9":"0.15943581","_rn_":"510","_row":"DVL1_INTRON2|PUSL1"},{"1":"DVL1_INTRON1","2":"TNFRSF18","3":"8.849763e-08","4":"9.309672e-07","5":"0.16537821","6":"3","7":"discordant","8":"medium","9":"1.13176772","_rn_":"380","_row":"DVL1_INTRON1|TNFRSF18"},{"1":"DVL1_INTRON1","2":"TNFRSF4","3":"3.559605e-05","4":"2.489860e-04","5":"0.23728164","6":"2","7":"up","8":"medium","9":"0.63577300","_rn_":"381","_row":"DVL1_INTRON1|TNFRSF4"},{"1":"DVL1_INTRON2","2":"LINC01409","3":"3.610774e-05","4":"2.512497e-04","5":"-0.03526504","6":"2","7":"down","8":"medium","9":"0.19519604","_rn_":"503","_row":"DVL1_INTRON2|LINC01409"},{"1":"DVL1_INTRON2","2":"MRPL20-AS1","3":"2.135059e-03","4":"1.015103e-02","5":"-0.03177858","6":"2","7":"down","8":"medium","9":"0.19862639","_rn_":"506","_row":"DVL1_INTRON2|MRPL20-AS1"},{"1":"DVL1_INTRON1","2":"DVL1","3":"5.323532e-03","4":"2.172669e-02","5":"-0.08022955","6":"2","7":"down","8":"medium","9":"0.08658387","_rn_":"360","_row":"DVL1_INTRON1|DVL1"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

Often the affected genes are hundreds of kilobases away from the target site. 


```r
summary(intron$distance)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##     553   45462  336691  325767  592460  839317
```

```r
ggplot(intron, aes(distance/1e3, colour=tier)) +
  geom_density(size=1) +
  scale_color_manual(values = GetColors(n=11, scheme = "sunset")[c(10,8)]) +
  xlab("gene-target distance (kb)") +
  th
```

![](06_hit_calling_files/figure-html/unnamed-chunk-16-1.png)<!-- -->

For example, downregulation of *ISG15* after `DVL1_INTRON` silencing, is over 300kb away and we observe consistent effects from 7 different gRNAs. But we cannot rule out the effect is due to *DVL1* downregulation, rather than silencing of the regulatory element.


```r
plots <- list()

target <- "DVL1_INTRON1"
gene <- "ISG15"
pair <- paste(target, gene, sep="|")
pvals <- mast[mast$pair_tgt==pair,]$p_value
names(pvals) <- mast[mast$pair_tgt==pair,]$gRNA_id

plots[[1]] <- plot_target_expr(sce = sce,
                 gene = gene,
                 guides = gRNA_ann[gRNA_ann$target==target,]$ID, 
                 target = target,
                 nt_cells = cells_neg,
                 per_guide = TRUE, 
                 fdr_thr = 0.05,
                 ann = pvals)

target <- "DVL1_INTRON2"
pair <- paste(target, gene, sep="|")
pvals <- mast[mast$pair_tgt==pair,]$p_value
names(pvals) <- mast[mast$pair_tgt==pair,]$gRNA_id

plots[[2]] <- plot_target_expr(sce = sce,
                 gene = gene,
                 guides = gRNA_ann[gRNA_ann$target==target,]$ID, 
                 target = target,
                 nt_cells = cells_neg,
                 per_guide = TRUE, 
                 fdr_thr = 0.05,
                 ann = pvals)
ggarrange(plotlist = plots)
```

![](06_hit_calling_files/figure-html/unnamed-chunk-17-1.png)<!-- -->

---

Finally, all 3 intergenic perturbations result in significant changes in expression within 1Mb.


```r
intergenic <- mast.sig[mast.sig$class == "INTERGENIC" & mast.sig$tier != "low",]
intergenic$tier <- droplevels(intergenic$tier)

table(intergenic$target, intergenic$tier)
```

```
##                     
##                      high medium
##   CXCR5_INTERGENIC      0      5
##   SMARCE1_INTERGENIC    1      1
##   TGFB1_INTERGENIC      2      3
```

DE genes tend to be far away from the perturbation loci, with 75% lying beyond 100kb.


```r
summary(intergenic$distance)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   22925  122867  385020  396316  561066  975560
```

```r
ggplot(intergenic, aes(distance/1e3, colour=tier)) +
  geom_density(size=1) +
  scale_color_manual(values = GetColors(n=11, scheme = "sunset")[c(10,8)]) +
  xlab("gene-target distance (kb)") +
  th
```

![](06_hit_calling_files/figure-html/unnamed-chunk-19-1.png)<!-- -->

Thus, only one of the DE genes is the nearest detected with respect to the targeting site (*BCL9L* for `CXCR5_INTERGENIC`).


```r
intergenic[order(intergenic$target, intergenic$tier, intergenic$FDR),c(1:2,6:9,12,16:17,19,20)]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["target"],"name":[1],"type":["chr"],"align":["left"]},{"label":["gene_id"],"name":[2],"type":["chr"],"align":["left"]},{"label":["pval_fisher"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["FDR"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["logFC_summary"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["n_gRNAs"],"name":[6],"type":["int"],"align":["right"]},{"label":["direction"],"name":[7],"type":["chr"],"align":["left"]},{"label":["tier"],"name":[8],"type":["fct"],"align":["left"]},{"label":["mean"],"name":[9],"type":["dbl"],"align":["right"]},{"label":["distance"],"name":[10],"type":["dbl"],"align":["right"]},{"label":["nearest_expr"],"name":[11],"type":["lgl"],"align":["right"]}],"data":[{"1":"CXCR5_INTERGENIC","2":"CD3D","3":"4.380505e-13","4":"7.051030e-12","5":"-0.23579327","6":"2","7":"down","8":"medium","9":"3.86876779","10":"529809","11":"FALSE","_rn_":"807","_row":"CXCR5_INTERGENIC|CD3D"},{"1":"CXCR5_INTERGENIC","2":"RPS25","3":"2.397841e-05","4":"1.713110e-04","5":"-0.13477394","6":"2","7":"down","8":"medium","9":"4.83031221","10":"146138","11":"FALSE","_rn_":"824","_row":"CXCR5_INTERGENIC|RPS25"},{"1":"CXCR5_INTERGENIC","2":"ATP5MG","3":"1.445679e-03","4":"7.316012e-03","5":"-0.10805248","6":"2","7":"down","8":"medium","9":"3.84347494","10":"470947","11":"FALSE","_rn_":"802","_row":"CXCR5_INTERGENIC|ATP5MG"},{"1":"CXCR5_INTERGENIC","2":"BCL9L","3":"2.445533e-03","4":"1.122760e-02","5":"0.13037690","6":"2","7":"up","8":"medium","9":"0.52707429","10":"53055","11":"TRUE","_rn_":"803","_row":"CXCR5_INTERGENIC|BCL9L"},{"1":"CXCR5_INTERGENIC","2":"JAML","3":"1.206371e-02","4":"4.344238e-02","5":"-0.08830887","6":"2","7":"down","8":"medium","9":"1.43039213","10":"647459","11":"FALSE","_rn_":"818","_row":"CXCR5_INTERGENIC|JAML"},{"1":"SMARCE1_INTERGENIC","2":"TOP2A","3":"2.032177e-04","4":"1.274643e-03","5":"0.26662170","6":"3","7":"up","8":"high","9":"1.00476256","10":"190377","11":"FALSE","_rn_":"3283","_row":"SMARCE1_INTERGENIC|TOP2A"},{"1":"SMARCE1_INTERGENIC","2":"RARA","3":"1.294658e-03","4":"6.652553e-03","5":"0.06248292","6":"2","7":"up","8":"medium","9":"0.22416216","10":"299093","11":"FALSE","_rn_":"3279","_row":"SMARCE1_INTERGENIC|RARA"},{"1":"TGFB1_INTERGENIC","2":"TGFB1","3":"8.005460e-121","4":"5.347647e-119","5":"-0.36483688","6":"3","7":"down","8":"high","9":"1.93175317","10":"25936","11":"FALSE","_rn_":"3362","_row":"TGFB1_INTERGENIC|TGFB1"},{"1":"TGFB1_INTERGENIC","2":"RPS19","3":"3.365968e-06","4":"2.692775e-05","5":"-0.07937312","6":"3","7":"down","8":"high","9":"6.37925668","10":"532269","11":"FALSE","_rn_":"3357","_row":"TGFB1_INTERGENIC|RPS19"},{"1":"TGFB1_INTERGENIC","2":"PAFAH1B3","3":"9.563737e-07","4":"8.461691e-06","5":"0.09308721","6":"2","7":"up","8":"medium","9":"0.27834287","10":"975560","11":"FALSE","_rn_":"3352","_row":"TGFB1_INTERGENIC|PAFAH1B3"},{"1":"TGFB1_INTERGENIC","2":"BLVRB","3":"1.314875e-03","4":"6.730548e-03","5":"-0.06442205","6":"2","7":"down","8":"medium","9":"0.32272250","10":"862222","11":"FALSE","_rn_":"3338","_row":"TGFB1_INTERGENIC|BLVRB"},{"1":"TGFB1_INTERGENIC","2":"TMEM91","3":"3.173163e-03","4":"1.427389e-02","5":"-0.04371367","6":"2","7":"down","8":"medium","9":"0.09608673","10":"22925","11":"FALSE","_rn_":"3363","_row":"TGFB1_INTERGENIC|TMEM91"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

The strongest effect is observed on *TGFB1* expression, after silencing a locus 25kb away.


```r
target <- "TGFB1_INTERGENIC"
gene <- "TGFB1"
pair <- paste(target, gene, sep="|")
pvals <- mast[mast$pair_tgt==pair,]$p_value
names(pvals) <- mast[mast$pair_tgt==pair,]$gRNA_id

plot_target_expr(sce = sce,
                 gene = gene,
                 guides = gRNA_ann[gRNA_ann$target==target,]$ID, 
                 target = target,
                 nt_cells = cells_neg,
                 per_guide = TRUE, 
                 fdr_thr = 0.05,
                 ann = pvals)
```

![](06_hit_calling_files/figure-html/unnamed-chunk-21-1.png)<!-- -->



```r
write.table(mast.tgt, paste0(dir, "results/06_DEresults_targetLevel.tsv"), quote = FALSE, sep="\t", row.names = FALSE)
write.table(mast.sig, paste0(dir, "results/06_DEgenes.tsv"), quote = FALSE, sep="\t", row.names = FALSE)
write.table(gene_distance_target, paste0(dir, "results/06_gene_target_distance.tsv"), quote = FALSE, sep="\t", row.names = FALSE)
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
##  [1] biomaRt_2.48.3              inlmisc_0.5.2              
##  [3] RColorBrewer_1.1-2          ggpubr_0.4.0               
##  [5] scales_1.1.1                dplyr_1.0.7                
##  [7] metapod_1.0.0               scran_1.20.1               
##  [9] scater_1.20.1               ggplot2_3.3.5              
## [11] scuttle_1.2.1               SingleCellExperiment_1.14.1
## [13] SummarizedExperiment_1.22.0 Biobase_2.52.0             
## [15] GenomicRanges_1.44.0        GenomeInfoDb_1.28.1        
## [17] IRanges_2.26.0              S4Vectors_0.30.0           
## [19] BiocGenerics_0.38.0         MatrixGenerics_1.4.2       
## [21] matrixStats_0.60.0         
## 
## loaded via a namespace (and not attached):
##   [1] readxl_1.3.1              backports_1.2.1          
##   [3] BiocFileCache_2.0.0       igraph_1.2.6             
##   [5] splines_4.1.1             sp_1.4-5                 
##   [7] BiocParallel_1.26.1       digest_0.6.27            
##   [9] htmltools_0.5.1.1         viridis_0.6.1            
##  [11] fansi_0.5.0               magrittr_2.0.1           
##  [13] checkmate_2.0.0           memoise_2.0.0            
##  [15] ScaledMatrix_1.0.0        cluster_2.1.2            
##  [17] openxlsx_4.2.4            limma_3.48.3             
##  [19] Biostrings_2.60.2         winch_0.0.6              
##  [21] prettyunits_1.1.1         colorspace_2.0-2         
##  [23] blob_1.2.2                rappdirs_0.3.3           
##  [25] haven_2.4.3               xfun_0.25                
##  [27] rgdal_1.5-23              crayon_1.4.1             
##  [29] RCurl_1.98-1.4            jsonlite_1.7.2           
##  [31] glue_1.4.2                gtable_0.3.0             
##  [33] zlibbioc_1.38.0           XVector_0.32.0           
##  [35] DelayedArray_0.18.0       car_3.0-11               
##  [37] BiocSingular_1.8.1        abind_1.4-5              
##  [39] DBI_1.1.1                 edgeR_3.34.0             
##  [41] rstatix_0.7.0             Rcpp_1.0.7               
##  [43] viridisLite_0.4.0         progress_1.2.2           
##  [45] dqrng_0.3.0               foreign_0.8-81           
##  [47] bit_4.0.4                 rsvd_1.0.5               
##  [49] httr_1.4.2                ellipsis_0.3.2           
##  [51] farver_2.1.0              pkgconfig_2.0.3          
##  [53] XML_3.99-0.7              sass_0.4.0               
##  [55] dbplyr_2.1.1              locfit_1.5-9.4           
##  [57] utf8_1.2.2                labeling_0.4.2           
##  [59] tidyselect_1.1.1          rlang_0.4.11             
##  [61] AnnotationDbi_1.54.1      munsell_0.5.0            
##  [63] cellranger_1.1.0          tools_4.1.1              
##  [65] cachem_1.0.5              generics_0.1.0           
##  [67] RSQLite_2.2.7             broom_0.7.9              
##  [69] evaluate_0.14             stringr_1.4.0            
##  [71] fastmap_1.1.0             yaml_2.2.1               
##  [73] knitr_1.33                bit64_4.0.5              
##  [75] zip_2.2.0                 purrr_0.3.4              
##  [77] KEGGREST_1.32.0           nlme_3.1-152             
##  [79] sparseMatrixStats_1.4.2   xml2_1.3.2               
##  [81] compiler_4.1.1            beeswarm_0.4.0           
##  [83] filelock_1.0.2            curl_4.3.2               
##  [85] png_0.1-7                 ggsignif_0.6.2           
##  [87] tibble_3.1.3              statmod_1.4.36           
##  [89] bslib_0.2.5.1             stringi_1.7.3            
##  [91] highr_0.9                 forcats_0.5.1            
##  [93] lattice_0.20-44           bluster_1.2.1            
##  [95] Matrix_1.3-4              vctrs_0.3.8              
##  [97] pillar_1.6.2              lifecycle_1.0.0          
##  [99] jquerylib_0.1.4           BiocNeighbors_1.10.0     
## [101] cowplot_1.1.1             data.table_1.14.0        
## [103] bitops_1.0-7              irlba_2.3.3              
## [105] raster_3.4-13             R6_2.5.1                 
## [107] gridExtra_2.3             rio_0.5.27               
## [109] vipor_0.4.5               codetools_0.2-18         
## [111] assertthat_0.2.1          withr_2.4.2              
## [113] GenomeInfoDbData_1.2.6    mgcv_1.8-36              
## [115] hms_1.1.0                 grid_4.1.1               
## [117] beachmat_2.8.1            tidyr_1.1.3              
## [119] rmarkdown_2.10            DelayedMatrixStats_1.14.2
## [121] carData_3.0-4             ggbeeswarm_0.6.0
```

