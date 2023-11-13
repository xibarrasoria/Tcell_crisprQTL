##### Figures
library(scran)
library(scuttle)
library(scales)
library(ggplot2)
library(ggExtra)
library(ggpubr)
library(ggrepel)
library(RColorBrewer)
library(inlmisc)
library(UpSetR)
library(ComplexHeatmap)
library(dplyr)

dir <- paste0(getwd(), "/")
data.dir <- paste0(dir, "data/")

th <- theme_bw() + theme(
  axis.text.x = element_text(size=12), 
  axis.title.x = element_text(size=12), 
  axis.text.y = element_text(size=12), 
  axis.title.y = element_text(size=12),
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(), 
  axis.line = element_line(colour = "black"), 
  panel.border = element_blank(), 
  plot.title = element_text(face="bold", hjust = 0.5, size=12),
  plot.subtitle = element_text(hjust = 0.5))

######## FUNCTIONS
plot_gene_gRNA <- function(sce=NULL, gene=NULL, guides=NULL, nt_cells=NULL, target=NULL,
                           per_guide=TRUE, ann=NULL, fdr_thr=0.05){
  ## get cells with gRNAs in `guides`
  cells.guides <- sapply(guides, function(x) names(which(counts(altExp(sce, 'gRNA_calls'))[x,])), simplify = FALSE)
  cells.pos <- rep(guides, times=unlist(lapply(cells.guides, length)))
  names(cells.pos) <- unname(do.call(c, cells.guides))
  
  ## get cells with NT controls
  cells.neg <- nt_cells
  
  ## retrieve expression of `gene`
  df <- data.frame(cell = union(names(cells.pos), cells.neg))
  df$expr <- logcounts(sce)[gene, df$cell]
  df$status = factor(ifelse(df$cell %in% names(cells.pos), "perturbed", "NT"),
                     levels=c("perturbed", "NT")) ## if a cell has both a targeting and NT gRNA, consider perturbed
  df <- cbind(df,
              sapply(unique(cells.pos), function(x) df$cell %in% names(cells.pos[cells.pos==x])))
  ## for targets with 2 or more gRNAs, account for possibly having more than one target gRNA per cell
  if(ncol(df)>4){
    df$guide <- ifelse(rowSums(df[,-c(1:3)])>1 , "multiple",
                       ifelse(rowSums(df[,-c(1:3)])==0, "NT", "guide"))
    df <- df[df$guide != "multiple",]
    for(g in unique(cells.pos)){
      df[df[,g]==TRUE & df$guide == "guide",]$guide <- g
    }
    df$guide <- factor(df$guide, levels=c(unique(cells.pos), "NT"))
  }else{
    df$guide <- ifelse(df[,4]==TRUE, colnames(df)[4], "NT")
    df$guide <- factor(df$guide, levels=c(colnames(df)[4], "NT"))
  }
  
  ## plot
  if(per_guide==TRUE){
    p <- ggplot(df, aes(guide, expr)) +
      geom_violin(aes(colour = status), lwd=1) + 
      # geom_jitter(stat="identity",  aes(colour = status),
      #             position=position_jitter(0.1), 
      #             # size=1, alpha=0.05) +
      #             size=1, alpha=0.2) +
      scale_color_manual(values = c("lightcoral", "grey")) +
      geom_boxplot(width=0.1, colour="black", fill=NA) +
      geom_hline(yintercept = median(df[df$status=="NT",]$expr), 
                 colour="grey40", lty=2) +
      xlab("") + ylab(expression('log'[2]*' expression')) +
      ggtitle(label = gene, 
              subtitle = paste("guides for:", target, "\n")) +
      th + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
                 panel.grid.major.y = element_line(),
                 legend.position = "none")
  }else{
    p <- ggplot(df, aes(status, expr)) +
      geom_violin(aes(colour = status), lwd=1) + 
      # geom_jitter(stat="identity",  aes(colour = status),
      #             position=position_jitter(0.1), 
      #             # size=1, alpha=0.05) +
      #             size=1, alpha=0.2) +
      scale_color_manual(values = c("lightcoral", "grey")) +
      geom_boxplot(width=0.1, colour="black", fill=NA) +
      geom_hline(yintercept = median(df[df$status=="NT",]$expr), 
                 colour="grey40", lty=2) +
      xlab("") + ylab(expression('log'[2]*' expression')) +
      ggtitle(label = gene, 
              subtitle = paste("guides for:", target, "\n")) +
      th + theme(axis.text.x = element_text(angle=90, hjust=1),
                 panel.grid.major.y = element_line(),
                 legend.position = "none")
  }
  if(per_guide==TRUE){
    if(is.null(ann)){
      p <- p + annotate("text", 
                        x = 1:length(levels(as.factor(df$guide))),
                        y = max(df$expr)+((max(df$expr)-min(df$expr))*0.1), 
                        label = paste0(
                          paste0("(", as.numeric(round(prop.table(table(df$guide, df$expr>0),1)*100,1)[,'TRUE']), ")"), "\n",
                          as.numeric(table(df$guide))),
                        size = 2, colour = "grey10")
    }else{
      p <- p + annotate("text", 
                        x = 1:length(levels(as.factor(df$guide))),
                        y = max(df$expr)+((max(df$expr)-min(df$expr))*0.2), 
                        label = paste0(
                          c(formatC(ann[setdiff(levels(as.factor(df$guide)),"NT")], 
                                    format = "e", digits = 1), ""), "\n",
                          paste0("(", as.numeric(round(prop.table(table(df$guide, df$expr>0),1)*100,1)[,'TRUE']), ")"), "\n",
                          as.numeric(table(df$guide))),
                        size = 2, 
                        colour = c(ifelse(ann[setdiff(levels(as.factor(df$guide)),"NT")] < fdr_thr, 
                                          "red", "grey10"), "grey10")) +
        coord_cartesian(clip = "off")
    }
  }else{
    if(is.null(ann)){
      p <- p + annotate("text", 
                        x = 1:length(levels(df$status)),
                        y = max(df$expr)+((max(df$expr)-min(df$expr))*0.1), 
                        label = paste0(
                          paste0("(", as.numeric(round(prop.table(table(df$status, df$expr>0),1)*100,1)[,'TRUE']), ")"), "\n",
                          as.numeric(table(df$status))),
                        size = 2, colour = "grey10") 
    }else{
      p <- p + annotate("text", 
                        x = 1:length(levels(df$status)),
                        y = max(df$expr)+((max(df$expr)-min(df$expr))*0.2), 
                        label = paste0(
                          c(formatC(ann, format = "e", digits = 1),
                            "NA"), 
                          "\n",
                          paste0("(", as.numeric(round(prop.table(table(df$status, df$expr>0),1)*100,1)[,'TRUE']), ")"), "\n",
                          as.numeric(table(df$status))),
                        size = 2, 
                        colour = c(ifelse(ann < fdr_thr, "red", "grey10"), "grey10")) +
        coord_cartesian(clip = "off")
    }
  }
  return(p)
}

plot_gene_target <- function(sce=NULL, nt_cells=NULL, gRNA_ann=NULL, hit=NULL, per_guide=TRUE){
  ## get relevant data
  gene <- hit$gene_id
  target <- hit$target
  guides <- unlist(strsplit(hit$gRNAs, ",", fixed=TRUE))
  pvals <- as.numeric(unlist(strsplit(hit$p_values, ",", fixed=TRUE)))
  names(pvals) <- guides
  fcs <- as.numeric(unlist(strsplit(hit$logFCs, ",", fixed=TRUE)))
  names(fcs) <- guides
  
  ## get cells with gRNAs in `guides`
  cells.guides <- sapply(guides, function(x) names(which(counts(altExp(sce, 'gRNA_calls'))[x,])), simplify = FALSE)
  cells.pos <- rep(guides, times=unlist(lapply(cells.guides, length)))
  names(cells.pos) <- unname(do.call(c, cells.guides))
  
  ## get cells with NT controls
  cells.neg <- nt_cells
  
  ## retrieve expression of `gene`
  df <- data.frame(cell = union(names(cells.pos), cells.neg))
  df$expr <- logcounts(sce[,df$cell])[gene,]
  # label perturbed vs NT
  df$status = factor(ifelse(df$cell %in% names(cells.pos), "perturbed", "NT"), levels=c("perturbed", "NT"))
  # annotate each gRNA
  df <- cbind(df, sapply(unique(cells.pos), function(x) df$cell %in% names(cells.pos[cells.pos==x])))
  # remove cells with several on-target gRNAs (to keep plots consistent and avoid having some with more violins than others)
  if(ncol(df)>4){
    df$guide <- ifelse(rowSums(df[,-c(1:3)])>1 , "multiple",
                       ifelse(rowSums(df[,-c(1:3)])==0, "NT", "guide"))
    df <- df[df$guide != "multiple",]
    for(g in unique(cells.pos)){
      df[df[,g]==TRUE & df$guide == "guide",]$guide <- g
    }
    df$guide <- factor(df$guide, levels=c(unique(cells.pos), "NT"))
  }else{
    df$guide <- ifelse(df[,4]==TRUE, colnames(df)[4], "NT")
    df$guide <- factor(df$guide, levels=c(colnames(df)[4], "NT"))
  }
  
  ## colour based on significance per-guide
  df$colour <- ifelse(df$status == "NT", "NT", "perturbed")
  
  ## plot expression levels
  if(per_guide){
    p1 <- ggplot(df, aes(guide, expr, colour = colour)) +
      geom_violin(lwd=1) +
      # geom_jitter(stat="identity",
      #             position=position_jitter(0.1), 
      #             size=1, alpha=0.2) +
      scale_color_manual(values = c(perturbed="lightcoral", NT="grey60")) +
      geom_boxplot(outlier.shape = NA, colour="grey20", width=0.1, fill=NA) +
      geom_hline(yintercept = median(df[df$status=="NT",]$expr), colour="grey20", lty=2) +
      xlab("") + 
      ylab(expression('log'[2]*' normalised counts')) +
      ggtitle(label = paste0(hit$target, " -> ", hit$gene_id), 
              subtitle = paste("FDR:", formatC(hit$FDR, format = "e", digits = 3), 
                               "\t\ttier:", toupper(hit$tier), "\t\tdistance:", hit$distance, "bp\n")) +
      th + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
                 legend.position = "none")
  
    ## annotate
    p1 <- p1 + annotate("text", 
                      x = 1:length(levels(df$guide)),
                      y = max(df$expr)+((max(df$expr)-min(df$expr))*0.3), 
                      label = paste0(
                        c(formatC(pvals[setdiff(levels(df$guide),"NT")], format = "e", digits = 2), ""), "\n",
                        c(round(fcs[setdiff(levels(df$guide),"NT")], 2), ""), "\n",
                        paste0("(", as.numeric(round(prop.table(table(df$guide, df$expr>0),1)*100,1)[,'TRUE']), ")"), "\n",
                        c(as.numeric(table(cells.pos)[setdiff(levels(df$guide),"NT")]), length(cells.neg))),
                      size = 2) +
      coord_cartesian(clip = "off")
  }else{
    p1 <- ggplot(df, aes(status, expr, colour = colour)) +
      geom_violin(lwd=1) +
      # geom_jitter(stat="identity",
      #             position=position_jitter(0.1), 
      #             size=1, alpha=0.2) +
      scale_color_manual(values = c(perturbed="lightcoral", NT="grey60")) +
      geom_boxplot(outlier.shape = NA, colour="grey20", width=0.1, fill=NA) +
      geom_hline(yintercept = median(df[df$status=="NT",]$expr), colour="grey20", lty=2) +
      xlab("") + 
      ylab(expression('log'[2]*' normalised counts')) +
      ggtitle(label = paste0(hit$target, " -> ", hit$gene_id), 
              subtitle = paste("FDR:", formatC(hit$FDR, format = "e", digits = 3), 
                               "\t\ttier:", toupper(hit$tier), "\t\tdistance:", hit$distance, "bp\n")) +
      th + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
                 legend.position = "none")
    
    ## annotate
    p1 <- p1 + annotate("text", 
                        x = 1:2,
                        y = max(df$expr)+((max(df$expr)-min(df$expr))*0.3), 
                        label = paste0(
                          paste0("(", as.numeric(round(prop.table(table(df$status, df$expr>0),1)*100,1)[,'TRUE']), ")"), "\n",
                          c(length(cells.pos), length(cells.neg))),
                        size = 2) +
      coord_cartesian(clip = "off")
  }
  
  ## plot proportion of zeroes
  df2 <- as.data.frame(round(prop.table(table(df$guide, df$expr>0),1)*100,1)[,'TRUE'])
  colnames(df2) <- "prop"
  df2$guide <- row.names(df2)
  df2$colour <- ifelse(df2$guide == "NT", "NT", "perturbed")
  
  p2 <- ggplot(df2, aes(guide, prop, colour=colour)) +
    geom_bar(stat="identity", fill=NA, lwd=1, width = 0.5) +
    scale_colour_manual(values = c(perturbed="lightcoral", NT="grey60")) +
    xlab("") +
    ylab("% cells with detected expression") +
    th + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
               legend.position = "none")
  ## add cell number
  p2 <- p2 + annotate("text", 
                      x = 1:length(levels(df$guide)),
                      y = max(df2$prop)+((max(df2$prop)-min(df2$prop))*0.2), 
                      label = c(as.numeric(table(cells.pos)[setdiff(levels(df$guide),"NT")]), length(cells.neg)),
                      size = 2)
  
  ## combine plots
  p <- ggarrange(plotlist = list(p1, p2), align="h", widths = c(0.6,0.4))
  return(p)
}

########

### Figure 1 ###############

## 1D ======
sce <- readRDS(paste0(data.dir, "data/flow_comparison_Tcells/sce_flow_comparison_cellecta.NORM.calls.Rds"))
sce <- sce[,-which(sce$call=="multiple")]
cells_nt <- names(which(counts(altExp(sce, 'gRNA_calls'))['gRNA_10',]))

plots <- list()
for(gene in c("CD4", "CD81", "BST2", "CD298")){
  guides <- rowData(altExp(sce))[grep(gene, rowData(altExp(sce))$Symbol),]$ID
  target <- gene
  if(gene == "CD298"){ gene <- 'ATP1B3' } # use Ensembl gene name
  plots[[gene]] <- plot_gene_gRNA(sce = sce,
                                  gene = gene, 
                                  guides = guides, 
                                  nt_cells = cells_nt, 
                                  target = target)
}
pdf(paste0(dir, "figures/Figure1D.pdf"), width = 12, height = 4.5, useDingbats = FALSE)
ggarrange(plotlist = plots, ncol=4)
dev.off()


### Figure S1 ###############

## S1E ======
## test for downregulation with a wilcoxon rank sum test
test_target <- function(sce=NULL, gene=NULL, guides=NULL, nt_cells=NULL, per_guide=TRUE, test="wilcox"){
  ## ensure NT controls have been provided
  if(is.null(nt_cells)){
    stop("Please provide IDs of cells with NT controls using the `nt_cell` argument.")
  }
  ## get cells with gRNAs in `guides`
  cells.guides <- sapply(guides, function(x) names(which(counts(altExp(sce, 'gRNA_calls'))[x,])),
                         simplify = FALSE)
  cells.pos <- rep(guides, times=unlist(lapply(cells.guides, length)))
  names(cells.pos) <- unname(do.call(c, cells.guides))
  
  ## get cells with NT controls
  cells.neg <- nt_cells
  
  ## check that the gene is expressed
  if(length(intersect(gene, row.names(sce)))==0){
    stop(paste("The gene", gene, "is not in the expression matrix."))
  }
  
  ## retrieve expression of `gene`
  df <- data.frame(cell = union(names(cells.pos), cells.neg))
  df$expr <- logcounts(sce)[gene, df$cell]
  df$status = factor(ifelse(df$cell %in% names(cells.pos), "perturbed", "NT"),
                     levels=c("perturbed", "NT")) ## if a cell has both a targeting and NT gRNA, consider perturbed
  df <- cbind(df,
              sapply(unique(cells.pos), function(x) df$cell %in% names(cells.pos[cells.pos==x])))
  ## for targets with 2 or more gRNAs, account for possibly having more than one target gRNA per cell
  if(ncol(df)>4){
    df$guide <- ifelse(rowSums(df[,-c(1:3)])>1 , "multiple",
                       ifelse(rowSums(df[,-c(1:3)])==0, "NT", "guide"))
    for(g in unique(cells.pos)){
      df[df[,g]==TRUE & df$guide == "guide",]$guide <- g
    }
  }else{
    df$guide <- ifelse(df[,4]==TRUE, colnames(df)[4], "NT")
  }
  
  if(per_guide==TRUE){
    if(test == "wilcox"){
      pvals <- sapply(setdiff(unique(df$guide), "NT"), function(x){
        wilcox.test(x = df[df$guide==x,]$expr,
                    y = df[df$guide=="NT",]$expr,
                    alternative = "less")$p.value
      })
    }else{
      pvals <- sapply(guides, function(x){
        ks.test(x = df[df$guide==x,]$expr,
                y = df[df$guide=="NT",]$expr,
                alternative = "greater")$p.value
      })
    }
  }else{
    if(test == "wilcox"){
      pvals <- wilcox.test(x = df[df$status=="perturbed",]$expr,
                           y = df[df$status=="NT",]$expr,
                           alternative = "less")$p.value
    }else{
      pvals <- ks.test(x = df[df$status=="perturbed",]$expr,
                       y = df[df$status=="NT",]$expr,
                       alternative = "greater")$p.value
    }
  }
  return(pvals)
}

# data from flow cytometry
flow <- c(65.2,71.3,66,28.4,41.9,18.9,13.2)
names(flow) <- c("CD4_1","CD81_1","CD81_2","BST2_1","BST2_2","CD298_1","CD298_2")

# gRNA annotation
gRNA_ann <- read.table(paste0(data.dir, "data/flow_comparison_Tcells/gRNA_annotation.tsv"), header = TRUE)
gRNA_ann$flow <- flow[match(gRNA_ann$Symbol, names(flow))]

# test 
target_pvals <- sapply(setdiff(unique(gRNA_ann$target), "NT"), function(target){
  test_target(sce = sce, 
              gene = unique(gRNA_ann[gRNA_ann$target==target,]$gene),
              guides = gRNA_ann[gRNA_ann$target==target,]$ID,
              nt_cells = cells_nt,
              test = "wilcox",
              per_guide = TRUE)
})
target_pvals <- as.data.frame(unlist(target_pvals))
colnames(target_pvals) <- "p_value"
target_pvals$target <- unlist(lapply(strsplit(row.names(target_pvals), ".", fixed=TRUE), '[[', 1))
target_pvals$guide <- unlist(lapply(strsplit(row.names(target_pvals), ".", fixed=TRUE), '[[', 2))
target_pvals <- target_pvals[,c(2,3,1)]
target_pvals$p_adj <- p.adjust(target_pvals$p_value)

## plot
target_pvals$flow <- gRNA_ann[match(target_pvals$guide, gRNA_ann$ID),]$flow

lm.flow <- lm(-log10(target_pvals$p_adj) ~ target_pvals$flow)
# summary(lm.flow)
# Multiple R-squared:  0.9831,	Adjusted R-squared:  0.9789 
# F-statistic: 233.1 on 1 and 4 DF,  p-value: 0.0001073

pdf(paste0(dir, "figures/FigureS1E.pdf"), width = 5, height = 3.5, useDingbats = FALSE)
ggplot(target_pvals, aes(flow, -log10(p_adj), label=row.names(target_pvals))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), lty=2, colour="grey") +
  geom_smooth(method=lm, se=FALSE, colour="lightblue") +
  geom_text_repel(size=3, colour="grey20") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.25, vjust = 1, 
           label = expression('r'^2*' = 0.9831'), size=3) +
  xlab("% knockdown cells by flow") +
  ylab(expression('-log'[10]*' FDR from CROP-seq')) +
  th
dev.off()



### Figure 2 ###############

## 2C ======
# data
sce <- readRDS(paste0(dir, "results/03_sce_POC_Tcells.goodQual.calls.NORM.Rds"))
# QC stats contains _all_ data; sce is only assigned cells
qc_metrics <- read.table(paste0(dir, "results/01_QCmetrics_POC_Tcells.tsv"))
# add calls to metadata - `sce` no longer contains unassigned
df <- as.data.frame(colData(sce))
qc_metrics$call <- df[match(row.names(qc_metrics), row.names(df)),]$call
qc_metrics[is.na(qc_metrics$call),]$call <- "unassigned"

# only keep good-quality
df <- qc_metrics[qc_metrics$passQC & qc_metrics$scDblFinder.class=="singlet" & !qc_metrics$nuclei,]

pdf(paste0(dir, "figures/Figure2C.pdf"), width = 5, height = 5, useDingbats = FALSE)
pie(prop.table(table(df$call))*100, col=c("indianred", "grey90", "grey40"),
    init.angle = -8, clockwise = TRUE, labels = "")
dev.off()

## 2D ======
ncells_per_gRNA <- rowSums(counts(altExp(sce, 'gRNA_calls')))

pdf(paste0(dir, "figures/Figure2D.pdf"), width = 3, height = 4, useDingbats = FALSE)
ggplot(as.data.frame(ncells_per_gRNA), aes(1, ncells_per_gRNA)) +
  scale_y_log10() +
  geom_violin() +
  geom_boxplot(width=0.1) +
  stat_summary(geom = "text", fun.y = quantile,
               aes(label=sprintf("%1.0f", 10^..y..)), 
               colour = "grey60",
               position = position_nudge(x=0.15), size=3.5) +
  annotation_logticks(side="l") +
  xlab("") +
  ylab("number of cells per gRNA") +
  th + theme(axis.text.x = element_blank(),
             axis.ticks.x = element_blank())
dev.off()


## 2E ======
gRNA_ann <- read.table(paste0(dir, "data/gRNA_library.tsv"), sep="\t", header = TRUE)
df <- data.frame(gRNA = names(ncells_per_gRNA),
                 n_cells = ncells_per_gRNA,
                 target = gRNA_ann[match(names(ncells_per_gRNA), gRNA_ann$ID),]$target)
df <- df %>% group_by(target) %>%
  summarise(ncells_per_target = sum(n_cells))

pdf(paste0(dir, "figures/Figure2E.pdf"), width = 3, height = 4, useDingbats = FALSE)
ggplot(df, aes(1, ncells_per_target)) +
  scale_y_log10() +
  geom_violin() +
  geom_boxplot(width=0.1) +
  stat_summary(geom = "text", fun = quantile,
               aes(label=sprintf("%1.0f", 10^..y..)), 
               colour = "grey60",
               position = position_nudge(x=0.15), size=3.5) +
  annotation_logticks(side="l") +
  xlab("") +
  ylab("number of cells per target") +
  th + theme(axis.text.x = element_blank(),
             axis.ticks.x = element_blank())
dev.off()



### Figure S2 ###############

## S2A ======
# define outliers and cut-offs
outliers.libSize <- isOutlier(qc_metrics$sum, type = "lower", log = TRUE, nmads = 3, batch = qc_metrics$sample)
outliers.ngenes <- isOutlier(qc_metrics$detected, type = "lower", log = TRUE, nmads = 3, batch = qc_metrics$sample)
outliers.mt <- isOutlier(qc_metrics$subsets_Mt_percent, type = "higher", log = FALSE, nmads = 3, batch = qc_metrics$sample)
outliers <- outliers.libSize | outliers.ngenes | outliers.mt
thresholds <- data.frame(libSize = attr(outliers.libSize, "thresholds")[1,],
                         nGenes = attr(outliers.ngenes, "thresholds")[1,],
                         mt = attr(outliers.mt, "thresholds")[2,])
thresholds$start <- 1:nrow(thresholds)-0.4
thresholds$end <- 1:nrow(thresholds)+0.4

# plot
qc_metrics$sample <- factor(qc_metrics$sample, levels=paste0("Tcell_CRISPRi_", c(1:4,6:8,10:40)))
plots <- list()
plots[['libSize']] <- ggplot(qc_metrics, aes(sample, sum)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(side="l") +
  geom_violin(aes(colour = colour)) + 
  # geom_boxplot(aes(colour = colour), width=0.05) + 
  scale_color_manual(values = as.character(GetColors(n=5, scheme = "pale"))) +
  geom_segment(data = thresholds, 
               aes(x=start, xend=end, y=libSize, yend=libSize),
               colour="grey20", lwd=0.5, lty=2) +
  ggtitle("library size") +
  xlab("") + 
  ylab("total UMIs") + 
  th +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
plots[['nGenes']] <- ggplot(qc_metrics, aes(sample, detected)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(side="l") +
  geom_violin(aes(colour = colour)) + 
  # geom_boxplot(aes(colour = colour), width=0.05) + 
  scale_color_manual(values = as.character(GetColors(n=5, scheme = "pale"))) +
  geom_segment(data = thresholds, 
               aes(x=start, xend=end, y=nGenes, yend=nGenes),
               colour="grey20", lwd=0.5, lty=2) +
  ggtitle("total genes per cell") +
  xlab("") + 
  ylab("number genes detected") +
  th + 
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
plots[['mt']] <- ggplot(qc_metrics, aes(sample, subsets_Mt_percent)) + 
  geom_violin(aes(colour = colour)) + 
  # geom_boxplot(aes(colour = colour), width=0.05) + 
  scale_color_manual(values = as.character(GetColors(n=5, scheme = "pale"))) +
  geom_segment(data = thresholds, 
               aes(x=start, xend=end, y=mt, yend=mt),
               colour="grey20", lwd=0.5, lty=2) +
  ggtitle("% reads in mitochondrial genes") +
  xlab("") + 
  ylab(expression('% reads')) +
  th + 
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

pdf(paste0(dir, "figures/FigureS2A.pdf"), width = 10, height = 7, useDingbats = FALSE)
ggarrange(plotlist = plots, ncol=1, nrow=length(plots), align = "v")
dev.off()

## S2B ======
# only keep good-quality
df <- as.data.frame(colData(sce))
df <- qc_metrics[qc_metrics$passQC & qc_metrics$scDblFinder.class=="singlet" & !qc_metrics$nuclei,]
df <- round(prop.table(table(df$sample, df$call),1)*100, 2)

pdf(paste0(dir, "figures/FigureS2B.pdf"), width = 7, height = 4, useDingbats = FALSE)
par(mar=c(6,4,2,2))
barplot(t(df[,c(3,1,2)]), col=c("grey40", "indianred", "grey90"), 
        ylab = "% of cells", las=2,
        names.arg = rep("", nrow(df)),
        width = 0.8, space = 0.2, xlim=c(0,42))
legend("topright", legend = c("unassigned", "multiple", "unique"),
       col = c("grey90", "indianred", "grey40"), pch = 15, cex=0.7)
dev.off()

## S2C ======
df <- data.frame(top_gRNA = apply(counts(altExp(sce, 'CRISPR')), 2, which.max),
                 top_cDNA = apply(counts(altExp(sce, 'gRNA_in_cDNA')), 2, which.max))
df$total_gRNA <- colSums(counts(altExp(sce, 'CRISPR')))
df$total_cDNA <- colSums(counts(altExp(sce, 'gRNA_in_cDNA')))
df$class <- ifelse(df$top_cDNA==1 & df$total_cDNA==0, "not_detected_cDNA",
                   ifelse(df$top_gRNA==1 & df$top_gRNA==0, "not_detected_gRNA",
                          ifelse(df$top_cDNA == df$top_gRNA, "concordant", "discordant")))
df[df$class=="not_detected_cDNA",]$top_cDNA <- -10

# only plot points that are in distinct locations
use <- iSEE::subsetPointsByGrid(df$top_gRNA, df$top_cDNA, resolution = 500)

plots <- list()
plots[[1]] <- ggplot(df[use,], aes(top_gRNA, top_cDNA, colour=class)) +
  geom_point(size=1, alpha = 0.5) +
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
  xlab("") + 
  ylab("number of cells") +
  th + theme(axis.text.x = element_blank(),
             axis.ticks.x = element_blank())

pdf(paste0(dir, "figures/FigureS2C.pdf"), width = 6, height = 5, useDingbats = FALSE)
ggarrange(plotlist = plots, ncol=2, widths = c(0.75, 0.25), 
          common.legend = TRUE, legend="bottom", align = "h")
dev.off()

## S2D ======
lib <- read.table(paste0(dir, "data/gRNA_count_plasmid_library.tsv"), header = TRUE)
ncells_per_gRNA <- rowSums(counts(altExp(sce, 'gRNA_calls')))
lib$n_cells <- ncells_per_gRNA[lib$gRNA_ID]

p <- ggplot(lib, aes(Count, n_cells)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  annotation_logticks() +
  geom_smooth(method="lm") +
  xlab("count in plasmid library") +
  ylab("count in recovered cells") +
  th

pdf(paste0(dir, "figures/FigureS2D.pdf"), width = 6, height = 4, useDingbats = FALSE)
ggMarginal(p, type="density")
dev.off()



### Figure 3 ###############

## 3A ======
mast <- read.table(paste0(dir, "results/05_DE_results.mast.tsv"), header = TRUE)
mast.exp <- mast[!is.na(mast$expected_gene),] # separate targets with expected effect
mast.tp <- mast.exp[mast.exp$gene_id == mast.exp$expected_gene,]

pvals <- matrix(-log10(mast.tp$FDR), ncol=4, nrow=length(unique(mast.tp$target)), byrow = TRUE)
colnames(pvals) <- paste0("gRNA", 1:4)
row.names(pvals) <- unlist(lapply(strsplit(unique(mast.tp$target), "_"), '[[', 1))
# some p-values are 0; assign smallest non-zero value to avoid Inf
pvals[pvals==Inf] <- -log10(min(mast.tp[mast.tp$FDR > 0,]$FDR))

# annotate by perturbation type
cols <- as.character(GetColors(n=5, scheme = "muted")[c(1,3,5)])
names(cols) <- unique(mast.tp$class)
# and count number of significnat gRNAs
# consider that if all the perturbed cells have 0 count a fold-change is not estimated, but it should be considered as a negative (downregulation) effect
ha  <- rowAnnotation(class = mast.tp[seq(1,nrow(mast.tp),4),]$class,
                     n_gRNAs = anno_barplot(sapply(seq(1,nrow(mast.tp),4), 
                                                function(x) sum(mast.tp[x:(x+3),]$sig & 
                                                                  (mast.tp[x:(x+3),]$logFC < 0 | is.na(mast.tp[x:(x+3),]$logFC))))),
                     col = list(class = cols),
                     show_legend = FALSE,
                     show_annotation_name = FALSE)

h <- Heatmap(pvals,
             col = circlize::colorRamp2(breaks = c(0,-log10(0.05),5,10,15,20),
                                        colors = c("white", brewer.pal(n=9, "BuPu")[3:7])),
             cluster_columns = FALSE, 
             right_annotation = ha,
             split = mast.tp[seq(1,nrow(mast.tp),4),]$class,
             show_column_names = FALSE,
             row_names_gp = gpar(fontsize = 6),
             show_row_dend = FALSE,
             heatmap_legend_param = list(title = "-log10(FDR)",
                                         direction = "horizontal"))
pdf(paste0(dir, "figures/Figure3A.pdf"), width = 2.5, height = 7, useDingbats = FALSE)
draw(h, heatmap_legend_side = "bottom")
dev.off()


## 3B ======
df <- mast.tp[mast.tp$FDR < 0.05,c('class', 'logFC')]

plots <- list()
plots[[1]] <- ggplot(df[df$class=="TSS",], aes(class, logFC, colour=class)) +
  geom_violin(size=1) +
  geom_boxplot(width=0.1, size=0.75) +
  geom_jitter(stat="identity",
              position=position_jitter(0.1), size=1, alpha=0.5) +
  geom_hline(yintercept = 0, lwd=1, col="grey60", lty=2) +
  scale_color_manual(values = cols) +
  ylim(c(-5, 0.5)) +
  ylab("log2 FC") +
  xlab("") +
  th + theme(legend.position = "none",
             axis.text.x = element_blank(),
             axis.ticks.x = element_blank())
plots[[2]] <- ggplot(df[df$class=="LCR",], aes(class, logFC, colour=class)) +
  geom_violin(size=1) +
  geom_boxplot(width=0.1, size=0.75) +
  geom_jitter(stat="identity",
              position=position_jitter(0.1), size=1, alpha=0.5) +
  geom_hline(yintercept = 0, lwd=1, col="grey60", lty=2) +
  scale_color_manual(values = cols) +
  scale_y_continuous(breaks = c(0,-1,-2), limits=c(-2, 0.5)) +
  ylab("log2 FC") +
  xlab("") +
  th + theme(legend.position = "none",
             axis.text.x = element_blank(),
             axis.ticks.x = element_blank())
plots[[3]] <- ggplot(df[df$class=="ENH",], aes(class, logFC, colour=class)) +
  geom_violin(size=1) +
  geom_boxplot(width=0.1, size=0.75) +
  geom_jitter(stat="identity",
              position=position_jitter(0.1), size=1, alpha=0.5) +
  geom_hline(yintercept = 0, lwd=1, col="grey60", lty=2) +
  scale_color_manual(values = cols) +
  scale_y_continuous(breaks = c(0,-1,-2), limits=c(-2, 0.5)) +
  ylab("log2 FC") +
  xlab("") +
  th + theme(legend.position = "none",
             axis.text.x = element_blank(),
             axis.ticks.x = element_blank())

pdf(paste0(dir, "figures/Figure3Ba.pdf"), width = 2, height = 5, useDingbats = FALSE)
plots[[1]]
dev.off()
pdf(paste0(dir, "figures/Figure3Bb.pdf"), width = 2, height = 2.5, useDingbats = FALSE)
plots[[2]]
dev.off()
pdf(paste0(dir, "figures/Figure3Bc.pdf"), width = 2, height = 2.5, useDingbats = FALSE)
plots[[3]]
dev.off()


## 3C ======
nt_guides <- gRNA_ann[gRNA_ann$class=="NT",]$ID
cells <- colnames(sce[,sce$call == "unique"])
cells_neg <- apply(counts(altExp(sce, 'gRNA_calls'))[nt_guides, cells], 1, function(x) names(which(x)))
cells_neg <- unlist(cells_neg)

# plot expression for each targeted TSS
plots <- list()
for(target in c("TMSB10_TSS", "CD2_TSS", "CD2_LCR1", "CTSC_ENH2")){
  # expected gene
  gene <- unique(gRNA_ann[gRNA_ann$target==target,]$expected_DE_gene)
  
  # collect p-values
  pvals <- mast[mast$target == target & mast$gene_id==gene,]$FDR
  names(pvals) <- mast[mast$target == target & mast$gene_id==gene,]$gRNA_id
  plots[[target]] <- plot_gene_gRNA(sce = sce,
                                    gene = gene,
                                    guides = gRNA_ann[gRNA_ann$target==target,]$ID, 
                                    target = target,
                                    nt_cells = cells_neg,
                                    per_guide = TRUE, 
                                    fdr_thr = 0.05,
                                    ann = pvals)
}

pdf(paste0(dir, "figures/Figure3C.pdf"), width = 2.5, height = 12, useDingbats = FALSE)
ggarrange(plotlist = plots, ncol=1, nrow=4, align="hv")
dev.off()


## 3D ======
means <- rowMeans(logcounts(sce))
genes_to_keep <- rowSums(counts(sce)>0) > ncol(sce)*0.05
genes_detected <- rownames(sce)[genes_to_keep]

degs <- read.table(paste0(dir, "results/06_DEgenes.tsv"), header=TRUE)

ranks <- means[order(means, decreasing = TRUE)]
ranks <- data.frame(mean=ranks,
                    rank=1:length(ranks))
ranks$detected <- row.names(ranks) %in% genes_detected
ranks$perturbed <- row.names(ranks) %in% gRNA_ann$expected_DE_gene
ranks$sig <- row.names(ranks) %in% degs$gene_id
ranks$class <- gRNA_ann[match(row.names(ranks), gRNA_ann$expected_DE_gene),]$class

ranks$dummy <- 1
ranks[ranks$sig,]$dummy <- 2.5
ranks[ranks$perturbed,]$dummy <- 6
ranks[which(ranks$sig & ranks$class != "TSS"),]$dummy <- 4
ranks[which(ranks$sig & ranks$class == "TSS"),]$dummy <- 5

ranks$status <- ifelse(ranks$sig & ranks$perturbed, "targeted + DE",
                       ifelse(ranks$sig, "non-targeted + DE",
                              ifelse(ranks$perturbed, "targeted + non-DE", 
                                     ifelse(ranks$detected, "tested", "not-tested"))))
ranks$group <- "not-targeted"
ranks[which(ranks$class == "TSS"),]$group <- "TSS"
ranks[which(ranks$class != "TSS"),]$group <- "INTERGENIC"
ranks$group <- factor(ranks$group, levels=c('not-targeted', 'TSS', 'INTERGENIC'))
ranks$col <- paste(ranks$status, ranks$group, sep=".")
ranks[which(ranks$class=="LCR"),]$col <- 'targeted + DE.LCR'

# for all non-targeted genes, subsample to only plot visible points
tmp <- ranks[ranks$col %in% c("not-tested.not-targeted", "tested.not-targeted"),]
tmp <- tmp[iSEE::subsetPointsByGrid(tmp$rank, rep(1, nrow(tmp)), resolution = 1e3),]
df <- rbind(tmp, ranks[!(ranks$col %in% c("not-tested.not-targeted", "tested.not-targeted")),])

pdf(paste0(dir, "figures/Figure3D.pdf"), width = 7, height = 3, useDingbats = FALSE)
ggplot(df[order(df$class),], aes(rank, dummy, shape=group, colour=col)) +
  geom_point(size=2) +
  scale_color_manual(values = c('tested.not-targeted'='grey50', 'not-tested.not-targeted'='grey80', 
                                'non-targeted + DE.not-targeted'='grey65', 
                                'targeted + DE.TSS'='#DDCC77', 'targeted + non-DE.TSS'='#DDCC7759',
                                'targeted + DE.LCR'='#88CCEE',
                                'targeted + DE.INTERGENIC'='#CC6677', 'targeted + non-DE.INTERGENIC'='#CC667759')) +
  ylim(0,7) +
  ylab("") +
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         shape=guide_legend(nrow=3, byrow=TRUE)) +
  labs(colour="", shape = "") +
  th + theme(axis.text.x = element_blank(),
             axis.text.y = element_blank(),
             axis.ticks.y = element_blank(),
             legend.position = "bottom")
dev.off()



### Figure S3 ###############

## S3A ======
mast <- read.table(paste0(dir, "results/05_DE_results.mast.tsv"), header = TRUE)
mast.exp <- mast[!is.na(mast$expected_gene),] # separate targets with expected effect
mast.exp$expected <- mast.exp$gene_id == mast.exp$expected_gene

# expected DEGs
tp <- mast.exp[mast.exp$expected & mast.exp$sig,]
# additional DEGs
other <- mast.exp[!mast.exp$expected & mast.exp$sig,]

# number sig gRNAs per target-gene pair
n_gRNAs <- table(tp$target)
tp$n_gRNAs <- as.numeric(n_gRNAs[match(tp$target, names(n_gRNAs))])

n_gRNAs <- table(paste0(other$target, other$gene_id))
other$n_gRNAs <- as.numeric(n_gRNAs[match(paste0(other$target, other$gene_id), names(n_gRNAs))])

tmp <- rbind(tp, other)
tmp <- unique(tmp[,c('expected', 'target', 'class', 'gene_id', 'n_gRNAs')])
tmp$n_gRNAs <- factor(tmp$n_gRNAs, levels=4:1)
tmp$hit <- ifelse(tmp$expected, "expected","other")
tmp$class <- factor(tmp$class, levels = c("TSS", "LCR", "ENH"))

pdf(paste0(dir, "figures/FigureS3A.pdf"), width = 8, height = 4, useDingbats = FALSE)
ggplot(tmp, aes(class, fill=n_gRNAs)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = brewer.pal(n=6, "PuBu")[4:1]) +
  facet_wrap(~hit) +
  xlab("") +
  ylab("fraction") +
  th
dev.off()


## S3B ======
mast.nt <- read.table(paste0(dir, "results/05_DE_results_NT.mast.tsv"), header = TRUE)

revlog_trans <- function(base = exp(10)){
  ## Define the desired transformation.
  trans <- function(x) {
    -log(x, base)
  }
  ## Define the reverse of the desired transformation
  inv <- function(x){
    base^(-x)
  }
  ## Creates the transformation
  trans_new(paste("revlog-", base, sep = ""),
            trans, ## The transformation function (can be defined using anonymous functions)
            inv,  ## The reverse of the transformation
            log_breaks(base = base), ## default way to define the scale breaks
            domain = c(1e-100, Inf) ## The domain over which the transformation is valued
  )
}

plot_qq_nt <- function(res=NULL, title=""){
  ci <- 0.95
  p_thresh <- 1e-8
  pvals_nt <- res[,c('gene_id', 'gRNA_id', 'p_value')]
  df <- pvals_nt %>%
    mutate(r = rank(p_value), expected = ppoints(n())[r],
           clower = qbeta(p=(1-ci)/2, shape1 = r, shape2 = n()+1-r),
           cupper = qbeta(p=(1+ci)/2, shape1 = r, shape2 = n()+1-r)) %>%
    ungroup() %>% mutate()
  df <- df[iSEE::subsetPointsByGrid(df$expected, df$p_value, resolution = 15e3),]
  p <- df %>% 
    mutate(p_value = ifelse(p_value < p_thresh, p_thresh, p_value)) %>%
    ggplot(aes(x = expected, y = p_value, ymin = clower, ymax = cupper)) +
    geom_point(size = 1, alpha = 0.5) +
    geom_ribbon(alpha = 0.2) +
    geom_abline(intercept = 0, slope = 1) +
    scale_x_continuous(trans = revlog_trans(base = 10)) + 
    scale_y_continuous(trans = revlog_trans(base = 10)) +
    xlab(expression(paste("Expected null p-value"))) +
    ylab(expression(paste("Observed p-value"))) +
    ggtitle(title) +
    th
  return(p)
}

pdf(paste0(dir, "figures/FigureS3B.pdf"), width = 5, height = 5, useDingbats = FALSE)
plot_qq_nt(res = mast.nt, title = "MAST")
dev.off()


## S3C ======
sceptre.nt <- read.table(paste0(dir, "results/05_DE_results_NT.sceptre.tsv"), header = TRUE)
limma.nt <- read.table(paste0(dir, "results/05_DE_results_NT.limma.tsv"), header = TRUE)
wilcox.nt <- read.table(paste0(dir, "results/05_DE_results_NT.wilcox.tsv"), header = TRUE)

plots <- list()
plots[[1]] <- plot_qq_nt(res = sceptre.nt, title = "sceptre")
plots[[2]] <- plot_qq_nt(res = limma.nt, title = "limma-voom")
plots[[3]] <- plot_qq_nt(res = wilcox.nt, title = "wilcoxson rank sum")

pdf(paste0(dir, "figures/FigureS3C.pdf"), width = 15, height = 5, useDingbats = FALSE)
ggarrange(plotlist = plots, ncol=3)
dev.off()


## S3D-E ======
compare <- read.table(paste0(dir, "results/05_DEGs_compareMethods.tsv"), header = TRUE)

tmp <- compare[compare$expected, c(seq(9,15,2),5,18)]
colnames(tmp)[1:4] <- c("sceptre", "MAST", "limma", "wilcoxon")

pdf(paste0(dir, "figures/FigureS3D_a.pdf"), width = 5, height = 4, useDingbats = FALSE)
upset(tmp, order.by = "freq", keep.order = TRUE,
      intersections = list(list("sceptre", "limma", "MAST", "wilcoxon"), list("sceptre", "limma", "MAST"), 
                           list("limma", "MAST", "wilcoxon"),
                           list("limma", "MAST"), list("sceptre", "limma"), list("limma", "wilcoxon"), 
                           list("wilcoxon"), list("sceptre"), list("limma")),
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
dev.off()

pdf(paste0(dir, "figures/FigureS3E_a.pdf"), width = 5, height = 4, useDingbats = FALSE)
upset(tmp, order.by = "freq", keep.order = TRUE,
      intersections = list(list("sceptre", "limma", "MAST", "wilcoxon"), list("sceptre", "limma", "MAST"), 
                           list("limma", "MAST", "wilcoxon"),
                           list("limma", "MAST"), list("sceptre", "limma"), list("limma", "wilcoxon"), 
                           list("wilcoxon"), list("sceptre"), list("limma")),
      queries = list(
        list(query = elements, 
             params = list("n_gRNAs", c(1:4)), query.name = "4",
             color = brewer.pal(n=6, "PuBu")[4], active = TRUE),
        list(query = elements, 
             params = list("n_gRNAs", c(1:3)), query.name = "3",
             color = brewer.pal(n=6, "PuBu")[3], active = TRUE),
        list(query = elements, 
             params = list("n_gRNAs", c(1:2)), query.name = "2",
             color = brewer.pal(n=6, "PuBu")[2], active = TRUE),
        list(query = elements, 
             params = list("n_gRNAs", c(1)), query.name = "1",
             color = brewer.pal(n=6, "PuBu")[1], active = TRUE)),
      query.legend="bottom")
dev.off()

tmp <- compare[!compare$expected, c(seq(9,15,2),5,18)]
colnames(tmp)[1:4] <- c("sceptre", "MAST", "limma", "wilcoxon")

pdf(paste0(dir, "figures/FigureS3D_b.pdf"), width = 7, height = 4, useDingbats = FALSE)
upset(tmp, order.by = "freq", keep.order = TRUE,
      intersections = list(list("sceptre", "limma", "MAST", "wilcoxon"), list("sceptre", "limma", "MAST"), 
                           list("limma", "MAST", "wilcoxon"), list("sceptre", "limma", "wilcoxon"), list("sceptre", "MAST", "wilcoxon"), 
                           list("limma", "MAST"), list("sceptre", "limma"), list("limma", "wilcoxon"), list("MAST", "wilcoxon"), 
                           list("sceptre", "wilcoxon"),  list("sceptre", "MAST"), 
                           list("wilcoxon"), list("sceptre"), list("limma"), list("MAST")),
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
dev.off()

pdf(paste0(dir, "figures/FigureS3E_b.pdf"), width = 7, height = 4, useDingbats = FALSE)
upset(tmp, order.by = "freq", keep.order = TRUE,
      intersections = list(list("sceptre", "limma", "MAST", "wilcoxon"), list("sceptre", "limma", "MAST"), 
                           list("limma", "MAST", "wilcoxon"), list("sceptre", "limma", "wilcoxon"), list("sceptre", "MAST", "wilcoxon"), 
                           list("limma", "MAST"), list("sceptre", "limma"), list("limma", "wilcoxon"), list("MAST", "wilcoxon"), 
                           list("sceptre", "wilcoxon"),  list("sceptre", "MAST"), 
                           list("wilcoxon"), list("sceptre"), list("limma"), list("MAST")),
      queries = list(
        list(query = elements, 
             params = list("n_gRNAs", c(1:4)), query.name = "4",
             color = brewer.pal(n=6, "PuBu")[4], active = TRUE),
        list(query = elements, 
             params = list("n_gRNAs", c(1:3)), query.name = "3",
             color = brewer.pal(n=6, "PuBu")[3], active = TRUE),
        list(query = elements, 
             params = list("n_gRNAs", c(1:2)), query.name = "2",
             color = brewer.pal(n=6, "PuBu")[2], active = TRUE),
        list(query = elements, 
             params = list("n_gRNAs", c(1)), query.name = "1",
             color = brewer.pal(n=6, "PuBu")[1], active = TRUE)),
      query.legend="bottom")
dev.off()






### Figure 4 ###############

## 4A ======
mast.tgt <- read.table(paste0(dir, "results/06_DEgenes.tsv"), header = TRUE)
mast.tgt$n_gRNAs_conc <- apply(mast.tgt, 1, function(x) max(x[c('n_down','n_up')]))
mast.tgt$n_gRNAs_conc <- factor(mast.tgt$n_gRNAs_conc, levels = 4:1)
mast.tgt$tier <- factor(mast.tgt$tier, levels = c("high","medium","low"))

pdf(paste0(dir, "figures/Figure4A.pdf"), width = 4, height = 4, useDingbats = FALSE)
ggplot(mast.tgt, aes(n_gRNAs_conc, FDR, colour=tier)) +
  geom_violin(scale="width") +
  geom_boxplot(width=0.1) +
  scale_color_manual(values = rev(brewer.pal(n=6, "Blues")[c(2,4,6)])) +
  xlab("number of gRNAs < 0.05 (raw gRNA-level p-values)") +
  ylab("target-level FDR") +
  th
dev.off()

## 4B ======
# only analyse itnergenic targets 
mast.tgt <- mast.tgt[mast.tgt$class %in% c("ENH", "INTRON", "INTERGENIC"),]
# remove low-confidence
mast.tgt <- mast.tgt[mast.tgt$tier != "low",]

# count total significant results - start with positive controls
n_DE <- as.data.frame(table(mast.tgt$target))
colnames(n_DE) <- c("target", "n_DEGs")
n_DE$class <- gRNA_ann[match(n_DE$target, gRNA_ann$target),'class']
n_DE$class <- factor(n_DE$class, levels=rev(c("ENH", "INTRON", "INTERGENIC")))

# plot
cols <- c(GetColors(n=1, scheme = "muted"), "#fca972", "#facb6e")
names(cols) <- c("ENH", "INTRON", "INTERGENIC")

pdf(paste0(dir, "figures/Figure4B.pdf"), width = 4, height = 4, useDingbats = FALSE)
ggplot(n_DE, aes(as.factor(n_DEGs), fill=class)) +
  geom_bar() +
  # coord_flip() +
  scale_fill_manual(values = cols) +
  xlab("number of DE genes within 1Mb") +
  ylab("number of targets") +
  th
dev.off()


## 4C ======
mast.tgt$class <- factor(mast.tgt$class, levels=c("ENH", "INTRON", "INTERGENIC"))

pdf(paste0(dir, "figures/Figure4C.pdf"), width = 5, height = 4, useDingbats = FALSE)
ggplot(mast.tgt, aes(distance/1e3, colour=class)) +
  geom_density(lwd=1) +
  scale_colour_manual(values = cols) +
  xlab("gene-target distance (kb)") +
  ylab("") +
  th
dev.off()

## 4D ======
mast.tgt$nearest_expr <- factor(mast.tgt$nearest_expr, levels=c(TRUE, FALSE))
pdf(paste0(dir, "figures/Figure4D.pdf"), width = 2, height = 4, useDingbats = FALSE)
ggplot(mast.tgt, aes(nearest_expr, distance/1e3)) +
  geom_boxplot() +
  ylab("gene-target distance (kb)") +
  xlab("") +
  th
dev.off()








### Figure 5 ###############

## 5A ======
plots <- list()
for(target in setdiff(unique(meta$target), "NT")){
  # expected gene
  gene <- unique(gRNA_ann[gRNA_ann$target==target,]$expected_DE_gene)
  
  hit <- mast.tgt[mast.tgt$pair_tgt == paste(target, gene, sep="|"),]
  p <- plot_gene_target(sce = sce, gRNA_ann = gRNA_ann, hit = hit, nt_cells = cells_neg, per_guide = FALSE)
  
  
  # collect p-values
  pvals <- mast[mast$target == target & mast$gene_id==gene,]$FDR
  names(pvals) <- mast[mast$target == target & mast$gene_id==gene,]$gRNA_id
  plots[[target]] <- plot_gene_gRNA(sce = sce,
                                    gene = gene,
                                    guides = gRNA_ann[gRNA_ann$target==target,]$ID, 
                                    target = target,
                                    nt_cells = cells_neg,
                                    per_guide = TRUE, 
                                    fdr_thr = 0.05,
                                    ann = pvals)
}

pdf(paste0(dir, "figures/Figure5A.pdf"), height = 3.5, width = 12, useDingbats = FALSE)
ggarrange(plotlist = plots[c(1,4,3,2)], ncol=4, nrow=1, align="hv")
dev.off()

## 5D ======
data.cpm <- read.table(paste0(dir, "results/08_validation_experiment_CPMexpr.tsv"), header=TRUE)

gene_ann <- read.table(paste0(dir, "data/validation_experiments/Homo_sapiens.GRCh38.96.ann"))[,c(1,3)]
colnames(gene_ann) <- c("id", "name")
meta <- read.table(paste0(dir, "data/validation_experiments/validation_metadata.tsv"), header = TRUE)
meta <- meta[meta$IRMS_Sample_ID %in% colnames(data.cpm),]
meta$target <- as.factor(meta$target)
meta$target <- relevel(meta$target, ref = "NT")

plots <- list()
for(target in levels(meta$target)[-1]){
  id <- gene_ann[gene_ann$name == meta[meta$target==target,'expected_gene'],'id']
  df <- meta[meta$target %in% c(target, "NT"),]
  df$expr <- as.numeric(data.cpm[id, df$IRMS_Sample_ID])
  df$target <- factor(df$target, levels=unique(df$target))
  
  plots[[target]] <- ggplot(df, aes(target, expr)) +
    geom_boxplot(colour="gray90") +
    geom_point(aes(colour=donor), size=2) +
    scale_color_manual(values = as.character(GetColors(n=2, scheme="bright"))) +
    xlab("") +
    ylab("log2 CPM") +
    ggtitle(target) +
    th
}

pdf(paste0(dir, "figures/Figure5D.pdf"), width = 10, height = 2.5, useDingbats = FALSE)
ggarrange(plotlist = plots[c(1,3,4,2)], ncol=4, nrow=1, common.legend = TRUE, legend = "bottom")
dev.off()





### Figure 6 ###############

## 6B ======
hit <- mast.tgt[mast.tgt$pair_tgt=="GIGYF1_ENH|GIGYF1",]
p <- plot_gene_target(sce = sce, gRNA_ann = gRNA_ann, hit = hit, nt_cells = cells_neg, per_guide = FALSE)

pdf(paste0(dir, "figures/Figure6B.pdf"), width = 6, height = 4, useDingbats = FALSE)
p
dev.off()


## 6F ======
hit <- mast.tgt[mast.tgt$pair_tgt=="PHF19_ENH|TRAF1",]
p <- plot_gene_target(sce = sce, gRNA_ann = gRNA_ann, hit = hit, nt_cells = cells_neg, per_guide = FALSE)

pdf(paste0(dir, "figures/Figure6F.pdf"), width = 6, height = 4, useDingbats = FALSE)
p
dev.off()

## 6J ======
hit <- mast.tgt[mast.tgt$pair_tgt=="CXCR5_INTERGENIC|CD3D",]
p <- plot_gene_target(sce = sce, gRNA_ann = gRNA_ann, hit = hit, nt_cells = cells_neg, per_guide = TRUE)

pdf(paste0(dir, "figures/Figure56.pdf"), width = 5, height = 4, useDingbats = FALSE)
p
dev.off()

