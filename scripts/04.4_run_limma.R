### limma-voom - T cell full pilot
library(scran)
library(limma)
dir <- paste0(getwd(), "/")

### DATA
sce <- readRDS(paste0(dir, "results/03_sce_POC_Tcells.goodQual.calls.NORM.Rds"))

## filter out lowly expressed genes
n_cells_detected <- rowSums(counts(sce)>0)
# keep only genes present in at least 5% of the cells
genes_detected <- read.table(paste0(dir, "results/04_genes_detected_5pct.tsv"))
genes_detected <- genes_detected$V1
sce <- sce[genes_detected,]


### LIMMA
## convert sce to dge to use with limma (this uses simple library size normalisation)
dge <- convertTo(sce)

## gRNA library
lib <- read.table(paste0(dir, "data/gRNA_library.tsv"), sep="\t", header=TRUE)
## to perform control test using NT gRNAs, we group them into 'targets' with 4 NT gRNAs each
lib[lib$class == "NT",]$target <- rep(paste0("NT_", 1:9), each=4)[1:35]

## gene-gRNA pairs in neighbourhood
gene_gRNA_pairs <- readRDS(paste0(dir, "results/04_gene_gRNA_pairs.Rds"))

## identify cells with NT controls
# only consider cells with a NT gRNA and nothing else (i.e. exclude cells with also a second gRNA)
cells <- colnames(sce[,sce$call == "unique"])
guides <- lib[lib$class == "NT",]$ID
cells_neg <- apply(counts(altExp(sce, 'gRNA_calls'))[guides, cells], 1, function(x) names(which(x)))
cells_neg <- unlist(cells_neg)
# select a random sample of 5K cells as NT group
set.seed(749)
cells_neg <- sample(cells_neg, size=5e3)


## fit model and test for each gRNA
res <- list()
res_1mb <- list()
start_time <- Sys.time()
for(target in unique(lib$target)){
  # select cells of interest
  guides <- lib[lib$target == target,]$ID
  cells_pos <- apply(counts(altExp(sce, 'gRNA_calls'))[guides,], 1, function(x) names(which(x)))
  cells_pos <- reshape2::melt(cells_pos)
  
  # subset dge
  tmp <- dge[,match(c(cells_pos$value, cells_neg), row.names(dge$samples))]

  # assign groups based on gRNA
  groups <- paste0("Perturbed", 1:length(guides))
  names(groups) <- unique(cells_pos$L1)
  guides <- unname(c(groups[cells_pos$L1], rep("NT", length(cells_neg))))
  guides <- factor(guides, levels=c("NT", paste0("Perturbed", 1:length(groups))))
  tmp$samples$condition <- guides

  # define model
  # we do not adjust for library size because size factors are already included
  design <- model.matrix(~batch + n_genes + condition, tmp$samples)
  
  # Standard voom-limma pipeline
  # we use all genes to fit the model
  v <- voom(tmp, design)
  fit <- lmFit(v, design)
  fit <- eBayes(fit, robust = TRUE)
  
  ## test the effect of each gRNA
  for(i in 1:length(groups)){
    test <- paste(target, names(groups)[i], sep=".")
    res[[test]] <- topTable(fit, n=Inf, coef=paste0('conditionPerturbed', i))[,-c(3)]
    ## add info on the test to result dataframe
    res[[test]]$target <- target
    res[[test]]$gRNA <- names(groups)[i]
    
    ## subset to +-1Mb neighbourhood and readjust p-values
    res_1mb[[test]] <- res[[test]][res[[test]]$Symbol %in% 
                                     gene_gRNA_pairs[gene_gRNA_pairs$gRNA_id == names(groups)[i],]$gene_id,]
    res_1mb[[test]]$adj.P.Val <- p.adjust(res_1mb[[test]]$P.Value, 'fdr')
  }
}
end_time <- Sys.time()

paste0("Runtime: ", round(as.numeric(end_time - start_time, units="mins"),2), " mins (",
       round(as.numeric(end_time - start_time, units="hours"),2), " hours)")


### SAVE
saveRDS(res, paste0(dir, "results/04.4_results_limma_genomeWide.Rds"))
saveRDS(res_1mb, paste0(dir, "results/04.4_results_limma_neighbourhood.Rds"))
