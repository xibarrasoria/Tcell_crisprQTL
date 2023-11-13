library(scran)
library(BiocParallel)
library(BiocSingular)
dir <- paste0(getwd(), "/")

## configure parallelisation options
bpp <- MulticoreParam(48)


### DATA
sce <- readRDS(paste0(dir, "results/03_sce_POC_Tcells.goodQual.calls.NORM.Rds"))


### TEST
## gRNA library
lib <- read.table(paste0(dir, "data/gRNA_library.tsv"), sep="\t", header=TRUE)
## to perform control test using NT gRNAs, we group them into 'targets' with 4 NT gRNAs each
lib[lib$class == "NT",]$target <- rep(paste0("NT_", 1:9), each=4)[1:35]

## gene-gRNA pairs in neighbourhood
gene_gRNA_pairs <- readRDS(paste0(dir, "results/04_gene_gRNA_pairs.Rds"))

## filter out lowly expressed genes
genes_detected <- read.table(paste0(dir, "results/04_genes_detected_5pct.tsv"))
genes_detected <- genes_detected$V1
sce <- sce[genes_detected,]

## identify cells with NT controls
# only consider cells with a NT gRNA and nothing else (i.e. exclude cells with also a second gRNA)
guides <- lib[lib$class == "NT",]$ID
cells <- colnames(sce[,sce$call == "unique"])
cells_neg <- apply(counts(altExp(sce, 'gRNA_calls'))[guides, cells], 1, function(x) names(which(x)))
cells_neg <- unlist(cells_neg)
# select a random sample of 5K cells as NT group
set.seed(749)
cells_neg <- sample(cells_neg, size=5e3)

### WILCOXON
res <- list()
res_1mb <- list()
start_time <- Sys.time()
for(guide in lib$ID){
  # select cells of interest
  cells_pos <- names(which(counts(altExp(sce, 'gRNA_calls'))[guide,]))

  # subset sce
  tmp <- sce[,match(c(cells_pos, cells_neg), colnames(sce))]
  
  # assign groups
  groups <- rep(c("Perturbed", "NT"), times=c(length(cells_pos), length(cells_neg)))
  groups <- factor(groups, levels=c("NT", "Perturbed"))
  tmp$condition <- groups
  
  # wilcoxon test within scran framework for parallelisation
  # we compute results genome-wide for later use
  res[[guide]] <- pairwiseWilcox(tmp, 
                                 groups = tmp$condition,
                                 block = tmp$batch,
                                 direction = "any",
                                 BPPARAM = bpp)$statistics[[2]] # select only one comparison
  res[[guide]]$gRNA <- guide
  
  ## subset to +-1Mb neighbourhood and test again
  res_1mb[[guide]] <- pairwiseWilcox(tmp[gene_gRNA_pairs[gene_gRNA_pairs$gRNA_id == guide,]$gene_id,], 
                                     groups = tmp$condition,
                                     block = tmp$batch,
                                     direction = "any",
                                     BPPARAM = bpp)$statistics[[2]] # select only one comparison
  res_1mb[[guide]]$gRNA <- guide
}
end_time <- Sys.time()
paste0("Runtime: ", round(as.numeric(end_time - start_time, units="mins"),2), " mins (",
       round(as.numeric(end_time - start_time, units="hours"),2), " hours)")


### SAVE
saveRDS(res, paste0(dir, "results/04.5_results_wilcoxon_genomeWide.Rds"))
saveRDS(res_1mb, paste0(dir, "results/04.5_results_wilcoxon_neighbourhood.Rds"))
