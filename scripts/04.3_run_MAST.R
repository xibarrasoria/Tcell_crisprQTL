### MAST - T cell full pilot
library(scran)
library(MAST)
library(data.table)
dir <- paste0(getwd(), "/")

## guides to test
args = commandArgs(trailingOnly=TRUE)
start <- args[1]
end <- args[2]
seed <- args[3]

## parallelisation options
options(mc.cores = 48)

### DATA
sce <- readRDS(paste0(dir, "results/03_sce_POC_Tcells.goodQual.calls.NORM.Rds"))

## compute the 'cellular detection rate'
cdr <- 10^sce$n_genes
sce$cdr <- scale(cdr)


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
cells <- colnames(sce[,sce$call == "unique"])
guides <- lib[lib$class == "NT",]$ID
cells_neg <- apply(counts(altExp(sce, 'gRNA_calls'))[guides, cells], 1, function(x) names(which(x)))
cells_neg <- unlist(cells_neg)
# select a random sample of 5K cells as NT group
set.seed(seed)
cells_neg <- sample(cells_neg, size=5e3)

## fit model and test for each gRNA
res <- list()
res_1mb <- list()
start_time <- Sys.time()
for(target in unique(lib$target)[start:end]){
  # select cells of interest
  guides <- lib[lib$target == target,]$ID
  cells_pos <- apply(counts(altExp(sce, 'gRNA_calls'))[guides,], 1, function(x) names(which(x)))
  cells_pos <- reshape2::melt(cells_pos)

  # subset sce and convert to sca
  tmp <- sce[,c(cells_pos$value, cells_neg)]
  tmp <- SceToSingleCellAssay(tmp)

  # assign groups based on gRNA
  groups <- paste0("Perturbed", 1:length(guides))
  names(groups) <- unique(cells_pos$L1)
  guides <- unname(c(groups[cells_pos$L1], rep("NT", length(cells_neg))))
  guides <- factor(guides, levels=c("NT", paste0("Perturbed", 1:length(groups))))
  tmp$condition <- guides
  
  # fit hurdle model using all data for the same target
  fit <- zlm(formula = ~condition + cdr + batch,
             sca = tmp,
             parallel = TRUE)
  
  ## test the effect of each gRNA
  for(i in 1:length(groups)){
    test <- paste(target, names(groups)[i], sep=".")
    res[[test]] <- summary(fit, doLRT=paste0('conditionPerturbed', i))
    res[[test]] <- res[[test]]$datatable
    res[[test]] <- merge(res[[test]][contrast==paste0('conditionPerturbed', i) & 
                                       component=='H',.(primerid, `Pr(>Chisq)`)],
                         res[[test]][contrast==paste0('conditionPerturbed', i) & 
                                       component=='logFC', .(primerid, coef, ci.hi, ci.lo)],
                         by='primerid')
    res[[test]]$FDR <- p.adjust(res[[test]]$`Pr(>Chisq)`, 'fdr') # adjust for multitesting
    ## add info on the test to result dataframe
    res[[test]]$target <- target
    res[[test]]$gRNA <- names(groups)[i]
    
    ## subset to +-1Mb neighbourhood and readjust p-values
    res_1mb[[test]] <- res[[test]][res[[test]]$primerid %in% 
                                     gene_gRNA_pairs[gene_gRNA_pairs$gRNA_id == names(groups)[i],]$gene_id,]
    res_1mb[[test]]$FDR <- p.adjust(res_1mb[[test]]$`Pr(>Chisq)`, 'fdr')
  }
  
  
}
end_time <- Sys.time()

paste0("Runtime: ", round(as.numeric(end_time - start_time, units="mins"),2), " mins (",
       round(as.numeric(end_time - start_time, units="hours"),2), " hours)")

### SAVE
saveRDS(res, paste0(dir, "results/04_MAST/04.3_results_MAST_genomeWide_seed", seed, "_guides", start, "-", end, ".Rds"))
saveRDS(res_1mb, paste0(dir, "results/04_MAST/04.3_results_MAST_neighbourhood_seed", seed, "_guides", start, "-", end, ".Rds"))

