## Define the gene-gRNA pairs to test for cis analysis
library(scran)
library(biomaRt)

dir <- paste0(getwd(), "/")

## gRNA annotation
gRNA_ann <- read.table(paste0(dir, "data/gRNA_library.tsv"), sep="\t", header = TRUE)

## data
sce <- readRDS(paste0(dir, "results/03_sce_POC_Tcells.goodQual.calls.NORM.Rds"))
# remove genes that are not expressed in at least 5% of the cells
genes_to_keep <- rowSums(counts(sce)>0) > ncol(sce)*0.05
genes_detected <- rownames(sce)[genes_to_keep]
# remove the Cas9 sequences added to the genome
genes_detected <- grep("Cas9", genes_detected, value = TRUE, invert = TRUE)


## connection with Ensembl to retrieve genes (matche version to annotation used by cellranger)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version=98)

## get the genes in the neighbourhood
distance <- 1e6
gene_gRNA_pairs <- sapply(setdiff(unique(gRNA_ann$target),"NT"), function(t){
  ## get coordinates of target
  coords <- unique(gRNA_ann[gRNA_ann$target==t, c('chr','start','end')])
  # remmove 'chr' prefix to match ensembl chr names
  coords[1] <- substr(coords[1], 4,5)
  
  ## retreive nearby genes
  genes <- getBM(attributes = c('external_gene_name'),
                 filters = c('chromosome_name','start','end'),
                 values = list(coords$chr, 
                               coords$start-distance, 
                               coords$end+distance),
                 mart = ensembl)
  
  ## restrict to genes expressed in enough cells
  genes <- intersect(genes$external_gene_name, genes_detected)
  
  ## select guides for target
  g <- gRNA_ann[which(gRNA_ann$target==t),]$ID
  
  ## pair each gene to the gRNAs for that target
  data.frame(gene_id = rep(genes, length(g)), gRNA_id = rep(g, each=length(genes)))
}, simplify = FALSE)
gene_gRNA_pairs <- do.call(rbind, t(gene_gRNA_pairs))

## also test each NT control versus all genes selected above
nt_guides <- gRNA_ann[gRNA_ann$class=="NT",]$ID
gene_NT_pairs <- do.call(rbind,
                         sapply(nt_guides, function(guide) 
                           data.frame(gene_id = unique(gene_gRNA_pairs$gene_id), gRNA_id = guide), 
                           simplify = FALSE))

## save
gene_gRNA_pairs <- rbind(gene_gRNA_pairs, gene_NT_pairs)
saveRDS(gene_gRNA_pairs, paste0(dir, "results/04_gene_gRNA_pairs.Rds"))
write.table(genes_detected, paste0(dir, "results/04_genes_detected_5pct.tsv"), quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
