### SCEPTRE
library(DropletUtils)
library(scran)
library(ondisc)
library(dplyr)

dir <- paste0(getwd(), "/")

#### Filter data
sce <- readRDS(paste0(dir, "results/03_sce_POC_Tcells.goodQual.calls.NORM.Rds"))

#### Save data in format required by `sceptre`
# directories to save outputs
out_dir <- paste0(dir, "results/04_sceptre/")
raw_data_dir <- paste0(out_dir, "raw_data/")
processed_data_dir <- paste0(out_dir, "processed_data/")
if (!dir.exists(raw_data_dir)) dir.create(path = raw_data_dir, recursive = TRUE)
if (!dir.exists(processed_data_dir)) dir.create(path = processed_data_dir)

## convert gene expression data
# write as mtx
write10xCounts(path = paste0(raw_data_dir, "gene_expression/"),
               x = counts(sce),
               barcodes = colnames(sce),
               gene.id = row.names(sce),
               overwrite = TRUE)
# save as on_disc matrix
sce_expr <- create_ondisc_matrix_from_mtx(
  mtx_fp = paste0(raw_data_dir, "gene_expression/matrix.mtx"),
  barcodes_fp = paste0(raw_data_dir, "gene_expression/barcodes.tsv"),
  features_fp = paste0(raw_data_dir, "gene_expression/genes.tsv"),
  n_lines_per_chunk = 35000000,
  on_disk_dir = processed_data_dir,
  file_name = "Tcells_POC_countData",
  return_metadata_ondisc_matrix = TRUE,
  progress = FALSE)
saveRDS(sce_expr, paste0(raw_data_dir, "Tcells_POC_countData.Rds"))

## convert gRNA calls
# write as mtx
write10xCounts(path = paste0(raw_data_dir, "gRNA_calls/"),
               x = counts(altExp(sce, 'gRNA_calls')),
               barcodes = colnames(sce),
               gene.id = row.names(altExp(sce, 'gRNA_calls')),
               overwrite = TRUE)
# save as on_disc matrix
sce_perturb <- create_ondisc_matrix_from_mtx(
  mtx_fp = paste0(raw_data_dir, "gRNA_calls/matrix.mtx"),
  barcodes_fp = paste0(raw_data_dir, "gRNA_calls/barcodes.tsv"),
  features_fp = paste0(raw_data_dir, "gRNA_calls/genes.tsv"),
  on_disk_dir = processed_data_dir,
  file_name = "Tcells_POC_gRNAcalls",
  return_metadata_ondisc_matrix = TRUE,
  progress = FALSE)
saveRDS(sce_perturb, paste0(raw_data_dir, "Tcells_POC_gRNAcalls.Rds"))

## create combined object
# keep genes in 5% of cells
genes_to_keep <- read.table(paste0(dir, "results/04_genes_detected_5pct.tsv"))
genes_to_keep <- genes_to_keep$V1
sce_expr <- sce_expr[genes_to_keep,]

crispr_experiment <- multimodal_ondisc_matrix(list(expressions = sce_expr,
                                                   perturbations = sce_perturb))
saveRDS(crispr_experiment@modalities$expressions@ondisc_matrix,
        paste0(processed_data_dir, "Tcells_POC_countData.Rds"))
saveRDS(crispr_experiment@modalities$perturbations@ondisc_matrix,
        paste0(processed_data_dir, "Tcells_POC_gRNAcalls.Rds"))


#### Save the covariates to use in the model
cell_covariate_matrix <- summarize(crispr_experiment@global_cell_covariates,
                                   lg_n_genes_expressed = log(expressions_n_nonzero))
cell_covariate_matrix$batch <- sce$batch

saveRDS(cell_covariate_matrix, paste0(processed_data_dir, "cell_covariate_matrix.Rds"))

#### Define gene-gRNA pairs
## defined in 04.1_gene_gRNA_pairs.R
gene_gRNA_pairs <- readRDS(paste0(dir, "results/04_gene_gRNA_pairs.Rds"))


#### Compute pod sizes to use in run
c(gene = ceiling(length(unique(gene_gRNA_pairs$gene_id))/48), 
  gRNA = ceiling(length(unique(gene_gRNA_pairs$gRNA_id))/48), 
  pair = ceiling(nrow(gene_gRNA_pairs)/48))


sessionInfo()
