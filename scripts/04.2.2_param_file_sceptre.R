library(ondisc)
dir <- paste0(getwd(), "/")

param_funct <- function(param) {
  switch(param,
         # modify the parameters below
         n_processors = 48,
         storage_dir = paste0(dir, "04_sceptre/"),
         expression_matrix = readRDS(paste0(dir, "04_sceptre/processed_data/Tcells_POC_countData.Rds")),
         perturbation_matrix = readRDS(paste0(dir, "04_sceptre/processed_data/Tcells_POC_gRNAcalls.Rds")),
         covariate_matrix = readRDS(paste0(dir, "04_sceptre/processed_data/cell_covariate_matrix.Rds")),
         gene_gRNA_pairs = readRDS(paste0(dir, "04_gene_gRNA_pairs.Rds")),
         side = "both",
         pod_sizes = c(gene = 24, gRNA = 8, pair = 950),
         regularization_amount = 1,
         seed = 749,
         B = 500)
}
