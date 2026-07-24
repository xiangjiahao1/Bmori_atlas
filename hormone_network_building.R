################################################################################
# GENIE3 Gene Regulatory Network Construction
# 
# This script constructs gene regulatory networks (GRNs) using GENIE3
# on single-nucleus transcriptomic data of the silkworm brain.
# 
# It focuses on specific cell subclusters (e.g., KC27) and extracts
# regulatory relationships targeting hormone receptor genes (EcR, Met).
#
# Input: Seurat object with cell type annotations
# Output: Weighted adjacency matrices, target gene statistics
#
# Author: Jiahao Xiang
# Date: 2026-07-24
################################################################################

# =============================================================================
# 0. Load Libraries
# =============================================================================
library(GENIE3)
library(Seurat)
library(dplyr)
library(ggplot2)

# =============================================================================
# 1. Configuration (USER EDIT HERE)
# =============================================================================

# Input paths
SEURAT_OBJECT_PATH <- "/share/home/xiangjiahao/01.sc_silkworm/01.integrated/Bmor_only_1replicate_per_sample_integrated.rds"
TF_LIST_PATH <- "/share/home/xiangjiahao/01.sc_silkworm/01.integrated/id_change_script/silkworm_tf/silkworm_tf_list_symbol"

# Output settings
OUTPUT_BASE_DIR <- "/share/home/xiangjiahao/01.sc_silkworm/01.integrated/26.network"

# Cluster to analyze
TARGET_CLUSTER <- 27
CLUSTER_NAME <- "kc27_network"

# GENIE3 parameters
GENE_EXPRESSION_THRESHOLD <- 0.01   # Minimum fraction of cells expressing a gene
N_CORES <- 8
SEED <- 1234

# Target genes of interest (for network extraction)
TARGET_GENES <- c("EcR", "gce")      # Note: 'gce' is the silkworm homolog of Met

# =============================================================================
# 2. Helper Functions
# =============================================================================

#' Filter genes by expression prevalence
#' 
#' @param assay_data Gene expression matrix (genes x cells)
#' @param min_cell_fraction Minimum fraction of cells expressing the gene
#' @return Character vector of gene names passing filter
filter_genes_by_prevalence <- function(assay_data, min_cell_fraction = 0.01) {
  n_cells <- ncol(assay_data)
  keep <- rowSums(assay_data != 0) / n_cells > min_cell_fraction
  genes_keep <- names(keep[keep])
  cat(sprintf("  Filtered genes: %d out of %d retained (%.1f%%)\n",
              length(genes_keep), nrow(assay_data),
              100 * length(genes_keep) / nrow(assay_data)))
  return(genes_keep)
}

#' Run GENIE3 network inference
#' 
#' @param expr_matrix Expression matrix (genes x cells)
#' @param regulators Character vector of candidate regulator genes
#' @param n_cores Number of CPU cores for parallelization
#' @param seed Random seed for reproducibility
#' @return Weighted adjacency matrix (regulators x target genes)
run_genie3_network <- function(expr_matrix, regulators, 
                               n_cores = 8, seed = 1234) {
  set.seed(seed)
  cat("Running GENIE3 network inference...\n")
  cat("  - Regulators:", length(regulators), "\n")
  cat("  - Target genes:", nrow(expr_matrix), "\n")
  
  weight_mat <- GENIE3(expr_matrix, 
                       regulators = regulators,
                       nCores = n_cores, 
                       verbose = TRUE)
  
  # Transpose: regulators as rows, targets as columns
  weight_mat <- t(as.data.frame(weight_mat))
  cat("  - Network complete. Matrix dimensions:", dim(weight_mat), "\n")
  
  return(weight_mat)
}

#' Extract regulatory links from adjacency matrix
#' 
#' @param weight_mat Weighted adjacency matrix
#' @param threshold Minimum weight threshold (default: 0.01)
#' @return Data frame with regulatoryGene, targetGene, weight
extract_links <- function(weight_mat, threshold = 0.01) {
  # Convert to list format
  links <- getLinkList(weight_mat, reportMax = NULL, threshold = threshold)
  links <- as.data.frame(links)
  colnames(links) <- c("regulatoryGene", "targetGene", "weight")
  cat("  - Extracted", nrow(links), "regulatory links (threshold =", threshold, ")\n")
  return(links)
}

#' Compute statistics for a target gene's regulators
#' 
#' @param links Data frame of regulatory links
#' @param target_gene Name of target gene
#' @param tf_list Character vector of transcription factors
#' @return List with total regulators and TF regulators
compute_target_stats <- function(links, target_gene, tf_list) {
  # Filter links targeting the gene of interest
  target_links <- links[which(links$targetGene == target_gene), ]
  
  # Total unique regulators
  total_regulators <- length(unique(target_links$regulatoryGene))
  
  # Regulators that are TFs
  tf_regulators <- intersect(unique(target_links$regulatoryGene), tf_list)
  tf_regulator_count <- length(tf_regulators)
  
  cat(sprintf("  %s: %d regulators, %d of which are TFs\n",
              target_gene, total_regulators, tf_regulator_count))
  
  return(list(
    target_gene = target_gene,
    total_regulators = total_regulators,
    tf_regulators = tf_regulator_count,
    regulator_list = unique(target_links$regulatoryGene),
    tf_regulator_list = tf_regulators,
    links = target_links
  ))
}

#' Save network summary statistics
#' 
#' @param stats_list List of statistics from compute_target_stats
#' @param cluster_name Name of the cluster
#' @param output_dir Output directory
#' @return Data frame of statistics
save_network_stats <- function(stats_list, cluster_name, output_dir) {
  stats_df <- data.frame(
    target_gene = character(),
    total_regulators = numeric(),
    tf_regulators = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (stat in stats_list) {
    stats_df <- rbind(stats_df, data.frame(
      target_gene = stat$target_gene,
      total_regulators = stat$total_regulators,
      tf_regulators = stat$tf_regulators,
      stringsAsFactors = FALSE
    ))
  }
  
  stats_df$cluster <- cluster_name
  
  write.csv(stats_df, 
            file.path(output_dir, "network_statistics.csv"),
            row.names = FALSE)
  
  return(stats_df)
}

# =============================================================================
# 3. Main Analysis Function
# =============================================================================

#' Run complete GENIE3 network analysis for a cell subcluster
#' 
#' @param seurat_path Path to Seurat object RDS file
#' @param tf_path Path to TF list file (one gene per line)
#' @param cluster_id Cluster identity to subset
#' @param cluster_name Name for output directory
#' @param output_base_dir Base output directory
#' @param min_cell_fraction Minimum fraction of cells expressing a gene
#' @param n_cores Number of cores for GENIE3
#' @param seed Random seed
#' @param target_genes Vector of target genes to extract statistics for
#' @return List containing network matrix, links, and statistics
run_network_analysis <- function(seurat_path, 
                                 tf_path,
                                 cluster_id, 
                                 cluster_name,
                                 output_base_dir,
                                 min_cell_fraction = 0.01,
                                 n_cores = 8,
                                 seed = 1234,
                                 target_genes = c("EcR", "gce")) {
  
  start_time <- Sys.time()
  cat("\n========================================\n")
  cat("Starting GENIE3 network analysis\n")
  cat("Cluster:", cluster_name, "(ID:", cluster_id, ")\n")
  cat("========================================\n\n")
  
  # ---- 3.1 Create output directory ----
  output_dir <- file.path(output_base_dir, cluster_name)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
  }
  setwd(output_dir)
  
  # ---- 3.2 Load data ----
  cat("Loading Seurat object...\n")
  ec <- readRDS(seurat_path)
  
  cat("Loading TF list...\n")
  tf <- read.table(tf_path, header = FALSE, stringsAsFactors = FALSE)
  tf_list <- tf$V1
  cat("  - Total TFs:", length(tf_list), "\n")
  
  # ---- 3.3 Subset cluster of interest ----
  Idents(ec) <- ec$integrated_snn_res.0.6
  ec_subset <- subset(ec, ident = cluster_id)
  cat("Subset cluster", cluster_id, ":", ncol(ec_subset), "cells\n")
  
  # ---- 3.4 Prepare expression matrix ----
  cat("Preparing expression matrix...\n")
  assay_data <- GetAssayData(object = ec_subset, 
                             slot = "counts", 
                             assay = "RNA")
  
  # Filter by expression prevalence
  keep_genes <- filter_genes_by_prevalence(assay_data, min_cell_fraction)
  expr_matrix <- as.matrix(assay_data[keep_genes, ])
  
  # Intersect with TFs for regulators
  regulators <- intersect(tf_list, keep_genes)
  cat("  - Regulators (TFs expressed in this cluster):", length(regulators), "\n")
  
  # ---- 3.5 Run GENIE3 ----
  weight_mat <- run_genie3_network(expr_matrix, 
                                   regulators = regulators,
                                   n_cores = n_cores, 
                                   seed = seed)
  
  # Save full network
  saveRDS(weight_mat, file.path(output_dir, "genie3_network_matrix.rds"))
  cat("  - Full network saved to:", file.path(output_dir, "genie3_network_matrix.rds"), "\n")
  
  # ---- 3.6 Extract links ----
  links <- extract_links(weight_mat, threshold = 0.01)
  write.csv(links, file.path(output_dir, "network_links_all.csv"), row.names = FALSE)
  
  # ---- 3.7 Extract target gene statistics ----
  cat("\nExtracting statistics for target genes...\n")
  stats_list <- list()
  
  for (target in target_genes) {
    # Check if target gene exists in the network
    if (target %in% colnames(weight_mat)) {
      stats <- compute_target_stats(links, target, tf_list)
      stats_list[[target]] <- stats
      
      # Save individual target gene network
      target_links <- stats$links
      write.csv(target_links, 
                file.path(output_dir, paste0(target, "_network.csv")),
                row.names = FALSE)
    } else {
      cat("  Warning: Target gene", target, "not found in network\n")
    }
  }
  
  # ---- 3.8 Save summary statistics ----
  stats_df <- save_network_stats(stats_list, cluster_name, output_dir)
  
  # ---- 3.9 Summary ----
  elapsed <- Sys.time() - start_time
  cat("\n========================================\n")
  cat("Analysis complete!\n")
  cat("  - Cluster:", cluster_name, "\n")
  cat("  - Cells:", ncol(ec_subset), "\n")
  cat("  - Genes in network:", nrow(expr_matrix), "\n")
  cat("  - Regulators:", length(regulators), "\n")
  cat("  - Regulatory links:", nrow(links), "\n")
  cat("  - Output directory:", output_dir, "\n")
  cat("  - Time elapsed:", round(elapsed, 2), units(elapsed), "\n")
  cat("========================================\n")
  
  # Return results
  return(list(
    seurat_subset = ec_subset,
    expression_matrix = expr_matrix,
    regulators = regulators,
    network_matrix = weight_mat,
    links = links,
    statistics = stats_list,
    stats_df = stats_df,
    output_dir = output_dir
  ))
}

# =============================================================================
# 4. Execute Analysis
# =============================================================================

# Run the full analysis
results <- run_network_analysis(
  seurat_path = SEURAT_OBJECT_PATH,
  tf_path = TF_LIST_PATH,
  cluster_id = TARGET_CLUSTER,
  cluster_name = CLUSTER_NAME,
  output_base_dir = OUTPUT_BASE_DIR,
  min_cell_fraction = GENE_EXPRESSION_THRESHOLD,
  n_cores = N_CORES,
  seed = SEED,
  target_genes = TARGET_GENES
)

# =============================================================================
# 5. Additional Analysis Functions
# =============================================================================

#' Visualize regulatory network for a target gene
#' 
#' @param target_stats Statistics list from compute_target_stats
#' @param top_n Number of top regulators to show
#' @param save_path Path to save plot
plot_target_network <- function(target_stats, top_n = 20, save_path = NULL) {
  if (is.null(target_stats) || nrow(target_stats$links) == 0) {
    cat("No data to plot\n")
    return(NULL)
  }
  
  library(ggplot2)
  
  # Get top regulators by weight
  top_regulators <- target_stats$links %>%
    arrange(desc(weight)) %>%
    head(top_n)
  
  p <- ggplot(top_regulators, aes(x = reorder(regulatoryGene, weight), 
                                   y = weight)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = paste("Top regulators of", target_stats$target_gene),
         x = "Regulator", y = "GENIE3 weight") +
    theme_minimal()
  
  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 8, height = 6)
    cat("Plot saved to:", save_path, "\n")
  }
  
  return(p)
}

#' Run network analysis on multiple clusters
#' 
#' @param cluster_ids Vector of cluster IDs to analyze
#' @param ... Additional arguments passed to run_network_analysis
run_batch_network_analysis <- function(cluster_ids, ...) {
  all_results <- list()
  
  for (cid in cluster_ids) {
    cluster_name <- paste0("cluster_", cid)
    cat("\n\n========== Processing cluster", cid, "==========\n")
    
    result <- run_network_analysis(
      cluster_id = cid,
      cluster_name = cluster_name,
      ...
    )
    
    all_results[[as.character(cid)]] <- result
  }
  
  return(all_results)
}

# =============================================================================
# 6. Session Info (for reproducibility)
# =============================================================================
print(sessionInfo())
