################################################################################
# Monocle2 Trajectory Analysis for Silkworm Brain snRNA-seq Data
# 
# This script performs pseudotime trajectory analysis using Monocle2
# on single-nucleus transcriptomic data of the silkworm brain across
# metamorphic stages.
#
# Input: Seurat object with annotated cell types and stage metadata
# Output: Pseudotime trajectories, branch analysis, gene expression heatmaps
#
# Author: Jiahao Xiang
# Date: 2026-07-24
################################################################################

# =============================================================================
# 0. Load Libraries
# =============================================================================
library(monocle)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggsci)
library(viridis)

# =============================================================================
# 1. Prepare Monocle Object from Seurat Object
# =============================================================================

#' Convert Seurat object to Monocle CellDataSet
#' 
#' @param seurat_obj A Seurat object
#' @param assay Assay to use (default: "RNA")
#' @param slot Slot to use (default: "counts")
#' @return A Monocle CellDataSet object
create_monocle_object <- function(seurat_obj, assay = "RNA", slot = "counts") {
  # Extract count matrix
  data <- GetAssayData(seurat_obj, assay = assay, slot = slot)
  
  # Extract metadata
  cell_metadata <- seurat_obj@meta.data
  
  # Create gene annotation (required for Monocle)
  gene_annotation <- data.frame(
    gene_short_name = row.names(seurat_obj),
    row.names = row.names(seurat_obj)
  )
  
  # Create AnnotatedDataFrame objects
  pd <- new("AnnotatedDataFrame", data = cell_metadata)
  fd <- new("AnnotatedDataFrame", data = gene_annotation)
  
  # Create CellDataSet
  sce <- newCellDataSet(data, phenoData = pd, featureData = fd)
  
  # Estimate size factors and dispersion
  sce <- estimateSizeFactors(sce)
  sce <- estimateDispersions(sce)
  
  return(sce)
}

# Usage:
# sce <- create_monocle_object(ec)

# =============================================================================
# 2. Select Ordering Genes
# =============================================================================

#' Select ordering genes for trajectory inference
#' 
#' Strategy: Combine marker genes from developmental stages (FindAllMarkers)
#' and optionally high-variance genes.
#' 
#' @param seurat_obj Seurat object with cell types annotated
#' @param sce Monocle CellDataSet object
#' @param p_val_adj_threshold Adjusted p-value threshold for marker selection
#' @param logfc_threshold Log fold change threshold for marker selection
#' @param min_pct Minimum percentage of cells expressing the gene
#' @param use_highvar Whether to also include Seurat variable features
#' @return Character vector of ordering genes
select_ordering_genes <- function(seurat_obj, sce, 
                                  p_val_adj_threshold = 0.05,
                                  logfc_threshold = 0.25,
                                  min_pct = 0.25,
                                  use_highvar = FALSE) {
  
  cat("Finding marker genes for each developmental stage...\n")
  
  # Find marker genes (stage-specific)
  markers <- FindAllMarkers(seurat_obj, 
                            only.pos = TRUE, 
                            min.pct = min_pct, 
                            logfc.threshold = logfc_threshold)
  
  # Filter significant markers
  sig_markers <- subset(markers, p_val_adj < p_val_adj_threshold)
  stage_marker_genes <- unique(sig_markers$gene)
  cat("  - Stage-specific marker genes:", length(stage_marker_genes), "\n")
  
  # Optionally add Seurat variable features
  if (use_highvar) {
    highvar_genes <- seurat_obj@assays[["RNA"]]@var.features
    cat("  - High-variable genes:", length(highvar_genes), "\n")
    ordering_genes <- unique(c(stage_marker_genes, highvar_genes))
  } else {
    ordering_genes <- stage_marker_genes
  }
  
  cat("  - Total ordering genes:", length(ordering_genes), "\n")
  
  # Also add specific genes of interest (e.g., EcR, Sox14, Utx)
  custom_genes <- c("EcR", "Sox14", "Utx")
  custom_genes <- custom_genes[custom_genes %in% row.names(seurat_obj)]
  if (length(custom_genes) > 0) {
    ordering_genes <- unique(c(ordering_genes, custom_genes))
    cat("  - Added custom genes:", paste(custom_genes, collapse = ", "), "\n")
  }
  
  return(ordering_genes)
}

# Usage:
# ordering_genes <- select_ordering_genes(ec, sce)
# sce <- setOrderingFilter(sce, ordering_genes)

# =============================================================================
# 3. Dimensionality Reduction and Cell Ordering
# =============================================================================

#' Perform dimensionality reduction and order cells
#' 
#' @param sce Monocle CellDataSet with ordering genes set
#' @param max_components Number of components for DDRTree (default: 2)
#' @param method Dimensionality reduction method (default: "DDRTree")
#' @return Monocle CellDataSet with cells ordered
run_trajectory <- function(sce, max_components = 2, method = "DDRTree") {
  sce <- reduceDimension(sce,
                         max_components = max_components,
                         method = method)
  sce <- orderCells(sce)
  return(sce)
}

# Usage:
# sce <- run_trajectory(sce)

# =============================================================================
# 4. Visualization Functions
# =============================================================================

#' Plot trajectory with specified coloring
#' 
#' @param sce Monocle CellDataSet
#' @param color_by Column name in pData to color by (e.g., "formal", "Pseudotime", "State")
#' @param show_branch_points Whether to show branch points
#' @param cell_size Size of cells in plot
#' @param save_path Optional file path to save plot
#' @return ggplot object
plot_trajectory <- function(sce, color_by, show_branch_points = FALSE, 
                            cell_size = 1, save_path = NULL) {
  p <- plot_cell_trajectory(sce, 
                            color_by = color_by, 
                            show_branch_points = show_branch_points,
                            cell_size = cell_size)
  
  if (color_by == "Pseudotime") {
    p <- p + scale_color_viridis_c(direction = -1)
  }
  
  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 8, height = 6)
    cat("Plot saved to:", save_path, "\n")
  }
  
  return(p)
}

#' Plot gene expression along pseudotime
#' 
#' @param sce Monocle CellDataSet
#' @param gene Gene name to plot
#' @param color_by Color by "Pseudotime" or "State" (default: "Pseudotime")
#' @param save_path Optional file path to save plot
#' @return ggplot object
plot_gene_pseudotime <- function(sce, gene, color_by = "Pseudotime", save_path = NULL) {
  # Add gene expression to pData
  pData(sce)[[gene]] <- log2(exprs(sce)[gene, ] + 1)
  
  p <- plot_cell_trajectory(sce, color_by = gene, cell_size = 1) +
    scale_color_viridis_c(direction = -1)
  
  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 6, height = 4)
    cat("Plot saved to:", save_path, "\n")
  }
  
  return(p)
}

#' Plot gene expression change along pseudotime (line plot)
#' 
#' @param sce Monocle CellDataSet
#' @param gene Gene name
#' @param save_path Optional file path to save plot
#' @return ggplot object
plot_gene_change <- function(sce, gene, save_path = NULL) {
  cds_subset <- sce[gene, ]
  p <- plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime") +
    scale_color_viridis_c(direction = -1)
  
  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 5, height = 2.5)
    cat("Plot saved to:", save_path, "\n")
  }
  
  return(p)
}

#' Generate pseudotime heatmap
#' 
#' @param sce Monocle CellDataSet
#' @param sig_genes Significant genes for heatmap
#' @param num_clusters Number of clusters for gene grouping
#' @param save_path Optional file path to save plot
#' @return List containing heatmap and cluster assignments
generate_pseudotime_heatmap <- function(sce, sig_genes, num_clusters = 3, 
                                        save_path = NULL) {
  # Filter sce to only significant genes
  sce_subset <- sce[sig_genes, ]
  
  p <- plot_pseudotime_heatmap(sce_subset,
                               num_clusters = num_clusters,
                               cores = 1,
                               show_rownames = FALSE,
                               return_heatmap = TRUE)
  
  if (!is.null(save_path)) {
    ggsave(save_path, p, height = 5, width = 4)
    cat("Heatmap saved to:", save_path, "\n")
  }
  
  return(p)
}

#' Perform branch expression analysis (BEAM)
#' 
#' @param sce Monocle CellDataSet with ordered cells
#' @param branch_point Branch point number
#' @param qval_threshold Q-value threshold for significance
#' @param num_clusters Number of gene clusters in heatmap
#' @param save_prefix Prefix for output files
#' @return List containing BEAM results and heatmap
run_branch_analysis <- function(sce, branch_point = 2, 
                                qval_threshold = 0.05,
                                num_clusters = 2,
                                save_prefix = "branch") {
  
  cat("Running BEAM analysis for branch point:", branch_point, "\n")
  
  # Run BEAM
  BEAM_res <- BEAM(sce, branch_point = branch_point, cores = 1)
  BEAM_res <- BEAM_res[order(BEAM_res$qval), ]
  BEAM_res <- BEAM_res[, c("gene_short_name", "pval", "qval")]
  
  # Save results
  saveRDS(BEAM_res, file = paste0(save_prefix, "_", branch_point, ".rds"))
  
  # Generate heatmap for significant genes
  sig_genes <- row.names(subset(BEAM_res, qval < qval_threshold))
  cat("  Significant genes:", length(sig_genes), "\n")
  
  if (length(sig_genes) > 0) {
    p <- plot_genes_branched_heatmap(
      sce[sig_genes, ],
      branch_point = branch_point,
      num_clusters = num_clusters,
      cores = 1,
      use_gene_short_name = TRUE,
      show_rownames = FALSE,
      return_heatmap = TRUE
    )
    
    ggsave(paste0(save_prefix, "_", branch_point, ".pdf"), 
           p$ph_res, width = 6.5, height = 10)
    
    # Extract cluster assignments
    cluster_assignments <- p$annotation_row
    write.csv(cluster_assignments, 
              file = paste0(save_prefix, "_", branch_point, "_clusters.csv"))
    
    return(list(BEAM_res = BEAM_res, heatmap = p, clusters = cluster_assignments))
  } else {
    warning("No significant genes found for branch point ", branch_point)
    return(list(BEAM_res = BEAM_res, heatmap = NULL, clusters = NULL))
  }
}

# =============================================================================
# 5. Main Analysis Pipeline
# =============================================================================

#' Run full Monocle2 trajectory analysis pipeline
#' 
#' @param seurat_obj Seurat object with cell type annotations
#' @param output_dir Directory to save outputs
#' @param stage_metadata Column name in metadata for developmental stages
#' @return Monocle CellDataSet object with trajectory analysis completed
run_full_analysis <- function(seurat_obj, output_dir = "./monocle_output/", 
                              stage_metadata = "formal") {
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  setwd(output_dir)
  
  # Step 1: Create Monocle object
  cat("Step 1: Creating Monocle object...\n")
  sce <- create_monocle_object(seurat_obj)
  saveRDS(sce, "monocle_object_init.rds")
  
  # Step 2: Select ordering genes
  cat("Step 2: Selecting ordering genes...\n")
  ordering_genes <- select_ordering_genes(seurat_obj, sce)
  sce <- setOrderingFilter(sce, ordering_genes)
  saveRDS(sce, "monocle_object_filtered.rds")
  
  # Step 3: Run trajectory
  cat("Step 3: Running trajectory inference...\n")
  sce <- run_trajectory(sce)
  saveRDS(sce, "monocle_object_trajectory.rds")
  
  # Step 4: Plot trajectory with different colorings
  cat("Step 4: Generating trajectory plots...\n")
  plot_trajectory(sce, "State", save_path = "trajectory_state.pdf")
  plot_trajectory(sce, stage_metadata, save_path = "trajectory_stage.pdf")
  plot_trajectory(sce, "Pseudotime", save_path = "trajectory_pseudotime.pdf")
  
  # Facet by stage
  p <- plot_cell_trajectory(sce, color_by = "Pseudotime", 
                            show_branch_points = FALSE) +
    facet_wrap(~formal, nrow = 1) +
    scale_color_viridis_c(direction = -1)
  ggsave("trajectory_stage_facet.pdf", p, height = 6, width = 20)
  
  # Step 5: Pseudotime heatmap
  cat("Step 5: Generating pseudotime heatmap...\n")
  diff_res <- differentialGeneTest(sce[ordering_genes, ],
                                   fullModelFormulaStr = "~sm.ns(Pseudotime)")
  sig_genes <- row.names(subset(diff_res, qval < 0.01))
  
  if (length(sig_genes) > 0) {
    generate_pseudotime_heatmap(sce, sig_genes, 
                                num_clusters = 3,
                                save_path = "pseudotime_heatmap.pdf")
  }
  
  # Step 6: Plot specific genes of interest
  cat("Step 6: Plotting genes of interest...\n")
  genes_of_interest <- c("EcR", "Sox14", "Utx")
  for (gene in genes_of_interest) {
    if (gene %in% row.names(sce)) {
      plot_gene_pseudotime(sce, gene, save_path = paste0(gene, "_trajectory.pdf"))
      plot_gene_change(sce, gene, save_path = paste0(gene, "_change.pdf"))
    }
  }
  
  # Step 7: Branch analysis (if branches exist)
  cat("Step 7: Running branch analysis...\n")
  # Note: You may need to manually set root_state based on your biology
  # sce <- orderCells(sce, root_state = 1)  # Uncomment and adjust as needed
  run_branch_analysis(sce, branch_point = 2, save_prefix = "branch")
  
  # Save final object
  saveRDS(sce, "monocle_object_final.rds")
  
  cat("Analysis complete! All outputs saved to:", output_dir, "\n")
  
  return(sce)
}

# =============================================================================
# 6. Utility Functions
# =============================================================================

#' Add custom metadata to Monocle object from Seurat
#' 
#' @param sce Monocle CellDataSet
#' @param seurat_obj Seurat object with additional metadata
#' @param meta_cols Columns to transfer
#' @return Updated Monocle CellDataSet
add_metadata <- function(sce, seurat_obj, meta_cols = NULL) {
  if (is.null(meta_cols)) {
    meta_cols <- colnames(seurat_obj@meta.data)
  }
  
  for (col in meta_cols) {
    if (col %in% colnames(seurat_obj@meta.data)) {
      sce@phenoData@data[[col]] <- seurat_obj@meta.data[rownames(sce@phenoData@data), col]
    }
  }
  
  return(sce)
}

#' Get cell counts per cluster/condition
#' 
#' @param seurat_obj Seurat object
#' @param group_by Column name for grouping (e.g., "stage_cluster")
#' @return Data frame with group and cell counts
get_cell_counts <- function(seurat_obj, group_by = "stage_cluster") {
  counts <- table(seurat_obj@meta.data[[group_by]])
  result <- data.frame(
    group = names(counts),
    count = as.numeric(counts)
  )
  return(result)
}

# =============================================================================
# 7. Execute Analysis
# =============================================================================

# Uncomment below to run the full analysis

# # Load your Seurat object
# ec <- readRDS("path_to_your_seurat_object.rds")
# 
# # Run full analysis pipeline
# sce <- run_full_analysis(ec, output_dir = "./monocle_output/")

# =============================================================================
# 8. Session Info (for reproducibility)
# =============================================================================

sessionInfo()
