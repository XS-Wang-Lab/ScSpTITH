#' ScSpTITH: A Rank-Correlation Framework for Robust Quantification of Multi-Dimensional Tumor Heterogeneity (ScSpTITHscore)
#'
#' @param data Gene expression matrix with rows as genes and columns as cells/samples
#' @param meta Metadata containing at minimum:
#' #'     - First column: cell/sample IDs (must match column names in data)
#'        - Patient ID column (required)
#'        - Condition/cluster column (optional, only required for cluster mode)
#' @param mode Calculation mode: "overall" - calculate for entire sample;
#'                               "cluster" - calculate separately for each cluster
#' @param patient_col Column name in metadata containing patient IDs, default "Patient"
#' @param condition_col Column name in metadata containing cluster/group information, default "condition"
#' @param top_n_genes Number of highly variable genes used for correlation calculation, default 5000
#' @param sample_fraction Sampling fraction, default 1.0 means use all cells
#'
#' @returns A data frame containing Patient,ScSpTITHscore, Cluster (mode:"cluster")
#' @importFrom stats sd cor median
#' @export
#'
#' @examples
#' # Mode 1: Calculate ITHscore for entire sample
#' #result1 <- ScSpTITHscore(data, meta, mode = "overall",patient_col = "Patient")
#' # Mode 2: Calculate ITHscore separately for each cluster
#' #result2 <- ScSpTITHscore(data_count, meta, mode = "cluster",patient_col = "Patient",condition_col = "condition")

ScSpTITHscore <- function(data,
                       meta,
                       mode = c("overall","cluster"),
                       patient_col = "Patient",
                       condition_col = "condition",
                       top_n_genes = 5000,
                       sample_fraction = 1.0) {# Sampling fraction, default 1.0 means use all cells

  # Parameter validation
  mode <- match.arg(mode)
  
  # set sample fraction
  if (!is.numeric(sample_fraction) || sample_fraction <= 0 || sample_fraction > 1) {
    stop("sample_fraction must be a number between 0 and 1")
  }

  if (!(is.matrix(data) || is.data.frame(data))) {
    stop("data must be a matrix or data.frame")
  }
  if (!is.data.frame(meta)) {
    stop("meta must be a data.frame")
  }
  if (!is.numeric(top_n_genes) || top_n_genes <= 0) {
    stop("top_n_genes must be a positive integer greater than 5")
  }

  # Check required columns
  required_cols <- c(patient_col, "cell_name", if(mode == "cluster") condition_col)
  missing_cols <- required_cols[!required_cols %in% colnames(meta)]

  if (length(missing_cols) > 0) {
    stop("Required columns missing in meta: ", paste(missing_cols, collapse = ", "))
  }

  # Main calculation function
  calculate_ITHscore <- function(data, meta, patient_id, group_id = NULL) {
    # Filter cells for current patient
    if (is.null(group_id)) {
      # Overall mode: select all cells from the patient
      meta_sub <- meta[meta[[patient_col]] == patient_id, ]
    } else {
      # Cluster mode: select cells from specific cluster of the patient
      meta_sub <- meta[meta[[patient_col]] == patient_id & meta[[condition_col]] == group_id, ]
    }

    # Extract cell IDs
    cell_ids <- meta_sub[["cell_name"]]
    count_sub <- data[, colnames(data) %in% cell_ids, drop = FALSE]

    # Skip if insufficient cells
    if (ncol(count_sub) < 2) {
      warning(paste0("Skip ", patient_id,
                     if(!is.null(group_id)) paste0("-", group_id),
                     ": <2 cells"))
      return(NULL)
    }
    
    # sample
    n_cells <- ncol(count_sub)
      n_sample <- round(n_cells * sample_fraction)
      n_sample <- max(2, min(n_sample, n_cells))
      
      set.seed(123)  
      sampled_cells <- sample(colnames(count_sub), size = n_sample, replace = FALSE)
      count_sub <- count_sub[, sampled_cells, drop = FALSE]
      
    #   message(paste0("Sampled ", n_sample, " cells from ", n_cells, 
    #                  " cells for ", patient_id,
    #                  if(!is.null(group_id)) paste0("-", group_id)))

    # Calculate standard deviation and select highly variable genes
    gene_sd <- apply(count_sub, 1, function(x) sd(x, na.rm = TRUE))
    high_var_genes <- names(sort(gene_sd, decreasing = TRUE))[1:min(top_n_genes, length(gene_sd))]

    # Extract expression matrix of highly variable genes
    exp_high_var <- count_sub[rownames(count_sub) %in% high_var_genes, , drop = FALSE]
    exp_high_var <- apply(exp_high_var, 2, as.numeric)
    rownames(exp_high_var) <- high_var_genes

    # Calculate correlation matrix
    spearman_corr_matrix <- cor(exp_high_var, method = "spearman")

    # Calculate median correlation for each cell with other cells
    cell_medians <- apply(spearman_corr_matrix, 1, function(row) {
      median(row[-which.max(row)], na.rm = TRUE)  # 
    })

    # Calculate overall median and ITHscore
    overall_median <- median(cell_medians, na.rm = TRUE)
    ScSpTITH_score <- 1 - overall_median

    result <- data.frame(
      Patient = patient_id,
      ScSpTITHscore = ScSpTITH_score 
    )

    if (!is.null(group_id)) {
      result$Cluster <- group_id
    }

    return(result)
  }

  # Get all patient IDs
  patient_ids <- unique(meta[[patient_col]])
  all_results <- data.frame()

  for (patient_id in patient_ids) {
    if (mode == "overall") {
      # Overall mode: calculate once per patient
      result <- calculate_ITHscore(data, meta, patient_id)
      if (!is.null(result)) {
        all_results <- rbind(all_results, result)
      }
    } else {
      # Cluster mode: calculate for each cluster per patient
      meta_patient <- meta[meta[[patient_col]] == patient_id, ]

      # Get all clusters and their cell counts for current patient
      cluster_counts <- table(meta_patient[[condition_col]])
      valid_clusters <- names(cluster_counts[cluster_counts > 2])

      if (length(valid_clusters) == 0) {
        message("Skipping patient ", patient_id, ": no valid clusters")
        next
      }

      for (cluster_id in valid_clusters) {
        result <- calculate_ITHscore(data, meta, patient_id, cluster_id)
        if (!is.null(result)) {
          all_results <- rbind(all_results, result)
        }
      }
    }
  }

  return(all_results)
}
