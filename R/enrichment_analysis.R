#' Perform driver gene enrichment analysis for chromoanagenesis events
#'
#' @param cohort_summary A list containing summarized data for each sample, 
#'                       including detected mechanisms and affected genes.
#' @param target_mechanism The mechanism to analyze (e.g., "CT_High", "CP_Likely").
#' @return A data frame containing Fisher's exact test results for each gene.
#' @export
calculate_gene_enrichment <- function(cohort_summary, target_mechanism = "CT_High") {
  
  # 1. Flatten results to sample-level data frame
  all_samples <- unlist(cohort_summary, recursive = FALSE)
  sample_ids <- names(all_samples)
  
  if (length(all_samples) == 0) stop("Cohort summary is empty")
  
  # Determine which samples have the mechanism
  has_mech <- sapply(all_samples, function(x) {
    val <- x$mechs[[target_mechanism]]
    if (is.null(val)) return(FALSE)
    return(val > 0)
  })
  
  mech_positive_samples <- sample_ids[has_mech]
  mech_negative_samples <- sample_ids[!has_mech]
  
  n_pos <- length(mech_positive_samples)
  n_neg <- length(mech_negative_samples)
  
  if (n_pos == 0) {
    warning(paste("No samples found with mechanism:", target_mechanism))
    return(NULL)
  }
  
  # 2. Get the list of all disrupted driver genes
  all_affected_genes <- unique(unlist(lapply(all_samples, function(x) x$genes)))
  
  if (length(all_affected_genes) == 0) return(NULL)
  
  # 3. Perform Fisher's Exact Test for each gene
  results <- lapply(all_affected_genes, function(gene) {
    # Count occurrences
    a <- sum(sapply(all_samples[mech_positive_samples], function(x) gene %in% x$genes))
    b <- sum(sapply(all_samples[mech_negative_samples], function(x) gene %in% x$genes))
    
    # Contingency Table:
    #            Mech+ | Mech-
    # Gene+ |    a    |   b
    # Gene- | n_pos-a | n_neg-b
    
    mat <- matrix(c(a, b, n_pos - a, n_neg - b), nrow = 2)
    ft <- fisher.test(mat, alternative = "greater")
    
    data.frame(
      Gene = gene,
      Mech_Pos_Hits = a,
      Mech_Neg_Hits = b,
      Mech_Pos_Freq = a / n_pos,
      Mech_Neg_Freq = b / n_neg,
      OddsRatio = ft$estimate,
      PValue = ft$p.value,
      stringsAsFactors = FALSE
    )
  })
  
  enrich_df <- rbindlist(results)
  if (nrow(enrich_df) > 0) {
    enrich_df$FDR <- p.adjust(enrich_df$PValue, method = "BH")
    return(enrich_df[order(enrich_df$PValue), ])
  }
  
  return(enrich_df)
}
