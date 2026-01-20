#' Plot Mechanism Prevalence in a Cohort
#'
#' Creates a publication-quality bar plot showing the percentage of samples
#' in the cohort affected by each chromoanagenesis mechanism.
#'
#' @param summary_df Data frame from summarize_cohort_results().
#' @param title Plot title.
#' @return A ggplot object.
#' @export
plot_mechanism_prevalence <- function(summary_df, title = "Prevalence of Chromoanagenesis Mechanisms") {
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
    if (!requireNamespace("tidyr", quietly = TRUE)) stop("tidyr is required.")
    
    n_samples <- nrow(summary_df)
    
    # Calculate counts
    counts <- data.frame(
        Mechanism = c("Chromothripsis (HC)", "Chromoplexy (Likely)", "Chromoanasynthesis (Likely)"),
        Count = c(
            sum(summary_df$ct_hc > 0, na.rm = TRUE),
            sum(summary_df$cp_likely > 0, na.rm = TRUE),
            sum(summary_df$cs_likely > 0, na.rm = TRUE)
        )
    )
    
    counts$Percentage <- (counts$Count / n_samples) * 100
    counts$Label <- sprintf("%d (%.1f%%)", counts$Count, counts$Percentage)
    
    p <- ggplot2::ggplot(counts, ggplot2::aes(x = Mechanism, y = Percentage, fill = Mechanism)) +
        ggplot2::geom_bar(stat = "identity", width = 0.6, show.legend = FALSE) +
        ggplot2::geom_text(ggplot2::aes(label = Label), vjust = -0.5, size = 4) +
        ggplot2::scale_fill_manual(values = c("#d32f2f", "#1976d2", "#388e3c")) +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = title, 
                      subtitle = paste("Total Cohort size: n =", n_samples),
                      y = "Samples Affected (%)", x = "") +
        ggplot2::ylim(0, 100) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(face = "bold"))
        
    return(p)
}

#' Plot Cohort Landscape (OncoPrint-style)
#'
#' Generates a comprehensive heatmap showing the distribution of CT, CP, and CS 
#' events across all samples, sorted by complexity.
#'
#' @param summary_df Data frame from summarize_cohort_results().
#' @return A ggplot object.
#' @export
plot_cohort_summary_landscape <- function(summary_df) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
    if (!requireNamespace("reshape2", quietly = TRUE)) stop("reshape2 is required.")
    
    # 1. Prepare data for heatmap
    # We create a binary matrix for mechanisms
    plot_data <- summary_df
    plot_data$CT <- ifelse(plot_data$ct_hc > 0, 1, 0)
    plot_data$CP <- ifelse(plot_data$cp_likely > 0, 1, 0)
    plot_data$CS <- ifelse(plot_data$cs_likely > 0, 1, 0)
    
    # Sort samples by number of mechanisms and then by total HC CT events
    plot_data$total_score <- (plot_data$CT + plot_data$CP + plot_data$CS) * 1000 + plot_data$ct_hc
    plot_data <- plot_data[order(plot_data$total_score, decreasing = TRUE), ]
    plot_data$sample_id <- factor(plot_data$sample_id, levels = plot_data$sample_id)
    
    # Melt for ggplot
    m_data <- reshape2::melt(plot_data[, c("sample_id", "CT", "CP", "CS")], id.vars = "sample_id")
    colnames(m_data) <- c("Sample", "Mechanism", "Present")
    
    # 2. Main Heatmap
    p1 <- ggplot2::ggplot(m_data, ggplot2::aes(x = Sample, y = Mechanism, fill = factor(Present))) +
        ggplot2::geom_tile(color = "white", linewidth = 0.5) +
        ggplot2::scale_fill_manual(values = c("0" = "#f5f5f5", "1" = "#424242"), 
                                  labels = c("Absent", "Detected"), name = "Status") +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = "Chromoanagenesis Landscape", x = "", y = "") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 8),
                      legend.position = "bottom")
    
    # 3. Barplot for Fusion Counts (Top annotation)
    p2 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = sample_id, y = n_fusions)) +
        ggplot2::geom_bar(stat = "identity", fill = "#FFC107") +
        ggplot2::theme_minimal() +
        ggplot2::labs(y = "Fusions", x = "") +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                      panel.grid.major.x = ggplot2::element_blank())
    
    if (requireNamespace("patchwork", quietly = TRUE)) {
        return(p2 / p1 + patchwork::plot_layout(heights = c(1, 2)))
    } else {
        return(p1) # Fallback to just heatmap
    }
}

#' Plot Cancer Gene Impact Landscape
#'
#' Analyzes which cancer driver genes are most frequently affected by 
#' chromoanagenesis mechanisms across the cohort.
#'
#' @param results_list A list of chromoanagenesis result objects.
#' @param top_n Number of top genes to show (default: 20).
#' @return A ggplot object.
#' @export
plot_gene_impact_landscape <- function(results_list, top_n = 20) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
    
    # 1. Extract all driver hits across samples
    all_hits_list <- lapply(names(results_list), function(sid) {
        res <- results_list[[sid]]
        # Check both direct annotation slot and sub-result
        hits <- if (!is.null(res$annotation$driver_hits)) res$annotation$driver_hits else NULL
        if (is.null(hits) && !is.null(res$results$annotation$driver_hits)) hits <- res$results$annotation$driver_hits
        
        if (!is.null(hits) && nrow(hits) > 0) {
            hits$sample_id <- sid
            return(hits[, c("sample_id", "Gene", "Mechanism")])
        }
        return(NULL)
    })
    
    all_hits <- do.call(rbind, all_hits_list)
    if (is.null(all_hits) || nrow(all_hits) == 0) {
        message("No driver gene hits found in the provided results.")
        return(NULL)
    }
    
    # 2. Calculate frequencies
    gene_stats <- as.data.frame(table(all_hits$Gene, all_hits$Mechanism))
    colnames(gene_stats) <- c("Gene", "Mechanism", "Count")
    
    # Calculate total hits per gene for sorting
    gene_totals <- aggregate(Count ~ Gene, data = gene_stats, sum)
    top_genes <- gene_totals$Gene[order(gene_totals$Count, decreasing = TRUE)][1:min(nrow(gene_totals), top_n)]
    
    plot_df <- gene_stats[gene_stats$Gene %in% top_genes, ]
    plot_df$Gene <- factor(plot_df$Gene, levels = rev(top_genes))
    
    # 3. Create Plot
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Gene, y = Count, fill = Mechanism)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::coord_flip() +
        ggplot2::scale_fill_manual(values = c("Chromothripsis" = "#d32f2f", 
                                             "Chromoplexy" = "#1976d2", 
                                             "Chromoanasynthesis" = "#388e3c")) +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = "Top Impacted Cancer Driver Genes",
                      subtitle = "Frequency of genes affected by complex rearrangements",
                      x = "", y = "Number of Samples Affected") +
        ggplot2::theme(axis.text.y = ggplot2::element_text(face = "bold"))
        
    return(p)
}
