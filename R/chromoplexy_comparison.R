################################################################################
# Chromoplexy Method Comparison Framework
#
# Compare ChainFinder-style (statistical adjacency) vs OncoImplexus (distance-based)
# chromoplexy detection methods
################################################################################

#' Compare chromoplexy detection methods
#'
#' Runs both ChainFinder-style and OncoImplexus detection methods on the same
#' data and provides a detailed comparison of results.
#'
#' @param SV.sample SV data
#' @param CNV.sample CNV data (optional)
#' @param min_chromosomes Minimum chromosomes (default: 3)
#' @param min_translocations Minimum translocations (default: 3)
#' @param genome Reference genome (default: "hg19")
#' @param verbose Print progress (default: TRUE)
#' @return Comparison results list
#' @export
compare_chromoplexy_methods <- function(SV.sample,
                                        CNV.sample = NULL,
                                        min_chromosomes = 3,
                                        min_translocations = 3,
                                        genome = "hg19",
                                        verbose = TRUE) {

    if (verbose) {
        cat("\n")
        cat("########################################################################\n")
        cat("#  CHROMOPLEXY METHOD COMPARISON                                      #\n")
        cat("#  ChainFinder-style (Statistical) vs OncoImplexus (Distance-based)   #\n")
        cat("########################################################################\n\n")
    }

    # Run ChainFinder-style detection
    if (verbose) cat("Running ChainFinder-style detection...\n")
    t1 <- Sys.time()
    result_chainfinder <- tryCatch({
        detect_chromoplexy_chainfinder(
            SV.sample = SV.sample,
            CNV.sample = CNV.sample,
            min_chromosomes = min_chromosomes,
            min_translocations = min_translocations,
            genome = genome,
            verbose = FALSE
        )
    }, error = function(e) {
        if (verbose) cat(sprintf("  Error: %s\n", e$message))
        return(NULL)
    })
    t1_elapsed <- difftime(Sys.time(), t1, units = "secs")

    # Run OncoImplexus detection
    if (verbose) cat("Running OncoImplexus detection...\n")
    t2 <- Sys.time()
    result_oncoimplexus <- tryCatch({
        detect_chromoplexy(
            SV.sample = SV.sample,
            CNV.sample = CNV.sample,
            min_chromosomes = min_chromosomes,
            min_translocations = min_translocations,
            genome = genome,
            verbose = FALSE
        )
    }, error = function(e) {
        if (verbose) cat(sprintf("  Error: %s\n", e$message))
        return(NULL)
    })
    t2_elapsed <- difftime(Sys.time(), t2, units = "secs")

    # Compare results
    comparison <- compare_results(
        result_chainfinder,
        result_oncoimplexus,
        verbose = verbose
    )

    comparison$timing <- list(
        chainfinder_seconds = as.numeric(t1_elapsed),
        oncoimplexus_seconds = as.numeric(t2_elapsed)
    )

    comparison$result_chainfinder <- result_chainfinder
    comparison$result_oncoimplexus <- result_oncoimplexus

    if (verbose) {
        cat("\n")
        cat("========================================================================\n")
        cat("COMPARISON SUMMARY\n")
        cat("========================================================================\n\n")
        print_comparison_summary(comparison)
    }

    class(comparison) <- c("chromoplexy_comparison", "list")
    return(comparison)
}


#' Compare two chromoplexy detection results
#' @keywords internal
compare_results <- function(result_cf, result_oi, verbose = TRUE) {

    comparison <- list()

    # Basic counts
    comparison$chainfinder <- list(
        total_chains = if (!is.null(result_cf)) result_cf$total_chains else 0,
        likely = if (!is.null(result_cf)) result_cf$likely_chromoplexy else 0,
        possible = if (!is.null(result_cf)) result_cf$possible_chromoplexy else 0
    )

    comparison$oncoimplexus <- list(
        total_chains = if (!is.null(result_oi)) result_oi$total_chains else 0,
        likely = if (!is.null(result_oi)) result_oi$likely_chromoplexy else 0,
        possible = if (!is.null(result_oi)) result_oi$possible_chromoplexy else 0
    )

    # Graph structure comparison
    if (!is.null(result_cf) && !is.null(result_cf$translocation_graph)) {
        cf_edges <- result_cf$translocation_graph$edges
        comparison$chainfinder$n_nodes <- length(result_cf$translocation_graph$nodes)
        comparison$chainfinder$n_edges <- nrow(cf_edges)
        comparison$chainfinder$n_translocation_edges <- sum(cf_edges$edge_type == "TRANSLOCATION")
        comparison$chainfinder$n_deletion_bridge_edges <- sum(cf_edges$edge_type == "DELETION_BRIDGE")
        comparison$chainfinder$n_adjacency_edges <- sum(cf_edges$edge_type == "STAT_ADJACENCY")
    }

    if (!is.null(result_oi) && !is.null(result_oi$translocation_graph)) {
        oi_edges <- result_oi$translocation_graph$edges
        comparison$oncoimplexus$n_nodes <- length(result_oi$translocation_graph$nodes)
        comparison$oncoimplexus$n_edges <- nrow(oi_edges)
        comparison$oncoimplexus$n_translocation_edges <- sum(oi_edges$edge_type == "TRANSLOCATION")
        comparison$oncoimplexus$n_deletion_bridge_edges <- sum(oi_edges$edge_type == "DELETION_BRIDGE")
        comparison$oncoimplexus$n_adjacency_edges <- sum(oi_edges$edge_type == "ADJACENCY")
    }

    # Chain overlap analysis
    if (!is.null(result_cf) && !is.null(result_oi) &&
        result_cf$total_chains > 0 && result_oi$total_chains > 0) {

        comparison$overlap <- analyze_chain_overlap(
            result_cf$chains,
            result_oi$chains
        )
    } else {
        comparison$overlap <- list(
            n_shared = 0,
            n_chainfinder_only = comparison$chainfinder$total_chains,
            n_oncoimplexus_only = comparison$oncoimplexus$total_chains,
            jaccard_index = 0
        )
    }

    # Classification concordance
    if (!is.null(result_cf) && !is.null(result_oi) &&
        nrow(result_cf$summary) > 0 && nrow(result_oi$summary) > 0) {

        comparison$concordance <- analyze_classification_concordance(
            result_cf$summary,
            result_oi$summary
        )
    }

    return(comparison)
}


#' Analyze overlap between chains detected by two methods
#' @keywords internal
analyze_chain_overlap <- function(chains_cf, chains_oi) {

    # Get SV indices for each chain
    cf_sv_sets <- lapply(chains_cf, function(x) sort(x$sv_indices))
    oi_sv_sets <- lapply(chains_oi, function(x) sort(x$sv_indices))

    # Find matching chains (same or >50% overlap in SVs)
    cf_matched <- rep(FALSE, length(cf_sv_sets))
    oi_matched <- rep(FALSE, length(oi_sv_sets))

    for (i in seq_along(cf_sv_sets)) {
        for (j in seq_along(oi_sv_sets)) {
            if (oi_matched[j]) next

            # Calculate overlap
            intersection <- length(intersect(cf_sv_sets[[i]], oi_sv_sets[[j]]))
            union_size <- length(union(cf_sv_sets[[i]], oi_sv_sets[[j]]))

            if (union_size > 0) {
                overlap_ratio <- intersection / union_size

                if (overlap_ratio > 0.5) {
                    cf_matched[i] <- TRUE
                    oi_matched[j] <- TRUE
                    break
                }
            }
        }
    }

    n_shared <- sum(cf_matched)
    n_cf_only <- sum(!cf_matched)
    n_oi_only <- sum(!oi_matched)

    # Jaccard index
    total_unique <- n_shared + n_cf_only + n_oi_only
    jaccard <- if (total_unique > 0) n_shared / total_unique else 0

    return(list(
        n_shared = n_shared,
        n_chainfinder_only = n_cf_only,
        n_oncoimplexus_only = n_oi_only,
        jaccard_index = jaccard,
        cf_matched_indices = which(cf_matched),
        oi_matched_indices = which(oi_matched)
    ))
}


#' Analyze classification concordance
#' @keywords internal
analyze_classification_concordance <- function(summary_cf, summary_oi) {

    # Count by classification
    cf_likely <- sum(summary_cf$classification == "Likely chromoplexy")
    cf_possible <- sum(summary_cf$classification == "Possible chromoplexy")

    oi_likely <- sum(summary_oi$classification == "Likely chromoplexy")
    oi_possible <- sum(summary_oi$classification == "Possible chromoplexy")

    return(list(
        chainfinder_likely = cf_likely,
        chainfinder_possible = cf_possible,
        oncoimplexus_likely = oi_likely,
        oncoimplexus_possible = oi_possible
    ))
}


#' Print comparison summary
#' @keywords internal
print_comparison_summary <- function(comparison) {

    cat("Detection Results:\n")
    cat("------------------------------------------------------------------------\n")
    cat(sprintf("                          ChainFinder    OncoImplexus\n"))
    cat(sprintf("Total chains:             %11d    %12d\n",
               comparison$chainfinder$total_chains,
               comparison$oncoimplexus$total_chains))
    cat(sprintf("Likely chromoplexy:       %11d    %12d\n",
               comparison$chainfinder$likely,
               comparison$oncoimplexus$likely))
    cat(sprintf("Possible chromoplexy:     %11d    %12d\n",
               comparison$chainfinder$possible,
               comparison$oncoimplexus$possible))

    if (!is.null(comparison$chainfinder$n_edges)) {
        cat("\nGraph Structure:\n")
        cat("------------------------------------------------------------------------\n")
        cat(sprintf("Nodes:                    %11d    %12d\n",
                   comparison$chainfinder$n_nodes,
                   comparison$oncoimplexus$n_nodes))
        cat(sprintf("Total edges:              %11d    %12d\n",
                   comparison$chainfinder$n_edges,
                   comparison$oncoimplexus$n_edges))
        cat(sprintf("Translocation edges:      %11d    %12d\n",
                   comparison$chainfinder$n_translocation_edges,
                   comparison$oncoimplexus$n_translocation_edges))
        cat(sprintf("Deletion bridge edges:    %11d    %12d\n",
                   comparison$chainfinder$n_deletion_bridge_edges,
                   comparison$oncoimplexus$n_deletion_bridge_edges))
        cat(sprintf("Adjacency edges:          %11d    %12d\n",
                   comparison$chainfinder$n_adjacency_edges,
                   comparison$oncoimplexus$n_adjacency_edges))
    }

    cat("\nChain Overlap:\n")
    cat("------------------------------------------------------------------------\n")
    cat(sprintf("Shared chains (>50%% SV overlap): %d\n", comparison$overlap$n_shared))
    cat(sprintf("ChainFinder only:                 %d\n", comparison$overlap$n_chainfinder_only))
    cat(sprintf("OncoImplexus only:                %d\n", comparison$overlap$n_oncoimplexus_only))
    cat(sprintf("Jaccard index:                    %.3f\n", comparison$overlap$jaccard_index))

    if (!is.null(comparison$timing)) {
        cat("\nExecution Time:\n")
        cat("------------------------------------------------------------------------\n")
        cat(sprintf("ChainFinder:    %.2f seconds\n", comparison$timing$chainfinder_seconds))
        cat(sprintf("OncoImplexus:   %.2f seconds\n", comparison$timing$oncoimplexus_seconds))
    }

    cat("\n")
}


#' Create comparison summary data frame for multiple samples
#'
#' @param comparison_results List of comparison results from multiple samples
#' @param sample_ids Vector of sample IDs
#' @return Data frame with comparison summary
#' @export
create_comparison_summary_df <- function(comparison_results, sample_ids) {

    summary_list <- list()

    for (i in seq_along(comparison_results)) {
        comp <- comparison_results[[i]]

        summary_list[[i]] <- data.frame(
            sample_id = sample_ids[i],
            # ChainFinder results
            cf_total_chains = comp$chainfinder$total_chains,
            cf_likely = comp$chainfinder$likely,
            cf_possible = comp$chainfinder$possible,
            cf_adjacency_edges = if (!is.null(comp$chainfinder$n_adjacency_edges))
                                    comp$chainfinder$n_adjacency_edges else NA,
            # OncoImplexus results
            oi_total_chains = comp$oncoimplexus$total_chains,
            oi_likely = comp$oncoimplexus$likely,
            oi_possible = comp$oncoimplexus$possible,
            oi_adjacency_edges = if (!is.null(comp$oncoimplexus$n_adjacency_edges))
                                    comp$oncoimplexus$n_adjacency_edges else NA,
            # Overlap
            n_shared = comp$overlap$n_shared,
            n_cf_only = comp$overlap$n_chainfinder_only,
            n_oi_only = comp$overlap$n_oncoimplexus_only,
            jaccard = comp$overlap$jaccard_index,
            # Timing
            cf_time = if (!is.null(comp$timing)) comp$timing$chainfinder_seconds else NA,
            oi_time = if (!is.null(comp$timing)) comp$timing$oncoimplexus_seconds else NA,
            stringsAsFactors = FALSE
        )
    }

    summary_df <- do.call(rbind, summary_list)
    return(summary_df)
}


#' Plot comparison of chromoplexy methods
#'
#' @param comparison_df Summary data frame from create_comparison_summary_df
#' @return ggplot object
#' @export
plot_mechanism_comparison <- function(comparison_df) {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' needed for this function.")
    }

    # Prepare data for plotting
    plot_data <- data.frame(
        sample = rep(comparison_df$sample_id, 2),
        method = rep(c("ChainFinder", "OncoImplexus"), each = nrow(comparison_df)),
        likely = c(comparison_df$cf_likely, comparison_df$oi_likely),
        possible = c(comparison_df$cf_possible, comparison_df$oi_possible)
    )

    plot_data$total <- plot_data$likely + plot_data$possible

    # Reshape for stacked bar
    plot_data_long <- rbind(
        data.frame(
            sample = plot_data$sample,
            method = plot_data$method,
            classification = "Likely",
            count = plot_data$likely
        ),
        data.frame(
            sample = plot_data$sample,
            method = plot_data$method,
            classification = "Possible",
            count = plot_data$possible
        )
    )

    p <- ggplot2::ggplot(plot_data_long,
                         ggplot2::aes(x = sample, y = count, fill = classification)) +
        ggplot2::geom_bar(stat = "identity", position = "stack") +
        ggplot2::facet_wrap(~method, ncol = 1) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::labs(
            title = "Chromoplexy Detection: ChainFinder vs OncoImplexus",
            x = "Sample",
            y = "Number of Events",
            fill = "Classification"
        ) +
        ggplot2::scale_fill_manual(values = c("Likely" = "#E41A1C", "Possible" = "#377EB8"))

    return(p)
}
