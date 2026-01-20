#' Diagnostic tool for chromothripsis classification failures
#'
#' Provides detailed breakdown of why chromosomes were not classified as chromothripsis,
#' showing which specific criteria were not met.
#'
#' @param chromoth_output Output from detect_chromothripsis function
#' @param show_all_chromosomes Show all chromosomes, not just those with clusters (default: FALSE)
#' @param pval_fragment_joins_threshold Significance threshold for fragment joins test (default: 0.05)
#' @param pval_clustering_threshold Significance threshold for clustering tests (default: 0.05)
#' @return A detailed diagnostic data frame
#' @export
diagnose_chromothripsis_classification <- function(chromoth_output,
                                                   show_all_chromosomes = FALSE,
                                                   pval_fragment_joins_threshold = 0.05,
                                                   pval_clustering_threshold = 0.05) {

    chromSummary <- chromoth_output@chromSummary

    # Filter chromosomes with clusters if requested
    if (!show_all_chromosomes) {
        chromSummary <- chromSummary[chromSummary$clusterSize > 0, ]
    }

    if (nrow(chromSummary) == 0) {
        message("No chromosomes with SV clusters found.")
        return(NULL)
    }

    # Initialize diagnostic data frame
    diagnostic <- data.frame(
        chrom = chromSummary$chrom,
        cluster_size_intra = chromSummary$clusterSize,
        cluster_size_total = chromSummary$clusterSize_including_TRA,
        num_TRA = chromSummary$number_TRA,
        cn_osc_2state = chromSummary$max_number_oscillating_CN_segments_2_states,
        cn_osc_3state = chromSummary$max_number_oscillating_CN_segments_3_states,
        pval_fragment_joins = chromSummary$pval_fragment_joins,
        pval_fragment_joins_inter = chromSummary$inter_pval_fragment_joins,
        pval_exp_cluster = chromSummary$pval_exp_cluster,
        pval_chr_enrich = chromSummary$chr_breakpoint_enrichment,
        stringsAsFactors = FALSE
    )

    # Check each criterion
    diagnostic$meets_intra_svs_high <- diagnostic$cluster_size_intra >= 6
    diagnostic$meets_intra_svs_low <- diagnostic$cluster_size_intra >= 6
    diagnostic$meets_intra_svs_multichrom <- diagnostic$cluster_size_intra >= 3
    diagnostic$meets_TRA_multichrom <- !is.na(diagnostic$num_TRA) & diagnostic$num_TRA >= 4

    diagnostic$meets_cn_osc_high <- !is.na(diagnostic$cn_osc_2state) & diagnostic$cn_osc_2state >= 7 |
                                     !is.na(diagnostic$cn_osc_3state) & diagnostic$cn_osc_3state >= 7

    diagnostic$meets_cn_osc_low <- (!is.na(diagnostic$cn_osc_2state) &
                                     diagnostic$cn_osc_2state >= 4 &
                                     diagnostic$cn_osc_2state <= 6) |
                                    (!is.na(diagnostic$cn_osc_3state) &
                                     diagnostic$cn_osc_3state >= 4 &
                                     diagnostic$cn_osc_3state <= 6)

    diagnostic$meets_fragment_joins <- !is.na(diagnostic$pval_fragment_joins) &
                                        diagnostic$pval_fragment_joins > pval_fragment_joins_threshold

    # For Type 2 (multi-chromosomal), use inter_pval if available
    diagnostic$meets_fragment_joins_inter <- (!is.na(diagnostic$pval_fragment_joins_inter) &
                                              diagnostic$pval_fragment_joins_inter > pval_fragment_joins_threshold) |
                                             (is.na(diagnostic$pval_fragment_joins_inter) &
                                              diagnostic$meets_fragment_joins)

    diagnostic$meets_clustering <- (!is.na(diagnostic$pval_exp_cluster) &
                                     diagnostic$pval_exp_cluster < pval_clustering_threshold) |
                                    (!is.na(diagnostic$pval_chr_enrich) &
                                     diagnostic$pval_chr_enrich < pval_clustering_threshold)

    # Determine potential classification
    diagnostic$could_be_high_type1 <- diagnostic$meets_intra_svs_high &
                                       diagnostic$meets_cn_osc_high &
                                       diagnostic$meets_fragment_joins &
                                       diagnostic$meets_clustering

    diagnostic$could_be_high_type2 <- diagnostic$meets_intra_svs_multichrom &
                                       diagnostic$meets_TRA_multichrom &
                                       diagnostic$meets_cn_osc_high &
                                       diagnostic$meets_fragment_joins_inter

    diagnostic$could_be_low <- diagnostic$meets_intra_svs_low &
                                diagnostic$meets_cn_osc_low &
                                diagnostic$meets_fragment_joins &
                                diagnostic$meets_clustering

    # Add blocking reasons
    diagnostic$blocking_reasons <- ""

    for (i in 1:nrow(diagnostic)) {
        reasons <- c()

        # Check fragment joins - use inter_pval for Type 2 if available
        is_type2_candidate <- !is.na(diagnostic$num_TRA[i]) && diagnostic$num_TRA[i] >= 4
        if (is_type2_candidate && !is.na(diagnostic$pval_fragment_joins_inter[i])) {
            if (!diagnostic$meets_fragment_joins_inter[i]) {
                reasons <- c(reasons, sprintf("Fragment joins FAILED (intra+TRA p=%.3f, need >%.2f)",
                                             diagnostic$pval_fragment_joins_inter[i],
                                             pval_fragment_joins_threshold))
            }
        } else if (!diagnostic$meets_fragment_joins[i]) {
            reasons <- c(reasons, sprintf("Fragment joins FAILED (p=%.3f, need >%.2f)",
                                         diagnostic$pval_fragment_joins[i],
                                         pval_fragment_joins_threshold))
        }

        if (!diagnostic$meets_cn_osc_high[i] && !diagnostic$meets_cn_osc_low[i]) {
            cn_val <- ifelse(!is.na(diagnostic$cn_osc_2state[i]),
                           diagnostic$cn_osc_2state[i],
                           diagnostic$cn_osc_3state[i])
            reasons <- c(reasons, sprintf("CN oscillations insufficient (%s, need >=4)",
                                         ifelse(is.na(cn_val), "NA", as.character(cn_val))))
        } else if (!diagnostic$meets_cn_osc_high[i]) {
            cn_val <- ifelse(!is.na(diagnostic$cn_osc_2state[i]),
                           diagnostic$cn_osc_2state[i],
                           diagnostic$cn_osc_3state[i])
            reasons <- c(reasons, sprintf("CN oscillations=%d (in low conf range 4-6, not high >=7)",
                                         cn_val))
        }

        if (!diagnostic$meets_intra_svs_high[i]) {
            reasons <- c(reasons, sprintf("Too few intrachromosomal SVs (%d, need >=6)",
                                         diagnostic$cluster_size_intra[i]))
        }

        if (!diagnostic$meets_clustering[i]) {
            reasons <- c(reasons, "Clustering tests FAILED (both p>=0.05)")
        }

        diagnostic$blocking_reasons[i] <- paste(reasons, collapse = "; ")
    }

    # Final classification
    diagnostic$final_classification <- "Not chromothripsis"
    diagnostic$final_classification[diagnostic$could_be_high_type1 |
                                     diagnostic$could_be_high_type2] <- "High confidence"
    diagnostic$final_classification[diagnostic$could_be_low &
                                     diagnostic$final_classification == "Not chromothripsis"] <- "Low confidence"

    # Store thresholds used for reference
    attr(diagnostic, "pval_fragment_joins_threshold") <- pval_fragment_joins_threshold
    attr(diagnostic, "pval_clustering_threshold") <- pval_clustering_threshold

    return(diagnostic)
}


#' Print diagnostic report in human-readable format
#'
#' @param diagnostic_result Output from diagnose_chromothripsis_classification
#' @export
print_diagnostic_report <- function(diagnostic_result) {

    if (is.null(diagnostic_result) || nrow(diagnostic_result) == 0) {
        cat("No data to report.\n")
        return(invisible(NULL))
    }

    # Get thresholds used (default to 0.05 if not set)
    pval_frag_thresh <- attr(diagnostic_result, "pval_fragment_joins_threshold")
    if (is.null(pval_frag_thresh)) pval_frag_thresh <- 0.05

    pval_clust_thresh <- attr(diagnostic_result, "pval_clustering_threshold")
    if (is.null(pval_clust_thresh)) pval_clust_thresh <- 0.05

    cat("\n")
    cat(rep("=", 80), "\n", sep="")
    cat("         CHROMOTHRIPSIS CLASSIFICATION DIAGNOSTIC REPORT\n")
    cat(rep("=", 80), "\n\n", sep="")

    # Show thresholds used
    cat(sprintf("P-value thresholds used:\n"))
    cat(sprintf("  - Fragment joins: p > %.2f (for random joining)\n", pval_frag_thresh))
    cat(sprintf("  - Clustering tests: p < %.2f (for significant clustering)\n\n", pval_clust_thresh))

    # Summary
    n_high <- sum(diagnostic_result$final_classification == "High confidence")
    n_low <- sum(diagnostic_result$final_classification == "Low confidence")
    n_not <- sum(diagnostic_result$final_classification == "Not chromothripsis")

    cat(sprintf("Total chromosomes analyzed: %d\n", nrow(diagnostic_result)))
    cat(sprintf("  - High confidence: %d\n", n_high))
    cat(sprintf("  - Low confidence:  %d\n", n_low))
    cat(sprintf("  - Not chromothripsis: %d\n\n", n_not))

    # Detailed breakdown for each chromosome
    for (i in 1:nrow(diagnostic_result)) {
        row <- diagnostic_result[i, ]

        cat(rep("-", 80), "\n", sep="")
        cat(sprintf("Chromosome %s: %s\n", row$chrom, row$final_classification))
        cat(rep("-", 80), "\n", sep="")

        # SV counts
        cat(sprintf("SV Counts:\n"))
        cat(sprintf("  - Intrachromosomal SVs: %d %s\n",
                   row$cluster_size_intra,
                   ifelse(row$meets_intra_svs_high, "[V] (>=6)", "[X] (<6)")))
        cat(sprintf("  - Translocations: %d %s\n",
                   ifelse(is.na(row$num_TRA), 0, row$num_TRA),
                   ifelse(row$meets_TRA_multichrom, "[V] (>=4)", "")))

        # CN oscillations
        cn_val <- ifelse(!is.na(row$cn_osc_2state), row$cn_osc_2state, row$cn_osc_3state)
        cn_status <- ifelse(row$meets_cn_osc_high, "[V] High (>=7)",
                           ifelse(row$meets_cn_osc_low, "~ Low (4-6)", "[X] (<4)"))
        cat(sprintf("  - CN oscillations: %s %s\n",
                   ifelse(is.na(cn_val), "NA", as.character(cn_val)),
                   cn_status))

        # Statistical tests
        cat(sprintf("\nStatistical Tests:\n"))

        # For Type 2 (>=4 TRA), show both fragment joins p-values if available
        is_type2 <- !is.na(row$num_TRA) && row$num_TRA >= 4
        if (is_type2 && !is.na(row$pval_fragment_joins_inter)) {
            cat(sprintf("  - Fragment joins (intra only): p=%.3f\n",
                       ifelse(is.na(row$pval_fragment_joins), NA, row$pval_fragment_joins)))
            cat(sprintf("  - Fragment joins (intra+TRA): p=%.3f %s [USED for Type 2]\n",
                       row$pval_fragment_joins_inter,
                       ifelse(row$meets_fragment_joins_inter,
                              sprintf("[V] (>%.2f = random)", pval_frag_thresh),
                              sprintf("[X] (<=%.2f = non-random)", pval_frag_thresh))))
        } else {
            cat(sprintf("  - Fragment joins: p=%.3f %s\n",
                       ifelse(is.na(row$pval_fragment_joins), NA, row$pval_fragment_joins),
                       ifelse(row$meets_fragment_joins,
                              sprintf("[V] (>%.2f = random)", pval_frag_thresh),
                              sprintf("[X] (<=%.2f = non-random)", pval_frag_thresh))))
        }
        cat(sprintf("  - Exp. distribution: p=%.2e %s\n",
                   ifelse(is.na(row$pval_exp_cluster), NA, row$pval_exp_cluster),
                   ifelse(!is.na(row$pval_exp_cluster) && row$pval_exp_cluster < pval_clust_thresh,
                          sprintf("[V] (<%.2f = clustered)", pval_clust_thresh), "")))
        cat(sprintf("  - Chr. enrichment: p=%.2e %s\n",
                   ifelse(is.na(row$pval_chr_enrich), NA, row$pval_chr_enrich),
                   ifelse(!is.na(row$pval_chr_enrich) && row$pval_chr_enrich < pval_clust_thresh,
                          sprintf("[V] (<%.2f = enriched)", pval_clust_thresh), "")))

        # Blocking reasons
        if (row$final_classification == "Not chromothripsis" &&
            nchar(row$blocking_reasons) > 0) {
            cat(sprintf("\nWARNING: BLOCKING ISSUES:\n"))
            reasons <- strsplit(row$blocking_reasons, "; ")[[1]]
            for (reason in reasons) {
                cat(sprintf("  [X] %s\n", reason))
            }
        }

        cat("\n")
    }

    cat(rep("=", 80), "\n", sep="")
    cat("\n")

    invisible(diagnostic_result)
}
