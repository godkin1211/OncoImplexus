################################################################################
# Diagnostic Tools for Chromoplexy Detection
# Provides detailed diagnostic information for troubleshooting and validation
################################################################################

#' Diagnose chromoplexy detection results
#'
#' Provides detailed diagnostic information about why chains were or weren't
#' classified as chromoplexy, including threshold comparisons and evidence breakdown.
#'
#' @param result Chromoplexy detection result (from detect_chromoplexy or detect_chromoplexy_v2)
#' @param show_all_chains Show all chains (default: FALSE, only show likely/possible)
#' @param verbose Print detailed information (default: TRUE)
#'
#' @return Invisibly returns diagnostic data frame
#'
#' @export
diagnose_chromoplexy_detection <- function(result,
                                           show_all_chains = FALSE,
                                           verbose = TRUE) {

    if (verbose) {
        cat("\n")
        cat(rep("=", 70), "\n", sep = "")
        cat("         CHROMOPLEXY DETECTION DIAGNOSTICS\n")
        cat(rep("=", 70), "\n\n", sep = "")
    }

    # Overall statistics
    if (verbose) {
        cat("Overall Detection Summary:\n")
        cat(rep("-", 70), "\n", sep = "")
        cat(sprintf("  Total chains detected: %d\n", result$total_chains))
        cat(sprintf("  Likely chromoplexy:    %d (%.1f%%)\n",
                    result$likely_chromoplexy,
                    if (result$total_chains > 0) 100 * result$likely_chromoplexy / result$total_chains else 0))
        cat(sprintf("  Possible chromoplexy:  %d (%.1f%%)\n",
                    result$possible_chromoplexy,
                    if (result$total_chains > 0) 100 * result$possible_chromoplexy / result$total_chains else 0))

        if (!is.null(result$deletion_bridges)) {
            cat(sprintf("  Deletion bridges found: %d\n", length(result$deletion_bridges)))
        }
        cat("\n")
    }

    if (result$total_chains == 0) {
        if (verbose) {
            cat("No chains detected. Possible reasons:\n")
            cat("  * Insufficient inter-chromosomal SVs (need >=3)\n")
            cat("  * SVs don't form connected chains\n")
            cat("  * Minimum chromosome threshold not met\n")
            cat("\n")
        }
        return(invisible(NULL))
    }

    # Per-chain diagnostics
    if (verbose) {
        cat("Per-Chain Diagnostics:\n")
        cat(rep("=", 70), "\n\n", sep = "")
    }

    summary_df <- result$summary

    # Filter chains if requested
    if (!show_all_chains) {
        summary_df <- summary_df[summary_df$classification %in%
                                c("Likely chromoplexy", "Possible chromoplexy"), ]
    }

    if (nrow(summary_df) == 0 && !show_all_chains) {
        if (verbose) {
            cat("No likely or possible chromoplexy chains detected.\n")
            cat("Use show_all_chains=TRUE to see all detected chains.\n\n")
        }
        return(invisible(NULL))
    }

    for (i in 1:nrow(summary_df)) {
        chain_row <- summary_df[i, ]

        if (verbose) {
            diagnose_single_chain(chain_row, result)
        }
    }

    if (verbose) {
        cat(rep("=", 70), "\n\n", sep = "")
    }

    invisible(summary_df)
}


#' Diagnose a single chain
#'
#' @keywords internal
diagnose_single_chain <- function(chain_row, result) {

    cat(sprintf("Chain %d: %s\n", chain_row$chain_id, chain_row$classification))
    cat(rep("-", 70), "\n", sep = "")

    # Basic info
    cat("Basic Information:\n")
    cat(sprintf("  Chromosomes: %s (n=%d)\n",
                chain_row$chromosomes_involved,
                chain_row$n_chromosomes))
    cat(sprintf("  Translocations: %d\n", chain_row$n_translocations))
    cat(sprintf("  Topology: %s\n", if (chain_row$is_cycle) "Circular (cycle)" else "Linear"))
    cat("\n")

    # Classification criteria
    cat("Classification Criteria:\n")

    # Criterion 1: Chromosomes
    chr_pass <- chain_row$n_chromosomes >= 3
    cat(sprintf("  [%s] Multiple chromosomes: %d (threshold: >=3)\n",
                if (chr_pass) "[V]" else "[X]",
                chain_row$n_chromosomes))

    # Criterion 2: Translocations
    tlx_pass <- chain_row$n_translocations >= 3
    cat(sprintf("  [%s] Multiple translocations: %d (threshold: >=3)\n",
                if (tlx_pass) "[V]" else "[X]",
                chain_row$n_translocations))

    # Criterion 3: CN stability
    cn_pass <- chain_row$cn_stability_score >= 0.7
    cat(sprintf("  [%s] CN stability: %.3f (threshold: >=0.7)\n",
                if (cn_pass) "[V]" else "[X]",
                chain_row$cn_stability_score))

    if ("cn_global_stability" %in% colnames(chain_row)) {
        cat(sprintf("      - Global stability: %.3f\n", chain_row$cn_global_stability))
        cat(sprintf("      - Local changes: %.3f\n", chain_row$cn_local_changes))
    }

    # Criterion 4: Complexity
    complexity_pass <- chain_row$complexity_score >= 0.3
    cat(sprintf("  [%s] Complexity: %.3f (threshold: >=0.3)\n",
                if (complexity_pass) "[V]" else "[X]",
                chain_row$complexity_score))

    # Criterion 5: Deletion bridges (if available)
    if ("deletion_bridge_score" %in% colnames(chain_row)) {
        del_pass <- chain_row$deletion_bridge_score >= 0.5
        cat(sprintf("  [%s] Deletion bridges: %.3f (threshold: >=0.5, n=%d)\n",
                    if (del_pass) "[V]" else "[X]",
                    chain_row$deletion_bridge_score,
                    chain_row$n_deletion_bridges))
    }

    # Criterion 6: Statistical significance (if available)
    if ("pvalue" %in% colnames(chain_row) && !is.na(chain_row$pvalue)) {
        stat_pass <- chain_row$is_statistically_significant
        cat(sprintf("  [%s] Statistical significance: p=%.4e, FDR=%.4e\n",
                    if (stat_pass) "[V]" else "[X]",
                    chain_row$pvalue,
                    chain_row$fdr))
    }

    cat("\n")

    # Overall score
    if ("combined_score" %in% colnames(chain_row)) {
        cat(sprintf("Combined Evidence Score: %.3f\n", chain_row$combined_score))
        cat("\n")
    }

    # Interpretation
    cat("Interpretation:\n")
    if (chain_row$classification == "Likely chromoplexy") {
        cat("  This chain shows strong evidence for chromoplexy:\n")
        cat("  - Meets most or all classification criteria\n")
        if ("is_statistically_significant" %in% colnames(chain_row) &&
            !is.na(chain_row$is_statistically_significant) &&
            chain_row$is_statistically_significant) {
            cat("  - Statistically significant (not random)\n")
        }
        if ("deletion_bridge_score" %in% colnames(chain_row) &&
            chain_row$deletion_bridge_score > 0.7) {
            cat("  - Strong deletion bridge signature\n")
        }
    } else if (chain_row$classification == "Possible chromoplexy") {
        cat("  This chain shows moderate evidence for chromoplexy:\n")
        cat("  - Meets some but not all criteria\n")
        cat("  - Consider:\n")
        if (!cn_pass) {
            cat("    -> CN instability may indicate chromothripsis instead\n")
        }
        if (!complexity_pass) {
            cat("    -> Low complexity may indicate simple rearrangement\n")
        }
        if ("is_statistically_significant" %in% colnames(chain_row) &&
            !is.na(chain_row$is_statistically_significant) &&
            !chain_row$is_statistically_significant) {
            cat("    -> Not statistically significant (may be random)\n")
        }
    } else {
        cat("  This chain shows weak evidence for chromoplexy:\n")
        cat("  - Fails multiple classification criteria\n")
        cat("  - May represent:\n")
        cat("    -> Random co-occurrence of translocations\n")
        cat("    -> Different type of rearrangement\n")
        cat("    -> Incomplete detection\n")
    }

    cat("\n\n")
}


#' Compare chromoplexy detection methods
#'
#' Compares results from original and enhanced detection methods.
#'
#' @param result_v1 Result from detect_chromoplexy()
#' @param result_v2 Result from detect_chromoplexy_v2()
#' @param verbose Print comparison (default: TRUE)
#'
#' @return Invisibly returns comparison data frame
#'
#' @export
compare_chromoplexy_methods <- function(result_v1, result_v2, verbose = TRUE) {

    if (verbose) {
        cat("\n")
        cat(rep("=", 70), "\n", sep = "")
        cat("         CHROMOPLEXY METHOD COMPARISON\n")
        cat(rep("=", 70), "\n\n", sep = "")

        cat("Detection Summary:\n")
        cat(rep("-", 70), "\n", sep = "")
        cat(sprintf("                      Original   Enhanced   Difference\n"))
        cat(sprintf("Total chains:         %8d   %8d   %+9d\n",
                    result_v1$total_chains,
                    result_v2$total_chains,
                    result_v2$total_chains - result_v1$total_chains))
        cat(sprintf("Likely chromoplexy:   %8d   %8d   %+9d\n",
                    result_v1$likely_chromoplexy,
                    result_v2$likely_chromoplexy,
                    result_v2$likely_chromoplexy - result_v1$likely_chromoplexy))
        cat(sprintf("Possible chromoplexy: %8d   %8d   %+9d\n",
                    result_v1$possible_chromoplexy,
                    result_v2$possible_chromoplexy,
                    result_v2$possible_chromoplexy - result_v1$possible_chromoplexy))
        cat("\n")

        # Enhanced features
        if (!is.null(result_v2$deletion_bridges)) {
            cat("Enhanced Features (v2 only):\n")
            cat(rep("-", 70), "\n", sep = "")
            cat(sprintf("  Deletion bridges identified: %d\n", length(result_v2$deletion_bridges)))

            if (result_v2$total_chains > 0) {
                n_sig <- sum(result_v2$summary$is_statistically_significant, na.rm = TRUE)
                cat(sprintf("  Statistically significant chains: %d / %d\n",
                            n_sig, result_v2$total_chains))

                if (n_sig > 0) {
                    sig_chains <- result_v2$summary[result_v2$summary$is_statistically_significant == TRUE, ]
                    cat(sprintf("  Mean p-value (significant): %.4e\n", mean(sig_chains$pvalue, na.rm = TRUE)))
                }
            }
            cat("\n")
        }

        # Algorithm improvements
        cat("Key Improvements:\n")
        cat(rep("-", 70), "\n", sep = "")
        if (result_v2$total_chains > result_v1$total_chains) {
            cat(sprintf("  [+] Found %d additional chains (comprehensive path search)\n",
                        result_v2$total_chains - result_v1$total_chains))
        }
        if (!is.null(result_v2$deletion_bridges) && length(result_v2$deletion_bridges) > 0) {
            cat(sprintf("  [+] Identified %d deletion bridges\n", length(result_v2$deletion_bridges)))
        }
        if (result_v2$total_chains > 0 && !is.na(result_v2$summary$pvalue[1])) {
            cat("  [+] Statistical significance testing implemented\n")
        }
        if ("cn_global_stability" %in% colnames(result_v2$summary)) {
            cat("  [+] Enhanced CN stability evaluation (multi-component)\n")
        }
        if ("combined_score" %in% colnames(result_v2$summary)) {
            cat("  [+] Combined evidence scoring for ranking\n")
        }
        cat("\n")

        cat(rep("=", 70), "\n\n", sep = "")
    }

    # Create comparison data frame
    comparison <- data.frame(
        metric = c("Total chains", "Likely", "Possible",
                   "Deletion bridges", "Statistical testing"),
        v1 = c(result_v1$total_chains,
               result_v1$likely_chromoplexy,
               result_v1$possible_chromoplexy,
               NA, NA),
        v2 = c(result_v2$total_chains,
               result_v2$likely_chromoplexy,
               result_v2$possible_chromoplexy,
               if (!is.null(result_v2$deletion_bridges)) length(result_v2$deletion_bridges) else 0,
               sum(result_v2$summary$is_statistically_significant, na.rm = TRUE)),
        stringsAsFactors = FALSE
    )
    comparison$difference <- comparison$v2 - comparison$v1

    invisible(comparison)
}


#' Print classification thresholds
#'
#' Shows the thresholds used for chromoplexy classification.
#'
#' @param use_statistical_testing Whether statistical testing is enabled
#' @param fdr_threshold FDR threshold if using statistical testing
#'
#' @export
print_chromoplexy_thresholds <- function(use_statistical_testing = TRUE,
                                         fdr_threshold = 0.01) {

    cat("\n")
    cat(rep("=", 70), "\n", sep = "")
    cat("         CHROMOPLEXY CLASSIFICATION THRESHOLDS\n")
    cat(rep("=", 70), "\n\n", sep = "")

    cat("Core Criteria:\n")
    cat("  1. Chromosomes:    >=3\n")
    cat("  2. Translocations: >=3\n")
    cat("  3. CN stability:   >=0.7\n")
    cat("  4. Complexity:     >=0.3\n")
    cat("\n")

    if (use_statistical_testing) {
        cat("Additional Criteria (Enhanced):\n")
        cat(sprintf("  5. Deletion bridges:  >=0.5 score\n"))
        cat(sprintf("  6. Statistical sig:   FDR < %.3f\n", fdr_threshold))
        cat("\n")

        cat("Classification Rules (Enhanced):\n")
        cat("  Likely chromoplexy:   >=5 criteria + statistically significant\n")
        cat("  Possible chromoplexy: >=4 criteria\n")
        cat("  Unlikely chromoplexy: >=3 criteria\n")
        cat("  Not chromoplexy:      <3 criteria\n")
    } else {
        cat("Classification Rules (Original):\n")
        cat("  Likely chromoplexy:   >=4 criteria\n")
        cat("  Possible chromoplexy: >=3 criteria\n")
        cat("  Unlikely chromoplexy: >=2 criteria\n")
        cat("  Not chromoplexy:      <2 criteria\n")
    }

    cat("\n")
    cat(rep("=", 70), "\n\n", sep = "")
}


#' Generate diagnostic report for chromoplexy detection
#'
#' Creates a comprehensive diagnostic report combining all diagnostic tools.
#'
#' @param result Chromoplexy detection result
#' @param output_file Optional file to write report (default: NULL, prints to console)
#'
#' @export
generate_chromoplexy_report <- function(result, output_file = NULL) {

    if (!is.null(output_file)) {
        sink(output_file)
        on.exit(sink())
    }

    cat("\n")
    cat(rep("=", 70), "\n", sep = "")
    cat("         CHROMOPLEXY DETECTION REPORT\n")
    cat(rep("=", 70), "\n\n", sep = "")

    cat(sprintf("Generated: %s\n", Sys.time()))
    cat("\n")

    # Detection parameters
    if (!is.null(result$parameters)) {
        cat("Detection Parameters:\n")
        cat(rep("-", 70), "\n", sep = "")
        params <- result$parameters
        cat(sprintf("  Minimum chromosomes:         %d\n", params$min_chromosomes))
        cat(sprintf("  Minimum translocations:      %d\n", params$min_translocations))
        cat(sprintf("  Statistical testing:         %s\n", params$use_statistical_testing))
        if (params$use_statistical_testing) {
            cat(sprintf("  FDR threshold:               %.3f\n", params$fdr_threshold))
        }
        cat(sprintf("  Deletion bridge detection:   %s\n", params$identify_deletion_bridges))
        cat("\n")
    }

    # Main diagnostics
    diagnose_chromoplexy_detection(result, show_all_chains = TRUE, verbose = TRUE)

    # Summary statistics
    if (result$total_chains > 0) {
        cat("Summary Statistics:\n")
        cat(rep("=", 70), "\n\n", sep = "")

        cat("Score Distributions:\n")
        cat(sprintf("  CN stability:     mean=%.3f, sd=%.3f, range=[%.3f, %.3f]\n",
                    mean(result$summary$cn_stability_score),
                    sd(result$summary$cn_stability_score),
                    min(result$summary$cn_stability_score),
                    max(result$summary$cn_stability_score)))

        cat(sprintf("  Complexity:       mean=%.3f, sd=%.3f, range=[%.3f, %.3f]\n",
                    mean(result$summary$complexity_score),
                    sd(result$summary$complexity_score),
                    min(result$summary$complexity_score),
                    max(result$summary$complexity_score)))

        if ("combined_score" %in% colnames(result$summary)) {
            cat(sprintf("  Combined score:   mean=%.3f, sd=%.3f, range=[%.3f, %.3f]\n",
                        mean(result$summary$combined_score),
                        sd(result$summary$combined_score),
                        min(result$summary$combined_score),
                        max(result$summary$combined_score)))
        }

        cat("\n")
    }

    cat(rep("=", 70), "\n", sep = "")
    cat("End of Report\n")
    cat(rep("=", 70), "\n\n", sep = "")

    if (!is.null(output_file)) {
        message(sprintf("Report written to: %s", output_file))
    }
}
