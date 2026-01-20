#' Calculate integrated chromothripsis confidence score
#'
#' This function integrates multiple statistical criteria to provide a unified
#' confidence score for chromothripsis detection. The score ranges from 0 to 1,
#' where higher values indicate stronger evidence for chromothripsis.
#'
#' @param chromSummary A data frame from ShatterSeek output containing statistical criteria
#' @param chr Optional: specific chromosome to score (default: all chromosomes with clusters)
#' @return A data frame with chromosome names and confidence scores
#' @details
#' The scoring system integrates five key criteria:
#' \itemize{
#'   \item CN oscillations (weight: 0.25)
#'   \item Fragment joins randomness (weight: 0.20)
#'   \item Breakpoint clustering (weight: 0.20)
#'   \item Chromosome enrichment (weight: 0.15)
#'   \item Cluster size (weight: 0.20)
#' }
#'
#' Based on Cortes-Ciriano et al. Nature Genetics 2020, a high-confidence
#' chromothripsis event typically has:
#' - Cluster size >= 10 SVs
#' - CN oscillations >= 4 segments
#' - p-value for fragment joins > 0.05
#' - p-value for exponential distribution < 0.05
#'
#' @export
calculate_confidence_score <- function(chromSummary, chr = NULL) {

    if (!is.null(chr)) {
        chromSummary <- chromSummary[chromSummary$chrom == chr, ]
    }

    # Filter chromosomes with clusters
    chromSummary <- chromSummary[chromSummary$clusterSize > 0, ]

    if (nrow(chromSummary) == 0) {
        warning("No chromosomes with SV clusters found")
        return(data.frame(chrom = character(0),
                         confidence_score = numeric(0),
                         confidence_category = character(0)))
    }

    scores <- data.frame(
        chrom = chromSummary$chrom,
        cluster_score = 0,
        cn_oscillation_score = 0,
        fragment_joins_score = 0,
        clustering_score = 0,
        enrichment_score = 0,
        total_score = 0,
        confidence_category = "",
        stringsAsFactors = FALSE
    )

    for (i in 1:nrow(chromSummary)) {
        row <- chromSummary[i, ]

        # 1. Cluster size score (weight: 0.20)
        # Based on literature: >=10 SVs is typical for chromothripsis
        cluster_size <- row$clusterSize_including_TRA
        if (is.na(cluster_size)) cluster_size <- row$clusterSize

        if (cluster_size >= 10) {
            scores$cluster_score[i] <- 1.0
        } else if (cluster_size >= 6) {
            scores$cluster_score[i] <- 0.6
        } else if (cluster_size >= 3) {
            scores$cluster_score[i] <- 0.3
        } else {
            scores$cluster_score[i] <- 0.1
        }

        # 2. CN oscillation score (weight: 0.25)
        # Chromothripsis shows oscillating CN patterns
        cn_osc_2state <- row$max_number_oscillating_CN_segments_2_states
        cn_osc_3state <- row$max_number_oscillating_CN_segments_3_states

        if (!is.na(cn_osc_2state) && cn_osc_2state >= 4) {
            scores$cn_oscillation_score[i] <- 1.0
        } else if (!is.na(cn_osc_3state) && cn_osc_3state >= 4) {
            scores$cn_oscillation_score[i] <- 0.8
        } else if (!is.na(cn_osc_2state) && cn_osc_2state >= 2) {
            scores$cn_oscillation_score[i] <- 0.4
        } else {
            scores$cn_oscillation_score[i] <- 0.0
        }

        # 3. Fragment joins randomness score (weight: 0.20)
        # Random rejoining: p-value should be > 0.05 (not significantly different from random)
        pval_joins <- row$pval_fragment_joins

        if (!is.na(pval_joins)) {
            if (pval_joins > 0.05) {
                scores$fragment_joins_score[i] <- 1.0
            } else if (pval_joins > 0.01) {
                scores$fragment_joins_score[i] <- 0.5
            } else {
                scores$fragment_joins_score[i] <- 0.0
            }
        }

        # 4. Breakpoint clustering score (weight: 0.20)
        # Clustered breakpoints: p-value should be < 0.05 (reject exponential distribution)
        pval_exp_cluster <- row$pval_exp_cluster
        pval_exp_chr <- row$pval_exp_chr

        if (!is.na(pval_exp_cluster) && pval_exp_cluster < 0.05) {
            scores$clustering_score[i] <- 1.0
        } else if (!is.na(pval_exp_chr) && pval_exp_chr < 0.05) {
            scores$clustering_score[i] <- 0.7
        } else if (!is.na(pval_exp_cluster) && pval_exp_cluster < 0.1) {
            scores$clustering_score[i] <- 0.5
        } else {
            scores$clustering_score[i] <- 0.0
        }

        # 5. Chromosome breakpoint enrichment score (weight: 0.15)
        # Enrichment: p-value should be < 0.05
        pval_enrich <- row$chr_breakpoint_enrichment

        if (!is.na(pval_enrich)) {
            if (pval_enrich < 0.001) {
                scores$enrichment_score[i] <- 1.0
            } else if (pval_enrich < 0.01) {
                scores$enrichment_score[i] <- 0.7
            } else if (pval_enrich < 0.05) {
                scores$enrichment_score[i] <- 0.4
            } else {
                scores$enrichment_score[i] <- 0.0
            }
        }

        # Calculate weighted total score
        scores$total_score[i] <- (
            scores$cluster_score[i] * 0.20 +
            scores$cn_oscillation_score[i] * 0.25 +
            scores$fragment_joins_score[i] * 0.20 +
            scores$clustering_score[i] * 0.20 +
            scores$enrichment_score[i] * 0.15
        )

        # Assign confidence category
        if (scores$total_score[i] >= 0.7) {
            scores$confidence_category[i] <- "High"
        } else if (scores$total_score[i] >= 0.4) {
            scores$confidence_category[i] <- "Moderate"
        } else {
            scores$confidence_category[i] <- "Low"
        }
    }

    # Rename for clarity
    names(scores)[names(scores) == "total_score"] <- "confidence_score"

    return(scores)
}


#' Classify chromothripsis events based on criteria from Cortes-Ciriano et al.
#'
#' Applies the classification criteria from the original ShatterSeek publication
#' (Cortes-Ciriano et al., Nature Genetics 2020) to classify chromosomes as
#' high confidence chromothripsis, low confidence chromothripsis, or not chromothripsis.
#'
#' @param chromoth_output Output from detect_chromothripsis function
#' @param pval_fragment_joins_threshold Significance threshold for fragment joins test (default: 0.05).
#'   Fragment joins p-value must be GREATER than this threshold to pass (indicating random joining).
#'   Use 0.01 for more stringent filtering (higher specificity, lower sensitivity).
#' @param pval_clustering_threshold Significance threshold for clustering tests (default: 0.05).
#'   Clustering p-values must be LESS than this threshold to pass (indicating significant clustering).
#'   Use 0.01 for more stringent filtering (higher specificity, lower sensitivity).
#' @return A data frame with classification results
#' @details
#' Classification criteria based on Korbel and Campbell (Cell 2013) and
#' Cortes-Ciriano et al. (Nature Genetics 2020):
#'
#' **High confidence chromothripsis (Type 1 - single chromosome)**:
#' - At least 6 interleaved intrachromosomal SVs
#' - 7 or more contiguous segments oscillating between 2 CN states
#' - Fragment joins test passed (p > pval_fragment_joins_threshold, default 0.05)
#' - Either chromosomal breakpoint enrichment OR exponential distribution test
#'   (p < pval_clustering_threshold, default 0.05)
#'
#' **High confidence chromothripsis (Type 2 - multi-chromosomal)**:
#' - At least 3 interleaved intrachromosomal SVs
#' - 4 or more interchromosomal SVs (translocations)
#' - 7 or more contiguous segments oscillating between 2 CN states
#' - Fragment joins test passed (p > pval_fragment_joins_threshold)
#'
#' **Low confidence chromothripsis**:
#' - At least 6 interleaved intrachromosomal SVs
#' - 4, 5, or 6 contiguous segments oscillating between 2 CN states
#' - Fragment joins test passed (p > pval_fragment_joins_threshold)
#' - Either chromosomal breakpoint enrichment OR exponential distribution test
#'   (p < pval_clustering_threshold)
#'
#' **Note on threshold selection**:
#' - Default (0.05): Standard significance level, balances sensitivity and specificity
#' - Stringent (0.01): Reduces false positives but may miss borderline cases
#' - The choice depends on your research goals (discovery vs. validation)
#'
#' @references
#' Korbel, J. O. & Campbell, P. J. Criteria for inference of chromothripsis in cancer genomes.
#' Cell 152, 1226-1236 (2013).
#'
#' Cortes-Ciriano, I. et al. Comprehensive analysis of chromothripsis in 2,658 human cancers
#' using whole-genome sequencing. Nat. Genet. 52, 331-341 (2020).
#'
#' @export
classify_chromothripsis <- function(chromoth_output,
                                    pval_fragment_joins_threshold = 0.05,
                                    pval_clustering_threshold = 0.05) {

    chromSummary <- chromoth_output@chromSummary

    # Filter chromosomes with clusters
    chromSummary <- chromSummary[chromSummary$clusterSize > 0, ]

    if (nrow(chromSummary) == 0) {
        return(data.frame(
            chrom = character(0),
            classification = character(0),
            cluster_size = numeric(0),
            cluster_size_with_TRA = numeric(0),
            cn_oscillations = numeric(0),
            pval_fragment_joins = numeric(0),
            pval_exp_cluster = numeric(0),
            chr_breakpoint_enrichment = numeric(0),
            confidence_score = numeric(0)
        ))
    }

    # Calculate confidence scores (for backward compatibility and ranking)
    scores <- calculate_confidence_score(chromSummary)

    # Merge scores with original summary
    result <- merge(chromSummary, scores[, c("chrom", "confidence_score", "confidence_category")],
                   by = "chrom", all.x = FALSE)

    # Apply classification criteria based on tutorial
    result$classification <- "Not chromothripsis"

    for (i in 1:nrow(result)) {
        # Extract key metrics
        cluster_size_intra <- result$clusterSize[i]
        cluster_size_total <- result$clusterSize_including_TRA[i]
        if (is.na(cluster_size_total)) cluster_size_total <- cluster_size_intra

        num_TRA <- result$number_TRA[i]
        if (is.na(num_TRA)) num_TRA <- 0

        cn_osc_2state <- result$max_number_oscillating_CN_segments_2_states[i]
        cn_osc_3state <- result$max_number_oscillating_CN_segments_3_states[i]

        # Use 2-state as primary, fallback to 3-state
        cn_osc <- cn_osc_2state
        if (is.na(cn_osc)) cn_osc <- cn_osc_3state

        pval_joins <- result$pval_fragment_joins[i]
        pval_joins_inter <- result$inter_pval_fragment_joins[i]
        pval_exp_cluster <- result$pval_exp_cluster[i]
        pval_chr_enrich <- result$chr_breakpoint_enrichment[i]

        # Check fragment joins test (CRITICAL: p > threshold means random joining)
        # For Type 1: use pval_fragment_joins (intrachromosomal only)
        # For Type 2: use inter_pval_fragment_joins (inter + intra) if available
        fragment_joins_passed <- !is.na(pval_joins) && pval_joins > pval_fragment_joins_threshold
        fragment_joins_passed_inter <- !is.na(pval_joins_inter) && pval_joins_inter > pval_fragment_joins_threshold

        # Check clustering tests (at least one should pass: p < threshold means clustering)
        clustering_passed <- (!is.na(pval_exp_cluster) && pval_exp_cluster < pval_clustering_threshold) ||
                            (!is.na(pval_chr_enrich) && pval_chr_enrich < pval_clustering_threshold)

        # HIGH CONFIDENCE TYPE 1: Single/few chromosome, many intrachromosomal SVs
        # - At least 6 intrachromosomal SVs
        # - 7+ oscillating CN segments
        # - Fragment joins test passed
        # - Clustering test passed
        if (!is.na(cluster_size_intra) && cluster_size_intra >= 6 &&
            !is.na(cn_osc) && cn_osc >= 7 &&
            fragment_joins_passed &&
            clustering_passed) {
            result$classification[i] <- "High confidence"
        }
        # HIGH CONFIDENCE TYPE 2: Multi-chromosomal
        # - At least 3 intrachromosomal SVs
        # - 4+ interchromosomal SVs
        # - 7+ oscillating CN segments
        # - Fragment joins test passed (prefer inter_pval_fragment_joins if available)
        else if (!is.na(cluster_size_intra) && cluster_size_intra >= 3 &&
                 num_TRA >= 4 &&
                 !is.na(cn_osc) && cn_osc >= 7 &&
                 (fragment_joins_passed_inter || (is.na(pval_joins_inter) && fragment_joins_passed))) {
            result$classification[i] <- "High confidence"
        }
        # LOW CONFIDENCE
        # - At least 6 intrachromosomal SVs
        # - 4-6 oscillating CN segments (note the upper bound!)
        # - Fragment joins test passed
        # - Clustering test passed
        else if (!is.na(cluster_size_intra) && cluster_size_intra >= 6 &&
                 !is.na(cn_osc) && cn_osc >= 4 && cn_osc <= 6 &&
                 fragment_joins_passed &&
                 clustering_passed) {
            result$classification[i] <- "Low confidence"
        }
    }

    # Select key columns for output
    output_cols <- c("chrom", "classification", "confidence_score", "confidence_category",
                    "clusterSize", "clusterSize_including_TRA", "number_TRA",
                    "max_number_oscillating_CN_segments_2_states",
                    "max_number_oscillating_CN_segments_3_states",
                    "pval_fragment_joins", "inter_pval_fragment_joins",
                    "pval_exp_cluster", "chr_breakpoint_enrichment")

    output_cols <- output_cols[output_cols %in% names(result)]
    result <- result[, output_cols]

    # Sort by classification priority, then confidence score
    result$class_priority <- ifelse(result$classification == "High confidence", 1,
                                   ifelse(result$classification == "Low confidence", 2, 3))
    result <- result[order(result$class_priority, -result$confidence_score), ]
    result$class_priority <- NULL

    return(result)
}


#' Generate a summary report of chromothripsis analysis
#'
#' Creates a human-readable summary of chromothripsis detection results,
#' highlighting high-confidence and low-confidence events based on the
#' criteria from Cortes-Ciriano et al. (Nature Genetics 2020).
#'
#' @param chromoth_output Output from detect_chromothripsis function
#' @param print_summary Logical, whether to print summary to console (default: TRUE)
#' @return A list containing summary statistics and classifications
#' @export
summarize_chromothripsis <- function(chromoth_output, print_summary = TRUE) {

    # Get classifications
    classifications <- classify_chromothripsis(chromoth_output)

    # Count classifications
    n_high <- sum(classifications$classification == "High confidence")
    n_low <- sum(classifications$classification == "Low confidence")
    n_not <- sum(classifications$classification == "Not chromothripsis")

    # Get high confidence chromosomes
    high_conf <- classifications[classifications$classification == "High confidence", ]
    low_conf <- classifications[classifications$classification == "Low confidence", ]

    # Prepare summary
    summary_list <- list(
        total_chromosomes_with_clusters = nrow(classifications),
        high_confidence_chromothripsis = n_high,
        low_confidence_chromothripsis = n_low,
        not_chromothripsis = n_not,
        high_confidence_chromosomes = if(nrow(high_conf) > 0) high_conf$chrom else character(0),
        low_confidence_chromosomes = if(nrow(low_conf) > 0) low_conf$chrom else character(0),
        detailed_results = classifications
    )

    if (print_summary) {
        cat("\n")
        cat(rep("=", 70), "\n", sep="")
        cat("         CHROMOTHRIPSIS DETECTION SUMMARY\n")
        cat(rep("=", 70), "\n\n", sep="")

        cat(sprintf("Total chromosomes with SV clusters: %d\n\n",
                   summary_list$total_chromosomes_with_clusters))

        cat("Classification Results (Cortes-Ciriano et al. 2020 criteria):\n")
        cat(sprintf("  - High confidence: %d chromosome(s)\n", n_high))
        cat(sprintf("  - Low confidence:  %d chromosome(s)\n", n_low))
        cat(sprintf("  - Not chromothripsis: %d chromosome(s)\n\n", n_not))

        if (length(summary_list$high_confidence_chromosomes) > 0) {
            cat("High-confidence chromothripsis chromosomes:\n")
            cat(sprintf("  %s\n\n", paste(summary_list$high_confidence_chromosomes, collapse=", ")))
        }

        if (length(summary_list$low_confidence_chromosomes) > 0) {
            cat("Low-confidence chromothripsis chromosomes:\n")
            cat(sprintf("  %s\n\n", paste(summary_list$low_confidence_chromosomes, collapse=", ")))
        }

        if (n_high == 0 && n_low == 0) {
            cat("No chromothripsis events detected.\n\n")
        }

        if (n_high > 0 || n_low > 0) {
            cat("Detailed Results:\n")
            cat(sprintf("%-8s %-20s %-10s %-10s %-8s %-12s\n",
                       "Chrom", "Classification", "CN Osc.", "Intra SVs", "TRA", "FragJoin p"))
            cat(rep("-", 70), "\n", sep="")

            for (i in 1:nrow(classifications)) {
                if (classifications$classification[i] != "Not chromothripsis") {
                    cn_osc <- classifications$max_number_oscillating_CN_segments_2_states[i]
                    if (is.na(cn_osc)) cn_osc <- classifications$max_number_oscillating_CN_segments_3_states[i]
                    if (is.na(cn_osc)) cn_osc <- "-"

                    cluster_size <- classifications$clusterSize[i]
                    num_tra <- classifications$number_TRA[i]
                    if (is.na(num_tra)) num_tra <- 0

                    pval_joins <- classifications$pval_fragment_joins[i]
                    pval_joins_str <- if (!is.na(pval_joins)) sprintf("%.3f", pval_joins) else "-"

                    cat(sprintf("%-8s %-20s %-10s %-10d %-8d %-12s\n",
                               classifications$chrom[i],
                               classifications$classification[i],
                               cn_osc,
                               cluster_size,
                               num_tra,
                               pval_joins_str))
                }
            }
        }

        cat("\n")
        cat("Note: Fragment joins p-value > 0.05 indicates random fragment joining\n")
        cat("      (characteristic of chromothripsis)\n")
        cat(rep("=", 70), "\n", sep="")
        cat("\n")
    }

    return(invisible(summary_list))
}


#' Check data quality and completeness
#'
#' Validates input data quality and provides warnings for potential issues
#' that might affect chromothripsis detection accuracy.
#'
#' @param SV.sample SVs object or data frame
#' @param CNV.sample CNVsegs object or data frame
#' @param verbose Logical, whether to print detailed messages (default: TRUE)
#' @return A list with data quality metrics and warnings
#' @export
check_data_quality <- function(SV.sample, CNV.sample = NULL, verbose = TRUE) {

    issues <- list()
    warnings_list <- character(0)

    # Convert to data frame if needed
    if (is(SV.sample, "SVs")) {
        SV.sample <- as(SV.sample, "data.frame")
    }

    if (!is.null(CNV.sample) && is(CNV.sample, "CNVsegs")) {
        CNV.sample <- as(CNV.sample, "data.frame")
    }

    # Check SV data
    issues$total_SVs <- nrow(SV.sample)
    issues$intrachromosomal_SVs <- sum(SV.sample$chrom1 == SV.sample$chrom2)
    issues$interchromosomal_SVs <- sum(SV.sample$chrom1 != SV.sample$chrom2)

    if (issues$intrachromosomal_SVs < 5) {
        warnings_list <- c(warnings_list,
                          "Very few intrachromosomal SVs (<5). Chromothripsis detection requires multiple clustered SVs.")
    }

    # Check SV types distribution
    if ("SVtype" %in% names(SV.sample)) {
        sv_types <- table(SV.sample$SVtype)
        issues$sv_type_distribution <- sv_types

        if (length(sv_types) < 2) {
            warnings_list <- c(warnings_list,
                             "Only one SV type detected. Chromothripsis typically involves multiple SV types.")
        }
    }

    # Check CNV data
    if (!is.null(CNV.sample)) {
        issues$total_CNV_segments <- nrow(CNV.sample)

        if (nrow(CNV.sample) == 0) {
            warnings_list <- c(warnings_list,
                             "No CNV data provided. CN oscillations cannot be evaluated.")
        } else if (nrow(CNV.sample) < 10) {
            warnings_list <- c(warnings_list,
                             "Very few CNV segments (<10). May affect CN oscillation detection.")
        }

        # Check for truly adjacent segments with same CN
        # Only check segments that are genomically adjacent (end[i] >= start[i+1])
        if (nrow(CNV.sample) > 1) {
            chrom_list <- unique(CNV.sample$chrom)
            duplicate_cn <- 0

            for (chr in chrom_list) {
                chr_cnv <- CNV.sample[CNV.sample$chrom == chr, ]
                if (nrow(chr_cnv) > 1) {
                    chr_cnv <- chr_cnv[order(chr_cnv$start), ]
                    for (i in 1:(nrow(chr_cnv)-1)) {
                        # Check if segments are truly adjacent (no gap) AND have same CN
                        if (chr_cnv$end[i] >= chr_cnv$start[i+1] &&
                            chr_cnv$total_cn[i] == chr_cnv$total_cn[i+1]) {
                            duplicate_cn <- duplicate_cn + 1
                        }
                    }
                }
            }

            if (duplicate_cn > 0) {
                warnings_list <- c(warnings_list,
                                 sprintf("Found %d truly adjacent CNV segments with identical copy numbers. Consider merging these segments.", duplicate_cn))
            }
        }
    } else {
        warnings_list <- c(warnings_list,
                         "No CNV data provided. CN oscillation analysis will be skipped.")
    }

    issues$warnings <- warnings_list
    issues$has_issues <- length(warnings_list) > 0

    if (verbose) {
        cat("\nDATA QUALITY CHECK\n")
        cat(rep("=", 50), "\n", sep="")
        cat(sprintf("Total SVs: %d\n", issues$total_SVs))
        cat(sprintf("  - Intrachromosomal: %d\n", issues$intrachromosomal_SVs))
        cat(sprintf("  - Interchromosomal: %d\n", issues$interchromosomal_SVs))

        if (!is.null(CNV.sample)) {
            cat(sprintf("Total CNV segments: %d\n", issues$total_CNV_segments))
        }

        if (issues$has_issues) {
            cat("\nWARNINGS:\n")
            for (i in 1:length(warnings_list)) {
                cat(sprintf("  [%d] %s\n", i, warnings_list[i]))
            }
        } else {
            cat("\nNo data quality issues detected.\n")
        }
        cat(rep("=", 50), "\n\n", sep="")
    }

    return(invisible(issues))
}
