################################################################################
# Statistical Significance Testing for Chromoplexy
# Based on ChainFinder methodology (Baca et al. 2013)
################################################################################

#' Calculate statistical significance for a chromoplexy chain
#'
#' Implements nearest-neighbor testing similar to ChainFinder.
#' Tests whether breakpoints in the chain are closer than expected by chance.
#'
#' @param chain Chain object
#' @param inter_chr_SVs All inter-chromosomal SVs
#' @param genome_size Total genome size
#' @param fdr_threshold FDR threshold for significance
#' @return List with p-values and significance assessment
#' @keywords internal
calculate_chain_significance <- function(chain,
                                            inter_chr_SVs,
                                            genome_size,
                                            fdr_threshold = 0.01) {

    # Calculate local rearrangement rate (μ)
    mu <- estimate_local_rearrangement_rate( # Updated call
        chromosomes = chain$chromosomes,
        total_SVs = nrow(inter_chr_SVs),
        genome_size = genome_size
    )

    # Calculate p-values for each consecutive pair of breakpoints
    pvalues <- c()
    distances <- c()

    if (length(chain$nodes) < 2) {
        return(list(
            pvalue = 1.0,
            fdr = 1.0,
            is_significant = FALSE,
            mean_distance = NA,
            mu = mu
        ))
    }

    for (i in 1:(length(chain$nodes) - 1)) {
        node1 <- chain$nodes[i]
        node2 <- chain$nodes[i + 1]

        # Calculate genomic distance between consecutive breakpoints
        dist <- calculate_genomic_distance(node1, node2, genome_size) # Updated call
        distances <- c(distances, dist)

        # Calculate p-value for this pair
        pval <- calculate_adjacency_pvalue( # Updated call
            distance = dist,
            mu = mu,
            genome_size = genome_size
        )

        pvalues <- c(pvalues, pval)
    }

    # Apply FDR correction (Benjamini-Hochberg)
    if (length(pvalues) > 0) {
        adjusted_pvalues <- p.adjust(pvalues, method = "BH")
        min_fdr <- min(adjusted_pvalues)

        # Combine p-values using Fisher's method
        if (all(!is.na(pvalues)) && all(pvalues > 0)) {
            fisher_stat <- -2 * sum(log(pvalues))
            df <- 2 * length(pvalues)
            chain_pvalue <- pchisq(fisher_stat, df, lower.tail = FALSE)
        } else {
            chain_pvalue <- 1.0
        }
    } else {
        chain_pvalue <- 1.0
        min_fdr <- 1.0
    }

    # Determine significance
    is_significant <- (chain_pvalue < fdr_threshold) && (min_fdr < fdr_threshold)

    return(list(
        pvalue = chain_pvalue,
        fdr = min_fdr,
        is_significant = is_significant,
        pvalues_per_edge = pvalues,
        adjusted_pvalues = if (exists("adjusted_pvalues")) adjusted_pvalues else NULL,
        mean_distance = mean(distances, na.rm = TRUE),
        mu = mu
    ))
}


#' Estimate local rearrangement rate
#'
#' Calculates the expected rate of rearrangements per megabase
#' for the chromosomes involved in the chain.
#'
#' @param chromosomes Chromosomes in chain
#' @param total_SVs Total number of SVs in sample
#' @param genome_size Genome size
#' @return Local rearrangement rate (μ)
#' @keywords internal
estimate_local_rearrangement_rate <- function(chromosomes,
                                                 total_SVs,
                                                 genome_size) {

    # Simple estimate: total SVs / genome size (in Mb)
    # This is a conservative baseline rate

    genome_size_mb <- genome_size / 1e6
    mu <- total_SVs / genome_size_mb

    # Could be refined to use chromosome-specific rates
    # or local window-based rates for better accuracy

    return(mu)
}


#' Calculate genomic distance between two breakpoints
#'
#' For inter-chromosomal breakpoints, uses a normalized distance metric.
#'
#' @param node1 First node (format: "chrom:pos")
#' @param node2 Second node
#' @param genome_size Total genome size
#' @return Genomic distance
#' @keywords internal
calculate_genomic_distance <- function(node1, node2, genome_size) {

    # Parse nodes
    parts1 <- strsplit(node1, ":")[[1]]
    chrom1 <- parts1[1]
    pos1 <- as.numeric(parts1[2])

    parts2 <- strsplit(node2, ":")[[1]]
    chrom2 <- parts2[1]
    pos2 <- as.numeric(parts2[2])

    if (chrom1 == chrom2) {
        # Same chromosome: direct distance
        distance <- abs(pos2 - pos1)
    } else {
        # Different chromosomes: use normalized inter-chromosomal distance
        # This is a simplified metric; could use more sophisticated 3D genome models
        distance <- genome_size / 10  # Arbitrary but reasonable default
    }

    return(distance)
}


#' Calculate p-value for breakpoint adjacency
#'
#' Calculates the probability that two independent breakpoints would occur
#' within the observed distance by chance.
#'
#' @param distance Observed distance between breakpoints
#' @param mu Local rearrangement rate
#' @param genome_size Genome size
#' @return P-value
#' @keywords internal
calculate_adjacency_pvalue <- function(distance, mu, genome_size) {

    # Model: Poisson process for breakpoint occurrence
    # P(at least one breakpoint in window) = 1 - exp(-λ)
    # where λ = mu * window_size (in Mb)

    window_size_mb <- distance / 1e6
    lambda <- mu * window_size_mb

    # Probability of observing another breakpoint within this distance
    # (this is the null model - independent breakpoints)
    pval <- 1 - exp(-lambda)

    # Ensure pval is in valid range
    pval <- max(pval, 1e-10)  # Avoid exactly 0
    pval <- min(pval, 1.0)

    return(pval)
}


#' Get genome size for reference genome
#'
#' @param genome Genome build ("hg19" or "hg38")
#' @return Total genome size in bp
#' @keywords internal
get_genome_size <- function(genome = "hg19") {

    if (genome == "hg38") {
        # hg38 genome size (excluding patches and alts)
        return(3088269832)
    } else {
        # hg19 genome size
        return(3137161264)
    }
}


################################################################################
# Enhanced Chain Evaluation
################################################################################

#' Evaluate chromoplexy chain with enhanced metrics
#'
#' Includes deletion bridge scoring and statistical significance.
#'
#' @param chain Chain to evaluate
#' @param inter_chr_SVs Inter-chromosomal SVs
#' @param CNV.sample CNV data
#' @param deletion_bridges Identified deletion bridges
#' @param max_cn_change Maximum allowed CN change
#' @param use_statistical_testing Enable statistical testing
#' @param genome_size Genome size
#' @param fdr_threshold FDR threshold
#' @return Enhanced evaluation results
#' @keywords internal
evaluate_chromoplexy_chain_v2 <- function(chain,
                                          inter_chr_SVs,
                                          CNV.sample,
                                          deletion_bridges = NULL,
                                          max_cn_change = 1,
                                          use_statistical_testing = TRUE,
                                          genome_size = NULL,
                                          fdr_threshold = 0.01) {

    # Get SVs in this chain
    chain_SVs <- inter_chr_SVs[chain$sv_indices, ]

    # 1. Evaluate copy number stability (enhanced)
    cn_eval <- evaluate_cn_stability_enhanced_v2(chain, chain_SVs, CNV.sample)

    # 2. Calculate chain complexity score
    complexity_score <- calculate_chain_complexity(chain, chain_SVs)

    # 3. Evaluate deletion bridge enrichment
    deletion_bridge_score <- 0
    n_bridges_in_chain <- 0

    if (!is.null(deletion_bridges) && length(deletion_bridges) > 0) {
        bridge_eval <- evaluate_deletion_bridge_enrichment_v2(
            chain = chain,
            deletion_bridges = deletion_bridges
        )
        deletion_bridge_score <- bridge_eval$enrichment_score
        n_bridges_in_chain <- bridge_eval$n_bridges
    }

    # 4. Statistical significance testing
    statistical_significance <- NULL
    if (use_statistical_testing && !is.null(genome_size)) {
        statistical_significance <- calculate_chain_significance_v2(
            chain = chain,
            inter_chr_SVs = inter_chr_SVs,
            genome_size = genome_size,
            fdr_threshold = fdr_threshold
        )
    }

    # 5. SV type diversity
    sv_types <- table(chain_SVs$SVtype)
    has_deletions <- "DEL" %in% names(sv_types)
    type_diversity <- length(sv_types)

    # 6. Calculate combined evidence score
    combined_score <- calculate_combined_evidence_score_v2(
        cn_stability_score = cn_eval$combined_score,
        complexity_score = complexity_score,
        deletion_bridge_score = deletion_bridge_score,
        statistical_significance = statistical_significance
    )

    # Create enhanced summary
    summary <- data.frame(
        chain_id = chain$id,
        n_chromosomes = chain$n_chromosomes,
        chromosomes_involved = paste(chain$chromosomes, collapse = ","),
        n_translocations = chain$n_translocations,
        is_cycle = chain$is_cycle,
        # CN metrics
        cn_stability_score = cn_eval$combined_score,
        cn_global_stability = cn_eval$components$global_stability,
        cn_local_changes = cn_eval$components$local_cn_changes,
        max_cn_deviation = cn_eval$max_deviation,
        # Complexity
        complexity_score = complexity_score,
        # Deletion bridges
        deletion_bridge_score = deletion_bridge_score,
        n_deletion_bridges = n_bridges_in_chain,
        # Statistical
        pvalue = if (!is.null(statistical_significance)) statistical_significance$pvalue else NA,
        fdr = if (!is.null(statistical_significance)) statistical_significance$fdr else NA,
        is_statistically_significant = if (!is.null(statistical_significance))
            statistical_significance$is_significant else NA,
        # Other
        has_deletions = has_deletions,
        sv_type_diversity = type_diversity,
        # Combined
        combined_score = combined_score,
        stringsAsFactors = FALSE
    )

    return(list(
        summary = summary,
        chain = chain,
        SVs = chain_SVs,
        cn_evaluation = cn_eval,
        statistical_significance = statistical_significance
    ))
}


#' Enhanced CN stability evaluation
#'
#' Evaluates both global and local CN stability patterns.
#'
#' @param chain Chain object
#' @param chain_SVs SVs in chain
#' @param CNV.sample CNV data
#' @return Enhanced CN stability metrics
#' @keywords internal
evaluate_cn_stability_enhanced_v2 <- function(chain, chain_SVs, CNV.sample) {

    components <- list()

    # 1. Global CN stability (original metric)
    cn_deviations <- c()
    for (chr in chain$chromosomes) {
        chr_cnv <- CNV.sample[CNV.sample$chrom == chr, ]
        if (nrow(chr_cnv) > 0) {
            deviations <- abs(chr_cnv$total_cn - 2)
            cn_deviations <- c(cn_deviations, deviations)
        }
    }

    max_deviation <- if (length(cn_deviations) > 0) max(cn_deviations, na.rm = TRUE) else 0
    mean_deviation <- if (length(cn_deviations) > 0) mean(cn_deviations, na.rm = TRUE) else 0
    components$global_stability <- exp(-mean_deviation / 2)

    # 2. Local CN changes near breakpoints
    local_cn_changes <- evaluate_local_cn_changes_v2(
        breakpoints = extract_breakpoints_v2(chain_SVs),
        CNV.sample = CNV.sample,
        window_size = 5e6  # 5Mb window
    )
    components$local_cn_changes <- local_cn_changes

    # 3. Deletion association
    components$deletion_association <- if ("DEL" %in% chain_SVs$SVtype) 0.8 else 0.5

    # Combined score (weighted harmonic mean for balance)
    weights <- c(0.5, 0.3, 0.2)
    values <- c(
        components$global_stability,
        1 - components$local_cn_changes,  # Convert to stability metric
        components$deletion_association
    )

    # Weighted harmonic mean
    combined_score <- 1 / sum(weights / values)

    return(list(
        combined_score = combined_score,
        components = components,
        max_deviation = max_deviation,
        mean_deviation = mean_deviation
    ))
}


#' Extract breakpoints from SVs
#'
#' @param chain_SVs SV data frame
#' @return List of breakpoints
#' @keywords internal
extract_breakpoints_v2 <- function(chain_SVs) {
    breakpoints <- list()

    for (i in 1:nrow(chain_SVs)) {
        sv <- chain_SVs[i, ]
        breakpoints[[length(breakpoints) + 1]] <- list(chrom = sv$chrom1, pos = sv$pos1)
        breakpoints[[length(breakpoints) + 1]] <- list(chrom = sv$chrom2, pos = sv$pos2)
    }

    return(breakpoints)
}


#' Evaluate local CN changes near breakpoints
#'
#' @param breakpoints List of breakpoints
#' @param CNV.sample CNV data
#' @param window_size Window size for local analysis
#' @return Local CN change score (0-1, higher = more changes)
#' @keywords internal
evaluate_local_cn_changes_v2 <- function(breakpoints, CNV.sample, window_size = 5e6) {

    if (is.null(CNV.sample) || nrow(CNV.sample) == 0) return(0)

    local_changes <- c()

    for (bp in breakpoints) {
        # Get CNV segments near this breakpoint
        chr_cnv <- CNV.sample[CNV.sample$chrom == bp$chrom, ]

        if (nrow(chr_cnv) > 0) {
            # Find segments within window
            nearby_segs <- chr_cnv[
                (chr_cnv$start >= bp$pos - window_size & chr_cnv$start <= bp$pos + window_size) |
                (chr_cnv$end >= bp$pos - window_size & chr_cnv$end <= bp$pos + window_size),
            ]

            if (nrow(nearby_segs) > 1) {
                # Calculate CN variation in this window
                cn_variation <- sd(nearby_segs$total_cn, na.rm = TRUE)
                local_changes <- c(local_changes, cn_variation)
            }
        }
    }

    # Return average local change (normalized)
    if (length(local_changes) > 0) {
        avg_change <- mean(local_changes, na.rm = TRUE)
        return(min(avg_change / 2, 1.0))  # Normalize to 0-1
    } else {
        return(0)
    }
}


#' Evaluate deletion bridge enrichment in chain
#'
#' @param chain Chain object
#' @param deletion_bridges All deletion bridges
#' @return Enrichment metrics
#' @keywords internal
evaluate_deletion_bridge_enrichment_v2 <- function(chain, deletion_bridges) {

    # Count bridges in this chain
    n_bridges <- 0
    total_confidence <- 0

    for (bridge in deletion_bridges) {
        if (bridge$sv_index %in% chain$sv_indices) {
            n_bridges <- n_bridges + 1
            total_confidence <- total_confidence + bridge$confidence
        }
    }

    # Calculate enrichment score
    if (n_bridges > 0) {
        # Proportion of translocations with bridges
        proportion_with_bridges <- n_bridges / chain$n_translocations

        # Average confidence
        avg_confidence <- total_confidence / n_bridges

        # Combined enrichment score
        enrichment_score <- (proportion_with_bridges * 0.6 + avg_confidence * 0.4)
    } else {
        enrichment_score <- 0
    }

    return(list(
        enrichment_score = enrichment_score,
        n_bridges = n_bridges,
        proportion_with_bridges = if (n_bridges > 0) n_bridges / chain$n_translocations else 0
    ))
}


#' Calculate combined evidence score
#'
#' Integrates multiple lines of evidence into overall confidence.
#'
#' @param cn_stability_score CN stability score
#' @param complexity_score Complexity score
#' @param deletion_bridge_score Deletion bridge score
#' @param statistical_significance Statistical test results
#' @return Combined score (0-1)
#' @keywords internal
calculate_combined_evidence_score_v2 <- function(cn_stability_score,
                                                 complexity_score,
                                                 deletion_bridge_score,
                                                 statistical_significance = NULL) {

    scores <- c(cn_stability_score, complexity_score, deletion_bridge_score)
    weights <- c(0.3, 0.3, 0.2)

    # Add statistical significance if available
    if (!is.null(statistical_significance) && !is.na(statistical_significance$is_significant)) {
        if (statistical_significance$is_significant) {
            # Significant: add high score
            stat_score <- 1 - min(statistical_significance$fdr * 10, 1.0)
            scores <- c(scores, stat_score)
            weights <- c(weights, 0.2)
        } else {
            # Not significant: penalize
            scores <- c(scores, 0.3)
            weights <- c(weights, 0.2)
        }
    }

    # Weighted average
    combined <- sum(scores * weights) / sum(weights)

    return(combined)
}


################################################################################
# Enhanced Classification
################################################################################

#' Classify chromoplexy events with enhanced criteria
#'
#' Uses statistical significance in addition to empirical criteria.
#'
#' @param summary_df Summary data frame
#' @param use_statistical_testing Whether statistical testing was used
#' @param fdr_threshold FDR threshold
#' @return Character vector of classifications
#' @keywords internal
classify_chromoplexy_event_v2 <- function(summary_df,
                                          use_statistical_testing = TRUE,
                                          fdr_threshold = 0.01) {

    classifications <- character(nrow(summary_df))

    for (i in 1:nrow(summary_df)) {
        row <- summary_df[i, ]

        # Core criteria (same as before)
        meets_chr_criteria <- row$n_chromosomes >= 3
        meets_tlx_criteria <- row$n_translocations >= 3
        meets_cn_criteria <- row$cn_stability_score >= 0.7
        meets_complexity <- row$complexity_score >= 0.3

        # NEW: Statistical significance criterion
        meets_statistical <- FALSE
        if (use_statistical_testing && !is.na(row$is_statistically_significant)) {
            meets_statistical <- row$is_statistically_significant
        }

        # NEW: Deletion bridge criterion
        meets_deletion_bridge <- FALSE
        if ("deletion_bridge_score" %in% colnames(summary_df)) {
            meets_deletion_bridge <- row$deletion_bridge_score >= 0.5
        }

        # Count criteria met
        all_criteria <- c(
            meets_chr_criteria,
            meets_tlx_criteria,
            meets_cn_criteria,
            meets_complexity,
            meets_statistical,
            meets_deletion_bridge
        )

        criteria_met <- sum(all_criteria)

        # Enhanced classification
        if (use_statistical_testing) {
            # Stricter classification when using statistical testing
            if (meets_statistical && criteria_met >= 5) {
                classifications[i] <- "Likely chromoplexy"
            } else if (criteria_met >= 4) {
                classifications[i] <- "Possible chromoplexy"
            } else if (criteria_met >= 3) {
                classifications[i] <- "Unlikely chromoplexy"
            } else {
                classifications[i] <- "Not chromoplexy"
            }
        } else {
            # Original classification scheme
            if (criteria_met >= 4) {
                classifications[i] <- "Likely chromoplexy"
            } else if (criteria_met >= 3) {
                classifications[i] <- "Possible chromoplexy"
            } else if (criteria_met >= 2) {
                classifications[i] <- "Unlikely chromoplexy"
            } else {
                classifications[i] <- "Not chromoplexy"
            }
        }
    }

    return(classifications)
}


#' Create empty chromoplexy result (v2)
#'
#' @return Empty result structure
#' @keywords internal
create_empty_chromoplexy_result_v2 <- function() {
    result <- list(
        chains = list(),
        chain_details = list(),
        summary = data.frame(
            chain_id = integer(0),
            n_chromosomes = integer(0),
            chromosomes_involved = character(0),
            n_translocations = integer(0),
            is_cycle = logical(0),
            cn_stability_score = numeric(0),
            max_cn_deviation = numeric(0),
            complexity_score = numeric(0),
            deletion_bridge_score = numeric(0),
            pvalue = numeric(0),
            fdr = numeric(0),
            has_deletions = logical(0),
            sv_type_diversity = integer(0),
            combined_score = numeric(0),
            classification = character(0)
        ),
        translocation_graph = NULL,
        deletion_bridges = NULL,
        total_chains = 0,
        likely_chromoplexy = 0,
        possible_chromoplexy = 0
    )

    class(result) <- c("chromoplexy_v2", "chromoplexy", "list")
    return(result)
}
