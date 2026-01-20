##############################################################################
## Chromoplexy Detection and Chromoanagenesis Analysis Examples
##############################################################################
## This file demonstrates how to use the new chromoplexy detection features
## and the integrated chromoanagenesis analysis in ShatterSeek.
##############################################################################

library(OncoImplexus)

# Load example data
data(DO17373)

##############################################################################
## 1. Data Preparation
##############################################################################

# Prepare SV data
SV_data <- SVs(chrom1=as.character(SV_DO17373$chrom1),
               pos1=as.numeric(SV_DO17373$start1),
               chrom2=as.character(SV_DO17373$chrom2),
               pos2=as.numeric(SV_DO17373$end2),
               SVtype=as.character(SV_DO17373$svclass),
               strand1=as.character(SV_DO17373$strand1),
               strand2=as.character(SV_DO17373$strand2))

# Prepare CNV data
CN_data <- CNVsegs(chrom=as.character(SCNA_DO17373$chromosome),
                   start=SCNA_DO17373$start,
                   end=SCNA_DO17373$end,
                   total_cn=SCNA_DO17373$total_cn)

##############################################################################
## 2. Chromoplexy Detection Only
##############################################################################
## Use this when you specifically want to detect chromoplexy patterns

cat("=== CHROMOPLEXY DETECTION ===\n\n")

# Basic chromoplexy detection
chromoplexy_results <- detect_chromoplexy(
    SV.sample = SV_data,
    CNV.sample = CN_data,
    min_chromosomes = 3,         # Minimum chromosomes in chain
    min_translocations = 3,      # Minimum translocations
    max_cn_change = 1,           # Maximum CN deviation
    allow_cycles = TRUE          # Allow circular chains
)

# View results
print(chromoplexy_results)

# Access summary
cat("\nChromoplexy Summary:\n")
print(chromoplexy_results$summary)

# Check for specific patterns
if (chromoplexy_results$likely_chromoplexy > 0) {
    cat("\nLikely chromoplexy events found!\n")

    likely_chains <- chromoplexy_results$summary[
        chromoplexy_results$summary$classification == "Likely chromoplexy",
    ]

    cat("\nDetails of likely chromoplexy chains:\n")
    print(likely_chains)

    # Check for cycles
    cycles <- likely_chains[likely_chains$is_cycle, ]
    if (nrow(cycles) > 0) {
        cat(sprintf("\nFound %d circular chromoplexy chain(s)!\n", nrow(cycles)))
    }
}

##############################################################################
## 3. Visualize Chromoplexy Chains
##############################################################################

if (chromoplexy_results$total_chains > 0) {
    cat("\nPlotting chromoplexy chains...\n")

    # Plot all chains
    chromoplexy_plots <- plot_chromoplexy(
        chromoplexy_result = chromoplexy_results,
        sample_name = "DO17373",
        genome = "hg19"
    )

    # If you have a specific chain ID to plot
    # plot_chromoplexy(chromoplexy_results, chain_id = 1, sample_name = "DO17373")

    # For circular chains, use circular visualization
    for (i in 1:nrow(chromoplexy_results$summary)) {
        if (chromoplexy_results$summary$is_cycle[i]) {
            cat(sprintf("\nPlotting circular chain %d...\n", i))

            circular_plot <- plot_chromoplexy_circular(
                chromoplexy_result = chromoplexy_results,
                chain_id = i,
                sample_name = "DO17373"
            )

            # Display or save the plot
            # print(circular_plot)
            # ggsave(sprintf("chromoplexy_cycle_%d.pdf", i), circular_plot)
        }
    }
}

##############################################################################
## 4. Comprehensive Chromoanagenesis Analysis
##############################################################################
## Use this for a complete analysis of both chromothripsis and chromoplexy

cat("\n\n=== COMPREHENSIVE CHROMOANAGENESIS ANALYSIS ===\n\n")

# Run integrated analysis
chromoanagenesis_results <- detect_chromoanagenesis(
    SV.sample = SV_data,
    CNV.sample = CN_data,
    genome = "hg19",
    detect_chromothripsis = TRUE,
    detect_chromoplexy = TRUE,
    min_chromothripsis_size = 1,
    min_chromoplexy_chromosomes = 3,
    verbose = TRUE
)

# View integrated results
print(chromoanagenesis_results)

# Detailed summary
summary(chromoanagenesis_results)

# Access specific components
cat("\n\nAccessing specific results:\n")
cat("Overall classification:", chromoanagenesis_results$integrated_summary$overall_classification, "\n")

##############################################################################
## 5. Comparing Chromothripsis and Chromoplexy Results
##############################################################################

compare_results <- function(results) {
    cat("\n")
    cat(rep("=", 70), "\n", sep = "")
    cat("         COMPARISON: CHROMOTHRIPSIS vs CHROMOPLEXY\n")
    cat(rep("=", 70), "\n\n", sep = "")

    # Chromothripsis
    if (!is.null(results$chromothripsis)) {
        cat("CHROMOTHRIPSIS:\n")
        cat(sprintf("  Likely events:      %d\n", results$chromothripsis$n_likely))
        cat(sprintf("  Affected chromsomes: %s\n",
                   results$integrated_summary$chromothripsis_chromosomes))

        if (results$chromothripsis$n_likely > 0) {
            likely_ct <- results$chromothripsis$classification[
                results$chromothripsis$classification$classification == "Likely chromothripsis",
            ]

            avg_cluster_size <- mean(likely_ct$clusterSize, na.rm = TRUE)
            avg_cn_osc <- mean(likely_ct$max_number_oscillating_CN_segments_2_states, na.rm = TRUE)

            cat(sprintf("  Avg cluster size:   %.1f SVs\n", avg_cluster_size))
            cat(sprintf("  Avg CN oscillations: %.1f segments\n", avg_cn_osc))
        }
        cat("\n")
    }

    # Chromoplexy
    if (!is.null(results$chromoplexy)) {
        cat("CHROMOPLEXY:\n")
        cat(sprintf("  Likely events:      %d\n", results$chromoplexy$likely_chromoplexy))
        cat(sprintf("  Total chains:       %d\n", results$chromoplexy$total_chains))
        cat(sprintf("  Involved chromosomes: %s\n",
                   results$integrated_summary$chromoplexy_chromosomes))

        if (results$chromoplexy$likely_chromoplexy > 0) {
            likely_cp <- results$chromoplexy$summary[
                results$chromoplexy$summary$classification == "Likely chromoplexy",
            ]

            avg_n_chr <- mean(likely_cp$n_chromosomes, na.rm = TRUE)
            avg_n_tlx <- mean(likely_cp$n_translocations, na.rm = TRUE)
            avg_cn_stability <- mean(likely_cp$cn_stability_score, na.rm = TRUE)

            cat(sprintf("  Avg chromosomes:    %.1f\n", avg_n_chr))
            cat(sprintf("  Avg translocations: %.1f\n", avg_n_tlx))
            cat(sprintf("  Avg CN stability:   %.2f\n", avg_cn_stability))
        }
        cat("\n")
    }

    # Interpretation
    cat("INTERPRETATION:\n")

    has_ct <- !is.null(results$chromothripsis) && results$chromothripsis$n_likely > 0
    has_cp <- !is.null(results$chromoplexy) && results$chromoplexy$likely_chromoplexy > 0

    if (has_ct && has_cp) {
        cat("  This sample shows evidence of BOTH chromothripsis and chromoplexy,\n")
        cat("  indicating multiple mechanisms of complex chromosomal rearrangement.\n")
    } else if (has_ct) {
        cat("  This sample shows primarily CHROMOTHRIPSIS:\n")
        cat("  - Localized shattering and reassembly\n")
        cat("  - Oscillating copy numbers\n")
        cat("  - Confined to specific chromosomal regions\n")
    } else if (has_cp) {
        cat("  This sample shows primarily CHROMOPLEXY:\n")
        cat("  - Chained translocations across multiple chromosomes\n")
        cat("  - Stable copy numbers\n")
        cat("  - Inter-chromosomal complexity\n")
    } else {
        cat("  No clear evidence of chromoanagenesis detected.\n")
        cat("  The sample may have:\n")
        cat("  - Simple/sporadic SVs\n")
        cat("  - Alternative complex rearrangement mechanisms\n")
        cat("  - Insufficient data for detection\n")
    }

    cat("\n")
    cat(rep("=", 70), "\n", sep = "")
    cat("\n")
}

# Run comparison
compare_results(chromoanagenesis_results)

##############################################################################
## 6. Export Results
##############################################################################

# Export chromoplexy results
if (chromoplexy_results$total_chains > 0) {
    # Save summary table
    write.csv(chromoplexy_results$summary,
             "chromoplexy_summary.csv",
             row.names = FALSE)

    cat("\nChromoplexy results exported to chromoplexy_summary.csv\n")
}

# Export chromoanagenesis integrated results
write.csv(chromoanagenesis_results$integrated_summary,
         "chromoanagenesis_summary.csv",
         row.names = FALSE)

if (!is.null(chromoanagenesis_results$chromothripsis)) {
    write.csv(chromoanagenesis_results$chromothripsis$classification,
             "chromothripsis_classification.csv",
             row.names = FALSE)
}

cat("Integrated results exported to chromoanagenesis_summary.csv\n")

##############################################################################
## 7. Advanced Usage: Custom Thresholds
##############################################################################

cat("\n\n=== ADVANCED: CUSTOM THRESHOLDS ===\n\n")

# Stringent chromoplexy detection (for high-confidence events only)
stringent_chromoplexy <- detect_chromoplexy(
    SV.sample = SV_data,
    CNV.sample = CN_data,
    min_chromosomes = 4,         # More chromosomes required
    min_translocations = 5,      # More translocations required
    max_cn_change = 0.5,         # Stricter CN stability
    allow_cycles = TRUE
)

cat("Stringent chromoplexy detection:\n")
print(stringent_chromoplexy)

# Lenient chromoplexy detection (for exploratory analysis)
lenient_chromoplexy <- detect_chromoplexy(
    SV.sample = SV_data,
    CNV.sample = CN_data,
    min_chromosomes = 2,         # Fewer chromosomes allowed
    min_translocations = 2,      # Fewer translocations allowed
    max_cn_change = 2,           # More lenient CN changes
    allow_cycles = TRUE
)

cat("\n\nLenient chromoplexy detection:\n")
print(lenient_chromoplexy)

##############################################################################
## End of Examples
##############################################################################

cat("\n\nAnalysis complete!\n\n")
