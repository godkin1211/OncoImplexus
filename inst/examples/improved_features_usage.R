##############################################################################
## Example Usage of ShatterSeek Improved Features
##############################################################################
## This file demonstrates how to use the new intermediate-level improvements
## added to ShatterSeek for better chromothripsis detection and analysis.
##############################################################################

library(OncoImplexus)

# Load example data
data(DO17373)

##############################################################################
## 1. Data Quality Check (NEW)
##############################################################################
## Before running ShatterSeek, check your data quality

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

# Check data quality (NEW FUNCTION)
quality_report <- check_data_quality(SV_data, CN_data, verbose = TRUE)

# Examine quality metrics
print(quality_report)

##############################################################################
## 2. Run ShatterSeek with Improved Validation
##############################################################################
## The detect_chromothripsis function now includes:
## - SVtype-strand consistency validation
## - CNV data completeness checks
## - Warnings for data quality issues

chromothripsis <- detect_chromothripsis(SV.sample=SV_data,
                                        seg.sample=CN_data,
                                        genome="hg19")

# View basic results
print(chromothripsis)

##############################################################################
## 3. Calculate Confidence Scores (NEW)
##############################################################################
## Get integrated confidence scores for chromothripsis events

confidence_scores <- calculate_confidence_score(chromothripsis@chromSummary)

# View confidence scores
print(confidence_scores)

# Filter high-confidence events
high_confidence <- confidence_scores[confidence_scores$confidence_category == "High", ]
print(high_confidence)

##############################################################################
## 4. Classify Chromothripsis Events (NEW)
##############################################################################
## Apply evidence-based criteria to classify events

classifications <- classify_chromothripsis(
    chromoth_output = chromothripsis,
    min_cluster_size = 6,          # Minimum SVs required
    min_cn_oscillations = 4,       # Minimum CN oscillating segments
    max_pval_clustering = 0.05     # Max p-value for clustering test
)

# View classifications
print(classifications)

# Get only likely chromothripsis events
likely_events <- classifications[
    classifications$classification == "Likely chromothripsis",
]
print(likely_events)

##############################################################################
## 5. Generate Summary Report (NEW)
##############################################################################
## Create a human-readable summary of results

summary_report <- summarize_chromothripsis(chromothripsis, print_summary = TRUE)

# Access detailed results programmatically
cat("\nNumber of likely chromothripsis events:",
    summary_report$likely_chromothripsis, "\n")

cat("High-confidence chromosomes:",
    paste(summary_report$high_confidence_chromosomes, collapse=", "), "\n")

##############################################################################
## 6. Example Workflow: Complete Analysis Pipeline
##############################################################################

complete_analysis <- function(SV_data, CN_data, sample_name = "Sample") {

    cat("\n", rep("=", 70), "\n", sep="")
    cat("  CHROMOTHRIPSIS ANALYSIS PIPELINE:", sample_name, "\n")
    cat(rep("=", 70), "\n\n", sep="")

    # Step 1: Quality check
    cat("Step 1: Checking data quality...\n")
    quality <- check_data_quality(SV_data, CN_data, verbose = FALSE)

    if (quality$has_issues) {
        cat("WARNING: Data quality issues detected!\n")
        cat("Number of warnings:", length(quality$warnings), "\n\n")
    } else {
        cat("Data quality: OK\n\n")
    }

    # Step 2: Run chromothripsis detection
    cat("Step 2: Running chromothripsis detection...\n")
    results <- detect_chromothripsis(SV.sample = SV_data,
                                     seg.sample = CN_data,
                                     genome = "hg19")
    cat("\n")

    # Step 3: Calculate confidence scores
    cat("Step 3: Calculating confidence scores...\n")
    scores <- calculate_confidence_score(results@chromSummary)
    cat(sprintf("Found %d chromosomes with SV clusters\n\n", nrow(scores)))

    # Step 4: Classify events
    cat("Step 4: Classifying chromothripsis events...\n")
    classifications <- classify_chromothripsis(results)

    n_likely <- sum(classifications$classification == "Likely chromothripsis")
    n_possible <- sum(classifications$classification == "Possible chromothripsis")

    cat(sprintf("  - Likely chromothripsis: %d\n", n_likely))
    cat(sprintf("  - Possible chromothripsis: %d\n\n", n_possible))

    # Step 5: Generate summary
    cat("Step 5: Generating final report...\n")
    summary <- summarize_chromothripsis(results, print_summary = TRUE)

    # Return all results
    return(list(
        detection_output = results,
        confidence_scores = scores,
        classifications = classifications,
        summary = summary,
        quality_report = quality
    ))
}

# Run complete analysis
analysis_results <- complete_analysis(SV_data, CN_data, "DO17373")

##############################################################################
## 7. Accessing Specific Results
##############################################################################

# Get chromosomes classified as likely chromothripsis
likely_chroms <- analysis_results$classifications[
    analysis_results$classifications$classification == "Likely chromothripsis",
    "chrom"
]

cat("\nChromosomes with likely chromothripsis:\n")
print(likely_chroms)

# Plot high-confidence chromothripsis events
if (length(likely_chroms) > 0) {
    for (chr in likely_chroms) {
        cat(sprintf("\nPlotting %s...\n", chr))
        plots <- plot_chromothripsis(
            ShatterSeek_output = analysis_results$detection_output,
            chr = chr,
            sample_name = "DO17373",
            genome = "hg19"
        )

        # Combine plots (requires gridExtra)
        # combined_plot <- arrangeGrob(plots[[1]], plots[[2]], plots[[3]], plots[[4]],
        #                             nrow=4, ncol=1, heights=c(0.2, 0.4, 0.4, 0.4))
        # grid.draw(combined_plot)
    }
}

##############################################################################
## 8. Custom Thresholds for Classification
##############################################################################

# Use more stringent criteria
stringent_classification <- classify_chromothripsis(
    chromoth_output = chromothripsis,
    min_cluster_size = 10,         # Higher threshold
    min_cn_oscillations = 6,       # More oscillations required
    max_pval_clustering = 0.01     # More stringent p-value
)

cat("\nStringent classification results:\n")
print(stringent_classification)

# Use more lenient criteria
lenient_classification <- classify_chromothripsis(
    chromoth_output = chromothripsis,
    min_cluster_size = 3,          # Lower threshold
    min_cn_oscillations = 2,       # Fewer oscillations required
    max_pval_clustering = 0.1      # More lenient p-value
)

cat("\nLenient classification results:\n")
print(lenient_classification)

##############################################################################
## End of Examples
##############################################################################
