#' Advanced Features Usage Examples
#'
#' This script demonstrates the usage of advanced ShatterSeek features:
#' 1. Integrated classifier for mixed chromoanagenesis mechanisms
#' 2. Breakpoint sequence analysis for DNA repair mechanism inference
#'
#' @author ShatterSeek Development Team
#' @date 2025

library(OncoImplexus)

# ============================================================================
# Example 1: Integrated Classifier for Mixed Mechanisms
# ============================================================================

# Load example data
data(DO17373)

# Prepare SV and CNV data
SV_data <- SVs(
    chrom1 = SV_DO17373$chrom1,
    pos1 = SV_DO17373$pos1,
    strand1 = SV_DO17373$strand1,
    chrom2 = SV_DO17373$chrom2,
    pos2 = SV_DO17373$pos2,
    strand2 = SV_DO17373$strand2,
    SVtype = SV_DO17373$SVtype
)

CN_data <- CNVsegs(
    chrom = SCNA_DO17373$chromosome,
    start = SCNA_DO17373$start,
    end = SCNA_DO17373$end,
    total_cn = SCNA_DO17373$total_cn
)

# ---- Step 1: Run comprehensive chromoanagenesis analysis ----
cat("\n=== Running Comprehensive Chromoanagenesis Analysis ===\n")

chromoanag_results <- detect_chromoanagenesis(
    SV.sample = SV_data,
    CNV.sample = CN_data,
    genome = "hg19",
    detect_chromothripsis = TRUE,
    detect_chromoplexy = TRUE,
    detect_chromoanasynthesis = TRUE,
    verbose = TRUE
)

# View overall results
print(chromoanag_results)
summary(chromoanag_results)


# ---- Step 2: Classify mixed mechanisms ----
cat("\n=== Classifying Mixed Mechanisms ===\n")

mixed_class <- classify_mixed_mechanisms(
    chromoanagenesis_result = chromoanag_results,
    overlap_threshold = 1e6,  # 1 Mb
    min_confidence = 0.3
)

# View classification results
print(mixed_class)
summary(mixed_class)

# Check if sample has mixed mechanisms
if (mixed_class$sample_classification$category == "mixed_overlapping") {
    cat("\nThis sample shows SPATIAL OVERLAP of different mechanisms!\n")
    cat("This suggests complex catastrophic rearrangements.\n")
} else if (mixed_class$sample_classification$category == "mixed_independent") {
    cat("\nThis sample has MULTIPLE INDEPENDENT mechanisms.\n")
} else {
    cat("\nThis sample shows a single dominant mechanism.\n")
}


# ---- Step 3: Visualize mixed mechanisms ----
cat("\n=== Creating Visualizations ===\n")

# 3.1 Mechanism landscape across chromosomes
p1 <- plot_mechanism_landscape(
    mixed_mechanisms_result = mixed_class,
    sample_name = "DO17373",
    show_confidence = TRUE
)
print(p1)

# 3.2 Chromosome-level mechanism heatmap
p2 <- plot_chromosome_mechanisms(
    mixed_mechanisms_result = mixed_class,
    sample_name = "DO17373"
)
print(p2)

# 3.3 Complexity breakdown
p3 <- plot_complexity_breakdown(
    mixed_mechanisms_result = mixed_class,
    sample_name = "DO17373"
)
print(p3)

# 3.4 Mechanism dominance pie chart
p4 <- plot_mechanism_dominance(
    mixed_mechanisms_result = mixed_class,
    sample_name = "DO17373"
)
print(p4)

# 3.5 Comprehensive report (combines all plots)
comprehensive_plot <- plot_mechanism_report(
    mixed_mechanisms_result = mixed_class,
    sample_name = "DO17373"
)
# Save to file
# ggsave("mechanism_report.pdf", comprehensive_plot, width = 14, height = 10)


# ---- Step 4: Analyze specific aspects ----

# Get chromosomes with mixed mechanisms
mixed_chroms <- mixed_class$chromosome_classification[
    mixed_class$chromosome_classification$is_mixed,
]

if (nrow(mixed_chroms) > 0) {
    cat("\n=== Chromosomes with Mixed Mechanisms ===\n")
    print(mixed_chroms)

    cat("\nThese chromosomes show evidence of multiple catastrophic events:\n")
    for (i in 1:nrow(mixed_chroms)) {
        cat(sprintf("  - %s: %s (dominant: %s, confidence: %.2f)\n",
                   mixed_chroms$chrom[i],
                   mixed_chroms$mechanisms[i],
                   mixed_chroms$dominant_mechanism[i],
                   mixed_chroms$max_confidence[i]))
    }
}

# Check complexity level
cat("\n=== Sample Complexity Assessment ===\n")
complexity <- mixed_class$complexity
cat(sprintf("Complexity Level: %s (score: %.3f)\n",
           complexity$complexity_level,
           complexity$complexity_score))
cat(sprintf("Total events: %d across %d chromosomes\n",
           complexity$total_events,
           complexity$n_chromosomes))

if (complexity$complexity_level %in% c("High", "Very High")) {
    cat("\nWARNING: This sample shows HIGH COMPLEXITY chromoanagenesis!\n")
    cat("Consider multiple catastrophic events or ongoing genomic instability.\n")
}


# ============================================================================
# Example 2: Breakpoint Sequence Analysis
# ============================================================================

cat("\n\n")
cat("============================================================================\n")
cat("=== Breakpoint Sequence Analysis ===\n")
cat("============================================================================\n\n")

# NOTE: Breakpoint sequence analysis requires a reference genome
# You can use BSgenome objects or FASTA files

# Example with BSgenome (uncomment if you have the package installed):
# library(BSgenome.Hsapiens.UCSC.hg19)
# genome_ref <- BSgenome.Hsapiens.UCSC.hg19

# For this example, we'll demonstrate the workflow
# (actual sequence extraction will return placeholders without genome)

# ---- Step 1: Analyze breakpoint sequences ----
cat("\n=== Analyzing Breakpoint Sequences ===\n")

# bp_analysis <- analyze_breakpoint_sequences(
#     SV.sample = SV_data,
#     genome = genome_ref,  # Provide BSgenome object
#     flank_size = 50,      # Extract 50 bp around each breakpoint
#     min_microhomology = 2,
#     max_microhomology = 25
# )

# Placeholder for demonstration (without actual genome)
cat("NOTE: To run breakpoint sequence analysis, you need:\n")
cat("  1. Install BSgenome package for your reference genome:\n")
cat("     BiocManager::install('BSgenome.Hsapiens.UCSC.hg19')\n")
cat("  2. Or provide a FASTA file path\n\n")

# The following shows the expected workflow:

# ---- Step 2: View results ----
# print(bp_analysis)
# summary(bp_analysis)

# Expected output:
# - Number of breakpoints analyzed
# - Proportion with microhomology
# - Mean/median microhomology length
# - Distribution of repair mechanisms


# ---- Step 3: Examine microhomology patterns ----
# Breakpoints with microhomology suggest MMEJ or alt-EJ repair
# mh_breakpoints <- bp_analysis$microhomology[
#     bp_analysis$microhomology$has_microhomology == TRUE,
# ]
#
# cat("\n=== Microhomology Pattern Analysis ===\n")
# cat(sprintf("Breakpoints with microhomology: %d\n", nrow(mh_breakpoints)))
# cat(sprintf("Mean MH length: %.1f bp\n",
#            mean(mh_breakpoints$microhomology_length)))
#
# if (nrow(mh_breakpoints) > 0) {
#     cat("\nTop microhomology sequences:\n")
#     print(head(mh_breakpoints[, c("sv_id", "microhomology_length",
#                                   "microhomology_seq")]))
# }


# ---- Step 4: Analyze repair mechanisms ----
# Understand which DNA repair pathways are active
# repair_dist <- table(bp_analysis$repair_mechanisms$repair_mechanism)
#
# cat("\n=== DNA Repair Mechanism Distribution ===\n")
# print(repair_dist)
#
# # Interpret results:
# # - High NHEJ: Classic double-strand break repair
# # - High MMEJ: Alternative end-joining, common in cancer
# # - FoSTeS/MMBIR: Serial replication, indicates chromoanasynthesis
# # - HR-like: Homologous recombination


# ---- Step 5: Visualize breakpoint features ----
# p_repair <- plot_repair_mechanisms(
#     breakpoint_result = bp_analysis,
#     sample_name = "DO17373"
# )
# print(p_repair)

# p_mh_dist <- plot_microhomology_distribution(
#     breakpoint_result = bp_analysis,
#     sample_name = "DO17373"
# )
# print(p_mh_dist)

# p_repair_by_sv <- plot_repair_by_svtype(
#     breakpoint_result = bp_analysis,
#     sv_data = SV_data,
#     sample_name = "DO17373"
# )
# print(p_repair_by_sv)

# Comprehensive breakpoint report
# bp_report <- plot_breakpoint_report(
#     breakpoint_result = bp_analysis,
#     sv_data = SV_data,
#     sample_name = "DO17373"
# )
# ggsave("breakpoint_report.pdf", bp_report, width = 12, height = 10)


# ============================================================================
# Example 3: Combined Analysis - Integrating Both Features
# ============================================================================

cat("\n\n")
cat("============================================================================\n")
cat("=== Integrated Analysis Workflow ===\n")
cat("============================================================================\n\n")

# The power of these features is in combining them:
# 1. Use chromoanagenesis detection to identify mechanism type
# 2. Use breakpoint analysis to understand repair pathways
# 3. Use mixed mechanism classifier to detect complex events

# Example interpretation workflow:
interpret_results <- function(chromoanag, mixed_class) {
    cat("\n=== Integrated Interpretation ===\n\n")

    # Check overall pattern
    cat("1. CHROMOANAGENESIS PATTERN:\n")
    cat(sprintf("   %s\n", mixed_class$sample_classification$classification))

    # Check dominant mechanism
    cat("\n2. DOMINANT MECHANISM:\n")
    cat(sprintf("   %s\n", mixed_class$dominance$dominant_mechanism))

    # Check complexity
    cat("\n3. COMPLEXITY:\n")
    cat(sprintf("   %s (score: %.2f)\n",
               mixed_class$complexity$complexity_level,
               mixed_class$complexity$complexity_score))

    # Clinical interpretation
    cat("\n4. CLINICAL INTERPRETATION:\n")

    if (mixed_class$sample_classification$category == "mixed_overlapping") {
        cat("   - Multiple overlapping catastrophic events detected\n")
        cat("   - Suggests chromothripsis crisis or ongoing instability\n")
        cat("   - High risk of therapy resistance\n")
    } else if (chromoanag$chromothripsis$n_likely > 0) {
        cat("   - Evidence of chromothripsis (catastrophic shattering)\n")
        cat("   - Single catastrophic event\n")
        cat("   - May respond to targeted therapy\n")
    } else if (chromoanag$chromoplexy$likely_chromoplexy > 0) {
        cat("   - Evidence of chromoplexy (chained rearrangements)\n")
        cat("   - Common in certain cancer types (e.g., prostate)\n")
    }

    # Repair mechanism insights (when available)
    cat("\n5. DNA REPAIR INSIGHTS:\n")
    cat("   [Run breakpoint sequence analysis for detailed repair mechanism data]\n")

    cat("\n")
}

# Run interpretation
interpret_results(chromoanag_results, mixed_class)


# ============================================================================
# Example 4: Batch Analysis
# ============================================================================

# For analyzing multiple samples:
batch_analysis <- function(sample_list) {
    results <- list()

    for (i in seq_along(sample_list)) {
        sample <- sample_list[[i]]

        cat(sprintf("\nAnalyzing sample %d/%d: %s\n",
                   i, length(sample_list), sample$name))

        # Run comprehensive analysis
        chromoanag <- detect_chromoanagenesis(
            SV.sample = sample$SV,
            CNV.sample = sample$CNV,
            verbose = FALSE
        )

        # Classify mechanisms
        mixed <- classify_mixed_mechanisms(
            chromoanagenesis_result = chromoanag,
            overlap_threshold = 1e6,
            min_confidence = 0.3
        )

        results[[sample$name]] <- list(
            chromoanagenesis = chromoanag,
            mixed_classification = mixed,
            complexity = mixed$complexity$complexity_level,
            dominant_mechanism = mixed$dominance$dominant_mechanism
        )
    }

    return(results)
}


# ============================================================================
# Summary
# ============================================================================

cat("\n")
cat("============================================================================\n")
cat("=== Analysis Complete ===\n")
cat("============================================================================\n\n")

cat("You have successfully demonstrated:\n")
cat("  ✓ Comprehensive chromoanagenesis detection\n")
cat("  ✓ Mixed mechanism classification\n")
cat("  ✓ Complexity scoring and visualization\n")
cat("  ✓ Breakpoint sequence analysis workflow\n")
cat("  ✓ Integrated interpretation approach\n\n")

cat("For more information, see:\n")
cat("  - ?detect_chromoanagenesis\n")
cat("  - ?classify_mixed_mechanisms\n")
cat("  - ?analyze_breakpoint_sequences\n\n")
