##############################################################################
## VCF Workflow Examples for ShatterSeek
##############################################################################
## This script demonstrates how to use ShatterSeek directly with VCF files
## from various SV and CNV callers, without manual data conversion.
##
## Author: ShatterSeek Development Team
## Date: 2025
##############################################################################

library(OncoImplexus)

##############################################################################
## Example 1: Basic VCF Workflow with Auto-detection
##############################################################################

# Read SV and CNV data from VCF files
# The caller type will be auto-detected
sv_data <- read_sv_vcf("sample.sv.vcf.gz")
cnv_data <- read_cnv_vcf("sample.cnv.vcf.gz")

# Run chromoanagenesis detection
results <- detect_chromoanagenesis(
    SV.sample = sv_data,
    CNV.sample = cnv_data,
    genome = "hg19"
)

# View results
print(results$chromothripsis$classification)
print(results$chromoplexy$summary)
print(results$chromoanasynthesis$summary)


##############################################################################
## Example 2: Manta SV Caller
##############################################################################

# Manta produces high-quality SV calls for Illumina data
# Typical output: sample.manta.vcf.gz

# Read Manta VCF (explicit caller specification)
sv_manta <- read_sv_vcf(
    vcf_file = "sample.manta.vcf.gz",
    caller = "manta",
    min_sv_size = 1000,      # Filter SVs < 1kb
    include_tra = TRUE        # Include translocations
)

# Read CNV data (e.g., from CNVkit)
cnv_cnvkit <- read_cnv_vcf(
    vcf_file = "sample.cnvkit.vcf.gz",
    caller = "cnvkit",
    merge_adjacent = TRUE     # Merge adjacent segments with same CN
)

# Detect chromothripsis
chromothripsis <- detect_chromothripsis(
    SV.sample = sv_manta,
    seg.sample = cnv_cnvkit,
    genome = "hg38"
)

# Calculate confidence scores
scores <- calculate_confidence_score(chromothripsis@chromSummary)
classification <- classify_chromothripsis(chromothripsis)

# Plot high-confidence events
likely_chroms <- classification$chrom[
    classification$classification == "Likely chromothripsis"
]

for (chr in likely_chroms) {
    plots <- plot_chromothripsis(
        ShatterSeek_output = chromothripsis,
        chr = chr,
        sample_name = "Manta_Sample",
        genome = "hg38"
    )
    print(plots)
}


##############################################################################
## Example 3: Delly SV Caller
##############################################################################

# Delly is another popular SV caller
# Uses CT (connection type) field for strand information

sv_delly <- read_sv_vcf(
    vcf_file = "sample.delly.vcf.gz",
    caller = "delly",
    min_sv_size = 500
)

# Can use GATK CNV caller output
cnv_gatk <- read_cnv_vcf(
    vcf_file = "sample.gatk.cnv.vcf.gz",
    caller = "gatk"
)

# Run integrated analysis
results_delly <- detect_chromoanagenesis(
    SV.sample = sv_delly,
    CNV.sample = cnv_gatk,
    genome = "hg19",
    min_chromothripsis_size = 6,
    verbose = TRUE
)


##############################################################################
## Example 4: Long-read Data (Sniffles)
##############################################################################

# Sniffles is designed for long-read sequencing (PacBio, ONT)
# Often produces larger SVs with different characteristics

sv_sniffles <- read_sv_vcf(
    vcf_file = "sample.sniffles.vcf.gz",
    caller = "sniffles",
    min_sv_size = 50    # Long-reads can detect smaller SVs accurately
)

# For long-read data, CNVs might come from specialized callers
cnv_lr <- read_cnv_vcf(
    vcf_file = "sample.longread.cnv.vcf.gz",
    cn_field = "CN"     # Explicitly specify CN field
)

# Detect chromoanagenesis
results_lr <- detect_chromoanagenesis(
    SV.sample = sv_sniffles,
    CNV.sample = cnv_lr,
    genome = "hg38"
)


##############################################################################
## Example 5: Multi-sample VCF
##############################################################################

# For multi-sample VCFs, specify which sample to extract

sv_multi <- read_sv_vcf(
    vcf_file = "cohort.sv.vcf.gz",
    sample_name = "SAMPLE_001",  # Extract specific sample
    min_sv_size = 1000
)

cnv_multi <- read_cnv_vcf(
    vcf_file = "cohort.cnv.vcf.gz",
    sample_name = "SAMPLE_001"
)

results_multi <- detect_chromoanagenesis(
    SV.sample = sv_multi,
    CNV.sample = cnv_multi
)


##############################################################################
## Example 6: Complete Workflow with VCF Input
##############################################################################

complete_vcf_analysis <- function(sv_vcf, cnv_vcf, sample_id,
                                  caller_sv = "auto",
                                  caller_cnv = "auto",
                                  genome = "hg19") {

    cat(sprintf("\n=== Analyzing %s ===\n\n", sample_id))

    # Step 1: Read VCF files
    cat("Step 1: Reading VCF files...\n")
    sv_data <- read_sv_vcf(
        vcf_file = sv_vcf,
        caller = caller_sv,
        min_sv_size = 1000,
        include_tra = TRUE
    )

    cnv_data <- read_cnv_vcf(
        vcf_file = cnv_vcf,
        caller = caller_cnv,
        merge_adjacent = TRUE
    )

    # Step 2: Quality check
    cat("\nStep 2: Quality check...\n")
    quality <- check_data_quality(sv_data, cnv_data)
    if (!is.null(quality$issues)) {
        cat("  Warning: Data quality issues detected:\n")
        print(quality$issues)
    }

    # Step 3: Detect all chromoanagenesis types
    cat("\nStep 3: Detecting chromoanagenesis events...\n")
    results <- detect_chromoanagenesis(
        SV.sample = sv_data,
        CNV.sample = cnv_data,
        genome = genome,
        detect_chromothripsis = TRUE,
        detect_chromoplexy = TRUE,
        detect_chromoanasynthesis = TRUE,
        verbose = TRUE
    )

    # Step 4: Mixed mechanism classification
    if (!is.null(results$chromothripsis) &&
        !is.null(results$chromoplexy) &&
        !is.null(results$chromoanasynthesis)) {

        cat("\nStep 4: Classifying mixed mechanisms...\n")
        mixed <- classify_mixed_mechanisms(
            chromoanagenesis_result = results,
            overlap_threshold = 1e6
        )

        results$mixed_mechanisms <- mixed

        # Plot mechanism landscape
        landscape_plot <- plot_mechanism_landscape(mixed, sample_id)
        print(landscape_plot)
    }

    # Step 5: Breakpoint sequence analysis (if genome available)
    cat("\nStep 5: Analyzing breakpoint sequences...\n")
    if (requireNamespace("BSgenome", quietly = TRUE)) {
        genome_pkg <- if (genome == "hg19") {
            "BSgenome.Hsapiens.UCSC.hg19"
        } else {
            "BSgenome.Hsapiens.UCSC.hg38"
        }

        if (requireNamespace(genome_pkg, quietly = TRUE)) {
            genome_obj <- get(genome_pkg, asNamespace(genome_pkg))

            bp_analysis <- analyze_breakpoint_sequences(
                SV.sample = sv_data,
                genome = genome_obj,
                flank_size = 50
            )

            results$breakpoint_analysis <- bp_analysis

            # Plot repair mechanisms
            repair_plot <- plot_repair_mechanisms(bp_analysis, sample_id)
            print(repair_plot)
        }
    }

    # Step 6: Generate summary report
    cat("\nStep 6: Generating summary...\n")
    if (!is.null(results$chromothripsis)) {
        summary <- summarize_chromothripsis(
            chromoth_output = results$chromothripsis$detection_output,
            print_summary = TRUE
        )
        results$summary <- summary
    }

    cat("\n=== Analysis Complete ===\n")

    return(results)
}


# Run complete analysis
# results <- complete_vcf_analysis(
#     sv_vcf = "patient_001.manta.vcf.gz",
#     cnv_vcf = "patient_001.cnvkit.vcf.gz",
#     sample_id = "Patient_001",
#     caller_sv = "manta",
#     caller_cnv = "cnvkit",
#     genome = "hg38"
# )


##############################################################################
## Example 7: Troubleshooting VCF Format Issues
##############################################################################

# If auto-detection fails or results look wrong, try:

# 1. Check VCF format
sv_test <- read_sv_vcf(
    vcf_file = "sample.vcf.gz",
    caller = "auto"  # Will print detected caller
)

# 2. If caller is misdetected, specify manually
sv_manual <- read_sv_vcf(
    vcf_file = "sample.vcf.gz",
    caller = "manta"  # or "delly", "gridss", "lumpy", "sniffles"
)

# 3. For CNV files, check available CN field
# You can inspect the VCF header to find the correct field name
cnv_custom <- read_cnv_vcf(
    vcf_file = "sample.cnv.vcf.gz",
    cn_field = "TOTAL_CN"  # Specify custom field name
)

# 4. If VCF is non-standard, use generic parser and inspect results
sv_generic <- read_sv_vcf(
    vcf_file = "custom.vcf.gz",
    caller = "generic"  # Uses basic SVTYPE/END fields
)

# Check the extracted data
head(as.data.frame(sv_generic))


##############################################################################
## Example 8: Batch Processing Multiple VCF Files
##############################################################################

# Process multiple samples
sample_list <- data.frame(
    sample_id = c("Sample_A", "Sample_B", "Sample_C"),
    sv_vcf = c("sampleA.sv.vcf.gz", "sampleB.sv.vcf.gz", "sampleC.sv.vcf.gz"),
    cnv_vcf = c("sampleA.cnv.vcf.gz", "sampleB.cnv.vcf.gz", "sampleC.cnv.vcf.gz"),
    stringsAsFactors = FALSE
)

batch_results <- list()

for (i in 1:nrow(sample_list)) {
    cat(sprintf("\n\nProcessing %s (%d/%d)...\n",
                sample_list$sample_id[i], i, nrow(sample_list)))

    tryCatch({
        batch_results[[sample_list$sample_id[i]]] <- complete_vcf_analysis(
            sv_vcf = sample_list$sv_vcf[i],
            cnv_vcf = sample_list$cnv_vcf[i],
            sample_id = sample_list$sample_id[i],
            genome = "hg38"
        )
    }, error = function(e) {
        cat(sprintf("  ERROR: %s\n", e$message))
        batch_results[[sample_list$sample_id[i]]] <<- NULL
    })
}

# Summarize batch results
cat("\n\n=== Batch Summary ===\n")
for (sample_id in names(batch_results)) {
    result <- batch_results[[sample_id]]
    if (!is.null(result)) {
        n_chromothripsis <- if (!is.null(result$chromothripsis)) {
            sum(result$chromothripsis$classification$classification == "Likely chromothripsis")
        } else {
            0
        }

        n_chromoplexy <- if (!is.null(result$chromoplexy)) {
            nrow(result$chromoplexy$chains)
        } else {
            0
        }

        cat(sprintf("%s: %d chromothripsis, %d chromoplexy events\n",
                    sample_id, n_chromothripsis, n_chromoplexy))
    }
}


##############################################################################
## Installation Instructions for Optional Dependencies
##############################################################################

# VCF reading requires either VariantAnnotation (Bioconductor) or vcfR:

# Option 1: VariantAnnotation (recommended for comprehensive VCF support)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("VariantAnnotation")

# Option 2: vcfR (lighter alternative, CRAN package)
# install.packages("vcfR")

# For breakpoint sequence analysis, you also need BSgenome:
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
