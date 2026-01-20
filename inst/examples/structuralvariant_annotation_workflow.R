#!/usr/bin/env Rscript
#
# Example workflow: Using StructuralVariantAnnotation with ShatterSeek
#
# This script demonstrates how to use StructuralVariantAnnotation package
# for robust VCF parsing and then feed the results to ShatterSeek for
# chromoanagenesis detection.
#
# When to use this approach:
# 1. Complex VCF formats with unusual BND records
# 2. Unsupported SV callers
# 3. Need additional SV QC and annotation
# 4. Already using Bioconductor ecosystem
#

library(OncoImplexus)

# Check if StructuralVariantAnnotation is available
if (!requireNamespace("StructuralVariantAnnotation", quietly = TRUE)) {
    message("Installing StructuralVariantAnnotation from Bioconductor...")
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("StructuralVariantAnnotation")
}

if (!requireNamespace("VariantAnnotation", quietly = TRUE)) {
    message("Installing VariantAnnotation from Bioconductor...")
    BiocManager::install("VariantAnnotation")
}

library(StructuralVariantAnnotation)
library(VariantAnnotation)

# ==============================================================================
# METHOD 1: Step-by-step (maximum control)
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("  METHOD 1: Step-by-Step Workflow with StructuralVariantAnnotation\n")
cat("================================================================================\n\n")

# Example VCF file path (replace with your actual file)
vcf_file <- "path/to/your/sample.vcf.gz"
genome <- "hg38"  # or "hg19"

if (file.exists(vcf_file)) {

    # Step 1: Read VCF with VariantAnnotation
    cat("Step 1: Reading VCF file...\n")
    vcf <- readVcf(vcf_file, genome)
    cat(sprintf("  ✓ Read %d variants\n\n", nrow(vcf)))

    # Step 2: Extract breakpoints with StructuralVariantAnnotation
    cat("Step 2: Extracting breakpoint ranges...\n")
    gr_bp <- breakpointRanges(vcf)
    cat(sprintf("  ✓ Found %d paired breakpoints\n", length(gr_bp)))

    cat("Step 3: Extracting breakend ranges...\n")
    gr_be <- breakendRanges(vcf)
    cat(sprintf("  ✓ Found %d unpaired breakends\n\n", length(gr_be)))

    # Step 4: Combine all breakpoints
    cat("Step 4: Combining breakpoints and breakends...\n")
    gr_all <- c(gr_bp, gr_be)
    cat(sprintf("  ✓ Total: %d breakpoint/breakend events\n\n", length(gr_all)))

    # Step 5: Convert to ShatterSeek SVs object
    cat("Step 5: Converting to ShatterSeek SVs format...\n")
    sv_data <- granges_to_svs(gr_all)
    cat("  ✓ Conversion complete\n\n")

    # Step 6: Load CNV data (use appropriate method for your data)
    cat("Step 6: Loading CNV data...\n")
    cn_data <- read_cnv_vcf("path/to/your/cnv.vcf.gz")
    cat("  ✓ CNV data loaded\n\n")

    # Step 7: Run ShatterSeek analysis
    cat("Step 7: Running chromoanagenesis detection...\n")
    results <- detect_chromoanagenesis(sv_data, cn_data, genome = genome)
    cat("  ✓ Detection complete\n\n")

    # Step 8: Classification and visualization
    cat("Step 8: Classifying and visualizing results...\n")
    classification <- classify_mixed_mechanisms(results)

    # Generate report
    plot_mechanism_report(classification, "Sample_StructuralVariantAnnotation")
    cat("  ✓ Analysis complete\n\n")

} else {
    cat("Note: Replace vcf_file path with your actual VCF file to run this example\n\n")
}


# ==============================================================================
# METHOD 2: One-step wrapper (convenience)
# ==============================================================================

cat("================================================================================\n")
cat("  METHOD 2: One-Step Wrapper (ShatterSeek v2.0.2+)\n")
cat("================================================================================\n\n")

if (file.exists(vcf_file)) {

    # Single function call that uses StructuralVariantAnnotation backend
    cat("Reading VCF with StructuralVariantAnnotation backend...\n")
    sv_data <- read_sv_vcf_structuralvariant(vcf_file, genome = genome)
    cat("  ✓ VCF reading complete\n\n")

    # Continue with standard ShatterSeek workflow
    cat("Running analysis...\n")
    results <- detect_chromoanagenesis(sv_data, cn_data, genome = genome)
    cat("  ✓ Analysis complete\n\n")

} else {
    cat("Note: Replace vcf_file path with your actual VCF file to run this example\n\n")
}


# ==============================================================================
# COMPARISON: Built-in vs StructuralVariantAnnotation
# ==============================================================================

cat("================================================================================\n")
cat("  COMPARISON: Built-in VCF Reader vs StructuralVariantAnnotation\n")
cat("================================================================================\n\n")

if (file.exists(vcf_file)) {

    # Approach A: Built-in ShatterSeek VCF reader
    cat("Approach A: Built-in ShatterSeek VCF reader\n")
    sv_builtin <- read_sv_vcf(vcf_file, caller = "auto")
    cat(sprintf("  ✓ Detected %d SVs\n\n", length(sv_builtin@chrom1)))

    # Approach B: StructuralVariantAnnotation
    cat("Approach B: StructuralVariantAnnotation\n")
    sv_structvar <- read_sv_vcf_structuralvariant(vcf_file, genome = genome)
    cat(sprintf("  ✓ Detected %d SVs\n\n", length(sv_structvar@chrom1)))

    # Compare SV counts
    cat("Comparison:\n")
    cat(sprintf("  Built-in reader:  %d SVs\n", length(sv_builtin@chrom1)))
    cat(sprintf("  StructuralVariant: %d SVs\n", length(sv_structvar@chrom1)))

    if (length(sv_builtin@chrom1) != length(sv_structvar@chrom1)) {
        cat("\n  ⚠ Different SV counts detected!\n")
        cat("    This may be due to:\n")
        cat("    - Different filtering criteria\n")
        cat("    - Different BND parsing strategies\n")
        cat("    - Unpaired breakends handling\n\n")
    } else {
        cat("  ✓ Same number of SVs detected\n\n")
    }

} else {
    cat("Note: Replace vcf_file path to run comparison\n\n")
}


# ==============================================================================
# DETAILED EXAMPLE: GRIDSS VCF (Complex BND Records)
# ==============================================================================

cat("================================================================================\n")
cat("  DETAILED EXAMPLE: GRIDSS VCF Processing\n")
cat("================================================================================\n\n")

cat("GRIDSS produces complex BND records with extensive metadata.\n")
cat("StructuralVariantAnnotation is particularly useful for GRIDSS.\n\n")

gridss_vcf <- "path/to/gridss.vcf.gz"

if (file.exists(gridss_vcf)) {

    cat("Step 1: Read GRIDSS VCF\n")
    vcf <- readVcf(gridss_vcf, genome)
    cat(sprintf("  ✓ Read %d GRIDSS variants\n\n", nrow(vcf)))

    cat("Step 2: Extract breakpoints\n")
    cat("  GRIDSS uses complex BND pairing - StructuralVariantAnnotation handles this\n")
    gr_bp <- breakpointRanges(vcf)
    cat(sprintf("  ✓ Paired %d breakpoints from %d BND records\n\n",
               length(gr_bp), nrow(vcf)))

    cat("Step 3: Quality filtering (optional)\n")
    cat("  You can filter GRanges before conversion to ShatterSeek\n")

    # Example: Filter by QUAL score
    if ("QUAL" %in% names(mcols(gr_bp))) {
        high_qual <- gr_bp[mcols(gr_bp)$QUAL > 100]
        cat(sprintf("  ✓ Filtered to %d high-quality breakpoints (QUAL > 100)\n\n",
                   length(high_qual)))
        gr_bp <- high_qual
    }

    cat("Step 4: Convert to ShatterSeek\n")
    sv_data <- granges_to_svs(gr_bp)
    cat("  ✓ Conversion complete\n\n")

    cat("Step 5: Analyze\n")
    results <- detect_chromoanagenesis(sv_data, cn_data, genome = genome)
    cat("  ✓ Analysis complete\n\n")

} else {
    cat("Note: This is a template for GRIDSS VCF processing\n")
    cat("      Replace gridss_vcf path with your actual file\n\n")
}


# ==============================================================================
# ADVANCED: Custom Filtering with GRanges
# ==============================================================================

cat("================================================================================\n")
cat("  ADVANCED: Custom Filtering with GRanges\n")
cat("================================================================================\n\n")

cat("One advantage of StructuralVariantAnnotation workflow is that you can\n")
cat("perform complex filtering on GRanges before converting to ShatterSeek.\n\n")

if (file.exists(vcf_file)) {

    vcf <- readVcf(vcf_file, genome)
    gr <- breakpointRanges(vcf)

    cat("Original breakpoints: ", length(gr), "\n")

    # Example filters:

    # 1. Filter by size (for intrachromosomal SVs)
    cat("\n1. Size filtering:\n")
    if ("svLen" %in% names(mcols(gr))) {
        gr_large <- gr[!is.na(mcols(gr)$svLen) & abs(mcols(gr)$svLen) > 10000]
        cat(sprintf("   Breakpoints > 10kb: %d\n", length(gr_large)))
    }

    # 2. Filter by chromosome (exclude chr Y, M)
    cat("\n2. Chromosome filtering:\n")
    valid_chroms <- paste0("chr", c(1:22, "X"))
    gr_autosomal <- gr[as.character(seqnames(gr)) %in% valid_chroms]
    cat(sprintf("   Autosomal + X: %d\n", length(gr_autosomal)))

    # 3. Filter by PASS status
    cat("\n3. PASS filtering:\n")
    if ("FILTER" %in% names(mcols(gr))) {
        gr_pass <- gr[mcols(gr)$FILTER == "PASS"]
        cat(sprintf("   PASS only: %d\n", length(gr_pass)))
    }

    # Convert filtered GRanges to ShatterSeek
    cat("\n4. Converting filtered data to ShatterSeek:\n")
    sv_data_filtered <- granges_to_svs(gr_autosomal)
    cat("   ✓ Conversion complete\n\n")

} else {
    cat("Note: Replace vcf_file path to run filtering examples\n\n")
}


# ==============================================================================
# RECOMMENDATIONS
# ==============================================================================

cat("================================================================================\n")
cat("  RECOMMENDATIONS\n")
cat("================================================================================\n\n")

cat("Use BUILT-IN VCF reader when:\n")
cat("  ✓ Standard SV callers (Manta, Delly, DRAGEN, LUMPY, Sniffles)\n")
cat("  ✓ You want the simplest workflow\n")
cat("  ✓ Speed is important\n\n")

cat("Use StructuralVariantAnnotation when:\n")
cat("  ✓ Complex VCF formats (GRIDSS, custom callers)\n")
cat("  ✓ Need additional SV QC and filtering\n")
cat("  ✓ Already using Bioconductor ecosystem\n")
cat("  ✓ Maximum robustness for BND parsing\n")
cat("  ✓ Unsupported callers\n\n")

cat("For more information, see:\n")
cat("  - VCF_READING_STRATEGIES.md\n")
cat("  - https://bioconductor.org/packages/StructuralVariantAnnotation/\n\n")

cat("================================================================================\n")
