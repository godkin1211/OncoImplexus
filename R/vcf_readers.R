#' Read structural variants from VCF file
#'
#' Reads SV data from VCF format and converts it to SVs object compatible
#' with ShatterSeek. Supports multiple SV callers including Manta, GRIDSS,
#' Delly, LUMPY, and Sniffles.
#'
#' @param vcf_file Path to VCF file (can be gzipped, .vcf or .vcf.gz)
#' @param caller SV caller name. One of: "auto", "manta", "dragen", "gridss", "delly",
#'   "lumpy", "sniffles". Default "auto" attempts automatic detection.
#' @param min_sv_size Minimum SV size to include (default: 50 bp)
#' @param include_tra Include translocations/BND events (default: TRUE)
#' @param sample_name Sample name to extract from multi-sample VCF (default: first sample)
#' @return An SVs object containing the structural variations
#' @note Chromosome names are automatically standardized to match ShatterSeek's
#'   expected format (1, 2, ..., 22, X). Any "chr" prefix will be removed.
#'   SVs/CNVs on unsupported chromosomes (Y, M/MT, contigs, decoys, etc.) will be
#'   automatically excluded.
#' @details
#' This function handles various VCF formats from different SV callers:
#'
#' - **Manta**: Uses SVTYPE (BND, DEL, DUP, INV) and standard fields
#' - **DRAGEN**: Illumina's accelerated pipeline, uses standard SVTYPE with BND/MATEID
#' - **GRIDSS**: Primarily BND records with MATEID linking
#' - **Delly**: Standard SVTYPE with CT (connection type) for strands
#' - **LUMPY**: Similar to Delly with STRANDS field
#' - **Sniffles**: Long-read caller with standard SVTYPE
#'
#' The function extracts:
#' - Chromosome and position for both breakpoints
#' - SV type (DEL, DUP, INV, TRA)
#' - Strand orientation for both breakpoints
#'
#' Strand encoding:
#' - DEL: +/- (deletion-like)
#' - DUP: -/+ (duplication-like)
#' - h2hINV: +/+ (head-to-head inversion)
#' - t2tINV: -/- (tail-to-tail inversion)
#' - TRA: determined from ALT field or INFO fields
#'
#' @examples
#' \dontrun{
#' # Read from Manta VCF
#' sv_data <- read_sv_vcf("sample.manta.vcf.gz", caller = "manta")
#'
#' # Auto-detect caller
#' sv_data <- read_sv_vcf("sample.sv.vcf.gz")
#'
#' # Filter small SVs and exclude translocations
#' sv_data <- read_sv_vcf("sample.vcf.gz",
#'                        min_sv_size = 10000,
#'                        include_tra = FALSE)
#' }
#' @param vcf_file Path to the SV VCF file
#' @param genome Reference genome (hg19 or hg38)
#' @param caller SV caller name (auto, manta, delly, lumpy, gridss, etc.)
#' @param min_sv_size Minimum size for intrachromosomal SVs
#' @param include_tra Whether to include interchromosomal translocations
#' @param sample_name Optional sample name
#' @export
read_sv_vcf <- function(vcf_file,
                        genome = "hg19",
                        caller = "auto",
                        min_sv_size = 50,
                        include_tra = TRUE,
                        sample_name = NULL) {

    # Check if file exists
    if (!file.exists(vcf_file)) {
        stop(sprintf("VCF file not found: %s", vcf_file))
    }

    # PRIORITY 1: Use StructuralVariantAnnotation if available (Most robust)
    if (requireNamespace("StructuralVariantAnnotation", quietly = TRUE) &&
        requireNamespace("VariantAnnotation", quietly = TRUE)) {
        
        message("Using StructuralVariantAnnotation for robust SV parsing...")
        
        # Use the logic from read_sv_vcf_structuralvariant (implemented in structural_variant_annotation_support.R)
        # We assume this function is available in the package namespace
        svs <- read_sv_vcf_structuralvariant(vcf_file, genome = genome, include_breakends = TRUE, caller = caller)
        
        # Apply filters (min_sv_size and include_tra) which are not native to read_sv_vcf_structuralvariant
        
        # Convert to data frame for filtering
        # Note: We need to access slots directly as there might not be a data.frame method
        chrom1 <- svs@chrom1
        pos1 <- svs@pos1
        chrom2 <- svs@chrom2
        pos2 <- svs@pos2
        svtype <- svs@SVtype
        strand1 <- svs@strand1
        strand2 <- svs@strand2
        sv_id <- svs@sv_id
        
        keep <- rep(TRUE, length(chrom1))
        
        # Filter by size
        if (!is.null(min_sv_size) && min_sv_size > 0) {
            is_intra <- chrom1 == chrom2
            sv_size <- abs(pos2 - pos1)
            # Keep if NOT intrachromosomal OR (is intra AND size >= min_sv_size)
            keep <- keep & (!is_intra | (is_intra & sv_size >= min_sv_size))
        }
        
        # Filter translocations
        if (!include_tra) {
            is_tra <- svtype == "TRA" | chrom1 != chrom2
            keep <- keep & (!is_tra)
        }
        
        # Apply filter if needed
        if (sum(!keep) > 0) {
            svs <- SVs(
                chrom1 = chrom1[keep],
                pos1 = pos1[keep],
                chrom2 = chrom2[keep],
                pos2 = pos2[keep],
                SVtype = svtype[keep],
                strand1 = strand1[keep],
                strand2 = strand2[keep],
                sv_id = sv_id[keep]
            )
            message(sprintf("Filtered %d SVs based on size/type criteria", sum(!keep)))
        }
        
        message(sprintf("Loaded %d SVs from %s", length(svs@chrom1), vcf_file))
        return(svs)
        
    } else {
        # Fallback to legacy methods
        message("Note: StructuralVariantAnnotation package not found. Using legacy parsers (less robust for complex events).")
        message("Recommend installing: BiocManager::install('StructuralVariantAnnotation')")
        
        if (requireNamespace("VariantAnnotation", quietly = TRUE)) {
            svs <- .read_sv_vcf_variantannotation(vcf_file, caller, min_sv_size,
                                                  include_tra, sample_name)
        } else if (requireNamespace("vcfR", quietly = TRUE)) {
            svs <- .read_sv_vcf_vcfr(vcf_file, caller, min_sv_size,
                                     include_tra, sample_name)
        } else {
            stop("Please install 'StructuralVariantAnnotation' (recommended), 'VariantAnnotation', or 'vcfR' package.")
        }
        
        return(svs)
    }
}


#' Read copy number variants from VCF file
#'
#' Reads CNV data from VCF format and converts it to CNVsegs object.
#' Supports various CNV callers including CNVkit, GATK, Control-FREEC,
#' and Canvas.
#'
#' @param vcf_file Path to VCF file (can be gzipped)
#' @param caller CNV caller name. One of: "auto", "dragen", "cnvkit", "gatk",
#'   "controlfreec", "canvas". Default "auto" attempts automatic detection.
#' @param cn_field INFO field containing copy number (default: auto-detect)
#' @param merge_adjacent Merge adjacent segments with same CN (default: TRUE)
#' @param sample_name Sample name to extract from multi-sample VCF (default: first sample)
#' @return A CNVsegs object containing copy number segments
#' @note Chromosome names are automatically standardized to match ShatterSeek's
#'   expected format (1, 2, ..., 22, X). Any "chr" prefix will be removed.
#'   SVs/CNVs on unsupported chromosomes (Y, M/MT, contigs, decoys, etc.) will be
#'   automatically excluded.
#' @details
#' This function handles CNV VCF formats from different callers:
#'
#' - **DRAGEN**: Uses CN field, SVTYPE=CNV with DEL/DUP ALT
#' - **CNVkit**: Uses CN field in INFO
#' - **GATK**: Uses CN field
#' - **Control-FREEC**: Uses CN or CNV field
#' - **Canvas**: Uses CN field
#'
#' The function extracts:
#' - Chromosome
#' - Start position
#' - End position (from INFO/END)
#' - Total copy number
#'
#' @examples
#' \dontrun{
#' # Read from CNVkit VCF
#' cnv_data <- read_cnv_vcf("sample.cnvkit.vcf.gz", caller = "cnvkit")
#'
#' # Auto-detect caller
#' cnv_data <- read_cnv_vcf("sample.cnv.vcf.gz")
#'
#' # Don't merge adjacent segments
#' cnv_data <- read_cnv_vcf("sample.vcf.gz", merge_adjacent = FALSE)
#' }
#' @export
read_cnv_vcf <- function(vcf_file,
                         caller = "auto",
                         cn_field = NULL,
                         merge_adjacent = TRUE,
                         sample_name = NULL) {

    # Check if file exists
    if (!file.exists(vcf_file)) {
        stop(sprintf("VCF file not found: %s", vcf_file))
    }

    # Load package
    # Note: Chromosome name standardization is done inside these functions
    # BEFORE creating CNVsegs object to avoid initialization errors
    if (requireNamespace("VariantAnnotation", quietly = TRUE)) {
        cnvs <- .read_cnv_vcf_variantannotation(vcf_file, caller, cn_field,
                                                merge_adjacent, sample_name)
    } else if (requireNamespace("vcfR", quietly = TRUE)) {
        cnvs <- .read_cnv_vcf_vcfr(vcf_file, caller, cn_field,
                                   merge_adjacent, sample_name)
    } else {
        stop("Please install either 'VariantAnnotation' (Bioconductor) or 'vcfR' package:\n",
             "  BiocManager::install('VariantAnnotation')\n",
             "  or\n",
             "  install.packages('vcfR')")
    }

    return(cnvs)
}


#' Internal function: Read SV VCF using VariantAnnotation
#' @keywords internal
.read_sv_vcf_variantannotation <- function(vcf_file, caller, min_sv_size,
                                           include_tra, sample_name) {

    # Read VCF
    vcf <- VariantAnnotation::readVcf(vcf_file)

    # Auto-detect caller if needed
    if (caller == "auto") {
        caller <- .detect_sv_caller(vcf)
        message(sprintf("Auto-detected SV caller: %s", caller))
    }

    # Extract variant info
    info <- VariantAnnotation::info(vcf)
    ranges <- GenomicRanges::granges(vcf)

    # Initialize vectors
    chrom1 <- as.character(GenomicRanges::seqnames(ranges))
    pos1 <- GenomicRanges::start(ranges)

    # Get SVTYPE
    svtype <- if ("SVTYPE" %in% colnames(info)) {
        as.character(info$SVTYPE)
    } else {
        rep(NA, length(ranges))
    }

    # Extract end position and chr2 based on caller
    if (caller == "manta") {
        result <- .parse_manta_sv(vcf, info, ranges)
    } else if (caller == "dragen") {
        # DRAGEN uses standard VCF format, use generic parser
        result <- .parse_generic_sv(vcf, info, ranges)
    } else if (caller == "gridss") {
        result <- .parse_gridss_sv(vcf, info, ranges)
    } else if (caller == "delly") {
        result <- .parse_delly_sv(vcf, info, ranges)
    } else if (caller == "lumpy") {
        result <- .parse_lumpy_sv(vcf, info, ranges)
    } else if (caller == "sniffles") {
        result <- .parse_sniffles_sv(vcf, info, ranges)
    } else {
        # Generic parser
        result <- .parse_generic_sv(vcf, info, ranges)
    }

    # Filter by size
    if (!is.null(min_sv_size) && min_sv_size > 0) {
        is_intra <- result$chrom1 == result$chrom2
        sv_size <- abs(result$pos2 - result$pos1)
        keep <- !is_intra | (is_intra & sv_size >= min_sv_size)

        result$chrom1 <- result$chrom1[keep]
        result$pos1 <- result$pos1[keep]
        result$chrom2 <- result$chrom2[keep]
        result$pos2 <- result$pos2[keep]
        result$svtype <- result$svtype[keep]
        result$strand1 <- result$strand1[keep]
        result$strand2 <- result$strand2[keep]
    }

    # Filter translocations if requested
    if (!include_tra) {
        is_tra <- result$svtype == "TRA" | result$chrom1 != result$chrom2
        keep <- !is_tra

        result$chrom1 <- result$chrom1[keep]
        result$pos1 <- result$pos1[keep]
        result$chrom2 <- result$chrom2[keep]
        result$pos2 <- result$pos2[keep]
        result$svtype <- result$svtype[keep]
        result$strand1 <- result$strand1[keep]
        result$strand2 <- result$strand2[keep]
    }

    # Standardize chromosome names BEFORE creating SVs object
    # Remove "chr" prefix to match ShatterSeek's expected format (1, 2, ..., 22, X)
    result$chrom1 <- gsub("^chr", "", result$chrom1, ignore.case = TRUE)
    result$chrom1 <- gsub("^Chr", "", result$chrom1)
    result$chrom2 <- gsub("^chr", "", result$chrom2, ignore.case = TRUE)
    result$chrom2 <- gsub("^Chr", "", result$chrom2)

    # Convert MT to M if present
    result$chrom1 <- gsub("^MT$", "M", result$chrom1)
    result$chrom2 <- gsub("^MT$", "M", result$chrom2)

    # Filter to only supported chromosomes (1-22, X)
    # ShatterSeek only supports autosomes and X chromosome
    supported_chroms <- c(as.character(1:22), "X")
    keep_chrom <- result$chrom1 %in% supported_chroms & result$chrom2 %in% supported_chroms

    if (sum(!keep_chrom) > 0) {
        excluded <- sum(!keep_chrom)
        message(sprintf("Excluded %d SVs on unsupported chromosomes (Y, M, contigs, etc.)", excluded))
    }

    result$chrom1 <- result$chrom1[keep_chrom]
    result$pos1 <- result$pos1[keep_chrom]
    result$chrom2 <- result$chrom2[keep_chrom]
    result$pos2 <- result$pos2[keep_chrom]
    result$svtype <- result$svtype[keep_chrom]
    result$strand1 <- result$strand1[keep_chrom]
    result$strand2 <- result$strand2[keep_chrom]
    if (!is.null(result$sv_id)) {
        result$sv_id <- result$sv_id[keep_chrom]
    }

    # Create SVs object
    svs <- SVs(
        chrom1 = result$chrom1,
        pos1 = result$pos1,
        chrom2 = result$chrom2,
        pos2 = result$pos2,
        SVtype = result$svtype,
        strand1 = result$strand1,
        strand2 = result$strand2,
        sv_id = if (!is.null(result$sv_id)) as.character(result$sv_id) else as.character(1:length(result$chrom1))
    )

    message(sprintf("Loaded %d SVs from %s", length(svs@chrom1), vcf_file))

    return(svs)
}


#' Internal function: Read SV VCF using vcfR
#' @keywords internal
.read_sv_vcf_vcfr <- function(vcf_file, caller, min_sv_size,
                              include_tra, sample_name) {

    # Read VCF
    vcf <- vcfR::read.vcfR(vcf_file, verbose = FALSE)

    # Auto-detect caller if needed
    if (caller == "auto") {
        caller <- .detect_sv_caller_vcfr(vcf)
        message(sprintf("Auto-detected SV caller: %s", caller))
    }

    # Extract fixed fields
    fix <- vcfR::getFIX(vcf)
    chrom1 <- fix[, "CHROM"]
    pos1 <- as.numeric(fix[, "POS"])

    # Extract INFO fields
    info_df <- vcfR::extract.info(vcf)

    # Parse based on caller
    if (caller == "manta") {
        result <- .parse_manta_sv_vcfr(vcf, fix, info_df)
    } else if (caller == "dragen") {
        # DRAGEN uses standard VCF format, use generic parser
        result <- .parse_generic_sv_vcfr(vcf, fix, info_df)
    } else if (caller == "gridss") {
        result <- .parse_gridss_sv_vcfr(vcf, fix, info_df)
    } else if (caller == "delly") {
        result <- .parse_delly_sv_vcfr(vcf, fix, info_df)
    } else if (caller == "lumpy") {
        result <- .parse_lumpy_sv_vcfr(vcf, fix, info_df)
    } else if (caller == "sniffles") {
        result <- .parse_sniffles_sv_vcfr(vcf, fix, info_df)
    } else {
        result <- .parse_generic_sv_vcfr(vcf, fix, info_df)
    }

    # Apply filters (similar to VariantAnnotation version)
    if (!is.null(min_sv_size) && min_sv_size > 0) {
        is_intra <- result$chrom1 == result$chrom2
        sv_size <- abs(result$pos2 - result$pos1)
        keep <- !is_intra | (is_intra & sv_size >= min_sv_size)

        for (field in names(result)) {
            result[[field]] <- result[[field]][keep]
        }
    }

    if (!include_tra) {
        is_tra <- result$svtype == "TRA" | result$chrom1 != result$chrom2
        keep <- !is_tra

        for (field in names(result)) {
            result[[field]] <- result[[field]][keep]
        }
    }

    # Standardize chromosome names BEFORE creating SVs object
    # Remove "chr" prefix to match ShatterSeek's expected format (1, 2, ..., 22, X)
    result$chrom1 <- gsub("^chr", "", result$chrom1, ignore.case = TRUE)
    result$chrom1 <- gsub("^Chr", "", result$chrom1)
    result$chrom2 <- gsub("^chr", "", result$chrom2, ignore.case = TRUE)
    result$chrom2 <- gsub("^Chr", "", result$chrom2)

    # Convert MT to M if present
    result$chrom1 <- gsub("^MT$", "M", result$chrom1)
    result$chrom2 <- gsub("^MT$", "M", result$chrom2)

    # Filter to only supported chromosomes (1-22, X)
    # ShatterSeek only supports autosomes and X chromosome
    supported_chroms <- c(as.character(1:22), "X")
    keep_chrom <- result$chrom1 %in% supported_chroms & result$chrom2 %in% supported_chroms

    if (sum(!keep_chrom) > 0) {
        excluded <- sum(!keep_chrom)
        message(sprintf("Excluded %d SVs on unsupported chromosomes (Y, M, contigs, etc.)", excluded))
    }

    for (field in names(result)) {
        result[[field]] <- result[[field]][keep_chrom]
    }

    # Create SVs object
    svs <- SVs(
        chrom1 = result$chrom1,
        pos1 = result$pos1,
        chrom2 = result$chrom2,
        pos2 = result$pos2,
        SVtype = result$svtype,
        strand1 = result$strand1,
        strand2 = result$strand2,
        sv_id = if (!is.null(result$sv_id)) as.character(result$sv_id) else as.character(1:length(result$chrom1))
    )

    message(sprintf("Loaded %d SVs from %s", length(svs@chrom1), vcf_file))

    return(svs)
}


#' Internal function: Read CNV VCF using VariantAnnotation
#' @keywords internal
.read_cnv_vcf_variantannotation <- function(vcf_file, caller, cn_field,
                                            merge_adjacent, sample_name) {

    # Read VCF
    vcf <- VariantAnnotation::readVcf(vcf_file)

    # Auto-detect caller if needed
    if (caller == "auto") {
        caller <- .detect_cnv_caller(vcf)
        message(sprintf("Auto-detected CNV caller: %s", caller))
    }

    # Auto-detect CN field if not specified
    if (is.null(cn_field)) {
        cn_field <- .detect_cn_field(vcf, caller)
    }

    # Extract variant info
    info <- VariantAnnotation::info(vcf)
    ranges <- GenomicRanges::granges(vcf)

    chrom <- as.character(GenomicRanges::seqnames(ranges))
    start <- GenomicRanges::start(ranges)

    # Get END position
    end <- if ("END" %in% colnames(info)) {
        as.numeric(info$END)
    } else {
        start  # Use start if END not available
    }

    # Get copy number - check if it's in FORMAT or INFO
    if (grepl("^FORMAT:", cn_field)) {
        # Extract from FORMAT field (sample-specific)
        field_name <- sub("^FORMAT:", "", cn_field)
        geno <- VariantAnnotation::geno(vcf)

        if (!field_name %in% names(geno)) {
            stop(sprintf("FORMAT field '%s' not found in VCF", field_name))
        }

        # Get CN values - if multi-sample, use first sample or specified sample
        cn_matrix <- geno[[field_name]]
        if (is.matrix(cn_matrix)) {
            if (!is.null(sample_name)) {
                if (!sample_name %in% colnames(cn_matrix)) {
                    stop(sprintf("Sample '%s' not found in VCF", sample_name))
                }
                cn <- as.numeric(cn_matrix[, sample_name])
            } else {
                # Use first sample
                cn <- as.numeric(cn_matrix[, 1])
            }
        } else {
            cn <- as.numeric(cn_matrix)
        }
    } else {
        # Extract from INFO field
        if (!cn_field %in% colnames(info)) {
            stop(sprintf("INFO field '%s' not found in VCF", cn_field))
        }
        cn <- as.numeric(info[[cn_field]])
    }

    # Remove NA values
    valid <- !is.na(chrom) & !is.na(start) & !is.na(end) & !is.na(cn)
    chrom <- chrom[valid]
    start <- start[valid]
    end <- end[valid]
    cn <- cn[valid]

    # Create data frame for merging
    df <- data.frame(
        chrom = chrom,
        start = start,
        end = end,
        CN = cn,
        stringsAsFactors = FALSE
    )

    # Sort by chromosome and position
    df <- df[order(df$chrom, df$start), ]

    # Merge adjacent segments if requested
    if (merge_adjacent) {
        df <- .merge_adjacent_segments(df)
    }

    # Standardize chromosome names BEFORE creating CNVsegs object
    # Remove "chr" prefix to match ShatterSeek's expected format (1, 2, ..., 22, X)
    df$chrom <- gsub("^chr", "", df$chrom, ignore.case = TRUE)
    df$chrom <- gsub("^Chr", "", df$chrom)

    # Convert MT to M if present
    df$chrom <- gsub("^MT$", "M", df$chrom)

    # Filter to only supported chromosomes (1-22, X)
    # ShatterSeek only supports autosomes and X chromosome
    supported_chroms <- c(as.character(1:22), "X")
    keep_chrom <- df$chrom %in% supported_chroms

    if (sum(!keep_chrom) > 0) {
        excluded <- sum(!keep_chrom)
        message(sprintf("Excluded %d CNV segments on unsupported chromosomes (Y, M, contigs, etc.)", excluded))
    }

    df <- df[keep_chrom, ]

    # Create CNVsegs object
    cnvs <- CNVsegs(
        chrom = df$chrom,
        start = df$start,
        end = df$end,
        total_cn = df$CN
    )

    message(sprintf("Loaded %d CNV segments from %s", length(cnvs@chrom), vcf_file))

    return(cnvs)
}


#' Internal function: Read CNV VCF using vcfR
#' @keywords internal
.read_cnv_vcf_vcfr <- function(vcf_file, caller, cn_field,
                               merge_adjacent, sample_name) {

    # Read VCF
    vcf <- vcfR::read.vcfR(vcf_file, verbose = FALSE)

    # Auto-detect caller if needed
    if (caller == "auto") {
        caller <- .detect_cnv_caller_vcfr(vcf)
        message(sprintf("Auto-detected CNV caller: %s", caller))
    }

    # Auto-detect CN field if not specified
    if (is.null(cn_field)) {
        cn_field <- .detect_cn_field_vcfr(vcf, caller)
    }

    # Extract fixed fields
    fix <- vcfR::getFIX(vcf)
    chrom <- fix[, "CHROM"]
    start <- as.numeric(fix[, "POS"])

    # Extract INFO (always need END)
    info_df <- vcfR::extract.info(vcf, element = "END")

    end <- if ("END" %in% colnames(info_df)) {
        as.numeric(info_df[, "END"])
    } else {
        start
    }

    # Get copy number - check if it's in FORMAT or INFO
    if (grepl("^FORMAT:", cn_field)) {
        # Extract from FORMAT field (sample-specific)
        field_name <- sub("^FORMAT:", "", cn_field)

        # Use extract.gt to get FORMAT field
        cn_matrix <- vcfR::extract.gt(vcf, element = field_name, as.numeric = TRUE)

        if (is.null(cn_matrix)) {
            stop(sprintf("FORMAT field '%s' not found in VCF", field_name))
        }

        # Get CN values - if multi-sample, use first sample or specified sample
        if (is.matrix(cn_matrix)) {
            if (!is.null(sample_name)) {
                if (!sample_name %in% colnames(cn_matrix)) {
                    stop(sprintf("Sample '%s' not found in VCF", sample_name))
                }
                cn <- as.numeric(cn_matrix[, sample_name])
            } else {
                # Use first sample
                cn <- as.numeric(cn_matrix[, 1])
            }
        } else {
            cn <- as.numeric(cn_matrix)
        }
    } else {
        # Extract from INFO field
        info_cn <- vcfR::extract.info(vcf, element = cn_field)

        if (is.null(info_cn) || !cn_field %in% colnames(info_cn)) {
            stop(sprintf("INFO field '%s' not found in VCF", cn_field))
        }
        cn <- as.numeric(info_cn[, cn_field])
    }

    # Remove NA values
    valid <- !is.na(chrom) & !is.na(start) & !is.na(end) & !is.na(cn)
    chrom <- chrom[valid]
    start <- start[valid]
    end <- end[valid]
    cn <- cn[valid]

    # Create data frame
    df <- data.frame(
        chrom = chrom,
        start = start,
        end = end,
        CN = cn,
        stringsAsFactors = FALSE
    )

    # Sort by chromosome and position
    df <- df[order(df$chrom, df$start), ]

    # Merge adjacent segments if requested
    if (merge_adjacent) {
        df <- .merge_adjacent_segments(df)
    }

    # Standardize chromosome names BEFORE creating CNVsegs object
    # Remove "chr" prefix to match ShatterSeek's expected format (1, 2, ..., 22, X)
    df$chrom <- gsub("^chr", "", df$chrom, ignore.case = TRUE)
    df$chrom <- gsub("^Chr", "", df$chrom)

    # Convert MT to M if present
    df$chrom <- gsub("^MT$", "M", df$chrom)

    # Filter to only supported chromosomes (1-22, X)
    # ShatterSeek only supports autosomes and X chromosome
    supported_chroms <- c(as.character(1:22), "X")
    keep_chrom <- df$chrom %in% supported_chroms

    if (sum(!keep_chrom) > 0) {
        excluded <- sum(!keep_chrom)
        message(sprintf("Excluded %d CNV segments on unsupported chromosomes (Y, M, contigs, etc.)", excluded))
    }

    df <- df[keep_chrom, ]

    # Create CNVsegs object
    cnvs <- CNVsegs(
        chrom = df$chrom,
        start = df$start,
        end = df$end,
        total_cn = df$CN
    )

    message(sprintf("Loaded %d CNV segments from %s", length(cnvs@chrom), vcf_file))

    return(cnvs)
}


#' Internal: Merge adjacent CNV segments with same CN
#' @keywords internal
.merge_adjacent_segments <- function(df) {
    if (nrow(df) == 0) return(df)

    merged <- list()
    current <- df[1, ]

    for (i in 2:nrow(df)) {
        row <- df[i, ]

        # Check if same chromosome, adjacent position, and same CN
        if (row$chrom == current$chrom &&
            row$start <= current$end + 1 &&
            row$CN == current$CN) {
            # Merge: extend current segment
            current$end <- max(current$end, row$end)
        } else {
            # Save current and start new
            merged[[length(merged) + 1]] <- current
            current <- row
        }
    }

    # Add last segment
    merged[[length(merged) + 1]] <- current

    # Convert back to data frame
    do.call(rbind, merged)
}


#' Internal: Add "chr" prefix to chromosome names for SVs object
#' @keywords internal
.add_chr_prefix_svs <- function(svs) {
    # Convert SVs object to data frame
    svs_df <- as.data.frame(svs)

    # Add "chr" prefix to numeric chromosomes and X, Y
    # Don't add to chromosomes that already have it
    svs_df$chrom1 <- ifelse(
        !grepl("^chr", svs_df$chrom1, ignore.case = TRUE),
        paste0("chr", svs_df$chrom1),
        svs_df$chrom1
    )

    svs_df$chrom2 <- ifelse(
        !grepl("^chr", svs_df$chrom2, ignore.case = TRUE),
        paste0("chr", svs_df$chrom2),
        svs_df$chrom2
    )

    # Create new SVs object - PRESERVE SV_ID
    svs_new <- SVs(
        chrom1 = svs_df$chrom1,
        pos1 = svs_df$pos1,
        chrom2 = svs_df$chrom2,
        pos2 = svs_df$pos2,
        SVtype = svs_df$SVtype,
        strand1 = svs_df$strand1,
        strand2 = svs_df$strand2,
        sv_id = as.character(svs_df$sv_id)
    )

    return(svs_new)
}


#' Internal: Add "chr" prefix to chromosome names for CNVsegs object
#' @keywords internal
.add_chr_prefix_cnvs <- function(cnvs) {
    # Convert CNVsegs object to data frame
    cnvs_df <- as.data.frame(cnvs)

    # Add "chr" prefix to numeric chromosomes and X, Y
    # Don't add to chromosomes that already have it
    cnvs_df$chrom <- ifelse(
        !grepl("^chr", cnvs_df$chrom, ignore.case = TRUE),
        paste0("chr", cnvs_df$chrom),
        cnvs_df$chrom
    )

    # Create new CNVsegs object
    cnvs_new <- CNVsegs(
        chrom = cnvs_df$chrom,
        start = cnvs_df$start,
        end = cnvs_df$end,
        total_cn = cnvs_df$total_cn
    )

    return(cnvs_new)
}


# Caller-specific parsers are implemented in vcf_parsers.R
# All R files in R/ directory are automatically loaded by the package
