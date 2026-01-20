#' VCF Parser Functions for Different SV Callers
#' Internal functions for parsing VCF files from various SV and CNV callers
#' @keywords internal


#' Detect SV caller from VCF (VariantAnnotation)
#' @keywords internal
.detect_sv_caller <- function(vcf) {
    header <- VariantAnnotation::header(vcf)
    meta <- VariantAnnotation::meta(header)

    # Check source field in header
    if ("source" %in% names(meta)) {
        source_obj <- meta$source
        # Handle different meta object types
        if (is(source_obj, "DataFrame") || is(source_obj, "data.frame")) {
            source <- tolower(paste(source_obj[, 1], collapse = " "))
        } else if (is(source_obj, "character")) {
            source <- tolower(paste(source_obj, collapse = " "))
        } else {
            source <- tolower(paste(as.character(source_obj), collapse = " "))
        }

        # Check for DRAGEN first (before GRIDSS, as they use similar format)
        if (grepl("dragen", source)) return("dragen")
        if (grepl("manta", source)) return("manta")
        if (grepl("gridss", source)) return("gridss")
        if (grepl("delly", source)) return("delly")
        if (grepl("lumpy", source)) return("lumpy")
        if (grepl("sniffles", source)) return("sniffles")
    }

    # Check INFO fields
    info <- VariantAnnotation::info(vcf)
    info_names <- colnames(info)

    if ("MATEID" %in% info_names && "SVTYPE" %in% info_names) {
        return("gridss")
    }
    if ("CT" %in% info_names) {
        return("delly")
    }
    if ("STRANDS" %in% info_names) {
        return("lumpy")
    }
    if ("RNAMES" %in% info_names) {
        return("sniffles")
    }

    # Default to generic
    return("generic")
}


#' Detect SV caller from VCF (vcfR)
#' @keywords internal
.detect_sv_caller_vcfr <- function(vcf) {
    # Check meta lines
    meta <- vcf@meta
    source_lines <- grep("^##source=", meta, value = TRUE, ignore.case = TRUE)

    if (length(source_lines) > 0) {
        source <- tolower(source_lines[1])
        # Check for DRAGEN first (before GRIDSS, as they use similar format)
        if (grepl("dragen", source)) return("dragen")
        if (grepl("manta", source)) return("manta")
        if (grepl("gridss", source)) return("gridss")
        if (grepl("delly", source)) return("delly")
        if (grepl("lumpy", source)) return("lumpy")
        if (grepl("sniffles", source)) return("sniffles")
    }

    # Check INFO field definitions
    info_lines <- grep("^##INFO=", meta, value = TRUE)

    if (any(grepl("MATEID", info_lines)) && any(grepl("SVTYPE", info_lines))) {
        return("gridss")
    }
    if (any(grepl("\\bCT=", info_lines))) {
        return("delly")
    }
    if (any(grepl("STRANDS", info_lines))) {
        return("lumpy")
    }
    if (any(grepl("RNAMES", info_lines))) {
        return("sniffles")
    }

    return("generic")
}


#' Detect CNV caller from VCF (VariantAnnotation)
#' @keywords internal
.detect_cnv_caller <- function(vcf) {
    header <- VariantAnnotation::header(vcf)
    meta <- VariantAnnotation::meta(header)

    if ("source" %in% names(meta)) {
        source_obj <- meta$source
        # Handle different meta object types
        if (is(source_obj, "DataFrame") || is(source_obj, "data.frame")) {
            source <- tolower(paste(source_obj[, 1], collapse = " "))
        } else if (is(source_obj, "character")) {
            source <- tolower(paste(source_obj, collapse = " "))
        } else {
            source <- tolower(paste(as.character(source_obj), collapse = " "))
        }

        if (grepl("dragen", source)) return("dragen")
        if (grepl("cnvkit", source)) return("cnvkit")
        if (grepl("gatk", source)) return("gatk")
        if (grepl("freec|control-freec", source)) return("controlfreec")
        if (grepl("canvas", source)) return("canvas")
    }

    # Check INFO fields
    info <- VariantAnnotation::info(vcf)
    if ("CN" %in% colnames(info)) return("generic")

    return("generic")
}


#' Detect CNV caller from VCF (vcfR)
#' @keywords internal
.detect_cnv_caller_vcfr <- function(vcf) {
    meta <- vcf@meta
    source_lines <- grep("^##source=", meta, value = TRUE, ignore.case = TRUE)

    if (length(source_lines) > 0) {
        source <- tolower(source_lines[1])
        if (grepl("dragen", source)) return("dragen")
        if (grepl("cnvkit", source)) return("cnvkit")
        if (grepl("gatk", source)) return("gatk")
        if (grepl("freec|control-freec", source)) return("controlfreec")
        if (grepl("canvas", source)) return("canvas")
    }

    return("generic")
}


#' Detect CN field in VCF
#' @keywords internal
.detect_cn_field <- function(vcf, caller) {
    info <- VariantAnnotation::info(vcf)
    info_names <- colnames(info)

    # Check FORMAT fields for DRAGEN and other callers that use sample-specific CN
    geno <- VariantAnnotation::geno(vcf)
    format_names <- names(geno)

    # DRAGEN uses FORMAT/CN instead of INFO/CN
    if (caller == "dragen" && "CN" %in% format_names) return("FORMAT:CN")

    # Try caller-specific INFO fields first
    if (caller == "cnvkit" && "CN" %in% info_names) return("CN")
    if (caller == "gatk" && "CN" %in% info_names) return("CN")
    if (caller == "controlfreec") {
        if ("CNV" %in% info_names) return("CNV")
        if ("CN" %in% info_names) return("CN")
    }
    if (caller == "canvas" && "CN" %in% info_names) return("CN")

    # Generic detection - check INFO first, then FORMAT
    if ("CN" %in% info_names) return("CN")
    if ("CNV" %in% info_names) return("CNV")
    if ("TOTAL_CN" %in% info_names) return("TOTAL_CN")
    if ("TCN" %in% info_names) return("TCN")
    if ("CN_TOTAL" %in% info_names) return("CN_TOTAL")

    # Check FORMAT fields as fallback
    if ("CN" %in% format_names) return("FORMAT:CN")
    if ("CNV" %in% format_names) return("FORMAT:CNV")
    if ("TCN" %in% format_names) return("FORMAT:TCN")
    if ("NG" %in% format_names) return("FORMAT:NG") # Sometimes used by GATK/others

    stop("Could not detect copy number field in VCF. Please specify 'cn_field' parameter.")
}


#' Detect CN field in VCF (vcfR)
#' @keywords internal
.detect_cn_field_vcfr <- function(vcf, caller) {
    # Extract INFO field names from meta
    info_lines <- grep("^##INFO=<ID=", vcf@meta, value = TRUE)
    info_names <- gsub("^##INFO=<ID=([^,]+),.*", "\\1", info_lines)

    # Extract FORMAT field names from meta
    format_lines <- grep("^##FORMAT=<ID=", vcf@meta, value = TRUE)
    format_names <- gsub("^##FORMAT=<ID=([^,]+),.*", "\\1", format_lines)

    # DRAGEN uses FORMAT/CN instead of INFO/CN
    if (caller == "dragen" && "CN" %in% format_names) return("FORMAT:CN")

    # Try caller-specific INFO fields first
    if (caller == "cnvkit" && "CN" %in% info_names) return("CN")
    if (caller == "gatk" && "CN" %in% info_names) return("CN")
    if (caller == "controlfreec") {
        if ("CNV" %in% info_names) return("CNV")
        if ("CN" %in% info_names) return("CN")
    }
    if (caller == "canvas" && "CN" %in% info_names) return("CN")

    # Generic detection - check INFO first, then FORMAT
    if ("CN" %in% info_names) return("CN")
    if ("CNV" %in% info_names) return("CNV")
    if ("TOTAL_CN" %in% info_names) return("TOTAL_CN")
    if ("TCN" %in% info_names) return("TCN")
    if ("CN_TOTAL" %in% info_names) return("CN_TOTAL")

    # Check FORMAT fields as fallback
    if ("CN" %in% format_names) return("FORMAT:CN")
    if ("CNV" %in% format_names) return("FORMAT:CNV")
    if ("TCN" %in% format_names) return("FORMAT:TCN")
    if ("NG" %in% format_names) return("FORMAT:NG") # Sometimes used by GATK/others

    stop("Could not detect copy number field in VCF. Please specify 'cn_field' parameter.")
}


# =============================================================================
# Manta SV Parser
# ============================================================================ 

#' Parse Manta SV VCF (VariantAnnotation)
#' @keywords internal
.parse_manta_sv <- function(vcf, info, ranges) {
    chrom1 <- as.character(GenomicRanges::seqnames(ranges))
    pos1 <- GenomicRanges::start(ranges)

    # Get SVTYPE
    svtype <- as.character(info$SVTYPE)

    # Get CHR2 and END
    chrom2 <- if ("CHR2" %in% colnames(info)) {
        as.character(info$CHR2)
    } else {
        chrom1
    }

    pos2 <- if ("END" %in% colnames(info)) {
        as.numeric(info$END)
    } else {
        pos1
    }

    # Determine strands based on SVTYPE
    strand1 <- rep("+", length(svtype))
    strand2 <- rep("-", length(svtype))

    # DEL: +/-
    strand1[svtype == "DEL"] <- "+"
    strand2[svtype == "DEL"] <- "-"

    # DUP: -/+
    strand1[svtype == "DUP"] <- "-"
    strand2[svtype == "DUP"] <- "+"

    # INV: check if h2h or t2t
    if ("INV3" %in% colnames(info) || "INV5" %in% colnames(info)) {
        # INV3 is h2h (+/+), INV5 is t2t (-/-)
        is_inv3 <- !is.na(info$INV3)
        is_inv5 <- !is.na(info$INV5)

        svtype[is_inv3] <- "h2hINV"
        svtype[is_inv5] <- "t2tINV"

        strand1[is_inv3] <- "+"
        strand2[is_inv3] <- "+"

        strand1[is_inv5] <- "-"
        strand2[is_inv5] <- "-"
    } else {
        # Generic INV - try to determine from other fields
        is_inv <- svtype == "INV"
        svtype[is_inv] <- "h2hINV"  # Default to h2h
        strand1[is_inv] <- "+"
        strand2[is_inv] <- "+"
    }

    # BND (translocations)
    is_bnd <- svtype == "BND"
    svtype[is_bnd] <- "TRA"

    # For BND, parse ALT field for strand info
    if (any(is_bnd)) {
        alt <- as.character(VariantAnnotation::alt(vcf))
        bnd_idx <- which(is_bnd)

        for (i in bnd_idx) {
            alt_str <- alt[i]
            # Parse bracket notation: ]chr:pos] or [chr:pos[
            if (grepl("]", alt_str)) {
                # ] means forward strand
                if (grepl("^]", alt_str)) {
                    strand1[i] <- "-"
                } else {
                    strand1[i] <- "+"
                }
                if (grepl("]$", alt_str)) {
                    strand2[i] <- "-"
                } else {
                    strand2[i] <- "+"
                }
            } else if (grepl("[", alt_str)) {
                # [ means reverse strand
                if (grepl("^[ ", alt_str)) {
                    strand1[i] <- "+"
                } else {
                    strand1[i] <- "-"
                }
                if (grepl("[$]", alt_str)) {
                    strand2[i] <- "+"
                } else {
                    strand2[i] <- "-"
                }
            }
        }
    }

    # Convert SVtype to ShatterSeek format
    svtype[svtype == "INS"] <- "DUP"  # Treat insertions as duplications

    list(
        chrom1 = chrom1,
        pos1 = pos1,
        chrom2 = chrom2,
        pos2 = pos2,
        svtype = svtype,
        strand1 = strand1,
        strand2 = strand2,
        sv_id = names(ranges)
    )
}


#' Parse Manta SV VCF (vcfR)
#' @keywords internal
.parse_manta_sv_vcfr <- function(vcf, fix, info_df) {
    chrom1 <- fix[, "CHROM"]
    pos1 <- as.numeric(fix[, "POS"])
    alt <- fix[, "ALT"]

    # Extract SVTYPE
    svtype <- if ("SVTYPE" %in% colnames(info_df)) {
        as.character(info_df[, "SVTYPE"])
    } else {
        rep(NA, nrow(fix))
    }

    # Extract CHR2 and END
    chrom2 <- if ("CHR2" %in% colnames(info_df)) {
        as.character(info_df[, "CHR2"])
    } else {
        chrom1
    }

    pos2 <- if ("END" %in% colnames(info_df)) {
        as.numeric(info_df[, "END"])
    } else {
        pos1
    }

    # Determine strands (similar logic to VariantAnnotation version)
    strand1 <- rep("+", length(svtype))
    strand2 <- rep("-", length(svtype))

    strand1[svtype == "DEL"] <- "+"
    strand2[svtype == "DEL"] <- "-"

    strand1[svtype == "DUP"] <- "-"
    strand2[svtype == "DUP"] <- "+"

    # Handle inversions
    if ("INV3" %in% colnames(info_df)) {
        is_inv3 <- !is.na(info_df[, "INV3"])
        svtype[is_inv3] <- "h2hINV"
        strand1[is_inv3] <- "+"
        strand2[is_inv3] <- "+"
    }

    if ("INV5" %in% colnames(info_df)) {
        is_inv5 <- !is.na(info_df[, "INV5"])
        svtype[is_inv5] <- "t2tINV"
        strand1[is_inv5] <- "-"
        strand2[is_inv5] <- "-"
    }

    # Handle BND
    is_bnd <- svtype == "BND"
    svtype[is_bnd] <- "TRA"

    list(
        chrom1 = chrom1,
        pos1 = pos1,
        chrom2 = chrom2,
        pos2 = pos2,
        svtype = svtype,
        strand1 = strand1,
        strand2 = strand2,
        sv_id = fix[, "ID"]
    )
}


# ============================================================================ 
# Delly SV Parser
# ============================================================================ 

#' Parse Delly SV VCF (VariantAnnotation)
#' @keywords internal
.parse_delly_sv <- function(vcf, info, ranges) {
    chrom1 <- as.character(GenomicRanges::seqnames(ranges))
    pos1 <- GenomicRanges::start(ranges)

    svtype <- as.character(info$SVTYPE)

    # Delly uses CHR2 for second chromosome
    chrom2 <- if ("CHR2" %in% colnames(info)) {
        as.character(info$CHR2)
    } else {
        chrom1
    }

    pos2 <- if ("END" %in% colnames(info)) {
        as.numeric(info$END)
    } else {
        pos1
    }

    # Delly uses CT field for connection type (strand info)
    # CT: 3to5 means +/-, 5to3 means -/+, 3to3 means +/+, 5to5 means -/-
    strand1 <- rep("+", length(svtype))
    strand2 <- rep("-", length(svtype))

    if ("CT" %in% colnames(info)) {
        ct <- as.character(info$CT)

        strand1[ct == "3to5"] <- "+"
        strand2[ct == "3to5"] <- "-"

        strand1[ct == "5to3"] <- "-"
        strand2[ct == "5to3"] <- "+"

        strand1[ct == "3to3"] <- "+"
        strand2[ct == "3to3"] <- "+"

        strand1[ct == "5to5"] <- "-"
        strand2[ct == "5to5"] <- "-"
    } else {
        # Fallback to SVTYPE-based inference
        strand1[svtype == "DEL"] <- "+"
        strand2[svtype == "DEL"] <- "-"

        strand1[svtype == "DUP"] <- "-"
        strand2[svtype == "DUP"] <- "+"

        strand1[svtype == "INV"] <- "+"
        strand2[svtype == "INV"] <- "+"
    }

    # Convert SVTYPE
    svtype[svtype == "TRA"] <- "TRA"
    svtype[svtype == "INV"] <- "h2hINV"  # Default to h2h

    list(
        chrom1 = chrom1,
        pos1 = pos1,
        chrom2 = chrom2,
        pos2 = pos2,
        svtype = svtype,
        strand1 = strand1,
        strand2 = strand2,
        sv_id = names(ranges)
    )
}


#' Parse Delly SV VCF (vcfR)
#' @keywords internal
.parse_delly_sv_vcfr <- function(vcf, fix, info_df) {
    # Similar implementation to VariantAnnotation version
    chrom1 <- fix[, "CHROM"]
    pos1 <- as.numeric(fix[, "POS"])

    svtype <- if ("SVTYPE" %in% colnames(info_df)) {
        as.character(info_df[, "SVTYPE"])
    } else {
        rep(NA, nrow(fix))
    }

    chrom2 <- if ("CHR2" %in% colnames(info_df)) {
        as.character(info_df[, "CHR2"])
    } else {
        chrom1
    }

    pos2 <- if ("END" %in% colnames(info_df)) {
        as.numeric(info_df[, "END"])
    } else {
        pos1
    }

    strand1 <- rep("+", length(svtype))
    strand2 <- rep("-", length(svtype))

    if ("CT" %in% colnames(info_df)) {
        ct <- as.character(info_df[, "CT"])

        strand1[ct == "3to5"] <- "+"
        strand2[ct == "3to5"] <- "-"

        strand1[ct == "5to3"] <- "-"
        strand2[ct == "5to3"] <- "+"

        strand1[ct == "3to3"] <- "+"
        strand2[ct == "3to3"] <- "+"

        strand1[ct == "5to5"] <- "-"
        strand2[ct == "5to5"] <- "-"
    }

    svtype[svtype == "TRA"] <- "TRA"
    svtype[svtype == "INV"] <- "h2hINV"

    list(
        chrom1 = chrom1,
        pos1 = pos1,
        chrom2 = chrom2,
        pos2 = pos2,
        svtype = svtype,
        strand1 = strand1,
        strand2 = strand2,
        sv_id = fix[, "ID"]
    )
}


# ============================================================================ 
# Generic SV Parser
# ============================================================================ 

#' Parse generic SV VCF (VariantAnnotation)
#' @keywords internal
.parse_generic_sv <- function(vcf, info, ranges) {
    chrom1 <- as.character(GenomicRanges::seqnames(ranges))
    pos1 <- GenomicRanges::start(ranges)

    svtype <- if ("SVTYPE" %in% colnames(info)) {
        as.character(info$SVTYPE)
    } else {
        rep("UNK", length(ranges))
    }

    # Try multiple field names for second chromosome
    chrom2 <- chrom1
    for (field in c("CHR2", "CHROM2", "chr2")) {
        if (field %in% colnames(info)) {
            chrom2 <- as.character(info[[field]])
            break
        }
    }

    # Try multiple field names for end position
    pos2 <- pos1
    for (field in c("END", "END2", "POS2")) {
        if (field %in% colnames(info)) {
            pos2 <- as.numeric(info[[field]])
            break
        }
    }

    # Parse BND ALT field if SVTYPE=BND and chrom2/pos2 not set
    is_bnd <- svtype == "BND"
    if (any(is_bnd)) {
        alt <- as.character(VariantAnnotation::alt(vcf))

        for (i in which(is_bnd)) {
            # Only parse ALT if we don't have chrom2/pos2 from INFO
            # Check for NA values and default values
            needs_parsing <- FALSE
            if (is.na(chrom2[i]) || is.na(pos2[i])) {
                needs_parsing <- TRUE
            } else if (as.character(chrom2[i]) == as.character(chrom1[i]) && pos2[i] == pos1[i]) {
                needs_parsing <- TRUE
            }

            if (needs_parsing) {
                alt_str <- alt[i]

                # Parse bracket notation: ]chr:pos] or [chr:pos[
                # Patterns: N[chr:pos[, [chr:pos[N, N]chr:pos], ]chr:pos]N
                pattern <- "([\\[\\]])([^:\\[\\]]+):([0-9]+)([\\[\\]])"
                match <- regexec(pattern, alt_str)
                if (match[[1]][1] != -1) {
                    matches <- regmatches(alt_str, match)[[1]]
                    if (length(matches) >= 4) {
                        mate_chr <- matches[3]
                        mate_pos <- as.numeric(matches[4])

                        # Remove "chr" prefix if present
                        mate_chr <- gsub("^chr", "", mate_chr, ignore.case = TRUE)

                        chrom2[i] <- mate_chr
                        pos2[i] <- mate_pos
                    }
                }
            }
        }
    }

    # Remove variants with NA coordinates or chromosomes (vital for robustness)
    valid_idx <- !is.na(chrom1) & !is.na(pos1) & !is.na(chrom2) & !is.na(pos2)
    if (sum(!valid_idx) > 0) {
        message(sprintf("Dropped %d SVs with NA coordinates/chromosomes", sum(!valid_idx)))
        chrom1 <- chrom1[valid_idx]
        pos1 <- pos1[valid_idx]
        chrom2 <- chrom2[valid_idx]
        pos2 <- pos2[valid_idx]
        svtype <- svtype[valid_idx]
    }

    # Infer strands from SVTYPE
    strand1 <- rep("+", length(svtype))
    strand2 <- rep("-", length(svtype))

    strand1[svtype == "DEL"] <- "+"
    strand2[svtype == "DEL"] <- "-"

    strand1[svtype == "DUP"] <- "-"
    strand2[svtype == "DUP"] <- "+"

    strand1[svtype == "INV"] <- "+"
    strand2[svtype == "INV"] <- "+"

    # INS should be treated as DUP (same strand orientation)
    strand1[svtype == "INS"] <- "-"
    strand2[svtype == "INS"] <- "+"

    # Convert to ShatterSeek format
    svtype[svtype == "BND"] <- "TRA"
    svtype[svtype == "INV"] <- "h2hINV"
    svtype[svtype == "INS"] <- "DUP"  # Treat insertions as duplications

    list(
        chrom1 = chrom1,
        pos1 = pos1,
        chrom2 = chrom2,
        pos2 = pos2,
        svtype = svtype,
        strand1 = strand1,
        strand2 = strand2,
        sv_id = names(ranges)
    )
}


#' Parse generic SV VCF (vcfR)
#' @keywords internal
.parse_generic_sv_vcfr <- function(vcf, fix, info_df) {
    chrom1 <- fix[, "CHROM"]
    pos1 <- as.numeric(fix[, "POS"])

    svtype <- if ("SVTYPE" %in% colnames(info_df)) {
        as.character(info_df[, "SVTYPE"])
    } else {
        rep("UNK", nrow(fix))
    }

    chrom2 <- chrom1
    for (field in c("CHR2", "CHROM2", "chr2")) {
        if (field %in% colnames(info_df)) {
            chrom2 <- as.character(info_df[, field])
            break
        }
    }

    pos2 <- pos1
    for (field in c("END", "END2", "POS2")) {
        if (field %in% colnames(info_df)) {
            pos2 <- as.numeric(info_df[, field])
            break
        }
    }

    # Parse BND ALT field if SVTYPE=BND and chrom2/pos2 not set
    is_bnd <- svtype == "BND"
    if (any(is_bnd)) {
        alt <- fix[, "ALT"]

        for (i in which(is_bnd)) {
            # Only parse ALT if we don't have chrom2/pos2 from INFO
            # Check for NA values and default values
            needs_parsing <- FALSE
            if (is.na(chrom2[i]) || is.na(pos2[i])) {
                needs_parsing <- TRUE
            } else if (as.character(chrom2[i]) == as.character(chrom1[i]) && pos2[i] == pos1[i]) {
                needs_parsing <- TRUE
            }

            if (needs_parsing) {
                alt_str <- alt[i]

                # Parse bracket notation: ]chr:pos] or [chr:pos[
                # Patterns: N[chr:pos[, [chr:pos[N, N]chr:pos], ]chr:pos]N
                pattern <- "([\\[\\]])([^:\\[\\]]+):([0-9]+)([\\[\\]])"
                match <- regexec(pattern, alt_str)
                if (match[[1]][1] != -1) {
                    matches <- regmatches(alt_str, match)[[1]]
                    if (length(matches) >= 4) {
                        mate_chr <- matches[3]
                        mate_pos <- as.numeric(matches[4])

                        # Remove "chr" prefix if present
                        mate_chr <- gsub("^chr", "", mate_chr, ignore.case = TRUE)

                        chrom2[i] <- mate_chr
                        pos2[i] <- mate_pos
                    }
                }
            }
        }
    }

    # Remove variants with NA coordinates or chromosomes (vital for robustness)
    valid_idx <- !is.na(chrom1) & !is.na(pos1) & !is.na(chrom2) & !is.na(pos2)
    if (sum(!valid_idx) > 0) {
        message(sprintf("Dropped %d SVs with NA coordinates/chromosomes", sum(!valid_idx)))
        chrom1 <- chrom1[valid_idx]
        pos1 <- pos1[valid_idx]
        chrom2 <- chrom2[valid_idx]
        pos2 <- pos2[valid_idx]
        svtype <- svtype[valid_idx]
    }

    strand1 <- rep("+", length(svtype))
    strand2 <- rep("-", length(svtype))

    strand1[svtype == "DEL"] <- "+"
    strand2[svtype == "DEL"] <- "-"

    strand1[svtype == "DUP"] <- "-"
    strand2[svtype == "DUP"] <- "+"

    strand1[svtype == "INV"] <- "+"
    strand2[svtype == "INV"] <- "+"

    # INS should be treated as DUP (same strand orientation)
    strand1[svtype == "INS"] <- "-"
    strand2[svtype == "INS"] <- "+"

    svtype[svtype == "BND"] <- "TRA"
    svtype[svtype == "INV"] <- "h2hINV"
    svtype[svtype == "INS"] <- "DUP"  # Treat insertions as duplications

    list(
        chrom1 = chrom1,
        pos1 = pos1,
        chrom2 = chrom2,
        pos2 = pos2,
        svtype = svtype,
        strand1 = strand1,
        strand2 = strand2,
        sv_id = fix[, "ID"]
    )
}


# Stub implementations for GRIDSS, LUMPY, and Sniffles
# These follow similar patterns to the above

#' Parse GRIDSS SV VCF (VariantAnnotation)
#' @keywords internal
.parse_gridss_sv <- function(vcf, info, ranges) {
    # GRIDSS uses primarily BND records
    # Implementation would handle MATEID pairing
    .parse_generic_sv(vcf, info, ranges)
}

#' Parse GRIDSS SV VCF (vcfR)
#' @keywords internal
.parse_gridss_sv_vcfr <- function(vcf, fix, info_df) {
    .parse_generic_sv_vcfr(vcf, fix, info_df)
}

#' Parse LUMPY SV VCF (VariantAnnotation)
#' @keywords internal
.parse_lumpy_sv <- function(vcf, info, ranges) {
    # LUMPY uses STRANDS field
    result <- .parse_generic_sv(vcf, info, ranges)

    # Override with STRANDS if available
    if ("STRANDS" %in% colnames(info)) {
        strands <- as.character(info$STRANDS)
        # STRANDS format: "++", "+- ", "-+", "--"
        result$strand1 <- substr(strands, 1, 1)
        result$strand2 <- substr(strands, 2, 2)
    }

    result
}

#' Parse LUMPY SV VCF (vcfR)
#' @keywords internal
.parse_lumpy_sv_vcfr <- function(vcf, fix, info_df) {
    result <- .parse_generic_sv_vcfr(vcf, fix, info_df)

    if ("STRANDS" %in% colnames(info_df)) {
        strands <- as.character(info_df[, "STRANDS"])
        result$strand1 <- substr(strands, 1, 1)
        result$strand2 <- substr(strands, 2, 2)
    }

    result
}

#' Parse Sniffles SV VCF (VariantAnnotation)
#' @keywords internal
.parse_sniffles_sv <- function(vcf, info, ranges) {
    .parse_generic_sv(vcf, info, ranges)
}

#' Parse Sniffles SV VCF (vcfR)
#' @keywords internal
.parse_sniffles_sv_vcfr <- function(vcf, fix, info_df) {
    .parse_generic_sv_vcfr(vcf, fix, info_df)
}