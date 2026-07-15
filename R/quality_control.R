#' @include chromoanagenesis_integrated.R
NULL

#' Assess input data quality for Chromoanagenesis analysis
#'
#' Calculates quality metrics for SV and CNV inputs, including counts,
#' distributions, and consistency checks.
#'
#' @param SV.sample Structural variant data (SVs object or data.frame)
#' @param CNV.sample Optional copy number data (CNVsegs object or data.frame)
#' @param genome Reference genome ("hg19" or "hg38")
#' @return A list containing quality metrics
#' @export
assess_data_quality <- function(SV.sample, CNV.sample = NULL, genome = "hg19") {
    
    # Standardize inputs
    if (is(SV.sample, "SVs")) {
        sv_data_df <- as(SV.sample, "data.frame")
    } else {
        sv_data_df <- SV.sample
    }

    if (is.null(CNV.sample)) {
        cnv_data_df <- data.frame(
            chrom = character(0),
            start = numeric(0),
            end = numeric(0),
            total_cn = numeric(0),
            stringsAsFactors = FALSE
        )
    } else if (is(CNV.sample, "CNVsegs")) {
        cnv_data_df <- as(CNV.sample, "data.frame")
    } else {
        cnv_data_df <- CNV.sample
    }

    metrics <- list()

    # 1. SV Statistics
    metrics$sv_count <- nrow(sv_data_df)
    if (metrics$sv_count > 0) {
        metrics$sv_types <- table(sv_data_df$SVtype)
        metrics$interchromosomal_rate <- sum(sv_data_df$chrom1 != sv_data_df$chrom2) / metrics$sv_count
    } else {
        metrics$sv_types <- c()
        metrics$interchromosomal_rate <- 0
    }

    # 2. CNV Statistics
    metrics$cnv_count <- nrow(cnv_data_df)
    if (metrics$cnv_count > 0) {
        metrics$avg_segment_size <- mean(cnv_data_df$end - cnv_data_df$start)
        # Estimate noise: number of very small segments (< 10kb)
        metrics$small_segments <- sum((cnv_data_df$end - cnv_data_df$start) < 10000)
        metrics$hyper_segmentation_risk <- metrics$small_segments / metrics$cnv_count
    } else {
        metrics$avg_segment_size <- NA_real_
        metrics$small_segments <- 0
        metrics$hyper_segmentation_risk <- NA_real_
    }

    # 3. Consistency Checks
    # Check chromosome naming
    sv_chroms <- unique(c(sv_data_df$chrom1, sv_data_df$chrom2))
    cnv_chroms <- unique(cnv_data_df$chrom)
    
    has_chr_prefix_sv <- any(grepl("^chr", sv_chroms))
    has_chr_prefix_cnv <- any(grepl("^chr", cnv_chroms))
    
    metrics$naming_consistent <- if (metrics$cnv_count > 0) {
        has_chr_prefix_sv == has_chr_prefix_cnv
    } else {
        NA
    }
    
    # 4. Overlap Analysis (Basic)
    # Check how many SV breakpoints fall within defined CNV segments
    if (metrics$sv_count > 0 && metrics$cnv_count > 0) {
        # Normalize chrom names for overlap check
        sv_c1 <- gsub("^chr", "", sv_data_df$chrom1)
        sv_p1 <- sv_data_df$pos1
        cnv_c <- gsub("^chr", "", cnv_data_df$chrom)
        
        # Simple check for first breakpoints (can be slow for large datasets, so we sample if large)
        n_check <- min(1000, nrow(sv_data_df))
        idx_check <- sample(nrow(sv_data_df), n_check)
        
        overlaps <- 0
        for (i in idx_check) {
            c <- sv_c1[i]
            p <- sv_p1[i]
            # Check if this point is in any segment on the same chromosome
            if (any(cnv_c == c & cnv_data_df$start <= p & cnv_data_df$end >= p)) {
                overlaps <- overlaps + 1
            }
        }
        metrics$sv_cnv_overlap_rate <- overlaps / n_check
    } else {
        metrics$sv_cnv_overlap_rate <- NA_real_
    }

    return(metrics)
}

#' Flag inter-chromosomal breakpoints that cluster unusually tightly within
#' a single sample
#'
#' Chromoplexy chain-finding builds its translocation graph from
#' inter-chromosomal (BND) breakpoints and assumes each one represents a
#' distinct, real rearrangement junction. When several independent BND
#' records land within a very small window (default 1000 bp) of each other,
#' that is more consistent with alignment ambiguity in a
#' repetitive/segmental-duplication region -- which produces multiple
#' slightly different candidate breakpoint calls for what is likely the same
#' underlying locus -- than with several independent translocation
#' junctions. Unlike \code{\link{flag_recurrent_chromoplexy_artifacts}}
#' (which needs multiple unrelated cohort samples to detect recurrence),
#' this check works from a single sample and is used by
#' \code{detect_chromoanagenesis(germline_mode = TRUE)} to pre-filter likely
#' artifacts before chain-finding.
#'
#' Only inter-chromosomal SVs (\code{chrom1 != chrom2}) are considered.
#' Intra-chromosomal clustering (small tandem DELs/DUPs/INVs sitting close
#' together) is deliberately left untouched: it is exactly the kind of local
#' signal chromothripsis and chromoanasynthesis detection legitimately rely
#' on, and both already have their own statistical/structural gates.
#'
#' @param SV.sample An SVs object.
#' @param window Window size in bp within which breakpoints on the same
#'   chromosome are considered part of the same cluster. Default 1000.
#' @param min_count Minimum number of distinct SV records in a cluster
#'   before it is flagged. Default 2.
#' @return A list with \code{flagged_sv_ids} (character vector of
#'   \code{sv_id} values with at least one breakend in a flagged cluster)
#'   and \code{clusters} (a data frame describing each flagged cluster:
#'   chrom, start, end, n_sv, sv_ids).
#' @export
flag_locally_clustered_breakpoints <- function(SV.sample, window = 1000, min_count = 2) {
    if (!is(SV.sample, "SVs")) stop("SV.sample must be an SVs object")

    inter_chr <- SV.sample@chrom1 != SV.sample@chrom2
    sides <- unique(rbind(
        data.frame(chrom = SV.sample@chrom1[inter_chr], pos = SV.sample@pos1[inter_chr],
                  sv_id = SV.sample@sv_id[inter_chr], stringsAsFactors = FALSE),
        data.frame(chrom = SV.sample@chrom2[inter_chr], pos = SV.sample@pos2[inter_chr],
                  sv_id = SV.sample@sv_id[inter_chr], stringsAsFactors = FALSE)
    ))

    cluster_rows <- list()
    for (ch in unique(sides$chrom)) {
        sub <- sides[sides$chrom == ch, , drop = FALSE]
        if (nrow(sub) < min_count) next
        sub <- sub[order(sub$pos), , drop = FALSE]
        # Single-linkage clustering by position gap; O(n log n) via one sort,
        # no per-breakpoint linear scans.
        cluster_id <- cumsum(c(Inf, diff(sub$pos)) > window)

        split_sub <- split(sub, cluster_id)
        for (grp in split_sub) {
            n_unique_sv <- length(unique(grp$sv_id))
            if (n_unique_sv >= min_count) {
                cluster_rows[[length(cluster_rows) + 1]] <- data.frame(
                    chrom = ch,
                    start = min(grp$pos),
                    end = max(grp$pos),
                    n_sv = n_unique_sv,
                    sv_ids = paste(sort(unique(grp$sv_id)), collapse = ","),
                    stringsAsFactors = FALSE
                )
            }
        }
    }

    if (length(cluster_rows) == 0) {
        return(list(flagged_sv_ids = character(0), clusters = data.frame()))
    }

    clusters <- do.call(rbind, cluster_rows)
    clusters <- clusters[order(-clusters$n_sv, clusters$chrom, clusters$start), , drop = FALSE]
    flagged_sv_ids <- unique(unlist(strsplit(clusters$sv_ids, ",", fixed = TRUE)))

    list(flagged_sv_ids = flagged_sv_ids, clusters = clusters)
}
