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
