#' Visualization Utility Functions
#'
#' Helper functions for genome-wide visualization and dashboard creation.
#'
#' @keywords internal


#' Get chromosome-level mechanism summary
#'
#' Summarizes chromoanagenesis mechanisms detected on each chromosome.
#'
#' @param chromoanagenesis_result Result from detect_chromoanagenesis()
#' @return Data frame with columns: chrom, dominant_mechanism, n_events, mechanisms_list
#' @keywords internal
get_chromosome_mechanism_summary <- function(chromoanagenesis_result) {

    all_chroms <- c(as.character(1:22), "X")
    summary_list <- list()

    for (chr in all_chroms) {
        mechanisms <- character()
        n_events <- 0

        # Check chromothripsis
        if (!is.null(chromoanagenesis_result$chromothripsis)) {
            ct_class <- chromoanagenesis_result$chromothripsis$classification
            if (!is.null(ct_class) && nrow(ct_class) > 0) {
                ct_chr <- ct_class[ct_class$chrom == chr &
                                  ct_class$classification %in%
                                  c("High confidence", "Low confidence"), ]
                if (nrow(ct_chr) > 0) {
                    mechanisms <- c(mechanisms, "chromothripsis")
                    n_events <- n_events + nrow(ct_chr)
                }
            }
        }

        # Check chromoplexy
        if (!is.null(chromoanagenesis_result$chromoplexy)) {
            cp_summary <- chromoanagenesis_result$chromoplexy$summary
            if (!is.null(cp_summary) && nrow(cp_summary) > 0) {
                # Chromoplexy involves multiple chromosomes
                for (i in seq_len(nrow(cp_summary))) {
                    if (cp_summary$classification[i] %in% c("Likely chromoplexy", "Possible chromoplexy")) {
                        chroms_str <- as.character(cp_summary$chromosomes_involved[i])
                        if (!is.na(chroms_str) && chroms_str != "") {
                            chroms <- unlist(strsplit(chroms_str, ","))
                            chroms <- trimws(chroms)
                            if (chr %in% chroms) {
                                if (!"chromoplexy" %in% mechanisms) {
                                    mechanisms <- c(mechanisms, "chromoplexy")
                                    n_events <- n_events + 1
                                }
                            }
                        }
                    }
                }
            }
        }

        # Check chromoanasynthesis
        if (!is.null(chromoanagenesis_result$chromoanasynthesis)) {
            cs_summary <- chromoanagenesis_result$chromoanasynthesis$summary
            if (!is.null(cs_summary) && nrow(cs_summary) > 0) {
                cs_chr <- cs_summary[cs_summary$chrom == chr &
                                    cs_summary$classification %in%
                                    c("Likely chromoanasynthesis", "Possible chromoanasynthesis"), ]
                if (nrow(cs_chr) > 0) {
                    mechanisms <- c(mechanisms, "chromoanasynthesis")
                    n_events <- n_events + nrow(cs_chr)
                }
            }
        }

        # Determine dominant mechanism
        if (length(mechanisms) == 0) {
            dominant <- "none"
        } else if (length(mechanisms) > 1) {
            dominant <- "mixed"
        } else {
            dominant <- mechanisms[1]
        }

        summary_list[[chr]] <- data.frame(
            chrom = chr,
            dominant_mechanism = dominant,
            n_events = n_events,
            mechanisms_list = paste(mechanisms, collapse = ", "),
            stringsAsFactors = FALSE
        )
    }

    summary_df <- do.call(rbind, summary_list)
    rownames(summary_df) <- NULL

    # Only return chromosomes with events
    summary_df <- summary_df[summary_df$n_events > 0, ]

    return(summary_df)
}


#' Calculate linear genomic positions for visualization
#'
#' Converts per-chromosome positions to linear genome coordinates
#' for genome-wide plots.
#'
#' @param cnv_data CNV data frame with chrom, start, end columns
#' @param genome Genome assembly ("hg19" or "hg38")
#' @return CNV data with added linear_start, linear_end columns
#' @keywords internal
calculate_linear_positions <- function(cnv_data, genome = "hg19") {

    # Approximate chromosome sizes (hg19/hg38)
    chr_sizes <- get_chromosome_sizes(genome)

    # Calculate cumulative offsets
    chr_order <- sort_chromosomes(unique(as.character(cnv_data$chrom)))
    cumulative_offset <- 0
    chr_offsets <- list()

    for (chr in chr_order) {
        chr_offsets[[chr]] <- cumulative_offset
        chr_size <- chr_sizes[[chr]]
        if (is.null(chr_size)) chr_size <- 150e6  # Default fallback
        cumulative_offset <- cumulative_offset + chr_size
    }

    # Add linear positions
    cnv_data$linear_start <- NA
    cnv_data$linear_end <- NA

    for (i in seq_len(nrow(cnv_data))) {
        chr <- as.character(cnv_data$chrom[i])
        offset <- chr_offsets[[chr]]
        if (!is.null(offset)) {
            cnv_data$linear_start[i] <- offset + cnv_data$start[i]
            cnv_data$linear_end[i] <- offset + cnv_data$end[i]
        }
    }

    return(cnv_data)
}


#' Get chromosome boundaries for linear plots
#'
#' Calculates start, end, and midpoint of each chromosome in
#' linear genomic coordinates.
#'
#' @param cnv_data_linear CNV data with linear positions
#' @return Data frame with chrom, start, end, mid columns
#' @keywords internal
get_chromosome_boundaries <- function(cnv_data_linear) {

    chr_bounds_list <- list()

    for (chr in unique(cnv_data_linear$chrom)) {
        chr_data <- cnv_data_linear[cnv_data_linear$chrom == chr, ]
        chr_bounds_list[[as.character(chr)]] <- data.frame(
            chrom = chr,
            start = min(chr_data$linear_start, na.rm = TRUE),
            end = max(chr_data$linear_end, na.rm = TRUE),
            stringsAsFactors = FALSE
        )
    }

    chr_bounds <- do.call(rbind, chr_bounds_list)
    chr_bounds$mid <- (chr_bounds$start + chr_bounds$end) / 2

    # Sort by chromosome order
    chr_order <- sort_chromosomes(unique(as.character(chr_bounds$chrom)))
    chr_bounds$chrom <- factor(chr_bounds$chrom, levels = chr_order)
    chr_bounds <- chr_bounds[order(chr_bounds$chrom), ]
    rownames(chr_bounds) <- NULL

    return(chr_bounds)
}


#' Annotate CNV segments with chromoanagenesis status
#'
#' Marks CNV segments that overlap with chromoanagenesis events.
#'
#' @param cnv_data CNV data frame
#' @param chromoanagenesis_result Chromoanagenesis results
#' @return CNV data with added has_chromoanagenesis column
#' @keywords internal
annotate_chromoanagenesis_regions <- function(cnv_data, chromoanagenesis_result) {

    cnv_data$has_chromoanagenesis <- FALSE

    # Mark chromothripsis regions
    if (!is.null(chromoanagenesis_result$chromothripsis)) {
        ct_class <- chromoanagenesis_result$chromothripsis$classification
        if (!is.null(ct_class) && nrow(ct_class) > 0) {
            ct_chroms <- ct_class$chrom[ct_class$classification %in%
                                       c("High confidence", "Low confidence")]
            cnv_data$has_chromoanagenesis[cnv_data$chrom %in% ct_chroms] <- TRUE
        }
    }

    # Mark chromoanasynthesis regions
    if (!is.null(chromoanagenesis_result$chromoanasynthesis)) {
        cs_summary <- chromoanagenesis_result$chromoanasynthesis$summary
        if (!is.null(cs_summary) && nrow(cs_summary) > 0) {
            cs_regions <- cs_summary[cs_summary$classification %in%
                                    c("Likely chromoanasynthesis", "Possible chromoanasynthesis"), ]
            for (i in seq_len(nrow(cs_regions))) {
                chr <- cs_regions$chrom[i]
                start <- cs_regions$start[i]
                end <- cs_regions$end[i]
                cnv_data$has_chromoanagenesis[
                    cnv_data$chrom == chr &
                    cnv_data$start >= start &
                    cnv_data$end <= end
                ] <- TRUE
            }
        }
    }

    # Note: Chromoplexy is multi-chromosomal, more complex to annotate
    # Could mark all chromosomes involved in a chain

    return(cnv_data)
}


#' Get approximate chromosome sizes
#'
#' Returns approximate chromosome sizes for hg19 or hg38.
#'
#' @param genome Genome assembly ("hg19" or "hg38")
#' @return Named list of chromosome sizes
#' @keywords internal
get_chromosome_sizes <- function(genome = "hg19") {
    if (genome == "hg38") {
        return(list(
            "1" = 248956422, "2" = 242193529, "3" = 198295559,
            "4" = 190214555, "5" = 181538259, "6" = 170805979,
            "7" = 159345973, "8" = 145138636, "9" = 138394717,
            "10" = 133797422, "11" = 135086622, "12" = 133275309,
            "13" = 114364328, "14" = 107349540, "15" = 101991189,
            "16" = 90338345, "17" = 83257441, "18" = 80373285,
            "19" = 58617616, "20" = 64444167, "21" = 46709983,
            "22" = 50818468, "X" = 156040895
        ))
    }
    # Default: hg19
    list(
        "1" = 249250621, "2" = 243199373, "3" = 198022430,
        "4" = 191154276, "5" = 180915260, "6" = 171115067,
        "7" = 159138663, "8" = 146364022, "9" = 141213431,
        "10" = 135534747, "11" = 135006516, "12" = 133851895,
        "13" = 115169878, "14" = 107349540, "15" = 102531392,
        "16" = 90354753, "17" = 81195210, "18" = 78077248,
        "19" = 59128983, "20" = 63025520, "21" = 48129895,
        "22" = 51304566, "X" = 155270560
    )
}


#' Get centromere positions
#'
#' Returns approximate centromere positions (hg19 or hg38) for karyogram visualization.
#'
#' @param genome Genome assembly ("hg19" or "hg38")
#' @return Data frame with chrom, size, centromere columns
#' @keywords internal
get_centromere_positions <- function(genome = "hg19") {
    if (genome == "hg38") {
        return(data.frame(
            chrom = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
                     "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"),
            size = c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
                    159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
                    114364328, 107349540, 101991189, 90338345, 83257441, 80373285,
                    58617616, 64444167, 46709983, 50818468, 156040895),
            centromere = c(123400000, 93900000, 90900000, 50000000, 48800000, 59800000,
                          60100000, 45200000, 43000000, 39800000, 53400000, 35500000,
                          17700000, 17200000, 17000000, 36800000, 25100000, 18500000,
                          26200000, 28100000, 12000000, 15000000, 61000000),
            stringsAsFactors = FALSE
        ))
    }
    # Default: hg19
    data.frame(
        chrom = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
                 "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"),
        size = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067,
                159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
                115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
                59128983, 63025520, 48129895, 51304566, 155270560),
        centromere = c(125000000, 93300000, 91000000, 50400000, 48400000, 61000000,
                      59900000, 45600000, 49000000, 40200000, 53700000, 35800000,
                      17900000, 17600000, 19000000, 36600000, 24000000, 17200000,
                      26500000, 27500000, 13200000, 14700000, 60600000),
        stringsAsFactors = FALSE
    )
}


#' Format genomic position for axis labels
#'
#' Converts genomic positions to human-readable format (Mb).
#'
#' @param pos Numeric vector of positions
#' @return Character vector of formatted positions
#' @keywords internal
format_genomic_position <- function(pos) {
    ifelse(pos >= 1e6,
           paste0(round(pos / 1e6, 1), " Mb"),
           paste0(round(pos / 1e3, 0), " Kb"))
}


#' Get mechanism color palette
#'
#' Returns consistent color scheme for chromoanagenesis mechanisms.
#'
#' @return Named vector of colors
#' @keywords internal
get_mechanism_colors <- function() {
    c(
        "chromothripsis" = "#E41A1C",
        "chromoplexy" = "#377EB8",
        "chromoanasynthesis" = "#4DAF4A",
        "mixed" = "#984EA3",
        "none" = "gray80"
    )
}
