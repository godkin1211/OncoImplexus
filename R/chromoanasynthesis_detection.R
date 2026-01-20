' Detect chromoanasynthesis events from SV and CNV data
#"
#' 
#' Chromoanasynthesis is characterized by replication-based mechanisms
#' (Fork Stalling and Template Switching - FoSTeS) resulting in
#' complex local rearrangements, linked SV chains, and copy number
#' gradients, often distinct from BFB cycles.
#'
#' @param SV.sample An instance of class SVs or data frame with SV data
#' @param CNV.sample An instance of class CNVsegs or data frame with CNV data
#' @param min_tandem_dups Minimum number of tandem duplications (default: 3)
#' @param min_cn_segments Minimum CN segments for gradient analysis (default: 5)
#' @param gradient_threshold Minimum gradient correlation (default: 0.4) - Relaxed to allow for "jumpy" profiles
#' @param max_region_size Maximum region size in bp (default: 10e6)
#' @param cn_segment_window_extension Integer specifying how many segments beyond `min_cn_segments` to consider for window size increments (default: 20).
#' @param cn_window_step_ratio Integer specifying the divisor for calculating the sliding window step size (default: 4, meaning step = win_size / 4).
#' @return A list containing chromoanasynthesis detection results
#' @details
#' Enhanced detection methodology v2.0:
#' 1. **CN Gradient Analysis**: Identifies regions with copy number instability and general increasing trend.
#' 2. **SV Topology Analysis**:
#'    - **Linked Hops**: Detects chains of SVs representing serial template switching (hallmark of FoSTeS).
#'    - **Fold-back Inversions**: Identifies short-span inversions characteristic of BFB cycles (to distinguish from chromoanasynthesis).
#'    - **Tandem Duplications**: Classic enrichment analysis.
#'
#' Classification now penalizes regions with high fold-back inversion content (likely BFB)
#' and rewards regions with clear linked SV chains.
#'
#' @references
#' Liu et al. (2011) Cell. Chromosome catastrophes involve replication mechanisms generating complex genomic rearrangements.
#' Lee et al. (2015) Nat Genet. Complex chromosomal rearrangements by single catastrophic pathways in cancer genomes.
#'
#' @export
detect_chromoanasynthesis <- function(SV.sample,
                                  CNV.sample,
                                  min_tandem_dups = 3,
                                  min_cn_segments = 5,
                                  gradient_threshold = 0.4,
                                  max_region_size = 10e6,
                                  cn_segment_window_extension = 20, # New parameter
                                  cn_window_step_ratio = 4) { # New parameter

    # Convert to data frames if needed
    if (is(SV.sample, "SVs")) {
        SV.sample <- as(SV.sample, "data.frame")
    }

    if (is.null(CNV.sample)) {
        stop("CNV data is required for chromoanasynthesis detection.")
    }

    if (is(CNV.sample, "CNVsegs")) {
        CNV.sample <- as(CNV.sample, "data.frame")
    }

    # Standardize SV column names
    # Support both pos1/pos2 and start1/start2 formats
    if (!"pos1" %in% names(SV.sample) && "start1" %in% names(SV.sample)) {
        SV.sample$pos1 <- SV.sample$start1
    }
    if (!"pos2" %in% names(SV.sample) && "start2" %in% names(SV.sample)) {
        SV.sample$pos2 <- SV.sample$start2
    }
    # Support both SVtype and svclass formats
    if (!"SVtype" %in% names(SV.sample) && "svclass" %in% names(SV.sample)) {
        SV.sample$SVtype <- SV.sample$svclass
    }

    # Standardize CNV column names
    # Support both chrom and chromosome formats
    if (!"chrom" %in% names(CNV.sample) && "chromosome" %in% names(CNV.sample)) {
        CNV.sample$chrom <- CNV.sample$chromosome
    }

    # Validate required columns exist
    sv_required <- c("chrom1", "chrom2", "pos1", "pos2", "SVtype")
    sv_missing <- setdiff(sv_required, names(SV.sample))
    if (length(sv_missing) > 0) {
        stop(sprintf("Missing required SV columns: %s\nAvailable columns: %s",
                    paste(sv_missing, collapse = ", "),
                    paste(names(SV.sample), collapse = ", ")))
    }

    cnv_required <- c("chrom", "start", "end", "total_cn")
    cnv_missing <- setdiff(cnv_required, names(CNV.sample))
    if (length(cnv_missing) > 0) {
        stop(sprintf("Missing required CNV columns: %s\nAvailable columns: %s",
                    paste(cnv_missing, collapse = ", "),
                    paste(names(CNV.sample), collapse = ", ")))
    }

    if (nrow(CNV.sample) < min_cn_segments) {
        warning(sprintf("Insufficient CNV segments (%d). Need at least %d for chromoanasynthesis detection.",
                       nrow(CNV.sample), min_cn_segments))
        return(create_empty_chromoanasynthesis_result())
    }

    cat("Detecting chromoanasynthesis patterns (Enhanced v2.0)...
")

    # Step 1: Identify regions with CN gradients / instability
    cat("Step 1: Analyzing copy number gradients and instability...
")
    gradient_regions <- detect_cn_gradients(
        CNV.sample,
        min_segments = min_cn_segments,
        gradient_threshold = gradient_threshold,
        max_region_size = max_region_size,
        cn_segment_window_extension = cn_segment_window_extension, # Pass new parameter
        cn_window_step_ratio = cn_window_step_ratio # Pass new parameter
    )

    if (length(gradient_regions) == 0) {
        cat("No significant CN gradient regions found.
")
        return(create_empty_chromoanasynthesis_result())
    }

    cat(sprintf("Found %d raw candidate window(s) with CN gradients.\n", length(gradient_regions)))

    # Step 1.5: Merge/Deduplicate overlapping gradient regions to save topology analysis time
    if (length(gradient_regions) > 1) {
        cat("Deduplicating overlapping regions...")
        gradient_regions <- remove_overlapping_windows(gradient_regions)
        cat(sprintf(" Reduced to %d distinct region(s).\n", length(gradient_regions)))
    }

    if (length(gradient_regions) == 0) {
        cat("No significant CN gradient regions after filtering.\n")
        return(create_empty_chromoanasynthesis_result())
    }

    # Step 2: Analyze SV Topology (Tandem Dups, Linked Chains, Fold-backs)
    # NOTE: Function moved to R/chromoanasynthesis_topology.R
    cat("Step 2: Analyzing SV topology (Linked Chains & Fold-backs)...
")
    sv_topology_analysis <- analyze_sv_topology(
        SV.sample,
        gradient_regions
    )

    # Step 3: Evaluate each region
    cat("Step 3: Evaluating chromoanasynthesis criteria...
")
    region_results <- list()

    for (i in 1:length(gradient_regions)) {
        region_results[[i]] <- evaluate_chromoanasynthesis_region(
            region = gradient_regions[[i]],
            SV.sample = SV.sample,
            CNV.sample = CNV.sample,
            topology = sv_topology_analysis[[i]],
            min_tandem_dups = min_tandem_dups
        )
    }

    # Step 4: Create summary
    summary_df <- do.call(rbind, lapply(region_results, function(x) x$summary))

    # Classify events
    summary_df$classification <- classify_chromoanasynthesis_event(summary_df)

    result <- list(
        regions = gradient_regions,
        region_details = region_results,
        summary = summary_df,
        total_regions = length(gradient_regions),
        likely_chromoanasynthesis = sum(summary_df$classification == "Likely chromoanasynthesis"),
        possible_chromoanasynthesis = sum(summary_df$classification == "Possible chromoanasynthesis")
    )

    class(result) <- c("chromoanasynthesis", "list")

    cat(sprintf("\nDetection complete: %d likely, %d possible chromoanasynthesis events.
",
               result$likely_chromoanasynthesis,
               result$possible_chromoanasynthesis))

    return(result)
}


#' Detect copy number gradients in chromosomes
#'
#' Identifies regions with gradual CN increases or high instability characteristic of chromoanasynthesis.
#'
#' @keywords internal
detect_cn_gradients <- function(CNV.sample,
                                min_segments = 5,
                                gradient_threshold = 0.4,
                                max_region_size = 10e6,
                                cn_segment_window_extension = 20, # New parameter
                                cn_window_step_ratio = 4) { # New parameter

    gradient_regions <- list()
    region_id <- 0

    chroms <- unique(CNV.sample$chrom)

    for (chr in chroms) {
        chr_cnv <- CNV.sample[CNV.sample$chrom == chr, ]

        if (nrow(chr_cnv) < min_segments) next

        # Sort by position
        chr_cnv <- chr_cnv[order(chr_cnv$start), ]

        # Use sliding window to detect gradient regions
        window_results <- find_gradient_windows(
            chr_cnv,
            min_segments = min_segments,
            gradient_threshold = gradient_threshold,
            max_size = max_region_size,
            cn_segment_window_extension = cn_segment_window_extension, # Pass new parameter
            cn_window_step_ratio = cn_window_step_ratio # Pass new parameter
        )

        if (length(window_results) > 0) {
            for (window in window_results) {
                region_id <- region_id + 1
                window$region_id <- region_id
                gradient_regions[[region_id]] <- window
            }
        }
    }

    return(gradient_regions)
}


#' Find windows with CN gradients using sliding window approach
#'
#' @keywords internal
find_gradient_windows <- function(chr_cnv,
                                  min_segments = 5,
                                  gradient_threshold = 0.4,
                                  max_size = 10e6,
                                  cn_segment_window_extension = 20, # New parameter
                                  cn_window_step_ratio = 4) { # New parameter

    windows <- list()
    n_segs <- nrow(chr_cnv)

    # Simplified sliding window strategy for efficiency
    # Step by 1, but window sizes in increments
    window_sizes <- unique(round(seq(min_segments, min(n_segs, min_segments + cn_segment_window_extension), length.out = 5)))
    
    for (win_size in window_sizes) {
        
        # Optimization: Don't slide by 1 for large windows, slide by win_size/2
        step <- max(1, floor(win_size / cn_window_step_ratio))
        
        for (start_idx in seq(1, n_segs - win_size + 1, by = step)) {
            end_idx <- start_idx + win_size - 1
            window_cnv <- chr_cnv[start_idx:end_idx, ]
            
            # Check window size constraint
            region_start <- window_cnv$start[1]
            region_end <- window_cnv$end[nrow(window_cnv)]
            region_size <- region_end - region_start

            if (region_size > max_size) next

            # Calculate gradient correlation
            segment_order <- 1:nrow(window_cnv)
            cn_values <- window_cnv$total_cn
            
            # Calculate Total Variation (Jumps)
            # Sum of absolute differences between adjacent segments
            total_variation <- sum(abs(diff(cn_values)))
            mean_cn <- mean(cn_values)
            variation_score <- if (!is.na(mean_cn) && mean_cn > 0) total_variation / mean_cn else 0

            # Spearman correlation with significance check
            gradient_cor <- NA
            gradient_pval <- 1.0
            n_segs_win <- nrow(window_cnv)
            
            if (length(unique(cn_values)) > 1) {
                test_res <- try(cor.test(segment_order, cn_values, method = "spearman", exact = FALSE), silent = TRUE)
                if (!inherits(test_res, "try-error")) {
                    gradient_cor <- test_res$estimate
                    gradient_pval <- test_res$p.value
                }
            }
            
            # Robust criteria for candidate selection
            is_candidate <- FALSE
            if (!is.na(gradient_cor)) {
                # 1. Significant positive correlation (for large windows)
                if (n_segs_win >= 6 && gradient_pval < 0.05 && gradient_cor > 0.3) {
                    is_candidate <- TRUE
                } 
                # 2. Strong correlation even if small window
                else if (gradient_cor >= 0.8) {
                    is_candidate <- TRUE
                }
                # 3. High variation/instability with moderate positive trend
                else if (gradient_cor > 0.3 && variation_score > 1.8) {
                    is_candidate <- TRUE
                }
            }
            
            if (is_candidate) {
                cn_range <- max(cn_values) - min(cn_values)
                cn_trend <- lm(cn_values ~ segment_order)
                slope <- coef(cn_trend)[2]

                # Only keep if there's actual increase or high activity
                if (slope > 0.05 || variation_score > 2.0) {
                    windows[[length(windows) + 1]] <- list(
                        chrom = window_cnv$chrom[1],
                        start = region_start,
                        end = region_end,
                        size = region_size,
                        n_segments = nrow(window_cnv),
                        gradient_correlation = if(is.na(gradient_cor)) 0 else gradient_cor,
                        cn_range = cn_range,
                        cn_mean = mean_cn,
                        slope = slope,
                        variation_score = variation_score,
                        segments = window_cnv
                    )
                }
            }
        }
    }

    # Remove overlapping windows, keep those with highest combined score
    if (length(windows) > 0) {
        windows <- remove_overlapping_windows(windows)
    }

    return(windows)
}


#' Remove redundant overlapping windows
#'
#' @keywords internal
remove_overlapping_windows <- function(windows, overlap_cutoff = 0.8) {
    if (length(windows) <= 1) return(windows)
    
    # Sort windows by chromosome and start position
    chroms <- sapply(windows, function(w) w$chrom)
    starts <- sapply(windows, function(w) w$start)
    ord <- order(chroms, starts)
    windows <- windows[ord]
    
    keep <- rep(TRUE, length(windows))
    
    for (i in 1:(length(windows)-1)) {
        if (!keep[i]) next
        
        for (j in (i+1):length(windows)) {
            if (windows[[i]]$chrom != windows[[j]]$chrom) break
            
            # Calculate overlap
            ovlp_start <- max(windows[[i]]$start, windows[[j]]$start)
            ovlp_end <- min(windows[[i]]$end, windows[[j]]$end)
            
            if (ovlp_end > ovlp_start) {
                ovlp_size <- ovlp_end - ovlp_start
                size_i <- windows[[i]]$end - windows[[i]]$start
                size_j <- windows[[j]]$end - windows[[j]]$start
                
                # If overlap is high relative to both, discard one
                if (ovlp_size / size_i > overlap_cutoff && ovlp_size / size_j > overlap_cutoff) {
                    # Keep the one with better correlation or more segments
                    if (windows[[i]]$gradient_correlation >= windows[[j]]$gradient_correlation) {
                        keep[j] <- FALSE
                    } else {
                        keep[i] <- FALSE
                        break
                    }
                }
            } else {
                # Windows don't overlap and are sorted, so no subsequent windows will overlap with i
                break
            }
        }
    }
    
    return(windows[keep])
}


#' Analyze SV Topology: Tandem Dups, Linked Chains, Fold-backs
#'
#' @keywords internal
analyze_sv_topology <- function(SV.sample, gradient_regions) {

    results_list <- list()

    for (i in 1:length(gradient_regions)) {
        region <- gradient_regions[[i]]

        # Extract SVs in this region (Intra-chromosomal only for topology)
        region_svs <- SV.sample[
            SV.sample$chrom1 == region$chrom &
            SV.sample$chrom2 == region$chrom &
            SV.sample$pos1 >= region$start &
            SV.sample$pos1 <= region$end &
            SV.sample$pos2 >= region$start &
            SV.sample$pos2 <= region$end,
        ]
        
        # 1. Tandem Duplications
        tandem_dups <- region_svs[
            region_svs$SVtype == "DUP" &
            abs(region_svs$pos2 - region_svs$pos1) < 1e6,
        ]

        # 2. Fold-back Inversions (BFB signature)
        # Definition: Inversions with very small span (< 20kb)
        fold_backs <- region_svs[
            region_svs$SVtype %in% c("h2hINV", "t2tINV") &
            abs(region_svs$pos2 - region_svs$pos1) < 20000,
        ]
        
        # 3. Linked SV Chains (FoSTeS signature)
        linked_chains <- detect_linked_sv_chains(region_svs, max_gap = 2000)

        results_list[[i]] <- list(
            n_tandem_dups = nrow(tandem_dups),
            n_fold_backs = nrow(fold_backs),
            n_linked_chains = length(linked_chains),
            max_chain_length = if(length(linked_chains)>0) max(sapply(linked_chains, length)) else 0,
            n_total_svs = nrow(region_svs),
            tandem_dups = tandem_dups,
            fold_backs = fold_backs,
            linked_chains = linked_chains,
            all_svs = region_svs
        )
    }

    return(results_list)
}


#' Detect chains of linked SVs (template switching events)
#'
#' @keywords internal
detect_linked_sv_chains <- function(sv_df, max_gap = 2000) {
    
    if (nrow(sv_df) < 2) return(list())

    # Ensure GenomicRanges and IRanges are available
    if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
        stop("Package 'GenomicRanges' needed for this function to work. Please install it.", call. = FALSE)
    }
    if (!requireNamespace("IRanges", quietly = TRUE)) {
        stop("Package 'IRanges' needed for this function to work. Please install it.", call. = FALSE)
    }
    if (!requireNamespace("S4Vectors", quietly = TRUE)) {
        stop("Package 'S4Vectors' needed for this function to work. Please install it.", call. = FALSE)
    }

    # Create GRanges for all breakpoints. Each SV has two breakpoints (pos1, pos2)
    # This approach correctly captures proximity between any breakpoint of one SV to any of another.
    bps_df <- rbind(
        data.frame(
            chrom = sv_df$chrom1,
            pos = sv_df$pos1,
            sv_idx = 1:nrow(sv_df),
            stringsAsFactors = FALSE
        ),
        data.frame(
            chrom = sv_df$chrom2, # assuming chrom1 == chrom2 for intra-chr SVs
            pos = sv_df$pos2,
            sv_idx = 1:nrow(sv_df),
            stringsAsFactors = FALSE
        )
    )

    # Create a GRanges object from all breakpoints
    bps_gr <- GenomicRanges::GRanges(
        seqnames = S4Vectors::Rle(bps_df$chrom),
        ranges = IRanges::IRanges(start = bps_df$pos, end = bps_df$pos),
        sv_idx = bps_df$sv_idx
    )

    # Find overlaps between *any* breakpoint and *any other* breakpoint within max_gap
    # maxgap parameter finds intervals separated by at most max_gap bases.
    hits <- GenomicRanges::findOverlaps(bps_gr, bps_gr, maxgap = max_gap, type = "any", ignore.strand = TRUE)

    # Convert hits to a data frame for easier processing
    hits_df <- as.data.frame(hits)

    # Filter out self-hits (a breakpoint with itself)
    # Also filter out connections between breakpoints belonging to the *same* SV
    hits_df <- hits_df[bps_gr$sv_idx[hits_df$queryHits] != bps_gr$sv_idx[hits_df$subjectHits], ]

    if (nrow(hits_df) == 0) return(list())

    # Construct adjacency list efficiently
    sv1_indices <- bps_gr$sv_idx[hits_df$queryHits]
    sv2_indices <- bps_gr$sv_idx[hits_df$subjectHits]
    
    # Bulk aggregate to avoid loop-based unique() calls
    all_edges <- data.frame(from = sv1_indices, to = sv2_indices)
    # Ensure undirected consistency and remove duplicates in one go
    all_edges <- unique(all_edges)
    
    # Fast split into list
    adj_list_sv <- split(all_edges$to, all_edges$from)
    
    # Find connected components (chains) using BFS/DFS on the SV-level adjacency list
    # Use names of adj_list_sv to check connectivity
    unvisited_ids <- as.integer(names(adj_list_sv))
    visited_sv_map <- rep(FALSE, nrow(sv_df))
    chains <- list()
    
    for (start_id in unvisited_ids) {
        if (!visited_sv_map[start_id]) {
            current_chain <- c()
            queue <- c(start_id)
            visited_sv_map[start_id] <- TRUE
            
            while(length(queue) > 0) {
                curr_sv_idx <- queue[1]
                queue <- queue[-1]
                current_chain <- c(current_chain, curr_sv_idx)
                
                # Fast lookup from named list
                neighbors <- adj_list_sv[[as.character(curr_sv_idx)]]
                
                if (!is.null(neighbors)) {
                    for (nbr_sv_idx in neighbors) {
                        if (nbr_sv_idx <= length(visited_sv_map) && !visited_sv_map[nbr_sv_idx]) {
                            visited_sv_map[nbr_sv_idx] <- TRUE
                            queue <- c(queue, nbr_sv_idx)
                        }
                    }
                }
            }
            
            if (length(current_chain) >= 3) { # Only keep chains of 3+ SVs
                chains[[length(chains) + 1]] <- sort(unique(current_chain))
            }
        }
    }
    
    return(chains)
}


#' Evaluate a chromoanasynthesis region
#'
#' @keywords internal
evaluate_chromoanasynthesis_region <- function(region,
                                           SV.sample,
                                           CNV.sample,
                                           topology,
                                           min_tandem_dups = 3) {

    # Calculate scores
    
    # 1. Fold-back penalty (High fold-backs suggest BFB)
    fold_back_ratio <- if (topology$n_total_svs > 0) topology$n_fold_backs / topology$n_total_svs else 0
    bfb_penalty <- min(fold_back_ratio * 2, 1.0) # Penalty up to 1.0 if >50% fold-backs
    
    # 2. Linked Chain bonus (Strong evidence of FoSTeS)
    chain_bonus <- min(topology$max_chain_length / 5, 1.0) # Bonus up to 1.0 for chains of length 5+

    # Calculate complexity score
    complexity_score <- calculate_chromoanasynthesis_complexity(
        region = region,
        topology = topology
    )

    # Check for inter-chromosomal events (should be minimal)
    inter_chr_svs <- SV.sample[
        SV.sample$chrom1 != SV.sample$chrom2 &
        ((SV.sample$chrom1 == region$chrom &
          SV.sample$pos1 >= region$start &
          SV.sample$pos1 <= region$end) |
         (SV.sample$chrom2 == region$chrom &
          SV.sample$pos2 >= region$start &
          SV.sample$pos2 <= region$end)),
    ]

    n_inter_chr <- nrow(inter_chr_svs)

    # Create summary
    summary <- data.frame(
        region_id = if (!is.null(region$region_id)) region$region_id else 1,
        chrom = region$chrom,
        start = region$start,
        end = region$end,
        size = region$size,
        n_cn_segments = region$n_segments,
        gradient_correlation = region$gradient_correlation,
        cn_variation = region$variation_score,
        
        # Topology stats
        n_tandem_dups = topology$n_tandem_dups,
        n_fold_backs = topology$n_fold_backs,
        n_linked_chains = topology$n_linked_chains,
        max_chain_length = topology$max_chain_length,
        
        n_inter_chr_svs = n_inter_chr,
        complexity_score = complexity_score,
        bfb_risk = fold_back_ratio,
        stringsAsFactors = FALSE
    )

    return(list(
        summary = summary,
        region = region,
        topology = topology
    ))
}


#' Calculate chromoanasynthesis complexity score
#'
#' @keywords internal
calculate_chromoanasynthesis_complexity <- function(region, topology) {

    # Complexity based on:
    # 1. Gradient strength / Variation
    # 2. Tandem duplication presence
    # 3. Linked Chains (New!)
    # 4. Fold-back presence (Negative impact on "Chromoanasynthesis" score, likely BFB)

    gradient_score <- region$gradient_correlation
    variation_score <- min(region$variation_score / 3, 1.0)
    
    tandem_score <- min(topology$n_tandem_dups / 5, 1.0)
    chain_score <- min(topology$max_chain_length / 5, 1.0)
    
    # Base score
    complexity <- (
        gradient_score * 0.2 +
        variation_score * 0.2 +
        tandem_score * 0.2 +
        chain_score * 0.4    # High weight for chains
    )
    
    return(complexity)
}


#' Classify chromoanasynthesis event
#'
#' @keywords internal
classify_chromoanasynthesis_event <- function(summary_df) {

    classifications <- character(nrow(summary_df))

    for (i in 1:nrow(summary_df)) {
        row <- summary_df[i, ]

        # Criteria for chromoanasynthesis (Refined):
        # 1. CN Profile: Gradient OR High Variation
        # 2. Topology: Linked Chains OR Multiple Tandem Dups
        # 3. Exclusion: Low Fold-back ratio (BFB exclusion)
        # 4. Locality: Low inter-chromosomal events

        has_cn_pattern <- (row$gradient_correlation >= 0.4) || (row$cn_variation >= 1.5)
        has_topology <- (row$max_chain_length >= 3) || (row$n_tandem_dups >= 3)
        not_bfb <- row$bfb_risk < 0.3  # Less than 30% fold-backs
        is_local <- row$n_inter_chr_svs < 3

        # Strong evidence: Chain length >= 4 AND Not BFB
        strong_evidence <- (row$max_chain_length >= 4) && not_bfb && is_local

        criteria_met <- sum(c(has_cn_pattern, has_topology, not_bfb, is_local))

        if (strong_evidence || criteria_met == 4) {
            classifications[i] <- "Likely chromoanasynthesis"
        } else if (criteria_met >= 3) {
            classifications[i] <- "Possible chromoanasynthesis"
        } else if (criteria_met >= 2) {
            classifications[i] <- "Unlikely chromoanasynthesis"
        } else {
            classifications[i] <- "Not chromoanasynthesis"
        }
    }

    return(classifications)
}


#' Create empty chromoanasynthesis result
#'
#' @return Empty result structure
#' @keywords internal
create_empty_chromoanasynthesis_result <- function() {
    result <- list(
        regions = list(),
        region_details = list(),
        summary = data.frame(
            region_id = integer(0),
            chrom = character(0),
            start = numeric(0),
            end = numeric(0),
            size = numeric(0),
            n_cn_segments = integer(0),
            gradient_correlation = numeric(0),
            cn_variation = numeric(0),
            n_tandem_dups = integer(0),
            n_fold_backs = integer(0),
            n_linked_chains = integer(0),
            max_chain_length = integer(0),
            n_inter_chr_svs = integer(0),
            complexity_score = numeric(0),
            bfb_risk = numeric(0),
            classification = character(0)
        ),
        total_regions = 0,
        likely_chromoanasynthesis = 0,
        possible_chromoanasynthesis = 0
    )

    class(result) <- c("chromoanasynthesis", "list")
    return(result)
}


#' Print method for chromoanasynthesis results
#'
#' @param x Chromoanasynthesis result object
#' @param ... Additional arguments
#' @export
print.chromoanasynthesis <- function(x, ...) {
    cat("\n")
    cat(rep("=", 70), "\n", sep = "")
    cat("         CHROMOANASYNTHESIS DETECTION RESULTS\n")
    cat(rep("=", 70), "\n\n", sep = "")

    cat(sprintf("Total regions detected: %d\n", x$total_regions))
    cat(sprintf("  - Likely chromoanasynthesis:   %d\n", x$likely_chromoanasynthesis))
    cat(sprintf("  - Possible chromoanasynthesis: %d\n", x$possible_chromoanasynthesis))
    cat("\n")

    if (x$total_regions > 0) {
        cat("Region Summary:\n")
        cat(rep("-", 70), "\n", sep = "")
        # Check if columns exist (for backward compatibility if run on old objects)
        cols <- c("region_id", "chrom", "gradient_correlation", "max_chain_length", "bfb_risk", "classification")
        cols <- intersect(cols, colnames(x$summary))
        print(x$summary[, cols])
    } else {
        cat("No chromoanasynthesis regions detected.\n")
    }

    cat("\n")
    cat(rep("=", 70), "\n", sep = "")
    cat("\n")

    invisible(x)
}


#' Summary method for chromoanasynthesis results
#'
#' @param object Chromoanasynthesis result object
#' @param ... Additional arguments
#' @export
summary.chromoanasynthesis <- function(object, ...) {
    print(object)
    invisible(object)
}
