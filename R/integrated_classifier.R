#' Integrated classifier for mixed chromoanagenesis mechanisms
#'
#' Analyzes chromoanagenesis results to identify and classify mixed mechanisms,
#' where multiple catastrophic events co-occur in the same sample or region.
#'
#' @param chromoanagenesis_result Result from detect_chromoanagenesis()
#' @param overlap_threshold Minimum overlap (in bp) to consider mechanisms co-occurring (default: 1e6)
#' @param min_confidence Minimum confidence score to include events (default: 0.3)
#' @return A list containing integrated classification results
#' @details
#' This function identifies complex patterns where multiple chromoanagenesis
#' mechanisms occur together, such as:
#' - Chromothripsis + Chromoplexy (catastrophic + chained rearrangements)
#' - Chromothripsis + Chromoanasynthesis (shattering + serial replication)
#' - Triple mechanism events (all three mechanisms)
#'
#' The classifier evaluates:
#' 1. Spatial overlap of different mechanisms
#' 2. Temporal relationships (primary vs secondary events)
#' 3. Mechanism dominance in each chromosome
#' 4. Overall sample complexity
#'
#' @export
classify_mixed_mechanisms <- function(chromoanagenesis_result,
                                     overlap_threshold = 1e6,
                                     min_confidence = 0.3) {

    if (!inherits(chromoanagenesis_result, "chromoanagenesis")) {
        stop("Input must be a chromoanagenesis result object from detect_chromoanagenesis()")
    }

    cat("\n")
    cat(rep("=", 70), "\n", sep = "")
    cat("     INTEGRATED MECHANISM CLASSIFICATION\n")
    cat(rep("=", 70), "\n\n", sep = "")

    results <- list()

    # 1. Extract mechanism locations
    cat("Step 1: Extracting mechanism locations...\n")
    mechanism_locations <- extract_mechanism_locations(
        chromoanagenesis_result,
        min_confidence
    )
    results$mechanism_locations <- mechanism_locations

    # 2. Detect overlapping mechanisms
    cat("Step 2: Detecting mechanism overlaps (optimized)...")
    overlaps <- detect_mechanism_overlaps(
        mechanism_locations,
        overlap_threshold
    )
    results$overlaps <- overlaps

    # 3. Chromosome-level classification
    cat("Step 3: Classifying mechanisms by chromosome...\n")
    chr_classification <- classify_by_chromosome(
        mechanism_locations,
        overlaps
    )
    results$chromosome_classification <- chr_classification

    # 4. Sample-level classification
    cat("Step 4: Sample-level integrated classification...\n")
    sample_classification <- classify_sample_level(
        chromoanagenesis_result,
        chr_classification,
        overlaps
    )
    results$sample_classification <- sample_classification

    # 5. Mechanism dominance analysis
    cat("Step 5: Analyzing mechanism dominance...\n")
    dominance <- analyze_mechanism_dominance(
        chromoanagenesis_result,
        chr_classification
    )
    results$dominance <- dominance

    # 6. Complexity scoring
    cat("Step 6: Calculating complexity scores...\n")
    complexity <- calculate_complexity_score(
        chromoanagenesis_result,
        overlaps,
        chr_classification
    )
    results$complexity <- complexity

    # 7. Global FDR Correction for P-values across mechanisms
    cat("Step 7: Applying global FDR correction for p-values...\n")
    fdr_results <- apply_global_fdr(chromoanagenesis_result)
    results$global_fdr <- fdr_results

    cat("\n")
    cat(rep("=", 70), "\n", sep = "")
    cat("CLASSIFICATION COMPLETE\n")
    cat(rep("=", 70), "\n\n", sep = "")

    class(results) <- c("mixed_mechanisms", "list")
    return(results)
}

#' Apply global FDR correction to p-values from different mechanisms
#'
#' This function collects p-values from chromothripsis and chromoplexy events,
#' applies a global Benjamini-Hochberg False Discovery Rate (FDR) correction,
#' and returns a list of original and adjusted p-values.
#' Chromoanasynthesis is excluded as it currently does not produce p-values.
#'
#' @param chromoanagenesis_result Result object from detect_chromoanagenesis()
#' @return A list containing original and adjusted p-values per event.
#' @keywords internal
apply_global_fdr <- function(chromoanagenesis_result) {
    all_pvalues <- list()
    pvalue_names <- character(0)
    original_pvalues_vec <- numeric(0)

    # --- 1. Extract Chromothripsis P-values ---
    if (!is.null(chromoanagenesis_result$chromothripsis) &&
        !is.null(chromoanagenesis_result$chromothripsis$classification) &&
        nrow(chromoanagenesis_result$chromothripsis$classification) > 0) {

        ct_summary <- chromoanagenesis_result$chromothripsis$classification
        ct_pvals <- numeric(0)
        ct_ids <- character(0)

        # For chromothripsis, we have multiple p-values per chromosome.
        # A pragmatic approach is to take the minimum p-value for each event (chromosome).
        # Relevant p-value columns from function_criteria.R are:
        # pval_fragment_joins, pval_exp_chr, pval_exp_cluster, chr_breakpoint_enrichment
        p_cols <- c("pval_fragment_joins", "pval_exp_chr", "pval_exp_cluster", "chr_breakpoint_enrichment")
        
        # Filter for existing p-value columns and ensure they are numeric
        existing_p_cols <- intersect(p_cols, colnames(ct_summary))
        if (length(existing_p_cols) > 0) {
            # Take the minimum p-value across relevant metrics for each chromosome event
            ct_pvals_matrix <- as.matrix(ct_summary[, existing_p_cols])
            # Handle NAs - treat as 1 (not significant) for min comparison, but keep original for reporting
            ct_pvals_matrix[is.na(ct_pvals_matrix)] <- 1
            ct_pvals <- apply(ct_pvals_matrix, 1, min)
            
            ct_ids <- paste0("CT_", ct_summary$chrom)
            
            original_pvalues_vec <- c(original_pvalues_vec, ct_pvals)
            pvalue_names <- c(pvalue_names, ct_ids)
            
            all_pvalues$chromothripsis <- data.frame(
                event_id = ct_ids,
                original_pvalue = ct_pvals,
                stringsAsFactors = FALSE
            )
        }
    }

    # --- 2. Extract Chromoplexy P-values ---
    if (!is.null(chromoanagenesis_result$chromoplexy) &&
        !is.null(chromoanagenesis_result$chromoplexy$summary) &&
        nrow(chromoanagenesis_result$chromoplexy$summary) > 0) {

        cp_summary <- chromoanagenesis_result$chromoplexy$summary
        if ("pvalue" %in% colnames(cp_summary)) {
            cp_pvals <- cp_summary$pvalue
            cp_ids <- paste0("CP_", cp_summary$chain_id)

            original_pvalues_vec <- c(original_pvalues_vec, cp_pvals)
            pvalue_names <- c(pvalue_names, cp_ids)

            all_pvalues$chromoplexy <- data.frame(
                event_id = cp_ids,
                original_pvalue = cp_pvals,
                stringsAsFactors = FALSE
            )
        }
    }

    # --- 3. Apply Global FDR Correction ---
    if (length(original_pvalues_vec) > 0) {
        # Remove NAs from p-values, if any (e.g., if a specific p-value was not calculated)
        valid_p_idx <- !is.na(original_pvalues_vec)
        valid_pvalues <- original_pvalues_vec[valid_p_idx]
        valid_p_names <- pvalue_names[valid_p_idx]

        if (length(valid_pvalues) > 0) {
            adjusted_p <- p.adjust(valid_pvalues, method = "BH")
            is_significant_fdr <- adjusted_p < 0.05
            
            # Create a lookup table for adjusted p-values and significance
            fdr_lookup <- data.frame(
                event_id = valid_p_names,
                adjusted_pvalue = adjusted_p,
                is_significant_fdr = is_significant_fdr,
                stringsAsFactors = FALSE
            )
            
            # Update chromothripsis results
            if (!is.null(all_pvalues$chromothripsis)) {
                all_pvalues$chromothripsis <- merge(all_pvalues$chromothripsis, fdr_lookup, by = "event_id", all.x = TRUE)
            }
            
            # Update chromoplexy results
            if (!is.null(all_pvalues$chromoplexy)) {
                all_pvalues$chromoplexy <- merge(all_pvalues$chromoplexy, fdr_lookup, by = "event_id", all.x = TRUE)
            }
        }
    } else {
        cat("  No p-values found from chromothripsis or chromoplexy for FDR correction.\n")
    }

    return(all_pvalues)
}


#' Extract genomic locations of each mechanism
#'
#' @param chromoanagenesis_result Chromoanagenesis result object
#' @param min_confidence Minimum confidence threshold
#' @return Data frame with mechanism locations
#' @keywords internal
extract_mechanism_locations <- function(chromoanagenesis_result, min_confidence) {

    # Initialize empty data frames for each mechanism
    df_ct <- NULL
    df_cp <- NULL
    df_cs <- NULL

    # --- 1. Chromothripsis Locations ---
    if (!is.null(chromoanagenesis_result$chromothripsis)) {
        ct_class <- chromoanagenesis_result$chromothripsis$classification
        
        cat(sprintf("  DEBUG: ct_class has %d rows\n", if(is.null(ct_class)) 0 else nrow(ct_class)))
        
        # Access the raw summary for coordinates if available
        # The detection_output usually contains 'chromSummary'
        ct_summary <- NULL
        if (!is.null(chromoanagenesis_result$chromothripsis$detection_output) && 
            !is.null(chromoanagenesis_result$chromothripsis$detection_output@chromSummary)) {
            ct_summary <- chromoanagenesis_result$chromothripsis$detection_output@chromSummary
        }

        if (!is.null(ct_class) && nrow(ct_class) > 0) {
            # Filter
            keep_idx <- which(
                !is.na(ct_class$classification) & 
                ct_class$classification %in% c("High confidence", "Low confidence") &
                (if ("confidence_score" %in% colnames(ct_class)) ct_class$confidence_score >= min_confidence else TRUE)
            )
            
            cat(sprintf("  DEBUG: CT keep_idx has %d elements\n", length(keep_idx)))
            
            if (length(keep_idx) > 0) {
                sel_ct <- ct_class[keep_idx, ]
                
                # Get start/end from summary if possible, matching by chromosome
                starts <- rep(1, length(keep_idx)) # Default to 1
                ends <- rep(250e6, length(keep_idx)) # Default to large number
                
                if (!is.null(ct_summary)) {
                    # Create a lookup
                    # Assuming chromSummary has 'chrom', 'start', 'end'
                    if (all(c("chrom", "start", "end") %in% colnames(ct_summary))) {
                        # Match rows
                        m <- match(sel_ct$chrom, ct_summary$chrom)
                        valid_m <- !is.na(m)
                        if (any(valid_m)) {
                            starts[valid_m] <- ct_summary$start[m[valid_m]]
                            ends[valid_m] <- ct_summary$end[m[valid_m]]
                        }
                    }
                }

                cat(sprintf("  DEBUG: sel_ct=%d, starts=%d, ends=%d\n", nrow(sel_ct), length(starts), length(ends)))

                df_ct <- data.frame(
                    mechanism = rep("chromothripsis", nrow(sel_ct)),
                    chrom = as.character(sel_ct$chrom),
                    start = as.numeric(starts),
                    end = as.numeric(ends),
                    confidence = if ("confidence_score" %in% colnames(sel_ct)) as.numeric(sel_ct$confidence_score) else rep(0.8, nrow(sel_ct)),
                    classification = as.character(sel_ct$classification),
                    event_id = paste0("CT_", sel_ct$chrom),
                    stringsAsFactors = FALSE
                )
            }
        }
    }

    # --- 2. Chromoplexy Locations ---
    if (!is.null(chromoanagenesis_result$chromoplexy)) {
        cp_summary <- chromoanagenesis_result$chromoplexy$summary
        cp_chains <- chromoanagenesis_result$chromoplexy$chains # Need this for coordinates

        if (!is.null(cp_summary) && nrow(cp_summary) > 0) {
            keep_idx <- which(
                cp_summary$classification %in% c("Likely chromoplexy", "Possible chromoplexy")
            )
            
            # Filter confidence if score exists
            if ("confidence_score" %in% colnames(cp_summary)) {
                keep_idx <- intersect(keep_idx, which(cp_summary$confidence_score >= min_confidence))
            } else if ("complexity_score" %in% colnames(cp_summary)) {
                # Map complexity to confidence roughly
                keep_idx <- intersect(keep_idx, which(cp_summary$complexity_score >= min_confidence))
            }

            if (length(keep_idx) > 0) {
                # We need to iterate here because one chain spans multiple chromosomes
                # and we want to record the region on EACH chromosome
                
                cp_list <- list()
                
                for (idx in keep_idx) {
                    chain_id <- cp_summary$chain_id[idx]
                    chain_data <- cp_chains[[chain_id]] # Assuming chain_id matches index or we search
                    
                    # If chain_id is not index, find it
                    if (is.null(chain_data) || chain_data$id != chain_id) {
                        # Search by ID
                        found <- FALSE
                        for (ch in cp_chains) {
                            if (ch$id == chain_id) {
                                chain_data <- ch
                                found <- TRUE
                                break
                            }
                        }
                        if (!found) next
                    }
                    
                    # Get nodes (chr:pos)
                    nodes <- chain_data$nodes
                    if (is.null(nodes)) next
                    
                    # Parse nodes to get chrom and pos
                    parts <- strsplit(nodes, ":")
                    chrs <- sapply(parts, `[`, 1)
                    pos <- as.numeric(sapply(parts, `[`, 2))
                    
                    # For each unique chromosome in the chain
                    uniq_chrs <- unique(chrs)
                    
                    for (uc in uniq_chrs) {
                        # Get range on this chromosome
                        chr_pos <- pos[chrs == uc]
                        
                        # Add a buffer? SVs are points, but the event affects the region between them (if >1)
                        # or a small region around them.
                        # Let's take min/max of the breakpoints on this chromosome
                        r_start <- min(chr_pos)
                        r_end <- max(chr_pos)
                        
                        # If single point, add buffer
                        if (r_start == r_end) {
                            r_start <- r_start - 10000
                            r_end <- r_end + 10000
                        }
                        
                        cp_list[[length(cp_list) + 1]] <- data.frame(
                            mechanism = "chromoplexy",
                            chrom = as.character(uc),
                            start = r_start,
                            end = r_end,
                            confidence = if ("confidence_score" %in% colnames(cp_summary)) cp_summary$confidence_score[idx] else 0.7,
                            classification = cp_summary$classification[idx],
                            event_id = paste0("CP_", chain_id),
                            stringsAsFactors = FALSE
                        )
                    }
                }
                if (length(cp_list) > 0) df_cp <- do.call(rbind, cp_list)
            }
        }
    }

    # --- 3. Chromoanasynthesis Locations ---
    if (!is.null(chromoanagenesis_result$chromoanasynthesis)) {
        cs_summary <- chromoanagenesis_result$chromoanasynthesis$summary
        
        if (!is.null(cs_summary) && nrow(cs_summary) > 0) {
            keep_idx <- which(
                !is.na(cs_summary$classification) &
                cs_summary$classification %in% c("Likely chromoanasynthesis", "Possible chromoanasynthesis")
            )
            
            # Confidence check
            conf_col <- if ("confidence_score" %in% colnames(cs_summary)) "confidence_score" else "complexity_score"
            if (conf_col %in% colnames(cs_summary)) {
                keep_idx <- intersect(keep_idx, which(cs_summary[[conf_col]] >= min_confidence))
            }
            
            if (length(keep_idx) > 0) {
                sel_cs <- cs_summary[keep_idx, ]
                
                # Check column names for start/end (might be region_start/region_end)
                s_col <- if ("start" %in% colnames(sel_cs)) "start" else "region_start"
                e_col <- if ("end" %in% colnames(sel_cs)) "end" else "region_end"
                
                df_cs <- data.frame(
                    mechanism = rep("chromoanasynthesis", nrow(sel_cs)),
                    chrom = as.character(sel_cs$chrom),
                    start = as.numeric(sel_cs[[s_col]]),
                    end = as.numeric(sel_cs[[e_col]]),
                    confidence = if (conf_col %in% colnames(sel_cs)) as.numeric(sel_cs[[conf_col]]) else rep(0.5, nrow(sel_cs)),
                    classification = as.character(sel_cs$classification),
                    event_id = paste0("CS_", sel_cs$region_id),
                    stringsAsFactors = FALSE
                )
            }
        }
    }

    # Combine all
    all_dfs <- list(df_ct, df_cp, df_cs)
    # Remove NULLs
    all_dfs <- all_dfs[!sapply(all_dfs, is.null)]
    
    if (length(all_dfs) == 0) {
        return(data.frame(
            mechanism = character(0),
            chrom = character(0),
            start = numeric(0),
            end = numeric(0),
            confidence = numeric(0),
            classification = character(0),
            event_id = character(0),
            stringsAsFactors = FALSE
        ))
    }
    
    final_df <- do.call(rbind, all_dfs)
    return(final_df)
}


#' Detect overlapping mechanisms using GenomicRanges
#'
#' optimized version using interval trees
#'
#' @param mechanism_locations Data frame from extract_mechanism_locations
#' @param overlap_threshold Minimum overlap in bp
#' @return List of overlapping mechanism pairs
#' @keywords internal
detect_mechanism_overlaps <- function(mechanism_locations, overlap_threshold) {

    # Early exit if insufficient data
    if (nrow(mechanism_locations) < 2) {
        return(list(n_overlaps = 0, overlaps = data.frame()))
    }
    
    # Check for GenomicRanges
    if (!requireNamespace("GenomicRanges", quietly = TRUE) || 
        !requireNamespace("IRanges", quietly = TRUE)) {
        warning("GenomicRanges/IRanges package not installed. Falling back to slow overlap detection.")
        return(.detect_mechanism_overlaps_slow(mechanism_locations, overlap_threshold))
    }

    # Remove rows with NA coordinates
    valid_locs <- mechanism_locations[!is.na(mechanism_locations$start) & 
                                      !is.na(mechanism_locations$end), ]
    
    if (nrow(valid_locs) < 2) {
        return(list(n_overlaps = 0, overlaps = data.frame()))
    }

    # Construct GRanges object
    gr <- GenomicRanges::GRanges(
        seqnames = valid_locs$chrom,
        ranges = IRanges::IRanges(start = valid_locs$start, end = valid_locs$end),
        mechanism = valid_locs$mechanism,
        event_id = valid_locs$event_id,
        confidence = valid_locs$confidence
    )

    # Find overlaps
    # minoverlap ensures significant spatial overlap
    hits <- GenomicRanges::findOverlaps(gr, gr, 
                                       minoverlap = overlap_threshold)
    
    if (length(hits) == 0) {
        return(list(n_overlaps = 0, overlaps = data.frame()))
    }
    
    # Convert to data frame indices
    q_idx <- S4Vectors::queryHits(hits)
    s_idx <- S4Vectors::subjectHits(hits)
    
    # Filter: 
    # 1. Keep only pairs where query index < subject index to avoid duplicates (A-B vs B-A)
    # 2. Keep only different mechanisms (don't care about overlapping CT-CT, usually merged anyway)
    
    keep <- (q_idx < s_idx) & 
            (valid_locs$mechanism[q_idx] != valid_locs$mechanism[s_idx])
            
    if (sum(keep) == 0) {
        return(list(n_overlaps = 0, overlaps = data.frame()))
    }
    
    q_idx <- q_idx[keep]
    s_idx <- s_idx[keep]
    
    # Calculate overlap size
    # GenomicRanges::pintersect gives intersection ranges
    intersections <- GenomicRanges::pintersect(gr[q_idx], gr[s_idx])
    overlap_sizes <- GenomicRanges::width(intersections)
    
    # Build result data frame
    overlaps_df <- data.frame(
        chrom = as.character(GenomicRanges::seqnames(gr)[q_idx]),
        mechanism1 = valid_locs$mechanism[q_idx],
        mechanism2 = valid_locs$mechanism[s_idx],
        event_id1 = valid_locs$event_id[q_idx],
        event_id2 = valid_locs$event_id[s_idx],
        overlap_size = overlap_sizes,
        confidence1 = valid_locs$confidence[q_idx],
        confidence2 = valid_locs$confidence[s_idx],
        stringsAsFactors = FALSE
    )
    
    # Add mechanism pair label
    overlaps_df$mechanism_pair <- apply(overlaps_df[, c("mechanism1", "mechanism2")], 1, function(x) {
        paste(sort(x), collapse = "+")
    })

    return(list(
        n_overlaps = nrow(overlaps_df),
        overlaps = overlaps_df
    ))
}

#' Fallback slow overlap detection (only if GenomicRanges missing)
#' @keywords internal
.detect_mechanism_overlaps_slow <- function(mechanism_locations, overlap_threshold) {
    # This is the original logic, kept as fallback
    if (nrow(mechanism_locations) == 0) {
        return(list(n_overlaps = 0, overlaps = data.frame()))
    }

    overlaps <- list()
    chroms <- unique(mechanism_locations$chrom)

    for (chr in chroms) {
        chr_mechs <- mechanism_locations[mechanism_locations$chrom == chr, ]
        if (nrow(chr_mechs) > 1) {
            for (i in 1:(nrow(chr_mechs) - 1)) {
                for (j in (i + 1):nrow(chr_mechs)) {
                    mech1 <- chr_mechs[i, ]
                    mech2 <- chr_mechs[j, ]
                    if (mech1$mechanism == mech2$mechanism) next
                    
                    overlap_detected <- FALSE
                    overlap_size <- NA
                    if (!is.na(mech1$start) && !is.na(mech2$start)) {
                        o_start <- max(mech1$start, mech2$start)
                        o_end <- min(mech1$end, mech2$end)
                        if (o_end >= o_start) {
                            sz <- o_end - o_start
                            if (sz >= overlap_threshold) {
                                overlap_detected <- TRUE
                                overlap_size <- sz
                            }
                        }
                    }
                    
                    if (overlap_detected) {
                        overlaps[[length(overlaps) + 1]] <- list(
                            chrom = chr,
                            mechanism1 = mech1$mechanism,
                            mechanism2 = mech2$mechanism,
                            event_id1 = mech1$event_id,
                            event_id2 = mech2$event_id,
                            overlap_size = overlap_size,
                            confidence1 = mech1$confidence,
                            confidence2 = mech2$confidence
                        )
                    }
                }
            }
        }
    }
    
    if (length(overlaps) == 0) return(list(n_overlaps=0, overlaps=data.frame()))
    
    overlaps_df <- do.call(rbind, lapply(overlaps, as.data.frame))
    overlaps_df$mechanism_pair <- paste(
        pmin(overlaps_df$mechanism1, overlaps_df$mechanism2),
        pmax(overlaps_df$mechanism1, overlaps_df$mechanism2),
        sep="+ "
    )
    
    return(list(n_overlaps = nrow(overlaps_df), overlaps = overlaps_df))
}


#' Classify mechanisms by chromosome
#'
#' @param mechanism_locations Data frame from extract_mechanism_locations
#' @param overlaps List from detect_mechanism_overlaps
#' @return Data frame with chromosome-level classification
#' @keywords internal
classify_by_chromosome <- function(mechanism_locations, overlaps) {

    if (nrow(mechanism_locations) == 0) {
        return(data.frame(
            chrom = character(0),
            mechanisms = character(0),
            dominant_mechanism = character(0),
            is_mixed = logical(0),
            n_mechanisms = numeric(0),
            stringsAsFactors = FALSE
        ))
    }

    chroms <- unique(mechanism_locations$chrom)
    chr_class <- list()

    for (chr in chroms) {
        chr_mechs <- mechanism_locations[mechanism_locations$chrom == chr, ]

        # Get unique mechanisms on this chromosome
        unique_mechs <- unique(chr_mechs$mechanism)
        n_mechs <- length(unique_mechs)
        is_mixed <- n_mechs > 1

        # Determine dominant mechanism (highest confidence)
        chr_mechs_agg <- aggregate(confidence ~ mechanism, data = chr_mechs, FUN = max)
        dominant_idx <- which.max(chr_mechs_agg$confidence)
        dominant_mech <- chr_mechs_agg$mechanism[dominant_idx]

        # Create classification label
        if (is_mixed) {
            mech_label <- paste(sort(unique_mechs), collapse = "+")
        } else {
            mech_label <- unique_mechs[1]
        }

        chr_class[[chr]] <- list(
            chrom = chr,
            mechanisms = mech_label,
            dominant_mechanism = dominant_mech,
            is_mixed = is_mixed,
            n_mechanisms = n_mechs,
            max_confidence = max(chr_mechs$confidence)
        )
    }

    # Convert to data frame
    chr_class_df <- do.call(rbind, lapply(chr_class, function(x) {
        data.frame(
            chrom = as.character(x$chrom),
            mechanisms = as.character(x$mechanisms),
            dominant_mechanism = as.character(x$dominant_mechanism),
            is_mixed = as.logical(x$is_mixed),
            n_mechanisms = as.numeric(x$n_mechanisms),
            max_confidence = as.numeric(x$max_confidence),
            stringsAsFactors = FALSE
        )
    }))

    # Sort by chromosome
    chr_class_df <- chr_class_df[order(chr_class_df$chrom), ]

    return(chr_class_df)
}


#' Sample-level integrated classification
#'
#' @param chromoanagenesis_result Chromoanagenesis result object
#' @param chr_classification Chromosome classification data frame
#' @param overlaps Overlap detection results
#' @return List with sample-level classification
#' @keywords internal
classify_sample_level <- function(chromoanagenesis_result, chr_classification, overlaps) {

    # Count mechanisms
    # Chromothripsis: check both high and low confidence
    has_chromothripsis <- FALSE
    if (!is.null(chromoanagenesis_result$chromothripsis)) {
        n_high <- chromoanagenesis_result$chromothripsis$n_high_confidence
        n_low <- chromoanagenesis_result$chromothripsis$n_low_confidence
        has_chromothripsis <- (!is.na(n_high) && n_high > 0) ||
                             (!is.na(n_low) && n_low > 0)
    }

    # Chromoplexy: check likely events
    has_chromoplexy <- FALSE
    if (!is.null(chromoanagenesis_result$chromoplexy)) {
        n_likely <- chromoanagenesis_result$chromoplexy$likely_chromoplexy
        has_chromoplexy <- !is.na(n_likely) && n_likely > 0
    }

    # Chromoanasynthesis: check likely events
    has_chromoanasynthesis <- FALSE
    if (!is.null(chromoanagenesis_result$chromoanasynthesis)) {
        n_likely <- chromoanagenesis_result$chromoanasynthesis$likely_chromoanasynthesis
        has_chromoanasynthesis <- !is.na(n_likely) && n_likely > 0
    }

    mechanisms_present <- c()
    if (has_chromothripsis) mechanisms_present <- c(mechanisms_present, "chromothripsis")
    if (has_chromoplexy) mechanisms_present <- c(mechanisms_present, "chromoplexy")
    if (has_chromoanasynthesis) mechanisms_present <- c(mechanisms_present, "chromoanasynthesis")

    n_mechanisms <- length(mechanisms_present)

    # Determine classification
    if (n_mechanisms == 0) {
        classification <- "No chromoanagenesis"
        category <- "normal"
    } else if (n_mechanisms == 1) {
        classification <- paste0("Pure ", mechanisms_present[1])
        category <- "single_mechanism"
    } else {
        # Mixed mechanism
        has_overlaps <- overlaps$n_overlaps > 0

        if (has_overlaps) {
            classification <- paste0("Mixed mechanisms with spatial overlap (",
                                   paste(mechanisms_present, collapse = "+"), ")")
            category <- "mixed_overlapping"
        } else {
            classification <- paste0("Multiple independent mechanisms (",
                                   paste(mechanisms_present, collapse = "+"), ")")
            category <- "mixed_independent"
        }
    }

    # Count mixed chromosomes
    n_mixed_chromosomes <- sum(chr_classification$is_mixed, na.rm = TRUE)

    return(list(
        classification = classification,
        category = category,
        n_mechanisms = n_mechanisms,
        mechanisms_present = mechanisms_present,
        n_mixed_chromosomes = n_mixed_chromosomes,
        has_spatial_overlap = overlaps$n_overlaps > 0,
        n_overlaps = overlaps$n_overlaps
    ))
}


#' Analyze mechanism dominance
#'
#' @param chromoanagenesis_result Chromoanagenesis result object
#' @param chr_classification Chromosome classification data frame
#' @return List with dominance analysis
#' @keywords internal
analyze_mechanism_dominance <- function(chromoanagenesis_result, chr_classification) {

    # Count events by mechanism
    n_chromothripsis <- 0
    n_chromoplexy <- 0
    n_chromoanasynthesis <- 0

    # Chromothripsis: count high + low confidence
    if (!is.null(chromoanagenesis_result$chromothripsis)) {
        n_high <- chromoanagenesis_result$chromothripsis$n_high_confidence
        n_low <- chromoanagenesis_result$chromothripsis$n_low_confidence
        n_chromothripsis <- sum(c(n_high, n_low), na.rm = TRUE)
    }

    # Chromoplexy: count likely + possible
    if (!is.null(chromoanagenesis_result$chromoplexy)) {
        n_likely <- chromoanagenesis_result$chromoplexy$likely_chromoplexy
        n_possible <- chromoanagenesis_result$chromoplexy$possible_chromoplexy
        n_chromoplexy <- sum(c(n_likely, n_possible), na.rm = TRUE)
    }

    # Chromoanasynthesis: count likely + possible
    if (!is.null(chromoanagenesis_result$chromoanasynthesis)) {
        n_likely <- chromoanagenesis_result$chromoanasynthesis$likely_chromoanasynthesis
        n_possible <- chromoanagenesis_result$chromoanasynthesis$possible_chromoanasynthesis
        n_chromoanasynthesis <- sum(c(n_likely, n_possible), na.rm = TRUE)
    }

    # Calculate proportions
    total_events <- n_chromothripsis + n_chromoplexy + n_chromoanasynthesis

    if (total_events == 0) {
        return(list(
            dominant_mechanism = "none",
            mechanism_proportions = data.frame(
                mechanism = c("chromothripsis", "chromoplexy", "chromoanasynthesis"),
                n_events = c(0, 0, 0),
                proportion = c(0, 0, 0)
            )
        ))
    }

    proportions <- data.frame(
        mechanism = c("chromothripsis", "chromoplexy", "chromoanasynthesis"),
        n_events = c(n_chromothripsis, n_chromoplexy, n_chromoanasynthesis),
        proportion = c(n_chromothripsis, n_chromoplexy, n_chromoanasynthesis) / total_events,
        stringsAsFactors = FALSE
    )

    # Determine dominant mechanism
    dominant_idx <- which.max(proportions$n_events)
    dominant_mechanism <- proportions$mechanism[dominant_idx]

    # Check if truly dominant (>50%)
    if (proportions$proportion[dominant_idx] <= 0.5) {
        dominant_mechanism <- "balanced"
    }

    return(list(
        dominant_mechanism = dominant_mechanism,
        mechanism_proportions = proportions,
        total_events = total_events
    ))
}


#' Calculate sample complexity score
#'
#' @param chromoanagenesis_result Chromoanagenesis result object
#' @param overlaps Overlap detection results
#' @param chr_classification Chromosome classification data frame
#' @return List with complexity metrics
#' @keywords internal
calculate_complexity_score <- function(chromoanagenesis_result, overlaps, chr_classification) {

    # Components of complexity:
    # 1. Number of different mechanisms (0-3)
    # 2. Number of chromosomes affected
    # 3. Presence of overlapping mechanisms
    # 4. Total number of events

    # Get basic counts
    # Count mechanisms present (with NA-safe checks)
    has_chromothripsis <- FALSE
    if (!is.null(chromoanagenesis_result$chromothripsis)) {
        n_high <- chromoanagenesis_result$chromothripsis$n_high_confidence
        n_low <- chromoanagenesis_result$chromothripsis$n_low_confidence
        has_chromothripsis <- (!is.na(n_high) && n_high > 0) ||
                             (!is.na(n_low) && n_low > 0)
    }

    has_chromoplexy <- FALSE
    if (!is.null(chromoanagenesis_result$chromoplexy)) {
        n_likely <- chromoanagenesis_result$chromoplexy$likely_chromoplexy
        has_chromoplexy <- !is.na(n_likely) && n_likely > 0
    }

    has_chromoanasynthesis <- FALSE
    if (!is.null(chromoanagenesis_result$chromoanasynthesis)) {
        n_likely <- chromoanagenesis_result$chromoanasynthesis$likely_chromoanasynthesis
        has_chromoanasynthesis <- !is.na(n_likely) && n_likely > 0
    }

    n_mechanisms <- sum(c(has_chromothripsis, has_chromoplexy, has_chromoanasynthesis))

    n_chromosomes <- nrow(chr_classification)
    n_mixed_chromosomes <- sum(chr_classification$is_mixed, na.rm = TRUE)
    n_overlaps <- overlaps$n_overlaps

    # Calculate total events
    total_events <- 0
    if (!is.null(chromoanagenesis_result$chromothripsis)) {
        n_high <- chromoanagenesis_result$chromothripsis$n_high_confidence
        n_low <- chromoanagenesis_result$chromothripsis$n_low_confidence
        total_events <- total_events + sum(c(n_high, n_low), na.rm = TRUE)
    }
    if (!is.null(chromoanagenesis_result$chromoplexy)) {
        n_likely <- chromoanagenesis_result$chromoplexy$likely_chromoplexy
        n_possible <- chromoanagenesis_result$chromoplexy$possible_chromoplexy
        total_events <- total_events + sum(c(n_likely, n_possible), na.rm = TRUE)
    }
    if (!is.null(chromoanagenesis_result$chromoanasynthesis)) {
        n_likely <- chromoanagenesis_result$chromoanasynthesis$likely_chromoanasynthesis
        n_possible <- chromoanagenesis_result$chromoanasynthesis$possible_chromoanasynthesis
        total_events <- total_events + sum(c(n_likely, n_possible), na.rm = TRUE)
    }

    # Complexity score (0-1 scale)
    # Weighted components:
    # - Mechanism diversity: 30%
    # - Spatial overlap: 25%
    # - Chromosome spread: 25%
    # - Event count: 20%

    mechanism_score <- n_mechanisms / 3  # Max 3 mechanisms
    overlap_score <- min(n_overlaps / 5, 1)  # Saturate at 5 overlaps
    chromosome_score <- min(n_chromosomes / 10, 1)  # Saturate at 10 chromosomes
    event_score <- min(total_events / 10, 1)  # Saturate at 10 events

    complexity_score <- (
        mechanism_score * 0.30 +
        overlap_score * 0.25 +
        chromosome_score * 0.25 +
        event_score * 0.20
    )

    # Classify complexity level
    if (complexity_score < 0.3) {
        complexity_level <- "Low"
    } else if (complexity_score < 0.6) {
        complexity_level <- "Moderate"
    } else if (complexity_score < 0.8) {
        complexity_level <- "High"
    } else {
        complexity_level <- "Very High"
    }

    return(list(
        complexity_score = complexity_score,
        complexity_level = complexity_level,
        n_mechanisms = n_mechanisms,
        n_chromosomes = n_chromosomes,
        n_mixed_chromosomes = n_mixed_chromosomes,
        n_overlaps = n_overlaps,
        total_events = total_events,
        components = data.frame(
            component = c("Mechanism diversity", "Spatial overlap",
                         "Chromosome spread", "Event count"),
            raw_score = c(mechanism_score, overlap_score,
                         chromosome_score, event_score),
            weight = c(0.30, 0.25, 0.25, 0.20),
            weighted_score = c(mechanism_score * 0.30, overlap_score * 0.25,
                              chromosome_score * 0.25, event_score * 0.20),
            stringsAsFactors = FALSE
        )
    ))
}


#' Print method for mixed mechanisms classification
#'
#' @param x Mixed mechanisms result object
#' @param ... Additional arguments
#' @export
print.mixed_mechanisms <- function(x, ...) {
    cat("\n")
    cat(rep("=", 70), "\n", sep = "")
    cat("         INTEGRATED MECHANISM CLASSIFICATION\n")
    cat(rep("=", 70), "\n\n", sep = "")

    # Sample-level classification
    cat("SAMPLE-LEVEL CLASSIFICATION:\n")
    cat(sprintf("  Classification: %s\n", x$sample_classification$classification))
    cat(sprintf("  Category: %s\n", x$sample_classification$category))
    cat(sprintf("  Number of mechanisms: %d\n", x$sample_classification$n_mechanisms))
    if (length(x$sample_classification$mechanisms_present) > 0) {
        cat(sprintf("  Mechanisms present: %s\n",
                   paste(x$sample_classification$mechanisms_present, collapse = ", ")))
    }
    cat("\n")

    # Complexity
    cat("COMPLEXITY ANALYSIS:\n")
    cat(sprintf("  Overall complexity: %s (score: %.3f)\n",
               x$complexity$complexity_level,
               x$complexity$complexity_score))
    cat(sprintf("  Total events: %d\n", x$complexity$total_events))
    cat(sprintf("  Chromosomes affected: %d\n", x$complexity$n_chromosomes))
    cat(sprintf("  Mixed chromosomes: %d\n", x$complexity$n_mixed_chromosomes))
    cat(sprintf("  Spatial overlaps: %d\n", x$complexity$n_overlaps))
    cat("\n")

    # Dominance
    cat("MECHANISM DOMINANCE:\n")
    cat(sprintf("  Dominant mechanism: %s\n", x$dominance$dominant_mechanism))
    cat("  Event distribution:\n")
    for (i in 1:nrow(x$dominance$mechanism_proportions)) {
        mech <- x$dominance$mechanism_proportions$mechanism[i]
        n <- x$dominance$mechanism_proportions$n_events[i]
        prop <- x$dominance$mechanism_proportions$proportion[i]
        cat(sprintf("    - %s: %d events (%.1f%%)\n", mech, n, prop * 100))
    }
    cat("\n")

    # Chromosome details
    if (nrow(x$chromosome_classification) > 0) {
        cat("CHROMOSOME-LEVEL DETAILS:\n")
        mixed_chrs <- x$chromosome_classification[x$chromosome_classification$is_mixed, ]
        if (nrow(mixed_chrs) > 0) {
            cat("  Mixed mechanism chromosomes:\n")
            for (i in 1:nrow(mixed_chrs)) {
                cat(sprintf("    - %s: %s (dominant: %s)\n",
                           mixed_chrs$chrom[i],
                           mixed_chrs$mechanisms[i],
                           mixed_chrs$dominant_mechanism[i]))
            }
        } else {
            cat("  No mixed mechanism chromosomes detected\n")
        }
    }

    cat("\n")
    cat(rep("=", 70), "\n", sep = "")
    cat("\n")

    invisible(x)
}


#' Summary method for mixed mechanisms classification
#'
#' @param object Mixed mechanisms result object
#' @param ... Additional arguments
#' @export
summary.mixed_mechanisms <- function(object, ...) {

    cat("\nMixed Mechanisms Summary:\n")
    cat("========================\n\n")

    # Print complexity components
    cat("Complexity Score Components:\n")
    print(object$complexity$components)
    cat("\n")

    # Print mechanism proportions
    cat("Mechanism Proportions:\n")
    print(object$dominance$mechanism_proportions)
    cat("\n")

    # Print chromosome classification
    cat("Chromosome Classification:\n")
    print(object$chromosome_classification)
    cat("\n")

    # Print overlaps if any
    if (object$overlaps$n_overlaps > 0) {
        cat("Mechanism Overlaps:\n")
        print(object$overlaps$overlaps[, c("chrom", "mechanism_pair",
                                          "overlap_size", "confidence1", "confidence2")])
        cat("\n")
    }

    invisible(object)
}