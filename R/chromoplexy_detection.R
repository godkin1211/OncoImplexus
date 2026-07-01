################################################################################
# Chromoplexy Detection - Version 3.0 (FIXED)
#
# CRITICAL FIX: Properly construct graph with THREE types of edges:
# 1. TRANSLOCATION edges (inter-chromosomal SVs)
# 2. DELETION_BRIDGE edges (deletions connecting nearby breakpoints)
# 3. GENOMIC_ADJACENCY edges (nearby breakpoints on same chromosome)
#
# This fixes the fatal flaw in v2.0 where the graph was disconnected.
################################################################################

#' Detect chromoplexy events with corrected graph construction
#'
#' This version fixes the critical bug in v2.0 where the translocation graph
#' was disconnected because it only included translocation edges. The new
#' implementation adds genomic adjacency edges and deletion bridge edges,
#' following the ChainFinder methodology.
#'
#' @param SV.sample An instance of class SVs or data frame with SV data
#' @param CNV.sample Optional: an instance of class CNVsegs or data frame with CNV data
#' @param min_chromosomes Minimum number of chromosomes involved (default: 3)
#' @param min_translocations Minimum number of translocations in chain (default: 3)
#' @param max_cn_change Maximum allowed copy number change (default: 1)
#' @param allow_cycles Allow circular chains (default: TRUE)
#' @param adjacency_distance Maximum distance for genomic adjacency edges (default: 10e6 = 10Mb)
#' @param deletion_bridge_distance Maximum distance to search for deletion bridges (default: 1e6 = 1Mb)
#' @param fdr_threshold FDR threshold for statistical significance (default: 0.01)
#' @param use_statistical_testing Enable statistical significance testing (default: TRUE)
#' @param likely_chromoplexy_threshold Minimum criteria met for "Likely chromoplexy" classification (default: 5)
#' @param possible_chromoplexy_threshold Minimum criteria met for "Possible chromoplexy" classification (default: 4)
#' @param max_path_search Maximum number of paths to search per node (default: 50).
#'   Reduced to prevent timeout on complex graphs.
#' @param max_search_depth Maximum chain length to explore (default: 10).
#'   Prevents exponential explosion in highly connected graphs.
#' @param genome Reference genome for calculating genomic distances (default: "hg19")
#' @return A list containing enhanced chromoplexy detection results
#'
#' @details
#' Key improvements in v3.0:
#' - Constructs complete graph with genomic adjacency edges
#' - Integrates deletion bridges as graph edges (not just annotations)
#' - Properly connects breakpoints on the same chromosome
#' - Enables detection of complex chromoplexy patterns
#'
#' IMPORTANT NOTE on Intra-chromosomal SVs:
#' This function filters to ONLY inter-chromosomal SVs (translocations) for
#' the primary graph construction, following the ChainFinder methodology.
#' Intra-chromosomal SVs (DEL, DUP, INV) are NOT explicitly included as edges.
#' However, their effects are captured through:
#' 1. GENOMIC ADJACENCY edges - connect nearby breakpoints regardless of SV type
#' 2. DELETION BRIDGE edges - deletions between translocation breakpoints
#' 3. CNV analysis - copy number changes are evaluated separately
#'
#' This approach is biologically justified because chromoplexy is characterized
#' by chains of TRANSLOCATIONS with minimal local rearrangement. Local
#' intra-chromosomal events are better captured as adjacency relationships
#' rather than explicit graph edges.
#'
#' @param SV.sample An SVs object
#' @param CNV.sample A CNVsegs object (optional)
#' @param min_chromosomes Minimum number of chromosomes in a chain (default: 3)
#' @param min_translocations Minimum number of translocations in a chain (default: 3)
#' @param max_cn_change Maximum allowed copy number change for stability (default: 1)
#' @param allow_cycles Whether to allow cycles in detected chains (default: TRUE)
#' @param adjacency_distance Maximum distance for genomic adjacency edges (default: 10Mb)
#' @param adjacency_decay_scale Distance scale for genomic adjacency weighting (default: 2e6)
#' @param deletion_bridge_distance Maximum distance for deletion bridges (default: 1Mb)
#' @param fdr_threshold FDR threshold for statistical significance (default: 0.01)
#' @param use_statistical_testing Whether to perform statistical testing (default: TRUE)
#' @param max_path_search Maximum number of paths to explore (default: 50)
#' @param max_neighbors Maximum number of neighbors to consider per node (default: 3)
#' @param max_adjacency_streak Maximum number of consecutive adjacency edges (default: 3)
#' @param genome Reference genome for calculating genomic distances (default: "hg19")
#' @param collapse_chains Collapse redundant chain enumerations into event-level
#'   connected components (default: TRUE)
#' @param collapse_classifications Chain classifications to include in collapsed
#'   events (default: likely chromoplexy only)
#' @param gene_granges Optional GRanges object for collapsed event gene annotation
#' @param breakpoint_padding Padding in bp around breakpoints for gene annotation
#' @param verbose Print progress messages (default: TRUE)
#' @return A list containing detected chromoplexy chains and summary statistics
#' @export
detect_chromoplexy <- function(SV.sample,
                               CNV.sample = NULL,
                               min_chromosomes = 3,
                               min_translocations = 3,
                               max_cn_change = 1,
                               allow_cycles = TRUE,
                               adjacency_distance = 10e6,
                               adjacency_decay_scale = 2e6,
                               deletion_bridge_distance = 1e6,
                               fdr_threshold = 0.01,
                               use_statistical_testing = TRUE,
                               likely_chromoplexy_threshold = 5,
                               possible_chromoplexy_threshold = 4,
                               max_path_search = 50,
                               max_search_depth = 10,
                               max_neighbors = 3,
                               max_adjacency_streak = 3, # New parameter
                               genome = "hg19",
                               collapse_chains = TRUE,
                               collapse_classifications = c("Likely chromoplexy"),
                               gene_granges = NULL,
                               breakpoint_padding = 1000,
                               verbose = TRUE) {
    # Convert to data frame if needed
    if (is(SV.sample, "SVs")) {
        SV.sample <- as(SV.sample, "data.frame")
    }

    if (!is.null(CNV.sample) && is(CNV.sample, "CNVsegs")) {
        CNV.sample <- as(CNV.sample, "data.frame")
    }
    cnv_available <- !is.null(CNV.sample) && !is.null(nrow(CNV.sample)) && nrow(CNV.sample) > 0
    sv_only <- !cnv_available
    if (!cnv_available) {
        CNV.sample <- NULL
    }

    # Extract inter-chromosomal SVs (translocations)
    # NOTE: Intra-chromosomal SVs (DEL, DUP, INV) are intentionally filtered out.
    # Their effects are captured through:
    #   - Genomic adjacency edges (connects nearby breakpoints)
    #   - Deletion bridge edges (deletions between translocation breakpoints)
    #   - CNV analysis (copy number evaluation)
    # This is consistent with ChainFinder methodology where chromoplexy is
    # defined by chains of TRANSLOCATIONS, not local intra-chromosomal events.
    inter_chr_SVs <- SV.sample[SV.sample$chrom1 != SV.sample$chrom2, ]

    if (nrow(inter_chr_SVs) < min_translocations) {
        warning(sprintf(
            "Only %d inter-chromosomal SVs found. Need at least %d for chromoplexy detection.",
            nrow(inter_chr_SVs), min_translocations
        ))
        return(create_empty_chromoplexy_result_v3(
            analysis_mode = if (sv_only) "SV-only chromoplexy" else "SV+CNV chromoplexy",
            limitations = if (sv_only) c(
                "CNV not provided: CN stability and deletion bridge support were not evaluated."
            ) else character(0)
        ))
    }

    if (verbose) {
        cat("========================================================================\n")
        cat("Chromoplexy Detection v3.0 - CORRECTED GRAPH CONSTRUCTION\n")
        if (sv_only) {
            cat("Mode: SV-only (CN stability and deletion bridges disabled)\n")
        }
        cat("========================================================================\n\n")
    }

    # STEP 1: Identify deletion bridges (before graph construction)
    deletion_bridges <- NULL
    if (!is.null(CNV.sample)) {
        if (verbose) cat("Step 1: Identifying deletion bridges...\n")
        deletion_bridges <- identify_deletion_bridges_v3(
            SV.sample = inter_chr_SVs,
            CNV.sample = CNV.sample,
            max_distance = deletion_bridge_distance
        )
        if (verbose) cat(sprintf("  -> Found %d deletion bridges\n\n", length(deletion_bridges)))
    }

    # STEP 2: Build COMPLETE translocation graph
    # This is the critical fix - graph now includes THREE types of edges
    if (verbose) cat("Step 2: Building complete translocation graph...\n")
    tlx_graph <- build_complete_translocation_graph_v3(
        inter_chr_SVs = inter_chr_SVs,
        deletion_bridges = deletion_bridges,
        adjacency_distance = adjacency_distance,
        adjacency_decay_scale = adjacency_decay_scale,
        max_neighbors = max_neighbors # Pass new parameter
    )

    if (verbose) {
        cat(sprintf(
            "  -> Graph has %d nodes and %d edges\n",
            length(tlx_graph$nodes), nrow(tlx_graph$edges)
        ))
        cat(sprintf(
            "     - Translocation edges: %d\n",
            sum(tlx_graph$edges$edge_type == "TRANSLOCATION")
        ))
        if (!is.null(deletion_bridges) && length(deletion_bridges) > 0) {
            cat(sprintf(
                "     - Deletion bridge edges: %d\n",
                sum(tlx_graph$edges$edge_type == "DELETION_BRIDGE")
            ))
        }
        cat(sprintf(
            "     - Genomic adjacency edges: %d\n\n",
            sum(tlx_graph$edges$edge_type == "ADJACENCY")
        ))
    }

    # STEP 3: Detect all possible chains using backtracking
    if (verbose) cat("Step 3: Detecting translocation chains (backtracking)...")
    all_chains <- detect_all_chains_v3(
        tlx_graph = tlx_graph,
        inter_chr_SVs = inter_chr_SVs,
        min_chromosomes = min_chromosomes,
        min_translocations = min_translocations,
        allow_cycles = allow_cycles,
        max_path_search = max_path_search,
        max_search_depth = max_search_depth,
        max_adjacency_streak = max_adjacency_streak # Pass new parameter
    )

    if (length(all_chains) == 0) {
        if (verbose) cat("  -> No chromoplexy chains detected.\n\n")
        return(create_empty_chromoplexy_result_v3(
            analysis_mode = if (sv_only) "SV-only chromoplexy" else "SV+CNV chromoplexy",
            limitations = if (sv_only) c(
                "CNV not provided: CN stability and deletion bridge support were not evaluated."
            ) else character(0)
        ))
    }

    if (verbose) cat(sprintf("  -> Found %d potential chromoplexy chain(s)\n\n", length(all_chains)))

    # STEP 4: Evaluate each chain
    if (verbose) cat("Step 4: Evaluating chains...\n")
    genome_size <- get_genome_size(genome)

    chain_results <- list()
    for (i in 1:length(all_chains)) {
        chain_results[[i]] <- evaluate_chromoplexy_chain_v3(
            chain = all_chains[[i]],
            inter_chr_SVs = inter_chr_SVs,
            CNV.sample = CNV.sample,
            tlx_graph = tlx_graph,
            max_cn_change = max_cn_change,
            use_statistical_testing = use_statistical_testing,
            genome_size = genome_size,
            fdr_threshold = fdr_threshold
        )
    }

    # Create summary
    summary_df <- do.call(rbind, lapply(chain_results, function(x) x$summary))

    # Classify events
    summary_df$classification <- classify_chromoplexy_event_v3(
        summary_df = summary_df,
        use_statistical_testing = use_statistical_testing,
        fdr_threshold = fdr_threshold,
        likely_threshold = likely_chromoplexy_threshold,
        possible_threshold = possible_chromoplexy_threshold,
        sv_only = sv_only
    )
    summary_df$evidence_mode <- if (sv_only) "SV-only" else "SV+CNV"

    # Sort chains by combined evidence score
    if ("combined_score" %in% colnames(summary_df)) {
        chain_order <- order(summary_df$combined_score, decreasing = TRUE)
        summary_df <- summary_df[chain_order, ]
        chain_results <- chain_results[chain_order]
        all_chains <- all_chains[chain_order]
    }

    # Update individual chain summaries with classification
    for (i in seq_along(chain_results)) {
        chain_results[[i]]$summary$classification <- summary_df$classification[i]
    }

    if (verbose) {
        cat(sprintf("  -> Likely chromoplexy: %d\n", sum(summary_df$classification == "Likely chromoplexy")))
        cat(sprintf("  -> Possible chromoplexy: %d\n\n", sum(summary_df$classification == "Possible chromoplexy")))
        cat("========================================================================\n")
        cat("Detection complete!\n")
        cat("========================================================================\n\n")
    }

    result <- list(
        chains = all_chains,
        chain_details = chain_results,
        summary = summary_df,
        translocation_graph = tlx_graph,
        deletion_bridges = deletion_bridges,
        total_chains = length(all_chains),
        likely_chromoplexy = sum(summary_df$classification == "Likely chromoplexy"),
        possible_chromoplexy = sum(summary_df$classification == "Possible chromoplexy"),
        analysis_mode = if (sv_only) "SV-only chromoplexy" else "SV+CNV chromoplexy",
        limitations = if (sv_only) c(
            "CNV not provided: CN stability and deletion bridge support were not evaluated."
        ) else character(0),
        parameters = list(
            min_chromosomes = min_chromosomes,
            min_translocations = min_translocations,
            adjacency_distance = adjacency_distance,
            adjacency_decay_scale = adjacency_decay_scale, # Added
            deletion_bridge_distance = deletion_bridge_distance,
            fdr_threshold = fdr_threshold,
            use_statistical_testing = use_statistical_testing,
            likely_chromoplexy_threshold = likely_chromoplexy_threshold, # Added
            possible_chromoplexy_threshold = possible_chromoplexy_threshold, # Added
            sv_only = sv_only
        ),
        version = "3.0"
    )

    if (collapse_chains) {
        result$collapsed_events <- collapse_chromoplexy_chains(
            chromoplexy_result = result,
            classifications = collapse_classifications,
            gene_granges = gene_granges,
            breakpoint_padding = breakpoint_padding
        )
    } else {
        result$collapsed_events <- empty_collapsed_chromoplexy_events()
    }

    class(result) <- c("chromoplexy_v3", "chromoplexy", "list")
    return(result)
}


################################################################################
# CRITICAL FIX: Complete Graph Construction with THREE Edge Types
################################################################################

#' Build complete translocation graph with all edge types
#'
#' This is the corrected version that builds a CONNECTED graph by including:
#' 1. Translocation edges (inter-chromosomal SVs)
#' 2. Deletion bridge edges (deletions connecting nearby breakpoints)
#' 3. Genomic adjacency edges (nearby breakpoints on same chromosome)
#'
#' @param inter_chr_SVs Inter-chromosomal SVs
#' @param deletion_bridges Identified deletion bridges
#' @param adjacency_distance Maximum distance for adjacency edges (default: 10Mb)
#' @param adjacency_decay_scale Decay factor for genomic adjacency evidence strength (default: 5e6)
#' @return Complete translocation graph
#' @keywords internal
build_complete_translocation_graph_v3 <- function(inter_chr_SVs,
                                                  deletion_bridges = NULL,
                                                  adjacency_distance = 10e6,
                                                  adjacency_decay_scale = 2e6,
                                                  max_neighbors = 3) { # New parameter

    # STEP 1: Collect all breakpoints
    all_breakpoints <- collect_all_breakpoints_v3(inter_chr_SVs)

    # Create nodes from breakpoints
    nodes <- unique(all_breakpoints$node_id)

    # STEP 2: Initialize edge list
    edge_list <- list()

    # STEP 3: Add TRANSLOCATION edges
    for (i in 1:nrow(inter_chr_SVs)) {
        sv <- inter_chr_SVs[i, ]
        node1 <- paste(sv$chrom1, sv$pos1, sep = ":")
        node2 <- paste(sv$chrom2, sv$pos2, sep = ":")

        edge_list[[length(edge_list) + 1]] <- data.frame(
            from = node1,
            to = node2,
            edge_type = "TRANSLOCATION",
            sv_index = i,
            distance = NA, # Inter-chromosomal, no genomic distance
            evidence_strength = 1.0, # Strong evidence
            stringsAsFactors = FALSE
        )
    }

    # STEP 4: Add DELETION BRIDGE edges
    if (!is.null(deletion_bridges) && length(deletion_bridges) > 0) {
        for (bridge in deletion_bridges) {
            # Bridge connects two breakpoints via deletion
            # Typically: the SV breakpoint and a nearby breakpoint

            # Get the SV breakpoint
            sv_bp_node <- paste(bridge$chrom, bridge$sv_breakpoint_pos, sep = ":")

            # Find nearby breakpoints that could be connected by this deletion
            nearby_bps <- find_breakpoints_in_region_v3(
                all_breakpoints = all_breakpoints,
                chrom = bridge$chrom,
                start = bridge$deletion_start,
                end = bridge$deletion_end
            )

            # Connect SV breakpoint to nearby breakpoints via deletion bridge
            for (nearby_node in nearby_bps) {
                if (nearby_node != sv_bp_node) {
                    edge_list[[length(edge_list) + 1]] <- data.frame(
                        from = sv_bp_node,
                        to = nearby_node,
                        edge_type = "DELETION_BRIDGE",
                        sv_index = bridge$sv_index,
                        distance = bridge$deletion_size,
                        evidence_strength = bridge$confidence,
                        stringsAsFactors = FALSE
                    )
                }
            }
        }
    }

    # STEP 5: Add GENOMIC ADJACENCY edges
    # Connect nearby breakpoints on the same chromosome
    adjacency_edges <- build_genomic_adjacency_edges_v3(
        all_breakpoints = all_breakpoints,
        adjacency_distance = adjacency_distance,
        adjacency_decay_scale = adjacency_decay_scale,
        max_neighbors = max_neighbors # Pass new parameter
    )

    if (length(adjacency_edges) > 0) {
        edge_list <- c(edge_list, adjacency_edges)
    }

    # STEP 6: Combine all edges
    if (length(edge_list) > 0) {
        edges <- do.call(rbind, edge_list)
    } else {
        edges <- data.frame(
            from = character(0),
            to = character(0),
            edge_type = character(0),
            sv_index = integer(0),
            distance = numeric(0),
            evidence_strength = numeric(0),
            stringsAsFactors = FALSE
        )
    }

    # STEP 7: Build adjacency list
    adj_list <- build_adjacency_list_v3(nodes, edges)

    return(list(
        nodes = nodes,
        edges = edges,
        adjacency_list = adj_list,
        breakpoints = all_breakpoints
    ))
}


#' Collect all breakpoints from SVs
#'
#' @param inter_chr_SVs Inter-chromosomal SVs
#' @return Data frame of breakpoints with metadata
#' @keywords internal
collect_all_breakpoints_v3 <- function(inter_chr_SVs) {
    breakpoints <- list()

    for (i in 1:nrow(inter_chr_SVs)) {
        sv <- inter_chr_SVs[i, ]

        # Breakpoint 1
        breakpoints[[length(breakpoints) + 1]] <- data.frame(
            node_id = paste(sv$chrom1, sv$pos1, sep = ":"),
            chrom = sv$chrom1,
            pos = sv$pos1,
            sv_index = i,
            breakpoint_num = 1,
            stringsAsFactors = FALSE
        )

        # Breakpoint 2
        breakpoints[[length(breakpoints) + 1]] <- data.frame(
            node_id = paste(sv$chrom2, sv$pos2, sep = ":"),
            chrom = sv$chrom2,
            pos = sv$pos2,
            sv_index = i,
            breakpoint_num = 2,
            stringsAsFactors = FALSE
        )
    }

    bp_df <- do.call(rbind, breakpoints)

    # Remove duplicates (same breakpoint used by multiple SVs)
    bp_df <- bp_df[!duplicated(bp_df$node_id), ]

    # Sort by chromosome and position
    bp_df <- bp_df[order(bp_df$chrom, bp_df$pos), ]

    return(bp_df)
}


#' Build genomic adjacency edges
#'
#' Connects nearby breakpoints on the same chromosome.
#' This is CRITICAL for detecting chromoplexy chains.
#'
#' @param all_breakpoints All breakpoints
#' @param adjacency_distance Maximum distance for adjacency (default: 10Mb)
#' @param adjacency_decay_scale Decay factor for genomic adjacency evidence strength (default: 5e6)
#' @return List of adjacency edges
#' @keywords internal
build_genomic_adjacency_edges_v3 <- function(all_breakpoints,
                                             adjacency_distance = 10e6,
                                             adjacency_decay_scale = 2e6,
                                             max_neighbors = 3) { # New parameter

    adjacency_edges <- list()
    chroms <- unique(all_breakpoints$chrom)

    # Ensure GenomicRanges is available
    if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
        stop("Package 'GenomicRanges' needed for this function to work. Please install it.",
            call. = FALSE
        )
    }

    for (chr in chroms) {
        chr_bps <- all_breakpoints[all_breakpoints$chrom == chr, ]

        if (nrow(chr_bps) < 2) next

        # Create GRanges object for breakpoints on the current chromosome
        bp_gr <- GenomicRanges::GRanges(
            seqnames = S4Vectors::Rle(chr_bps$chrom),
            ranges = IRanges::IRanges(start = chr_bps$pos, end = chr_bps$pos),
            node_id = chr_bps$node_id,
            sv_index = chr_bps$sv_index
        )

        # Find overlaps (i.e., breakpoints within adjacency_distance of each other)
        hits <- GenomicRanges::findOverlaps(bp_gr, bp_gr,
            maxgap = adjacency_distance,
            type = "any",
            ignore.strand = TRUE
        )

        # Filter out self-hits and duplicate pairs (i, j) and (j, i)
        # Keep only i < j (downstream neighbors)
        hits <- hits[S4Vectors::queryHits(hits) < S4Vectors::subjectHits(hits)]

        if (length(hits) == 0) next

        # Optimize: Keep only k-nearest neighbors for each breakpoint
        # Since chr_bps is sorted by position, the hits for a query `i` will be `i+1, i+2...`
        # which are already sorted by distance. We just need to take the first k hits for each query.

        # Convert hits to data frame for filtering
        hits_df <- data.frame(
            query = S4Vectors::queryHits(hits),
            subject = S4Vectors::subjectHits(hits)
        )

        # Group by query and limit to max_neighbors
        # We can do this efficiently without heavy dependencies
        hits_df$rank <- ave(hits_df$subject, hits_df$query, FUN = seq_along)
        hits_df <- hits_df[hits_df$rank <= max_neighbors, ]

        if (nrow(hits_df) == 0) next

        for (k in 1:nrow(hits_df)) {
            query_idx <- hits_df$query[k]
            subject_idx <- hits_df$subject[k]

            bp1 <- chr_bps[query_idx, ]
            bp2 <- chr_bps[subject_idx, ]

            distance <- bp2$pos - bp1$pos

            adjacency_edges[[length(adjacency_edges) + 1]] <- data.frame(
                from = bp1$node_id,
                to = bp2$node_id,
                edge_type = "ADJACENCY",
                sv_index = NA,
                distance = distance,
                evidence_strength = exp(-distance / adjacency_decay_scale),
                stringsAsFactors = FALSE
            )
        }
    }

    return(adjacency_edges)
}


#' Find breakpoints in a genomic region
#'
#' @param all_breakpoints All breakpoints
#' @param chrom Chromosome
#' @param start Start position
#' @param end End position
#' @param padding Buffer zone for caller coordinate differences (default: 1000bp)
#' @return Vector of node IDs in region
#' @keywords internal
find_breakpoints_in_region_v3 <- function(all_breakpoints, chrom, start, end,
                                          padding = 1000) {
    chr_bps <- all_breakpoints[all_breakpoints$chrom == chrom, ]

    # Add padding to account for coordinate system differences between callers
    # Different SV/CNV callers may report slightly different positions due to:
    # - Microhomology handling
    # - 0-based vs 1-based coordinates
    # - Caller-specific algorithms
    in_region <- chr_bps[chr_bps$pos >= (start - padding) &
        chr_bps$pos <= (end + padding), ]

    return(in_region$node_id)
}


#' Build adjacency list from edges
#'
#' @param nodes All nodes
#' @param edges All edges
#' @return Adjacency list
#' @keywords internal
build_adjacency_list_v3 <- function(nodes, edges) {
    # Initialize adjacency list
    adj_list <- list()
    for (node in nodes) {
        adj_list[[node]] <- list(
            neighbors = character(0),
            edge_types = character(0),
            sv_indices = integer(0),
            distances = numeric(0),
            evidence_strengths = numeric(0)
        )
    }

    # Add edges to adjacency list
    for (i in 1:nrow(edges)) {
        edge <- edges[i, ]

        from_node <- edge$from
        to_node <- edge$to

        # Add to 'from' node's neighbors
        adj_list[[from_node]]$neighbors <- c(adj_list[[from_node]]$neighbors, to_node)
        adj_list[[from_node]]$edge_types <- c(adj_list[[from_node]]$edge_types, edge$edge_type)
        adj_list[[from_node]]$sv_indices <- c(
            adj_list[[from_node]]$sv_indices,
            if (is.na(edge$sv_index)) NA else edge$sv_index
        )
        adj_list[[from_node]]$distances <- c(adj_list[[from_node]]$distances, edge$distance)
        adj_list[[from_node]]$evidence_strengths <- c(
            adj_list[[from_node]]$evidence_strengths,
            edge$evidence_strength
        )

        # Add to 'to' node's neighbors (undirected graph)
        adj_list[[to_node]]$neighbors <- c(adj_list[[to_node]]$neighbors, from_node)
        adj_list[[to_node]]$edge_types <- c(adj_list[[to_node]]$edge_types, edge$edge_type)
        adj_list[[to_node]]$sv_indices <- c(
            adj_list[[to_node]]$sv_indices,
            if (is.na(edge$sv_index)) NA else edge$sv_index
        )
        adj_list[[to_node]]$distances <- c(adj_list[[to_node]]$distances, edge$distance)
        adj_list[[to_node]]$evidence_strengths <- c(
            adj_list[[to_node]]$evidence_strengths,
            edge$evidence_strength
        )
    }

    return(adj_list)
}


################################################################################
# Chromoplexy Detection v3.0 - Part 2
# Deletion Bridges and Chain Detection
################################################################################

#' Identify deletion bridges (v3 - enhanced)
#'
#' Identifies deletions that connect nearby translocation breakpoints.
#' This version returns bridges in a format suitable for graph construction.
#'
#' @param SV.sample Inter-chromosomal SVs
#' @param CNV.sample CNV data
#' @param max_distance Maximum distance to search (default: 1Mb)
#' @return List of deletion bridges with detailed information
#' @keywords internal
identify_deletion_bridges_v3 <- function(SV.sample,
                                         CNV.sample,
                                         max_distance = 1e6) {
    bridges <- list()

    # For each translocation breakpoint
    for (i in 1:nrow(SV.sample)) {
        sv <- SV.sample[i, ]

        # Check both breakpoints
        for (bp_num in 1:2) {
            chrom <- if (bp_num == 1) sv$chrom1 else sv$chrom2
            pos <- if (bp_num == 1) sv$pos1 else sv$pos2

            # Find deletions near this breakpoint
            nearby_deletions <- find_nearby_deletions_v3(
                chrom = chrom,
                pos = pos,
                CNV.sample = CNV.sample,
                max_distance = max_distance
            )

            if (length(nearby_deletions) > 0) {
                for (del in nearby_deletions) {
                    bridges[[length(bridges) + 1]] <- list(
                        sv_index = i,
                        chrom = chrom,
                        sv_breakpoint_pos = pos,
                        breakpoint_num = bp_num,
                        deletion_start = del$start,
                        deletion_end = del$end,
                        deletion_size = del$size,
                        deletion_cn = del$cn,
                        distance_to_breakpoint = del$distance_to_breakpoint,
                        confidence = calculate_bridge_confidence_v3(
                            distance = del$distance_to_breakpoint,
                            deletion_size = del$size,
                            cn = del$cn
                        )
                    )
                }
            }
        }
    }

    return(bridges)
}


#' Find deletions near a breakpoint (v3)
#'
#' @param chrom Chromosome
#' @param pos Position
#' @param CNV.sample CNV data
#' @param max_distance Maximum search distance
#' @return List of nearby deletions
#' @keywords internal
find_nearby_deletions_v3 <- function(chrom, pos, CNV.sample, max_distance) {
    # Get CNV segments on this chromosome
    chr_cnv <- CNV.sample[CNV.sample$chrom == chrom, ]

    if (nrow(chr_cnv) == 0) {
        return(list())
    }

    deletions <- list()

    for (i in 1:nrow(chr_cnv)) {
        seg <- chr_cnv[i, ]

        # Check if this is a deletion (CN < 2). Handle potential NAs.
        if (is.na(seg$total_cn) || seg$total_cn >= 2) next

        # Check if deletion overlaps or is near the breakpoint
        distance_to_start <- abs(seg$start - pos)
        distance_to_end <- abs(seg$end - pos)
        min_distance <- min(distance_to_start, distance_to_end)

        # Also check if breakpoint is inside deletion
        if (pos >= seg$start && pos <= seg$end) {
            min_distance <- 0
        }

        if (min_distance <= max_distance) {
            deletions[[length(deletions) + 1]] <- list(
                start = seg$start,
                end = seg$end,
                cn = seg$total_cn,
                size = seg$end - seg$start,
                distance_to_breakpoint = min_distance
            )
        }
    }

    return(deletions)
}


#' Calculate deletion bridge confidence (v3)
#'
#' @param distance Distance to breakpoint
#' @param deletion_size Deletion size
#' @param cn Copy number
#' @return Confidence score (0-1)
#' @keywords internal
calculate_bridge_confidence_v3 <- function(distance, deletion_size, cn) {
    # Component 1: Distance score (closer = higher confidence)
    distance_score <- exp(-distance / 2e5) # 200kb scale

    # Component 2: Deletion size score
    # Typical deletion bridges are 100bp - 10kb
    size_score <- if (deletion_size >= 100 && deletion_size <= 1e4) {
        1.0
    } else if (deletion_size < 100) {
        deletion_size / 100
    } else {
        exp(-(deletion_size - 1e4) / 1e4)
    }

    # Component 3: CN score (lower CN = more confident deletion)
    cn_score <- exp(-(cn - 0) / 1) # CN=0 is highest confidence

    # Combine (weighted average)
    confidence <- (distance_score * 0.5 + size_score * 0.3 + cn_score * 0.2)

    return(confidence)
}


################################################################################
# Chain Detection (adapted for new graph structure)
################################################################################

#' Detect all possible chains in the graph (v3)
#'
#' Adapted to work with the new graph structure that includes
#' genomic adjacency edges.
#'
#' @param tlx_graph Complete translocation graph
#' @param inter_chr_SVs Inter-chromosomal SVs
#' @param min_chromosomes Minimum chromosomes
#' @param min_translocations Minimum translocations
#' @param allow_cycles Allow cycles
#' @param max_path_search Maximum paths to search (default: 1000).
#'   Increased to handle complex chromoplexy patterns with many translocations.
#' @param max_search_depth Maximum chain length to explore (default: 10).
#' @return List of detected chains
#' @keywords internal
detect_all_chains_v3 <- function(tlx_graph,
                                 inter_chr_SVs,
                                 min_chromosomes = 3,
                                 min_translocations = 3,
                                 allow_cycles = TRUE,
                                 max_path_search = 1000,
                                 max_search_depth = 10,
                                 max_adjacency_streak = 3) {
    adj_list <- tlx_graph$adjacency_list

    # Prune only nodes that have ABSOLUTELY no translocations
    relevant_nodes <- character(0)
    for (node in tlx_graph$nodes) {
        if (any(adj_list[[node]]$edge_types == "TRANSLOCATION")) {
            relevant_nodes <- c(relevant_nodes, node)
        }
    }

    node_scores <- sapply(relevant_nodes, function(node) {
        sum(adj_list[[node]]$edge_types == "TRANSLOCATION")
    })
    sorted_start_nodes <- relevant_nodes[order(node_scores, decreasing = TRUE)]

    all_chains <- list()
    chain_id <- 0
    # Keep track of unique signatures to avoid redundant evaluation later
    found_signatures <- new.env(hash = TRUE)

    for (start_node in sorted_start_nodes) {
        paths_from_node <- backtrack_all_paths_v3(
            start_node = start_node,
            adj_list = adj_list,
            max_paths_total = max_path_search,
            min_length_req = min_translocations,
            max_depth_req = max_search_depth,
            max_adjacency_streak = max_adjacency_streak
        )

        if (length(paths_from_node) == 0) next

        for (path_result in paths_from_node) {
            path <- path_result$path
            edge_types <- path_result$edge_types
            sv_indices <- unique(path_result$sv_indices[!is.na(path_result$sv_indices)])

            # Check criteria
            if (length(path) >= min_translocations) {
                chroms <- unique(sapply(path, function(x) strsplit(x, ":")[[1]][1]))
                n_translocation_edges <- sum(edge_types == "TRANSLOCATION", na.rm = TRUE)

                if (length(chroms) >= min_chromosomes && n_translocation_edges >= min_translocations) {
                    # Quick signature check to avoid obviously redundant chains in the list
                    sig <- paste(sort(sv_indices), collapse = "_")
                    if (!exists(sig, envir = found_signatures)) {
                        chain_id <- chain_id + 1
                        assign(sig, TRUE, envir = found_signatures)

                        is_cycle <- FALSE
                        if (allow_cycles && length(path) > 1) {
                            first_node <- path[1]
                            last_node <- path[length(path)]
                            if (last_node %in% adj_list[[first_node]]$neighbors) {
                                is_cycle <- TRUE
                            }
                        }

                        all_chains[[chain_id]] <- list(
                            id = chain_id,
                            nodes = path,
                            edge_types = edge_types,
                            sv_indices = sv_indices,
                            chromosomes = chroms,
                            n_chromosomes = length(chroms),
                            n_translocations = n_translocation_edges,
                            n_deletion_bridges = sum(edge_types == "DELETION_BRIDGE", na.rm = TRUE),
                            n_adjacency_edges = sum(edge_types == "ADJACENCY", na.rm = TRUE),
                            is_cycle = is_cycle
                        )
                    }
                }
            }
        }

        # Adjust dynamic limit for very dense graphs
        if (length(all_chains) > 20000) break
    }

    # Final deduplication
    if (length(all_chains) > 0) {
        all_chains <- deduplicate_chains_v3(all_chains)
    }

    return(all_chains)
}


#' Backtrack to find all possible paths (v3)
#'
#' Enhanced version that tracks edge types.
#'
#' NOTE: Future optimization could implement greedy heuristic to prioritize
#' TRANSLOCATION edges over ADJACENCY edges when exploring paths. This would
#' reduce search space while maintaining biological relevance.
#'
#' @keywords internal
backtrack_all_paths_v3 <- function(start_node,
                                   adj_list,
                                   max_paths_total,
                                   min_length_req,
                                   max_depth_req,
                                   max_adjacency_streak) {
    paths_found <- list()

    # Stack stores: list(current_node, path, edge_types, sv_indices, visited, adj_streak)
    stack <- list(list(
        node = start_node,
        path = start_node,
        edge_types = character(0),
        sv_indices = integer(0),
        visited = start_node,
        adj_streak = 0
    ))

    while (length(stack) > 0) {
        # Pop from stack (LIFO for DFS)
        curr <- stack[[length(stack)]]
        stack <- stack[-length(stack)]

        # 1. Check if current path meets criteria
        if (length(curr$path) >= min_length_req) {
            if (sum(curr$edge_types == "TRANSLOCATION") >= 2) {
                paths_found[[length(paths_found) + 1]] <- list(
                    path = curr$path,
                    edge_types = curr$edge_types,
                    sv_indices = curr$sv_indices
                )
            }
        }

        # Stop if we hit global limit
        if (length(paths_found) >= max_paths_total) break

        # Stop if we hit depth limit
        if (length(curr$path) >= max_depth_req) next

        # 2. Get neighbors
        neighbors_data <- adj_list[[as.character(curr$node)]]
        if (is.null(neighbors_data)) next

        neighbors <- neighbors_data$neighbors
        edge_types <- neighbors_data$edge_types
        sv_indices <- neighbors_data$sv_indices

        # Prioritize translocations by sorting (reverse order for stack push)
        # We want to process Translocations FIRST, so we push them LAST onto the stack
        ord <- order(edge_types == "TRANSLOCATION", edge_types == "DELETION_BRIDGE", decreasing = FALSE)

        for (i in ord) {
            neighbor <- neighbors[i]
            edge_type <- edge_types[i]

            # Pruning
            if (neighbor %in% curr$visited) next

            new_streak <- if (edge_type == "ADJACENCY") curr$adj_streak + 1 else 0
            if (new_streak > max_adjacency_streak) next

            # Push to stack
            stack[[length(stack) + 1]] <- list(
                node = neighbor,
                path = c(curr$path, neighbor),
                edge_types = c(curr$edge_types, edge_type),
                sv_indices = c(curr$sv_indices, sv_indices[i]),
                visited = c(curr$visited, neighbor),
                adj_streak = new_streak
            )
        }
    }

    return(paths_found)
}


#' Deduplicate chains (v3)
#'
#' @param chains List of chains
#' @return Deduplicated chains
#' @keywords internal
deduplicate_chains_v3 <- function(chains) {
    if (length(chains) <= 1) {
        return(chains)
    }
    sigs <- sapply(chains, function(x) paste(sort(x$sv_indices), collapse = "_"))
    unique_chains <- chains[which(!duplicated(sigs))]
    for (i in seq_along(unique_chains)) unique_chains[[i]]$id <- i

    cat(sprintf("  Deduplicated to %d unique chains\n", length(unique_chains)))
    return(unique_chains)
}


################################################################################
# Chain Evaluation (v3)
################################################################################

#' Evaluate chromoplexy chain (v3)
#'
#' Enhanced to consider edge types and graph structure.
#'
#' @keywords internal
evaluate_chromoplexy_chain_v3 <- function(chain,
                                          inter_chr_SVs,
                                          CNV.sample,
                                          tlx_graph,
                                          max_cn_change = 1,
                                          use_statistical_testing = TRUE,
                                          genome_size = NULL,
                                          fdr_threshold = 0.01) {
    # Get SVs in this chain
    chain_SVs <- inter_chr_SVs[chain$sv_indices, ]

    # 1. Evaluate copy number stability
    cn_eval <- evaluate_cn_stability_enhanced_v3(chain, chain_SVs, CNV.sample)

    # 2. Calculate chain complexity score
    complexity_score <- calculate_chain_complexity_v3(chain, chain_SVs)

    # 3. Evaluate edge type composition
    edge_composition <- evaluate_edge_composition_v3(chain)

    # 4. Statistical significance testing
    statistical_significance <- NULL
    if (use_statistical_testing && !is.null(genome_size)) {
        # Source from chromoplexy_statistics.R
        if (exists("calculate_chain_significance")) {
            statistical_significance <- calculate_chain_significance(
                chain = chain,
                inter_chr_SVs = inter_chr_SVs,
                genome_size = genome_size,
                fdr_threshold = fdr_threshold
            )
        }
    }

    # 5. Calculate combined evidence score
    combined_score <- calculate_combined_evidence_score_v3(
        cn_stability_score = cn_eval$combined_score,
        complexity_score = complexity_score,
        edge_composition = edge_composition,
        statistical_significance = statistical_significance
    )

    # Create summary
    summary <- data.frame(
        chain_id = chain$id,
        n_chromosomes = chain$n_chromosomes,
        chromosomes_involved = paste(chain$chromosomes, collapse = ","),
        n_translocations = chain$n_translocations,
        n_deletion_bridges = chain$n_deletion_bridges,
        n_adjacency_edges = chain$n_adjacency_edges,
        is_cycle = chain$is_cycle,
        # CN metrics
        cn_stability_score = cn_eval$combined_score,
        max_cn_deviation = cn_eval$max_deviation,
        # Complexity
        complexity_score = complexity_score,
        # Edge composition
        deletion_bridge_fraction = edge_composition$deletion_bridge_fraction,
        adjacency_fraction = edge_composition$adjacency_fraction,
        # Statistical
        pvalue = if (!is.null(statistical_significance)) statistical_significance$pvalue else NA,
        fdr = if (!is.null(statistical_significance)) statistical_significance$fdr else NA,
        is_statistically_significant = if (!is.null(statistical_significance)) {
            statistical_significance$is_significant
        } else {
            NA
        },
        # Combined
        combined_score = combined_score,
        stringsAsFactors = FALSE
    )

    return(list(
        summary = summary,
        chain = chain,
        SVs = chain_SVs,
        cn_evaluation = cn_eval,
        edge_composition = edge_composition,
        statistical_significance = statistical_significance
    ))
}


#' Evaluate edge type composition (v3)
#'
#' @keywords internal
evaluate_edge_composition_v3 <- function(chain) {
    total_edges <- length(chain$edge_types)

    if (total_edges == 0) {
        return(list(
            deletion_bridge_fraction = 0,
            adjacency_fraction = 0,
            translocation_fraction = 0
        ))
    }

    n_del_bridge <- sum(chain$edge_types == "DELETION_BRIDGE", na.rm = TRUE)
    n_adjacency <- sum(chain$edge_types == "ADJACENCY", na.rm = TRUE)
    n_translocation <- sum(chain$edge_types == "TRANSLOCATION", na.rm = TRUE)

    return(list(
        deletion_bridge_fraction = n_del_bridge / total_edges,
        adjacency_fraction = n_adjacency / total_edges,
        translocation_fraction = n_translocation / total_edges,
        n_deletion_bridges = n_del_bridge,
        n_adjacency_edges = n_adjacency,
        n_translocations = n_translocation
    ))
}


#' Enhanced CN stability evaluation (v3)
#'
#' @keywords internal
evaluate_cn_stability_enhanced_v3 <- function(chain, chain_SVs, CNV.sample) {
    # Similar to v2, but simplified for now
    if (is.null(CNV.sample) || is.null(nrow(CNV.sample)) || nrow(CNV.sample) == 0) {
        return(list(
            combined_score = NA_real_,
            max_deviation = NA_real_,
            mean_deviation = NA_real_,
            evaluated = FALSE
        ))
    }

    cn_deviations <- c()

    for (chr in chain$chromosomes) {
        chr_cnv <- CNV.sample[CNV.sample$chrom == chr, ]
        if (nrow(chr_cnv) > 0) {
            deviations <- abs(chr_cnv$total_cn - 2)
            cn_deviations <- c(cn_deviations, deviations)
        }
    }

    max_deviation <- if (length(cn_deviations) > 0) max(cn_deviations, na.rm = TRUE) else 0
    mean_deviation <- if (length(cn_deviations) > 0) mean(cn_deviations, na.rm = TRUE) else 0

    combined_score <- exp(-mean_deviation / 2)

    return(list(
        combined_score = combined_score,
        max_deviation = max_deviation,
        mean_deviation = mean_deviation,
        evaluated = TRUE
    ))
}


#' Calculate chain complexity (v3)
#'
#' @keywords internal
calculate_chain_complexity_v3 <- function(chain, chain_SVs) {
    chr_score <- min(chain$n_chromosomes / 5, 1.0)
    tlx_score <- min(chain$n_translocations / 10, 1.0)
    cycle_bonus <- if (chain$is_cycle) 0.2 else 0.0

    complexity_score <- (chr_score * 0.4 + tlx_score * 0.4 + cycle_bonus)

    return(complexity_score)
}


#' Calculate combined evidence score (v3)
#'
#' @keywords internal
calculate_combined_evidence_score_v3 <- function(cn_stability_score,
                                                 complexity_score,
                                                 edge_composition,
                                                 statistical_significance = NULL) {
    scores <- c(cn_stability_score, complexity_score)
    weights <- c(0.3, 0.3)

    # Add deletion bridge score
    del_bridge_score <- edge_composition$deletion_bridge_fraction
    scores <- c(scores, del_bridge_score)
    weights <- c(weights, 0.2)

    # Add statistical significance if available
    if (!is.null(statistical_significance) && !is.na(statistical_significance$is_significant)) {
        if (statistical_significance$is_significant) {
            stat_score <- 1 - min(statistical_significance$fdr * 10, 1.0)
            scores <- c(scores, stat_score)
            weights <- c(weights, 0.2)
        } else {
            scores <- c(scores, 0.3)
            weights <- c(weights, 0.2)
        }
    }

    valid <- !is.na(scores) & !is.na(weights)
    scores <- scores[valid]
    weights <- weights[valid]
    if (length(scores) == 0 || sum(weights) == 0) {
        return(NA_real_)
    }

    # Weighted average
    combined <- sum(scores * weights) / sum(weights)

    return(combined)
}


#' Classify chromoplexy event (v3)
#'
#' @param summary_df Chain-level summary data frame
#' @param use_statistical_testing Whether to use statistical significance in classification
#' @param fdr_threshold FDR threshold for statistical significance (default: 0.01)
#' @param likely_threshold Minimum criteria met for "Likely chromoplexy" (default: 5)
#' @param possible_threshold Minimum criteria met for "Possible chromoplexy" (default: 4)
#' @param sv_only Whether CNV evidence is unavailable and classification should use SV-only criteria
#' @keywords internal
classify_chromoplexy_event_v3 <- function(summary_df,
                                          use_statistical_testing = TRUE,
                                          fdr_threshold = 0.01,
                                          likely_threshold = 5, # New parameter
                                          possible_threshold = 4, # New parameter
                                          sv_only = FALSE) {

    classifications <- character(nrow(summary_df))

    for (i in 1:nrow(summary_df)) {
        row <- summary_df[i, ]

        # Core criteria
        meets_chr_criteria <- row$n_chromosomes >= 3
        meets_tlx_criteria <- row$n_translocations >= 3
        meets_cn_criteria <- !is.na(row$cn_stability_score) && row$cn_stability_score >= 0.7
        meets_complexity <- row$complexity_score >= 0.3

        # Enhanced criteria
        meets_deletion_bridge <- row$deletion_bridge_fraction >= 0.2 # At least 20% of edges
        meets_statistical <- FALSE
        if (use_statistical_testing && !is.na(row$is_statistically_significant)) {
            meets_statistical <- row$is_statistically_significant
        }

        all_criteria <- c(
            meets_chr_criteria,
            meets_tlx_criteria,
            meets_cn_criteria,
            meets_complexity,
            meets_deletion_bridge,
            meets_statistical
        )

        criteria_met <- sum(all_criteria)

        # Classification
        if (sv_only) {
            core_criteria_met <- sum(c(meets_chr_criteria, meets_tlx_criteria, meets_complexity))
            if (core_criteria_met == 3 && (!use_statistical_testing || meets_statistical)) {
                classifications[i] <- "Likely chromoplexy"
            } else if (core_criteria_met == 3) {
                classifications[i] <- "Possible chromoplexy"
            } else if (core_criteria_met >= 2) {
                classifications[i] <- "Unlikely chromoplexy"
            } else {
                classifications[i] <- "Not chromoplexy"
            }
        } else if (use_statistical_testing) {
            if (meets_statistical && criteria_met >= likely_threshold) {
                classifications[i] <- "Likely chromoplexy"
            } else if (criteria_met >= possible_threshold) {
                classifications[i] <- "Possible chromoplexy"
            } else if (criteria_met >= 3) {
                classifications[i] <- "Unlikely chromoplexy"
            } else {
                classifications[i] <- "Not chromoplexy"
            }
        } else {
            # When statistical testing is disabled, max score is reduced by 1
            # Adjust thresholds accordingly
            adjusted_likely <- max(3, likely_threshold - 1)
            adjusted_possible <- max(2, possible_threshold - 1)

            if (criteria_met >= adjusted_likely) {
                classifications[i] <- "Likely chromoplexy"
            } else if (criteria_met >= adjusted_possible) {
                classifications[i] <- "Possible chromoplexy"
            } else if (criteria_met >= 2) {
                classifications[i] <- "Unlikely chromoplexy"
            } else {
                classifications[i] <- "Not chromoplexy"
            }
        }
    }

    return(classifications)
}


#' Create empty result (v3)
#'
#' @param analysis_mode Short description of the analysis mode
#' @param limitations Character vector describing limitations for this result
#' @return Empty chromoplexy result object
#' @keywords internal
create_empty_chromoplexy_result_v3 <- function(analysis_mode = "SV+CNV chromoplexy",
                                               limitations = character(0)) {
    result <- list(
        chains = list(),
        chain_details = list(),
        summary = data.frame(),
        translocation_graph = NULL,
        deletion_bridges = NULL,
        total_chains = 0,
        likely_chromoplexy = 0,
        possible_chromoplexy = 0,
        collapsed_events = empty_collapsed_chromoplexy_events(),
        analysis_mode = analysis_mode,
        limitations = limitations,
        version = "3.0"
    )

    class(result) <- c("chromoplexy_v3", "chromoplexy", "list")
    return(result)
}
