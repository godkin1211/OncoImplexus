################################################################################
# ChainFinder-style Chromoplexy Detection
#
# A standalone implementation following the original ChainFinder methodology
# (Baca et al. 2013, Cell 153(3):666-677)
#
# Key features:
# 1. Statistical null model for adjacency edges (not distance-based)
# 2. Local background rearrangement rate (mu) calculation
# 3. P-value threshold for edge creation
#
# This is independent of the OncoImplexus default chromoplexy detection
################################################################################

#' Detect chromoplexy using ChainFinder-style statistical model
#'
#' This function implements chromoplexy detection following the original
#' ChainFinder methodology (Baca et al. 2013), which uses statistical testing
#' to determine whether breakpoint adjacency is significant given the local
#' background rearrangement rate.
#'
#' @param SV.sample SV data (data frame or SVs object)
#' @param CNV.sample CNV data (data frame or CNVsegs object), optional
#' @param min_chromosomes Minimum chromosomes in chain (default: 3)
#' @param min_translocations Minimum translocations in chain (default: 3)
#' @param adjacency_pvalue_threshold P-value threshold for adjacency (default: 0.01)
#' @param adjacency_max_distance Maximum distance to consider (default: 10Mb)
#' @param deletion_bridge_distance Distance for deletion bridge search (default: 1Mb)
#' @param use_chromosome_specific_mu Use per-chromosome mu (default: TRUE)
#' @param max_path_search Maximum paths to search (default: 100)
#' @param max_search_depth Maximum chain length (default: 10)
#' @param genome Reference genome (default: "hg19")
#' @param verbose Print progress (default: TRUE)
#' @return A list containing chromoplexy detection results
#' @export
detect_chromoplexy_chainfinder <- function(SV.sample,
                                           CNV.sample = NULL,
                                           min_chromosomes = 3,
                                           min_translocations = 3,
                                           adjacency_pvalue_threshold = 0.01,
                                           adjacency_max_distance = 10e6,
                                           deletion_bridge_distance = 1e6,
                                           use_chromosome_specific_mu = TRUE,
                                           max_path_search = 100,
                                           max_search_depth = 10,
                                           genome = "hg19",
                                           verbose = TRUE) {

    # Convert inputs to data frames
    sv_df <- .cf_to_dataframe(SV.sample, "SV")
    cnv_df <- if (!is.null(CNV.sample)) .cf_to_dataframe(CNV.sample, "CNV") else NULL

    # Extract inter-chromosomal SVs only
    inter_chr_svs <- sv_df[sv_df$chrom1 != sv_df$chrom2, ]

    if (nrow(inter_chr_svs) < min_translocations) {
        if (verbose) {
            message(sprintf("Only %d inter-chromosomal SVs found. Need >= %d.",
                          nrow(inter_chr_svs), min_translocations))
        }
        return(.cf_empty_result())
    }

    if (verbose) {
        cat("========================================================================\n")
        cat("ChainFinder-style Chromoplexy Detection\n")
        cat("(Statistical Adjacency Model - Baca et al. 2013)\n")
        cat("========================================================================\n\n")
        cat(sprintf("Input: %d inter-chromosomal SVs (translocations)\n", nrow(inter_chr_svs)))
        cat(sprintf("Adjacency P-value threshold: %.4f\n\n", adjacency_pvalue_threshold))
    }

    # Get chromosome sizes
    chr_sizes <- .cf_get_chromosome_sizes(genome)

    # Step 1: Calculate local rearrangement rates (mu)
    if (verbose) cat("Step 1: Calculating local rearrangement rates (mu)...\n")
    mu_estimates <- .cf_calculate_mu(sv_df, chr_sizes, use_chromosome_specific_mu)

    if (verbose) {
        cat(sprintf("  Global mu: %.6f breakpoints/Mb\n", mu_estimates$global_mu))
    }

    # Step 2: Identify deletion bridges
    deletion_bridges <- list()
    if (!is.null(cnv_df)) {
        if (verbose) cat("\nStep 2: Identifying deletion bridges...\n")
        deletion_bridges <- .cf_find_deletion_bridges(
            inter_chr_svs, cnv_df, deletion_bridge_distance
        )
        if (verbose) cat(sprintf("  Found %d deletion bridges\n", length(deletion_bridges)))
    }

    # Step 3: Build graph with statistical adjacency model
    if (verbose) cat("\nStep 3: Building translocation graph (statistical model)...\n")
    graph <- .cf_build_graph(
        inter_chr_svs = inter_chr_svs,
        deletion_bridges = deletion_bridges,
        mu_estimates = mu_estimates,
        pvalue_threshold = adjacency_pvalue_threshold,
        max_distance = adjacency_max_distance
    )

    if (verbose) {
        cat(sprintf("  Nodes: %d, Edges: %d\n", length(graph$nodes), nrow(graph$edges)))
        cat(sprintf("    - Translocation edges: %d\n",
                   sum(graph$edges$edge_type == "TRANSLOCATION")))
        cat(sprintf("    - Deletion bridge edges: %d\n",
                   sum(graph$edges$edge_type == "DELETION_BRIDGE")))
        cat(sprintf("    - Statistical adjacency edges: %d\n",
                   sum(graph$edges$edge_type == "STAT_ADJACENCY")))
    }

    # Step 4: Find chains
    if (verbose) cat("\nStep 4: Detecting translocation chains...\n")
    chains <- .cf_find_chains(
        graph = graph,
        min_chromosomes = min_chromosomes,
        min_translocations = min_translocations,
        max_path_search = max_path_search,
        max_search_depth = max_search_depth
    )

    if (length(chains) == 0) {
        if (verbose) cat("  No chromoplexy chains detected.\n")
        return(.cf_empty_result())
    }

    if (verbose) cat(sprintf("  Found %d potential chains\n", length(chains)))

    # Step 5: Evaluate and classify chains
    if (verbose) cat("\nStep 5: Evaluating and classifying chains...\n")
    evaluated_chains <- lapply(chains, function(chain) {
        .cf_evaluate_chain(chain, inter_chr_svs, cnv_df, graph, mu_estimates)
    })

    # Create summary
    summary_df <- do.call(rbind, lapply(evaluated_chains, function(x) x$summary))
    summary_df$classification <- .cf_classify_chains(summary_df)

    n_likely <- sum(summary_df$classification == "Likely chromoplexy")
    n_possible <- sum(summary_df$classification == "Possible chromoplexy")

    if (verbose) {
        cat(sprintf("  Likely chromoplexy: %d\n", n_likely))
        cat(sprintf("  Possible chromoplexy: %d\n", n_possible))
        cat("\n========================================================================\n")
    }

    result <- list(
        chains = chains,
        chain_details = evaluated_chains,
        summary = summary_df,
        graph = graph,
        deletion_bridges = deletion_bridges,
        mu_estimates = mu_estimates,
        total_chains = length(chains),
        likely_chromoplexy = n_likely,
        possible_chromoplexy = n_possible,
        parameters = list(
            min_chromosomes = min_chromosomes,
            min_translocations = min_translocations,
            adjacency_pvalue_threshold = adjacency_pvalue_threshold,
            adjacency_max_distance = adjacency_max_distance,
            use_chromosome_specific_mu = use_chromosome_specific_mu
        ),
        method = "ChainFinder",
        version = "1.0"
    )

    class(result) <- c("chromoplexy_chainfinder", "list")
    return(result)
}


################################################################################
# Internal Helper Functions (all prefixed with .cf_)
################################################################################

#' Convert input to data frame
#' @keywords internal
.cf_to_dataframe <- function(x, type = "SV") {
    if (is.data.frame(x)) {
        return(x)
    } else if (inherits(x, "SVs")) {
        return(data.frame(
            chrom1 = x@chrom1,
            pos1 = x@pos1,
            strand1 = x@strand1,
            chrom2 = x@chrom2,
            pos2 = x@pos2,
            strand2 = x@strand2,
            SVtype = x@SVtype,
            stringsAsFactors = FALSE
        ))
    } else if (inherits(x, "CNVsegs")) {
        return(data.frame(
            chrom = x@chrom,
            start = x@start,
            end = x@end,
            total_cn = x@total_cn,
            stringsAsFactors = FALSE
        ))
    } else {
        stop("Unknown input type")
    }
}

#' Get chromosome sizes
#' @keywords internal
.cf_get_chromosome_sizes <- function(genome = "hg19") {
    if (genome == "hg38") {
        sizes <- c(
            "1" = 248956422, "2" = 242193529, "3" = 198295559, "4" = 190214555,
            "5" = 181538259, "6" = 170805979, "7" = 159345973, "8" = 145138636,
            "9" = 138394717, "10" = 133797422, "11" = 135086622, "12" = 133275309,
            "13" = 114364328, "14" = 107043718, "15" = 101991189, "16" = 90338345,
            "17" = 83257441, "18" = 80373285, "19" = 58617616, "20" = 64444167,
            "21" = 46709983, "22" = 50818468, "X" = 156040895, "Y" = 57227415
        )
    } else {
        # hg19
        sizes <- c(
            "1" = 249250621, "2" = 243199373, "3" = 198022430, "4" = 191154276,
            "5" = 180915260, "6" = 171115067, "7" = 159138663, "8" = 146364022,
            "9" = 141213431, "10" = 135534747, "11" = 135006516, "12" = 133851895,
            "13" = 115169878, "14" = 107349540, "15" = 102531392, "16" = 90354753,
            "17" = 81195210, "18" = 78077248, "19" = 59128983, "20" = 63025520,
            "21" = 48129895, "22" = 51304566, "X" = 155270560, "Y" = 59373566
        )
    }
    return(sizes)
}

#' Calculate local rearrangement rate (mu)
#' @keywords internal
.cf_calculate_mu <- function(sv_df, chr_sizes, use_chromosome_specific = TRUE) {

    # Count breakpoints per chromosome
    bp_counts <- table(c(as.character(sv_df$chrom1), as.character(sv_df$chrom2)))

    # Global mu
    total_bp <- sum(bp_counts)
    total_size_mb <- sum(chr_sizes) / 1e6
    global_mu <- total_bp / total_size_mb

    # Chromosome-specific mu with shrinkage toward global
    chr_mu <- numeric()
    if (use_chromosome_specific) {
        for (chr in names(chr_sizes)) {
            n_bp <- if (chr %in% names(bp_counts)) as.numeric(bp_counts[chr]) else 0
            size_mb <- chr_sizes[chr] / 1e6

            # Empirical Bayes shrinkage
            chr_specific <- n_bp / size_mb
            weight <- min(n_bp / 10, 1)  # More data = more trust in local estimate
            chr_mu[chr] <- weight * chr_specific + (1 - weight) * global_mu
        }
    }

    list(
        global_mu = global_mu,
        chr_mu = chr_mu,
        total_breakpoints = total_bp
    )
}

#' Find deletion bridges
#' @keywords internal
.cf_find_deletion_bridges <- function(inter_chr_svs, cnv_df, max_distance) {

    bridges <- list()

    # Find deletions (CN < 2)
    deletions <- cnv_df[!is.na(cnv_df$total_cn) & cnv_df$total_cn < 2, ]
    if (nrow(deletions) == 0) return(bridges)

    for (i in 1:nrow(inter_chr_svs)) {
        sv <- inter_chr_svs[i, ]

        # Check both breakpoints
        for (bp_num in 1:2) {
            chr <- if (bp_num == 1) sv$chrom1 else sv$chrom2
            pos <- if (bp_num == 1) sv$pos1 else sv$pos2

            # Find nearby deletions
            chr_del <- deletions[deletions$chrom == chr, ]
            for (j in seq_len(nrow(chr_del))) {
                del <- chr_del[j, ]

                # Check if deletion is near breakpoint
                dist_start <- abs(del$start - pos)
                dist_end <- abs(del$end - pos)
                min_dist <- min(dist_start, dist_end)

                # Or if breakpoint is inside deletion
                if (pos >= del$start && pos <= del$end) {
                    min_dist <- 0
                }

                if (min_dist <= max_distance) {
                    bridges[[length(bridges) + 1]] <- list(
                        sv_index = i,
                        chrom = chr,
                        sv_pos = pos,
                        del_start = del$start,
                        del_end = del$end,
                        del_cn = del$total_cn,
                        distance = min_dist
                    )
                }
            }
        }
    }

    return(bridges)
}

#' Build graph with statistical adjacency model
#' @keywords internal
.cf_build_graph <- function(inter_chr_svs, deletion_bridges, mu_estimates,
                            pvalue_threshold, max_distance) {

    # Collect all breakpoints
    breakpoints <- data.frame(
        node_id = character(),
        chrom = character(),
        pos = numeric(),
        sv_index = integer(),
        stringsAsFactors = FALSE
    )

    for (i in 1:nrow(inter_chr_svs)) {
        sv <- inter_chr_svs[i, ]
        breakpoints <- rbind(breakpoints, data.frame(
            node_id = paste(sv$chrom1, sv$pos1, sep = ":"),
            chrom = as.character(sv$chrom1),
            pos = as.numeric(sv$pos1),
            sv_index = i,
            stringsAsFactors = FALSE
        ))
        breakpoints <- rbind(breakpoints, data.frame(
            node_id = paste(sv$chrom2, sv$pos2, sep = ":"),
            chrom = as.character(sv$chrom2),
            pos = as.numeric(sv$pos2),
            sv_index = i,
            stringsAsFactors = FALSE
        ))
    }

    breakpoints <- breakpoints[!duplicated(breakpoints$node_id), ]
    nodes <- breakpoints$node_id

    edges <- data.frame(
        from = character(),
        to = character(),
        edge_type = character(),
        sv_index = integer(),
        distance = numeric(),
        pvalue = numeric(),
        stringsAsFactors = FALSE
    )

    # 1. Translocation edges
    for (i in 1:nrow(inter_chr_svs)) {
        sv <- inter_chr_svs[i, ]
        edges <- rbind(edges, data.frame(
            from = paste(sv$chrom1, sv$pos1, sep = ":"),
            to = paste(sv$chrom2, sv$pos2, sep = ":"),
            edge_type = "TRANSLOCATION",
            sv_index = i,
            distance = NA,
            pvalue = NA,
            stringsAsFactors = FALSE
        ))
    }

    # 2. Deletion bridge edges
    for (bridge in deletion_bridges) {
        sv_node <- paste(bridge$chrom, bridge$sv_pos, sep = ":")

        # Find other breakpoints in the deletion region
        in_region <- breakpoints[
            breakpoints$chrom == bridge$chrom &
            breakpoints$pos >= bridge$del_start &
            breakpoints$pos <= bridge$del_end &
            breakpoints$node_id != sv_node,
        ]

        for (j in seq_len(nrow(in_region))) {
            edges <- rbind(edges, data.frame(
                from = sv_node,
                to = in_region$node_id[j],
                edge_type = "DELETION_BRIDGE",
                sv_index = bridge$sv_index,
                distance = bridge$del_end - bridge$del_start,
                pvalue = NA,
                stringsAsFactors = FALSE
            ))
        }
    }

    # 3. Statistical adjacency edges (KEY ChainFinder feature)
    for (chr in unique(breakpoints$chrom)) {
        chr_bps <- breakpoints[breakpoints$chrom == chr, ]
        if (nrow(chr_bps) < 2) next

        chr_bps <- chr_bps[order(chr_bps$pos), ]

        # Get local mu
        mu <- if (chr %in% names(mu_estimates$chr_mu)) {
            mu_estimates$chr_mu[chr]
        } else {
            mu_estimates$global_mu
        }
        mu_per_bp <- mu / 1e6  # Convert to per-bp

        # Test all pairs
        for (i in 1:(nrow(chr_bps) - 1)) {
            for (j in (i + 1):nrow(chr_bps)) {
                dist <- chr_bps$pos[j] - chr_bps$pos[i]

                if (dist > max_distance) break

                # ChainFinder null model: P = 1 - (1 - 2*mu)^L
                # Probability of NO breakpoint in distance L
                log_p_none <- dist * log(1 - 2 * mu_per_bp)
                pvalue <- 1 - exp(log_p_none)

                if (pvalue < pvalue_threshold) {
                    edges <- rbind(edges, data.frame(
                        from = chr_bps$node_id[i],
                        to = chr_bps$node_id[j],
                        edge_type = "STAT_ADJACENCY",
                        sv_index = NA,
                        distance = dist,
                        pvalue = pvalue,
                        stringsAsFactors = FALSE
                    ))
                }
            }
        }
    }

    # Build adjacency list
    adj_list <- list()
    for (node in nodes) {
        adj_list[[node]] <- list(neighbors = character(), edge_types = character())
    }

    for (i in seq_len(nrow(edges))) {
        e <- edges[i, ]
        adj_list[[e$from]]$neighbors <- c(adj_list[[e$from]]$neighbors, e$to)
        adj_list[[e$from]]$edge_types <- c(adj_list[[e$from]]$edge_types, e$edge_type)
        adj_list[[e$to]]$neighbors <- c(adj_list[[e$to]]$neighbors, e$from)
        adj_list[[e$to]]$edge_types <- c(adj_list[[e$to]]$edge_types, e$edge_type)
    }

    list(nodes = nodes, edges = edges, adj_list = adj_list, breakpoints = breakpoints)
}

#' Find chains using DFS
#' @keywords internal
.cf_find_chains <- function(graph, min_chromosomes, min_translocations,
                            max_path_search, max_search_depth) {

    adj_list <- graph$adj_list
    all_chains <- list()
    found_sigs <- new.env(hash = TRUE)

    # Start from nodes with translocation edges
    start_nodes <- names(adj_list)[sapply(names(adj_list), function(n) {
        any(adj_list[[n]]$edge_types == "TRANSLOCATION")
    })]

    for (start in start_nodes) {
        # DFS using stack
        stack <- list(list(
            node = start,
            path = start,
            edge_types = character(),
            visited = start
        ))

        paths_found <- 0

        while (length(stack) > 0 && paths_found < max_path_search) {
            curr <- stack[[length(stack)]]
            stack <- stack[-length(stack)]

            # Check if valid chain
            if (length(curr$path) >= min_translocations) {
                chroms <- unique(sapply(curr$path, function(x) strsplit(x, ":")[[1]][1]))
                n_tra <- sum(curr$edge_types == "TRANSLOCATION")

                if (length(chroms) >= min_chromosomes && n_tra >= min_translocations) {
                    # Signature for deduplication
                    sig <- paste(sort(curr$path), collapse = "_")

                    if (!exists(sig, envir = found_sigs)) {
                        assign(sig, TRUE, envir = found_sigs)
                        paths_found <- paths_found + 1

                        all_chains[[length(all_chains) + 1]] <- list(
                            id = length(all_chains) + 1,
                            nodes = curr$path,
                            edge_types = curr$edge_types,
                            chromosomes = chroms,
                            n_chromosomes = length(chroms),
                            n_translocations = n_tra
                        )
                    }
                }
            }

            # Continue search if not too deep
            if (length(curr$path) < max_search_depth) {
                neighbors <- adj_list[[curr$node]]

                for (i in seq_along(neighbors$neighbors)) {
                    nbr <- neighbors$neighbors[i]
                    etype <- neighbors$edge_types[i]

                    if (!(nbr %in% curr$visited)) {
                        stack[[length(stack) + 1]] <- list(
                            node = nbr,
                            path = c(curr$path, nbr),
                            edge_types = c(curr$edge_types, etype),
                            visited = c(curr$visited, nbr)
                        )
                    }
                }
            }
        }
    }

    return(all_chains)
}

#' Evaluate a chain
#' @keywords internal
.cf_evaluate_chain <- function(chain, inter_chr_svs, cnv_df, graph, mu_estimates) {

    # CN stability
    cn_score <- 1.0
    max_cn_dev <- 0

    if (!is.null(cnv_df)) {
        cn_devs <- c()
        for (chr in chain$chromosomes) {
            chr_cn <- cnv_df[cnv_df$chrom == chr, ]
            if (nrow(chr_cn) > 0) {
                devs <- abs(chr_cn$total_cn - 2)
                cn_devs <- c(cn_devs, devs)
            }
        }
        if (length(cn_devs) > 0) {
            max_cn_dev <- max(cn_devs, na.rm = TRUE)
            cn_score <- exp(-mean(cn_devs, na.rm = TRUE) / 2)
        }
    }

    # Chain p-value (Fisher's method on adjacency p-values)
    adj_edges <- graph$edges[
        graph$edges$edge_type == "STAT_ADJACENCY" &
        !is.na(graph$edges$pvalue),
    ]

    pvalues <- c()
    for (i in 1:(length(chain$nodes) - 1)) {
        n1 <- chain$nodes[i]
        n2 <- chain$nodes[i + 1]

        edge_pv <- adj_edges$pvalue[
            (adj_edges$from == n1 & adj_edges$to == n2) |
            (adj_edges$from == n2 & adj_edges$to == n1)
        ]

        if (length(edge_pv) > 0) {
            pvalues <- c(pvalues, edge_pv[1])
        }
    }

    chain_pvalue <- NA
    if (length(pvalues) > 0) {
        pvalues <- pvalues[pvalues > 0 & pvalues < 1]
        if (length(pvalues) > 0) {
            fisher_stat <- -2 * sum(log(pvalues))
            chain_pvalue <- pchisq(fisher_stat, 2 * length(pvalues), lower.tail = FALSE)
        }
    }

    # Edge composition
    n_del_bridge <- sum(chain$edge_types == "DELETION_BRIDGE")
    n_stat_adj <- sum(chain$edge_types == "STAT_ADJACENCY")

    summary <- data.frame(
        chain_id = chain$id,
        n_chromosomes = chain$n_chromosomes,
        chromosomes = paste(chain$chromosomes, collapse = ","),
        n_translocations = chain$n_translocations,
        n_deletion_bridges = n_del_bridge,
        n_stat_adjacency = n_stat_adj,
        cn_stability_score = cn_score,
        max_cn_deviation = max_cn_dev,
        chain_pvalue = chain_pvalue,
        stringsAsFactors = FALSE
    )

    list(summary = summary, chain = chain)
}

#' Classify chains
#' @keywords internal
.cf_classify_chains <- function(summary_df) {

    sapply(1:nrow(summary_df), function(i) {
        row <- summary_df[i, ]

        meets_chr <- row$n_chromosomes >= 3
        meets_tra <- row$n_translocations >= 3
        meets_cn <- row$cn_stability_score >= 0.7
        meets_stat <- !is.na(row$chain_pvalue) && row$chain_pvalue < 0.01
        has_bridges <- row$n_deletion_bridges > 0

        n_met <- sum(c(meets_chr, meets_tra, meets_cn, meets_stat, has_bridges))

        if (meets_stat && n_met >= 4) {
            "Likely chromoplexy"
        } else if (n_met >= 3) {
            "Possible chromoplexy"
        } else if (n_met >= 2) {
            "Unlikely chromoplexy"
        } else {
            "Not chromoplexy"
        }
    })
}

#' Create empty result
#' @keywords internal
.cf_empty_result <- function() {
    result <- list(
        chains = list(),
        chain_details = list(),
        summary = data.frame(),
        graph = NULL,
        deletion_bridges = list(),
        mu_estimates = NULL,
        total_chains = 0,
        likely_chromoplexy = 0,
        possible_chromoplexy = 0,
        method = "ChainFinder",
        version = "1.0"
    )
    class(result) <- c("chromoplexy_chainfinder", "list")
    result
}

#' Print method
#' @export
print.chromoplexy_chainfinder <- function(x, ...) {
    cat("\n========================================================================\n")
    cat("ChainFinder-style Chromoplexy Detection Results\n")
    cat("========================================================================\n\n")

    cat(sprintf("Total chains: %d\n", x$total_chains))
    cat(sprintf("  Likely chromoplexy:   %d\n", x$likely_chromoplexy))
    cat(sprintf("  Possible chromoplexy: %d\n", x$possible_chromoplexy))

    if (!is.null(x$mu_estimates)) {
        cat(sprintf("\nBackground rate (mu): %.4f breakpoints/Mb\n", x$mu_estimates$global_mu))
    }

    if (!is.null(x$graph)) {
        cat(sprintf("\nGraph: %d nodes, %d edges\n",
                   length(x$graph$nodes), nrow(x$graph$edges)))
        cat(sprintf("  Statistical adjacency edges: %d\n",
                   sum(x$graph$edges$edge_type == "STAT_ADJACENCY")))
    }

    if (x$total_chains > 0 && nrow(x$summary) > 0) {
        cat("\nChain Summary:\n")
        print(x$summary[, c("chain_id", "n_chromosomes", "n_translocations",
                           "chain_pvalue", "classification")])
    }

    cat("\n========================================================================\n")
    invisible(x)
}
