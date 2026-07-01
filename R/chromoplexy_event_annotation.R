#' Collapse redundant chromoplexy chains into event-level components
#'
#' Chromoplexy path search can enumerate multiple overlapping chains from the
#' same underlying rearrangement. This helper collapses selected chains into
#' connected event components based on shared breakpoint nodes and reports
#' event-level QC, SV support, breakpoints, and optional gene impacts.
#'
#' @param chromoplexy_result A result from \code{detect_chromoplexy()} or a
#'   \code{chromoanagenesis} result containing \code{$chromoplexy}.
#' @param classifications Character vector of chain classifications to collapse.
#'   Defaults to likely chromoplexy chains.
#' @param gene_granges Optional GRanges object with gene annotations.
#' @param breakpoint_padding Padding in bp around each breakpoint for gene
#'   overlap annotation.
#' @param driver_genes Optional vector of driver gene symbols. If NULL, bundled
#'   driver genes are used when available.
#' @return A list with event summaries, chain-to-event mapping, breakpoints, and
#'   optional gene impact tables.
#' @export
collapse_chromoplexy_chains <- function(chromoplexy_result,
                                        classifications = c("Likely chromoplexy"),
                                        gene_granges = NULL,
                                        breakpoint_padding = 1000,
                                        driver_genes = NULL) {
    if (!is.null(chromoplexy_result$chromoplexy)) {
        chromoplexy_result <- chromoplexy_result$chromoplexy
    }

    empty <- empty_collapsed_chromoplexy_events()
    summary_df <- chromoplexy_result$summary
    details <- chromoplexy_result$chain_details

    if (is.null(summary_df) || nrow(summary_df) == 0 ||
        is.null(details) || length(details) == 0) {
        return(empty)
    }
    if (!"chain_id" %in% colnames(summary_df) ||
        !"classification" %in% colnames(summary_df)) {
        return(empty)
    }

    keep_ids <- summary_df$chain_id[summary_df$classification %in% classifications]
    if (length(keep_ids) == 0) {
        return(empty)
    }

    detail_ids <- vapply(details, function(x) {
        if (!is.null(x$summary$chain_id)) as.integer(x$summary$chain_id[1]) else NA_integer_
    }, integer(1))
    details <- details[detail_ids %in% keep_ids]
    detail_ids <- detail_ids[detail_ids %in% keep_ids]

    if (length(details) == 0) {
        return(empty)
    }

    all_nodes <- sort(unique(unlist(lapply(details, function(x) {
        as.character(x$chain$nodes)
    }), use.names = FALSE)))
    if (length(all_nodes) == 0) {
        return(empty)
    }

    parent <- stats::setNames(all_nodes, all_nodes)
    find_root <- function(x) {
        while (!identical(parent[[x]], x)) {
            parent[[x]] <<- parent[[parent[[x]]]]
            x <- parent[[x]]
        }
        x
    }
    union_nodes <- function(a, b) {
        ra <- find_root(a)
        rb <- find_root(b)
        if (!identical(ra, rb)) parent[[rb]] <<- ra
    }

    for (detail in details) {
        nodes <- as.character(detail$chain$nodes)
        if (length(nodes) > 1) {
            for (i in seq_len(length(nodes) - 1)) {
                union_nodes(nodes[i], nodes[i + 1])
            }
        }
    }

    root_key <- vapply(all_nodes, find_root, character(1))
    names(root_key) <- all_nodes
    chain_root <- vapply(details, function(x) {
        roots <- unique(root_key[as.character(x$chain$nodes)])
        roots[1]
    }, character(1))

    raw_roots <- sort(unique(chain_root))
    event_rows <- lapply(raw_roots, function(root) {
        cids <- detail_ids[chain_root == root]
        edetails <- details[chain_root == root]
        s <- summary_df[match(cids, summary_df$chain_id), , drop = FALSE]

        event_nodes <- sort(unique(unlist(lapply(edetails, function(x) {
            as.character(x$chain$nodes)
        }), use.names = FALSE)))
        sv_indices <- sort(unique(unlist(lapply(edetails, function(x) {
            as.integer(x$chain$sv_indices)
        }), use.names = FALSE)))
        sv_df <- unique(do.call(rbind, lapply(edetails, function(x) x$SVs)))
        sv_ids <- if (!is.null(sv_df) && nrow(sv_df) > 0 && "sv_id" %in% colnames(sv_df)) {
            sort(unique(as.character(sv_df$sv_id)))
        } else {
            as.character(sv_indices)
        }
        chroms <- sort(unique(vapply(strsplit(event_nodes, ":", fixed = TRUE), `[`, character(1), 1)))

        representative_idx <- which.max(replace_na_numeric(s$combined_score, -Inf))
        score <- score_collapsed_chromoplexy_event(s, event_nodes, sv_indices)
        data.frame(
            raw_component_id = root,
            n_chains = length(cids),
            representative_chain_id = s$chain_id[representative_idx],
            representative_combined_score = max_or_na(s$combined_score),
            mean_combined_score = mean_or_na(s$combined_score),
            n_unique_svs = length(sv_indices),
            n_breakpoints = length(event_nodes),
            n_chromosomes = length(chroms),
            chromosomes_involved = paste(chroms, collapse = ","),
            n_translocations = sum(s$n_translocations, na.rm = TRUE),
            has_cycle = any(s$is_cycle, na.rm = TRUE),
            evidence_mode = paste(sort(unique(s$evidence_mode)), collapse = ","),
            event_qc_score = score$event_qc_score,
            event_confidence = score$event_confidence,
            chain_ids = paste(sort(cids), collapse = ","),
            sv_ids = paste(sv_ids, collapse = ","),
            stringsAsFactors = FALSE
        )
    })

    event_summary <- do.call(rbind, event_rows)
    event_summary <- event_summary[order(
        -event_summary$event_qc_score,
        -event_summary$n_unique_svs,
        event_summary$representative_chain_id
    ), , drop = FALSE]
    event_summary$collapsed_event_id <- sprintf("CE%03d", seq_len(nrow(event_summary)))
    event_summary <- event_summary[, c(
        "collapsed_event_id",
        setdiff(colnames(event_summary), "collapsed_event_id")
    )]

    root_to_event <- stats::setNames(event_summary$collapsed_event_id, event_summary$raw_component_id)
    chain_to_event <- data.frame(
        chain_id = detail_ids,
        raw_component_id = chain_root,
        collapsed_event_id = unname(root_to_event[chain_root]),
        stringsAsFactors = FALSE
    )
    chain_to_event <- merge(
        chain_to_event,
        summary_df[, intersect(c(
            "chain_id", "classification", "combined_score", "pvalue", "fdr",
            "n_chromosomes", "n_translocations", "chromosomes_involved",
            "evidence_mode"
        ), colnames(summary_df)), drop = FALSE],
        by = "chain_id",
        all.x = TRUE
    )
    chain_to_event <- chain_to_event[order(chain_to_event$collapsed_event_id,
                                           -chain_to_event$combined_score,
                                           chain_to_event$chain_id), , drop = FALSE]

    event_breakpoints <- build_collapsed_event_breakpoints(event_summary, details,
                                                           detail_ids, chain_to_event)

    gene_tables <- annotate_collapsed_event_genes(
        event_summary = event_summary,
        event_breakpoints = event_breakpoints,
        gene_granges = gene_granges,
        breakpoint_padding = breakpoint_padding,
        driver_genes = driver_genes
    )

    if (nrow(gene_tables$event_gene_summary) > 0) {
        event_summary <- merge(event_summary, gene_tables$event_gene_summary,
                               by = "collapsed_event_id", all.x = TRUE)
        event_summary <- event_summary[order(event_summary$collapsed_event_id), , drop = FALSE]
    } else {
        event_summary$n_genes <- 0L
        event_summary$genes <- ""
        event_summary$n_driver_genes <- 0L
        event_summary$driver_genes <- ""
    }
    event_summary$n_genes[is.na(event_summary$n_genes)] <- 0L
    event_summary$genes[is.na(event_summary$genes)] <- ""
    event_summary$n_driver_genes[is.na(event_summary$n_driver_genes)] <- 0L
    event_summary$driver_genes[is.na(event_summary$driver_genes)] <- ""

    list(
        event_summary = event_summary,
        chain_to_event = chain_to_event,
        event_breakpoints = event_breakpoints,
        gene_detail = gene_tables$gene_detail,
        event_gene_summary = gene_tables$event_gene_summary,
        gene_event_summary = gene_tables$gene_event_summary,
        parameters = list(
            classifications = classifications,
            breakpoint_padding = breakpoint_padding,
            gene_annotation = !is.null(gene_granges)
        )
    )
}

#' Annotate collapsed chromoplexy events with gene overlaps
#'
#' @param chromoplexy_result A chromoplexy or chromoanagenesis result object.
#' @param gene_granges GRanges object with gene annotations.
#' @param classifications Chain classifications to collapse before annotation.
#' @param breakpoint_padding Padding in bp around each breakpoint.
#' @param driver_genes Optional vector of driver gene symbols.
#' @return The same list returned by \code{collapse_chromoplexy_chains()}, with
#'   gene tables populated.
#' @export
annotate_chromoplexy_events <- function(chromoplexy_result,
                                        gene_granges,
                                        classifications = c("Likely chromoplexy"),
                                        breakpoint_padding = 1000,
                                        driver_genes = NULL) {
    collapse_chromoplexy_chains(
        chromoplexy_result = chromoplexy_result,
        classifications = classifications,
        gene_granges = gene_granges,
        breakpoint_padding = breakpoint_padding,
        driver_genes = driver_genes
    )
}

empty_collapsed_chromoplexy_events <- function() {
    list(
        event_summary = data.frame(),
        chain_to_event = data.frame(),
        event_breakpoints = data.frame(),
        gene_detail = data.frame(),
        event_gene_summary = data.frame(),
        gene_event_summary = data.frame(),
        parameters = list()
    )
}

score_collapsed_chromoplexy_event <- function(chain_summary, event_nodes, sv_indices) {
    combined <- chain_summary$combined_score
    combined <- combined[!is.na(combined)]
    evidence_score <- if (length(combined) > 0) max(combined) else 0
    chain_support_score <- min(nrow(chain_summary) / 3, 1)
    sv_support_score <- min(length(unique(sv_indices)) / 6, 1)
    chrom_support_score <- min(length(unique(vapply(strsplit(event_nodes, ":", fixed = TRUE),
                                                   `[`, character(1), 1))) / 5, 1)
    cycle_score <- if (any(chain_summary$is_cycle, na.rm = TRUE)) 1 else 0

    event_qc_score <- (
        evidence_score * 0.45 +
        sv_support_score * 0.20 +
        chrom_support_score * 0.20 +
        chain_support_score * 0.10 +
        cycle_score * 0.05
    )
    event_qc_score <- max(0, min(1, event_qc_score))

    event_confidence <- if (event_qc_score >= 0.75) {
        "High"
    } else if (event_qc_score >= 0.50) {
        "Moderate"
    } else {
        "Low"
    }

    list(event_qc_score = event_qc_score, event_confidence = event_confidence)
}

build_collapsed_event_breakpoints <- function(event_summary, details, detail_ids,
                                             chain_to_event) {
    rows <- list()
    for (i in seq_along(details)) {
        detail <- details[[i]]
        chain_id <- detail_ids[i]
        event_id <- chain_to_event$collapsed_event_id[match(chain_id, chain_to_event$chain_id)]
        nodes <- as.character(detail$chain$nodes)
        if (length(nodes) == 0 || is.na(event_id)) next
        for (node in nodes) {
            parsed <- strsplit(node, ":", fixed = TRUE)[[1]]
            if (length(parsed) < 2) next
            rows[[length(rows) + 1]] <- data.frame(
                collapsed_event_id = event_id,
                chain_id = chain_id,
                breakpoint_id = node,
                chrom = parsed[1],
                pos = suppressWarnings(as.integer(parsed[2])),
                stringsAsFactors = FALSE
            )
        }
    }
    if (length(rows) == 0) {
        return(data.frame())
    }
    bp <- unique(do.call(rbind, rows))
    bp <- bp[order(bp$collapsed_event_id, bp$chrom, bp$pos, bp$chain_id), , drop = FALSE]
    bp
}

annotate_collapsed_event_genes <- function(event_summary, event_breakpoints,
                                           gene_granges = NULL,
                                           breakpoint_padding = 1000,
                                           driver_genes = NULL) {
    empty <- list(
        gene_detail = data.frame(),
        event_gene_summary = data.frame(),
        gene_event_summary = data.frame()
    )
    if (is.null(gene_granges) || nrow(event_breakpoints) == 0) {
        return(empty)
    }
    if (!inherits(gene_granges, "GRanges")) {
        stop("gene_granges must be a GRanges object")
    }

    gene_cols <- colnames(GenomicRanges::mcols(gene_granges))
    symbol_col <- first_existing_column(gene_cols, c("symbol", "gene_name", "hgnc_symbol", "gene_id"))
    gene_id_col <- first_existing_column(gene_cols, c("gene_id", "ensembl_gene_id", symbol_col))
    if (is.na(symbol_col)) {
        stop("gene_granges must contain a symbol, gene_name, hgnc_symbol, or gene_id column")
    }

    if (is.null(driver_genes)) {
        driver_genes <- get_default_drivers()
    }
    cancer_genes <- load_cancer_genes()
    gene_type <- function(symbol) {
        if (is.null(cancer_genes) || !"Hugo_Symbol" %in% colnames(cancer_genes)) {
            return("Unknown")
        }
        idx <- match(symbol, cancer_genes$Hugo_Symbol)
        if (is.na(idx) || !"Gene_Type" %in% colnames(cancer_genes)) "Unknown" else cancer_genes$Gene_Type[idx]
    }

    target_has_chr <- any(grepl("^chr", as.character(GenomeInfoDb::seqlevels(gene_granges))))
    if (!target_has_chr && length(GenomeInfoDb::seqlevels(gene_granges)) == 0) {
        target_has_chr <- any(grepl("^chr", as.character(GenomicRanges::seqnames(gene_granges))))
    }
    norm_chrom <- function(x) {
        if (target_has_chr) {
            ifelse(grepl("^chr", x), x, paste0("chr", x))
        } else {
            gsub("^chr", "", x)
        }
    }

    bp <- event_breakpoints
    bp$annot_chrom <- norm_chrom(bp$chrom)
    bp$start <- pmax(1L, as.integer(bp$pos) - as.integer(breakpoint_padding))
    bp$end <- as.integer(bp$pos) + as.integer(breakpoint_padding)

    bp_gr <- GenomicRanges::GRanges(
        seqnames = bp$annot_chrom,
        ranges = IRanges::IRanges(start = bp$start, end = bp$end)
    )
    hits <- GenomicRanges::findOverlaps(bp_gr, gene_granges, ignore.strand = TRUE)
    if (length(hits) == 0) {
        return(empty)
    }

    qh <- S4Vectors::queryHits(hits)
    gh <- S4Vectors::subjectHits(hits)
    gm <- GenomicRanges::mcols(gene_granges)
    symbols <- as.character(gm[[symbol_col]][gh])
    gene_ids <- as.character(gm[[gene_id_col]][gh])
    distances <- pmax(0L, pmax(GenomicRanges::start(gene_granges)[gh] - bp$pos[qh],
                               bp$pos[qh] - GenomicRanges::end(gene_granges)[gh]))

    gene_detail <- data.frame(
        collapsed_event_id = bp$collapsed_event_id[qh],
        chain_id = bp$chain_id[qh],
        breakpoint_id = bp$breakpoint_id[qh],
        chrom = bp$chrom[qh],
        pos = bp$pos[qh],
        gene_id = gene_ids,
        symbol = symbols,
        gene_type = vapply(symbols, gene_type, character(1)),
        is_driver = symbols %in% driver_genes,
        gene_chrom = as.character(GenomicRanges::seqnames(gene_granges)[gh]),
        gene_start = GenomicRanges::start(gene_granges)[gh],
        gene_end = GenomicRanges::end(gene_granges)[gh],
        distance_to_breakpoint = distances,
        stringsAsFactors = FALSE
    )
    gene_detail <- unique(gene_detail)
    gene_detail <- gene_detail[order(gene_detail$collapsed_event_id,
                                     gene_detail$symbol,
                                     gene_detail$breakpoint_id), , drop = FALSE]

    event_gene_summary <- summarize_event_genes(gene_detail)
    gene_event_summary <- summarize_gene_events(gene_detail)

    list(
        gene_detail = gene_detail,
        event_gene_summary = event_gene_summary,
        gene_event_summary = gene_event_summary
    )
}

summarize_event_genes <- function(gene_detail) {
    if (nrow(gene_detail) == 0) return(data.frame())
    rows <- lapply(sort(unique(gene_detail$collapsed_event_id)), function(event_id) {
        sub <- gene_detail[gene_detail$collapsed_event_id == event_id, , drop = FALSE]
        drivers <- sort(unique(sub$symbol[sub$is_driver]))
        data.frame(
            collapsed_event_id = event_id,
            n_genes = length(unique(sub$gene_id)),
            genes = paste(sort(unique(sub$symbol)), collapse = ","),
            n_driver_genes = length(drivers),
            driver_genes = paste(drivers, collapse = ","),
            n_gene_overlapping_breakpoints = length(unique(sub$breakpoint_id)),
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, rows)
}

summarize_gene_events <- function(gene_detail) {
    if (nrow(gene_detail) == 0) return(data.frame())
    rows <- lapply(sort(unique(gene_detail$gene_id)), function(gid) {
        sub <- gene_detail[gene_detail$gene_id == gid, , drop = FALSE]
        data.frame(
            gene_id = gid,
            symbol = sub$symbol[1],
            gene_type = sub$gene_type[1],
            is_driver = any(sub$is_driver),
            n_collapsed_events = length(unique(sub$collapsed_event_id)),
            collapsed_event_ids = paste(sort(unique(sub$collapsed_event_id)), collapse = ","),
            n_breakpoints = length(unique(sub$breakpoint_id)),
            chromosomes = paste(sort(unique(sub$chrom)), collapse = ","),
            stringsAsFactors = FALSE
        )
    })
    out <- do.call(rbind, rows)
    out[order(-out$n_collapsed_events, -out$is_driver, out$symbol), , drop = FALSE]
}

first_existing_column <- function(columns, candidates) {
    hit <- candidates[candidates %in% columns]
    if (length(hit) == 0) NA_character_ else hit[1]
}

replace_na_numeric <- function(x, replacement) {
    x[is.na(x)] <- replacement
    x
}

max_or_na <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) NA_real_ else max(x)
}

mean_or_na <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) NA_real_ else mean(x)
}
