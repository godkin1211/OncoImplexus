#' Plot chromoplexy chains
#'
#' Visualizes chromoplexy translocation chains across chromosomes.
#'
#' @param chromoplexy_result Result from detect_chromoplexy()
#' @param chain_id Optional: specific chain ID to plot (default: plot all)
#' @param sample_name Sample name for plot title
#' @param genome Reference genome ("hg19" or "hg38")
#' @param show_cn Whether to show copy number tracks (default: FALSE)
#' @param CNV.sample CNV data for copy number tracks
#' @return A ggplot object or list of ggplot objects
#' @export
plot_chromoplexy <- function(chromoplexy_result,
                            chain_id = NULL,
                            sample_name = "",
                            genome = "hg19",
                            show_cn = FALSE,
                            CNV.sample = NULL) {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required for plotting chromoplexy.")
    }

    if (chromoplexy_result$total_chains == 0) {
        message("No chromoplexy chains to plot.")
        return(NULL)
    }

    # Get chromosome lengths
    if (genome == "hg38") {
        chr_lengths <- info_mappa_hg38
    } else {
        chr_lengths <- info_mappa
    }

    # Select chains to plot
    if (!is.null(chain_id)) {
        chains_to_plot <- list(chromoplexy_result$chain_details[[chain_id]])
    } else {
        chains_to_plot <- chromoplexy_result$chain_details
    }

    plots <- list()

    for (i in 1:length(chains_to_plot)) {
        chain_detail <- chains_to_plot[[i]]
        chain <- chain_detail$chain
        chain_SVs <- chain_detail$SVs

        # Create plot
        p <- create_chromoplexy_chain_plot(
            chain = chain,
            chain_SVs = chain_SVs,
            chr_lengths = chr_lengths,
            sample_name = sample_name,
            show_cn = show_cn,
            CNV.sample = CNV.sample
        )

        plots[[i]] <- p
    }

    if (length(plots) == 1) {
        return(plots[[1]])
    } else {
        return(plots)
    }
}

#' Plot a collapsed chromoplexy event graph
#'
#' Visualizes one event-level chromoplexy component from
#' \code{collapse_chromoplexy_chains()}. Nodes are breakpoints, edges are the
#' chain graph edge types, and optional labels show genes overlapping each
#' breakpoint.
#'
#' @param chromoplexy_result Result from \code{detect_chromoplexy()} or
#'   \code{detect_chromoanagenesis()}.
#' @param event_id Collapsed event ID such as "CE001". If NULL, the first event
#'   in the summary table is plotted.
#' @param sample_name Optional sample label for the plot title.
#' @param genome Reference genome ("hg19" or "hg38").
#' @param show_genes Whether to show gene labels at breakpoints.
#' @param max_gene_labels Maximum number of breakpoint gene labels to draw.
#' @return A ggplot object.
#' @export
plot_collapsed_chromoplexy_event <- function(chromoplexy_result,
                                             event_id = NULL,
                                             sample_name = "",
                                             genome = "hg19",
                                             show_genes = TRUE,
                                             max_gene_labels = 20) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required for plotting.")
    }
    if (!is.null(chromoplexy_result$chromoplexy)) {
        chromoplexy_result <- chromoplexy_result$chromoplexy
    }
    ce <- chromoplexy_result$collapsed_events
    if (is.null(ce) || is.null(ce$event_summary) || nrow(ce$event_summary) == 0) {
        stop("No collapsed chromoplexy events available")
    }

    events <- ce$event_summary
    if (is.null(event_id)) {
        event_id <- events$collapsed_event_id[1]
    } else if (is.numeric(event_id)) {
        event_id <- events$collapsed_event_id[event_id[1]]
    }
    if (!event_id %in% events$collapsed_event_id) {
        stop("event_id not found in collapsed event summary")
    }

    bp <- ce$event_breakpoints[ce$event_breakpoints$collapsed_event_id == event_id, , drop = FALSE]
    bp <- unique(bp[, c("breakpoint_id", "chrom", "pos"), drop = FALSE])
    if (nrow(bp) == 0) {
        stop("Selected collapsed event has no breakpoint table")
    }

    edge_df <- build_collapsed_event_edge_table(chromoplexy_result, ce, event_id)
    chroms <- sort(unique(c(bp$chrom, edge_df$chrom1, edge_df$chrom2)))
    chr_lengths <- if (genome == "hg38") info_mappa_hg38 else info_mappa
    chr_info <- data.frame(
        chrom = chroms,
        length = vapply(chroms, function(chr) {
            idx <- which(chr_lengths$V1 == chr)
            if (length(idx) > 0) chr_lengths$tot[idx[1]] else max(bp$pos[bp$chrom == chr], na.rm = TRUE)
        }, numeric(1)),
        y = seq(length(chroms), 1),
        stringsAsFactors = FALSE
    )

    bp <- merge(bp, chr_info[, c("chrom", "y")], by = "chrom", all.x = TRUE)
    if (nrow(edge_df) > 0) {
        edge_df <- merge(edge_df, chr_info[, c("chrom", "y")], by.x = "chrom1", by.y = "chrom", all.x = TRUE)
        colnames(edge_df)[colnames(edge_df) == "y"] <- "y1"
        edge_df <- merge(edge_df, chr_info[, c("chrom", "y")], by.x = "chrom2", by.y = "chrom", all.x = TRUE)
        colnames(edge_df)[colnames(edge_df) == "y"] <- "y2"
    }

    p <- ggplot2::ggplot()
    p <- p + ggplot2::geom_segment(
        data = chr_info,
        ggplot2::aes(x = 0, xend = length, y = y, yend = y),
        color = "gray75",
        linewidth = 3,
        lineend = "round"
    )
    p <- p + ggplot2::geom_text(
        data = chr_info,
        ggplot2::aes(x = 0, y = y, label = chrom),
        hjust = 1.2,
        size = 3.5
    )

    if (nrow(edge_df) > 0) {
        same_chr <- edge_df[edge_df$chrom1 == edge_df$chrom2, , drop = FALSE]
        inter_chr <- edge_df[edge_df$chrom1 != edge_df$chrom2, , drop = FALSE]
        if (nrow(same_chr) > 0) {
            p <- p + ggplot2::geom_segment(
                data = same_chr,
                ggplot2::aes(x = pos1, xend = pos2, y = y1 + 0.08, yend = y2 + 0.08, color = edge_type),
                linewidth = 0.8,
                alpha = 0.8
            )
        }
        if (nrow(inter_chr) > 0) {
            p <- p + ggplot2::geom_curve(
                data = inter_chr,
                ggplot2::aes(x = pos1, xend = pos2, y = y1, yend = y2, color = edge_type),
                curvature = 0.18,
                linewidth = 0.8,
                alpha = 0.75
            )
        }
    }

    p <- p + ggplot2::geom_point(
        data = bp,
        ggplot2::aes(x = pos, y = y),
        color = "black",
        fill = "#F4D35E",
        shape = 21,
        size = 2.6,
        stroke = 0.6
    )

    if (show_genes && !is.null(ce$gene_detail) && nrow(ce$gene_detail) > 0) {
        labels <- collapsed_event_gene_labels(ce$gene_detail, event_id, max_gene_labels)
        labels <- merge(labels, bp, by = "breakpoint_id", all.x = TRUE)
        if (nrow(labels) > 0) {
            p <- p + ggplot2::geom_text(
                data = labels,
                ggplot2::aes(x = pos, y = y + 0.25, label = label),
                size = 3,
                angle = 25,
                hjust = 0
            )
        }
    }

    event_row <- events[events$collapsed_event_id == event_id, , drop = FALSE]
    subtitle <- sprintf(
        "%s | %s | %d SVs, %d breakpoints, %d chromosomes",
        event_row$event_confidence[1],
        round(event_row$event_qc_score[1], 3),
        event_row$n_unique_svs[1],
        event_row$n_breakpoints[1],
        event_row$n_chromosomes[1]
    )

    p + ggplot2::scale_color_manual(
        values = c(
            TRANSLOCATION = "#C44536",
            ADJACENCY = "#1976D2",
            DELETION_BRIDGE = "#6A4C93"
        ),
        drop = FALSE
    ) +
        ggplot2::labs(
            title = paste(trimws(sample_name), "Collapsed Chromoplexy Event", event_id),
            subtitle = subtitle,
            x = "Genomic position (bp)",
            y = "",
            color = "Edge type"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            panel.grid.major.y = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(face = "bold")
        )
}

build_collapsed_event_edge_table <- function(chromoplexy_result, collapsed_events, event_id) {
    chain_map <- collapsed_events$chain_to_event
    if (is.null(chain_map) || nrow(chain_map) == 0) return(data.frame())
    chain_ids <- chain_map$chain_id[chain_map$collapsed_event_id == event_id]
    rows <- list()
    for (detail in chromoplexy_result$chain_details) {
        if (!detail$summary$chain_id %in% chain_ids) next
        nodes <- as.character(detail$chain$nodes)
        edge_types <- as.character(detail$chain$edge_types)
        if (length(nodes) < 2 || length(edge_types) == 0) next
        for (i in seq_len(min(length(edge_types), length(nodes) - 1))) {
            from <- parse_breakpoint_node(nodes[i])
            to <- parse_breakpoint_node(nodes[i + 1])
            if (is.null(from) || is.null(to)) next
            rows[[length(rows) + 1]] <- data.frame(
                from = nodes[i],
                to = nodes[i + 1],
                chrom1 = from$chrom,
                pos1 = from$pos,
                chrom2 = to$chrom,
                pos2 = to$pos,
                edge_type = edge_types[i],
                stringsAsFactors = FALSE
            )
        }
    }
    if (length(rows) == 0) return(data.frame())
    unique(do.call(rbind, rows))
}

parse_breakpoint_node <- function(node) {
    parsed <- strsplit(as.character(node), ":", fixed = TRUE)[[1]]
    if (length(parsed) < 2) return(NULL)
    list(chrom = parsed[1], pos = suppressWarnings(as.integer(parsed[2])))
}

collapsed_event_gene_labels <- function(gene_detail, event_id, max_gene_labels) {
    sub <- gene_detail[gene_detail$collapsed_event_id == event_id, , drop = FALSE]
    if (nrow(sub) == 0) return(data.frame())
    rows <- lapply(sort(unique(sub$breakpoint_id)), function(bp) {
        genes <- sort(unique(sub$symbol[sub$breakpoint_id == bp]))
        genes <- genes[nzchar(genes)]
        if (length(genes) == 0) return(NULL)
        data.frame(
            breakpoint_id = bp,
            label = paste(head(genes, 3), collapse = ","),
            stringsAsFactors = FALSE
        )
    })
    out <- bind_rows_or_empty(rows)
    if (nrow(out) > max_gene_labels) out <- out[seq_len(max_gene_labels), , drop = FALSE]
    out
}


#' Create a single chromoplexy chain plot
#'
#' @param chain Chain object
#' @param chain_SVs SVs in chain
#' @param chr_lengths Chromosome length data
#' @param sample_name Sample name
#' @param show_cn Show copy number
#' @param CNV.sample CNV data
#' @return ggplot object
#' @keywords internal
create_chromoplexy_chain_plot <- function(chain,
                                         chain_SVs,
                                         chr_lengths,
                                         sample_name,
                                         show_cn,
                                         CNV.sample) {

    # Prepare chromosome layout
    chroms <- chain$chromosomes
    n_chroms <- length(chroms)

    # Get chromosome lengths and positions
    chr_info <- data.frame(
        chrom = chroms,
        length = sapply(chroms, function(chr) {
            idx <- which(chr_lengths$V1 == chr)
            if (length(idx) > 0) chr_lengths$tot[idx] else 1e8
        }),
        y_pos = seq(n_chroms, 1, -1),
        stringsAsFactors = FALSE
    )

    # Prepare translocation arcs
    arc_data <- prepare_translocation_arcs(chain_SVs, chr_info)

    # Create base plot
    p <- ggplot2::ggplot()

    # Draw chromosome bars
    for (i in 1:nrow(chr_info)) {
        chr_name <- chr_info$chrom[i]
        chr_len <- chr_info$length[i]
        y_pos <- chr_info$y_pos[i]

        p <- p + ggplot2::geom_rect(
            ggplot2::aes(xmin = 0, xmax = chr_len,
                        ymin = y_pos - 0.3, ymax = y_pos + 0.3),
            fill = "gray80", color = "black", size = 0.5
        )

        # Add chromosome label
        p <- p + ggplot2::annotate(
            "text",
            x = -chr_len * 0.05,
            y = y_pos,
            label = chr_name,
            hjust = 1,
            size = 4
        )
    }

    # Draw translocation arcs
    if (nrow(arc_data) > 0) {
        p <- p + ggplot2::geom_curve(
            data = arc_data,
            ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
            curvature = 0.3,
            arrow = ggplot2::arrow(length = ggplot2::unit(0.02, "npc")),
            color = "red",
            size = 1.2,
            alpha = 0.7
        )

        # Mark breakpoints
        breakpoint_data <- rbind(
            data.frame(x = arc_data$x1, y = arc_data$y1),
            data.frame(x = arc_data$x2, y = arc_data$y2)
        )

        p <- p + ggplot2::geom_point(
            data = breakpoint_data,
            ggplot2::aes(x = x, y = y),
            color = "darkred",
            size = 2
        )
    }

    # Add title and labels
    title_text <- sprintf("%s - Chromoplexy Chain %d",
                         sample_name, chain$id)
    subtitle_text <- sprintf("%d chromosomes, %d translocations%s",
                            chain$n_chromosomes,
                            chain$n_translocations,
                            if (chain$is_cycle) " (cycle)" else "")

    p <- p + ggplot2::labs(
        title = title_text,
        subtitle = subtitle_text,
        x = "Genomic Position (bp)",
        y = ""
    )

    # Theme
    p <- p + ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            panel.grid.major.y = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(size = 14, face = "bold"),
            plot.subtitle = ggplot2::element_text(size = 11)
        )

    # Adjust plot limits
    max_chr_len <- max(chr_info$length)
    p <- p + ggplot2::xlim(-max_chr_len * 0.15, max_chr_len * 1.05)

    return(p)
}


#' Prepare translocation arc data for plotting
#'
#' @param chain_SVs SVs in chain
#' @param chr_info Chromosome position info
#' @return Data frame with arc coordinates
#' @keywords internal
prepare_translocation_arcs <- function(chain_SVs, chr_info) {

    if (nrow(chain_SVs) == 0) {
        return(data.frame(x1 = numeric(0), y1 = numeric(0),
                         x2 = numeric(0), y2 = numeric(0)))
    }

    arc_list <- list()

    for (i in 1:nrow(chain_SVs)) {
        sv <- chain_SVs[i, ]

        # Get y positions for chromosomes
        y1_idx <- which(chr_info$chrom == sv$chrom1)
        y2_idx <- which(chr_info$chrom == sv$chrom2)

        if (length(y1_idx) == 0 || length(y2_idx) == 0) next

        arc_list[[i]] <- data.frame(
            x1 = sv$pos1,
            y1 = chr_info$y_pos[y1_idx],
            x2 = sv$pos2,
            y2 = chr_info$y_pos[y2_idx]
        )
    }

    if (length(arc_list) > 0) {
        arc_data <- do.call(rbind, arc_list)
    } else {
        arc_data <- data.frame(x1 = numeric(0), y1 = numeric(0),
                              x2 = numeric(0), y2 = numeric(0))
    }

    return(arc_data)
}


#' Create circular plot of chromoplexy chain
#'
#' For chains that form cycles, create a circular visualization.
#'
#' @param chromoplexy_result Result from detect_chromoplexy()
#' @param chain_id Chain ID to plot
#' @param sample_name Sample name
#' @return A ggplot object
#' @export
plot_chromoplexy_circular <- function(chromoplexy_result,
                                     chain_id,
                                     sample_name = "") {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required for plotting.")
    }

    chain_detail <- chromoplexy_result$chain_details[[chain_id]]
    chain <- chain_detail$chain

    if (!chain$is_cycle) {
        warning("Chain is not a cycle. Use plot_chromoplexy() instead.")
        return(plot_chromoplexy(chromoplexy_result, chain_id, sample_name))
    }

    # Create circular layout
    chroms <- chain$chromosomes
    n_chroms <- length(chroms)

    # Assign angles to chromosomes
    angles <- seq(0, 2 * pi, length.out = n_chroms + 1)[1:n_chroms]

    chr_positions <- data.frame(
        chrom = chroms,
        angle = angles,
        x = cos(angles),
        y = sin(angles),
        stringsAsFactors = FALSE
    )

    # Prepare arc data
    chain_SVs <- chain_detail$SVs
    arc_data <- prepare_circular_arcs(chain_SVs, chr_positions)

    # Create plot
    p <- ggplot2::ggplot()

    # Draw chromosome points
    p <- p + ggplot2::geom_point(
        data = chr_positions,
        ggplot2::aes(x = x, y = y),
        size = 10,
        color = "steelblue"
    )

    # Add chromosome labels
    p <- p + ggplot2::geom_text(
        data = chr_positions,
        ggplot2::aes(x = x * 1.2, y = y * 1.2, label = chrom),
        size = 5
    )

    # Draw translocation connections
    if (nrow(arc_data) > 0) {
        p <- p + ggplot2::geom_segment(
            data = arc_data,
            ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
            arrow = ggplot2::arrow(length = ggplot2::unit(0.02, "npc")),
            color = "red",
            size = 1,
            alpha = 0.6
        )
    }

    # Add title
    title_text <- sprintf("%s - Chromoplexy Cycle Chain %d",
                         sample_name, chain$id)

    p <- p + ggplot2::labs(
        title = title_text,
        subtitle = sprintf("%d chromosomes in cycle", n_chroms)
    )

    # Theme
    p <- p + ggplot2::theme_void() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
            plot.subtitle = ggplot2::element_text(size = 11, hjust = 0.5)
        )

    # Set equal aspect ratio
    p <- p + ggplot2::coord_fixed()

    return(p)
}


#' Prepare arc data for circular plot
#'
#' @param chain_SVs SVs in chain
#' @param chr_positions Chromosome positions in circular layout
#' @return Data frame with arc coordinates
#' @keywords internal
prepare_circular_arcs <- function(chain_SVs, chr_positions) {

    if (nrow(chain_SVs) == 0) {
        return(data.frame(x1 = numeric(0), y1 = numeric(0),
                         x2 = numeric(0), y2 = numeric(0)))
    }

    arc_list <- list()

    for (i in 1:nrow(chain_SVs)) {
        sv <- chain_SVs[i, ]

        idx1 <- which(chr_positions$chrom == sv$chrom1)
        idx2 <- which(chr_positions$chrom == sv$chrom2)

        if (length(idx1) == 0 || length(idx2) == 0) next

        arc_list[[i]] <- data.frame(
            x1 = chr_positions$x[idx1],
            y1 = chr_positions$y[idx1],
            x2 = chr_positions$x[idx2],
            y2 = chr_positions$y[idx2]
        )
    }

    if (length(arc_list) > 0) {
        return(do.call(rbind, arc_list))
    } else {
        return(data.frame(x1 = numeric(0), y1 = numeric(0),
                         x2 = numeric(0), y2 = numeric(0)))
    }
}
