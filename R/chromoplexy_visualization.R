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
