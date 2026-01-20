#' gGnome-style Regional Visualization
#'
#' Functions for creating detailed regional views of chromoanagenesis events
#' with gGnome-style aesthetics, showing structural complexity through
#' SV arcs and copy number profiles.
#'
#' @name ggnome_style_viz
#' @keywords internal


#' Plot regional chromoanagenesis detail with SV arcs
#'
#' Creates a detailed view of a chromoanagenesis region showing:
#' - Chromosome ideogram (optional)
#' - SV arcs (colored by type)
#' - Copy number profile
#'
#' @param chromoanagenesis_result Result from detect_chromoanagenesis()
#' @param SV.sample Original SV data
#' @param CNV.sample Original CNV data
#' @param chrom Chromosome to plot
#' @param region_idx Region index (if multiple regions on same chromosome)
#' @param mechanism Type of mechanism to plot ("chromothripsis", "chromoplexy", "chromoanasynthesis")
#' @param show_ideogram Include chromosome ideogram (default: TRUE)
#' @param sample_name Sample name for title
#' @return List of ggplot objects
#' @export
plot_chromoanagenesis_region <- function(chromoanagenesis_result,
                                         SV.sample,
                                         CNV.sample,
                                         chrom,
                                         region_idx = 1,
                                         mechanism = "chromothripsis",
                                         show_ideogram = TRUE,
                                         sample_name = "") {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    # Convert S4 objects if needed
    if (is(SV.sample, "SVs")) {
        sv_data <- as(SV.sample, "data.frame")
    } else {
        sv_data <- SV.sample
    }

    if (is(CNV.sample, "CNVsegs")) {
        cnv_data <- as(CNV.sample, "data.frame")
    } else {
        cnv_data <- CNV.sample
    }

    # Standardize chromosome name
    chrom <- gsub("^chr", "", chrom)

    # Get the region coordinates based on mechanism
    region_coords <- get_region_coordinates(chromoanagenesis_result, chrom, mechanism, region_idx)

    if (is.null(region_coords)) {
        message(sprintf("No %s region found on chromosome %s", mechanism, chrom))
        return(NULL)
    }

    # Filter SVs and CNVs to region
    sv_region <- filter_sv_to_region(sv_data, chrom, region_coords$start, region_coords$end)
    cnv_region <- filter_cnv_to_region(cnv_data, chrom, region_coords$start, region_coords$end)

    # Create plots
    plots <- list()

    # 1. Ideogram (optional)
    if (show_ideogram) {
        plots$ideogram <- plot_region_ideogram(chrom, region_coords$start, region_coords$end)
    }

    # 2. SV arcs plot
    plots$sv_arcs <- plot_sv_arcs(sv_region, region_coords$start, region_coords$end,
        mechanism = mechanism
    )

    # 3. Copy number profile
    plots$cn_profile <- plot_region_cn(cnv_region, region_coords$start, region_coords$end)

    # Add title information
    title_text <- sprintf(
        "%s Region - Chr%s:%s-%s",
        tools::toTitleCase(mechanism),
        chrom,
        format_genomic_position(region_coords$start),
        format_genomic_position(region_coords$end)
    )

    if (sample_name != "") {
        title_text <- paste0(title_text, " - ", sample_name)
    }

    attr(plots, "title") <- title_text

    return(plots)
}


#' Get region coordinates for a specific mechanism
#'
#' @keywords internal
get_region_coordinates <- function(chromoanagenesis_result, chrom, mechanism, region_idx = 1) {
    if (mechanism == "chromothripsis") {
        if (!is.null(chromoanagenesis_result$chromothripsis)) {
            ct_class <- chromoanagenesis_result$chromothripsis$classification
            if (!is.null(ct_class) && nrow(ct_class) > 0) {
                regions <- ct_class[ct_class$chrom == chrom &
                    ct_class$classification %in% c("High confidence", "Low confidence"), ]
                if (nrow(regions) >= region_idx) {
                    # Get start/end from detection_output chromSummary
                    ct_summary <- chromoanagenesis_result$chromothripsis$detection_output@chromSummary
                    coords <- ct_summary[ct_summary$chrom == chrom, c("start", "end")]

                    if (nrow(coords) >= region_idx) {
                        return(list(
                            start = coords$start[region_idx],
                            end = coords$end[region_idx]
                        ))
                    }
                }
            }
        }
    } else if (mechanism == "chromoplexy") {
        if (!is.null(chromoanagenesis_result$chromoplexy)) {
            cp_summary <- chromoanagenesis_result$chromoplexy$summary
            if (!is.null(cp_summary) && nrow(cp_summary) > 0) {
                regions <- cp_summary[cp_summary$classification %in% c("Likely chromoplexy", "Possible chromoplexy"), ]
                # For chromoplexy, extract region from involved chromosomes
                if (nrow(regions) >= region_idx) {
                    # Get chains involving this chromosome
                    chains <- chromoanagenesis_result$chromoplexy$chains
                    if (!is.null(chains) && length(chains) >= region_idx) {
                        chain <- chains[[region_idx]]
                        chr_svs <- chain[chain$chrom1 == chrom | chain$chrom2 == chrom, ]
                        if (nrow(chr_svs) > 0) {
                            all_pos <- c(
                                chr_svs$pos1[chr_svs$chrom1 == chrom],
                                chr_svs$pos2[chr_svs$chrom2 == chrom]
                            )
                            return(list(
                                start = min(all_pos),
                                end = max(all_pos)
                            ))
                        }
                    }
                }
            }
        }
    } else if (mechanism == "chromoanasynthesis") {
        if (!is.null(chromoanagenesis_result$chromoanasynthesis)) {
            cs_summary <- chromoanagenesis_result$chromoanasynthesis$summary
            if (!is.null(cs_summary) && nrow(cs_summary) > 0) {
                regions <- cs_summary[cs_summary$chrom == chrom &
                    cs_summary$classification %in% c("Likely chromoanasynthesis", "Possible chromoanasynthesis"), ]
                if (nrow(regions) >= region_idx) {
                    return(list(
                        start = regions$start[region_idx],
                        end = regions$end[region_idx]
                    ))
                }
            }
        }
    }

    return(NULL)
}


#' Filter SVs to a specific region
#'
#' @keywords internal
filter_sv_to_region <- function(sv_data, chrom, start_pos, end_pos) {
    # Include SVs that overlap the region
    sv_data[sv_data$chrom1 == chrom &
        sv_data$pos1 >= start_pos &
        sv_data$pos2 <= end_pos, ]
}


#' Filter CNVs to a specific region
#'
#' @keywords internal
filter_cnv_to_region <- function(cnv_data, chrom, start_pos, end_pos) {
    # Include CNV segments that overlap the region
    cnv_overlap <- cnv_data[cnv_data$chrom == chrom &
        ((cnv_data$start >= start_pos & cnv_data$start <= end_pos) |
            (cnv_data$end >= start_pos & cnv_data$end <= end_pos) |
            (cnv_data$start <= start_pos & cnv_data$end >= end_pos)), ]

    # Trim to region boundaries
    if (nrow(cnv_overlap) > 0) {
        cnv_overlap$start[cnv_overlap$start < start_pos] <- start_pos
        cnv_overlap$end[cnv_overlap$end > end_pos] <- end_pos
    }

    return(cnv_overlap)
}


#' Plot region ideogram
#'
#' @keywords internal
plot_region_ideogram <- function(chrom, start_pos, end_pos) {
    # Simplified ideogram - just a chromosome bar with highlighted region
    # Create data frame for geom_rect
    ideogram_data <- data.frame(
        xmin = start_pos,
        xmax = end_pos,
        ymin = 0,
        ymax = 1
    )

    p <- ggplot2::ggplot(ideogram_data) +
        ggplot2::geom_rect(
            ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray80", color = "black", size = 0.5
        ) +
        ggplot2::scale_x_continuous(
            expand = c(0.01, 0.01),
            labels = function(x) {
                format_genomic_position(x)
            }
        ) +
        ggplot2::labs(x = NULL, y = NULL) +
        ggplot2::theme_void() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(size = 8),
            plot.margin = ggplot2::margin(5, 5, 0, 5)
        )

    return(p)
}


#' Plot SV arcs
#'
#' Creates arcs representing structural variants with gGnome-style aesthetics.
#'
#' @keywords internal
plot_sv_arcs <- function(sv_data, start_pos, end_pos, mechanism = "chromothripsis") {
    if (nrow(sv_data) == 0) {
        p <- ggplot2::ggplot() +
            ggplot2::annotate("text",
                x = (start_pos + end_pos) / 2, y = 8,
                label = "No SVs in this region",
                size = 5, color = "gray50"
            ) +
            ggplot2::theme_void()
        return(p)
    }

    # Calculate arc curvatures based on SV length
    sv_data$length <- abs(sv_data$pos2 - sv_data$pos1)
    max_length <- max(sv_data$length)
    sv_data$curvature <- 1 - (sv_data$length / max_length)

    # Adjust curvature for better visualization
    sv_data$curvature[sv_data$length / max_length > 0.8] <- 0.08
    sv_data$curvature[sv_data$length / max_length > 0.2 & sv_data$length / max_length <= 0.8] <- 0.15
    sv_data$curvature[sv_data$length / max_length <= 0.2] <- 1

    # gGnome-style colors (more muted, professional)
    sv_colors <- c(
        "DEL" = "#f4a582", # Coral
        "DUP" = "#92c5de", # Sky blue
        "t2tINV" = "#66c2a5", # Teal
        "h2hINV" = "#fc8d62" # Orange
    )

    # Y positions for different SV types
    y_del_dup <- 4
    y_inv <- 12

    # Create base plot
    p <- ggplot2::ggplot() +
        ggplot2::geom_hline(yintercept = y_del_dup, color = "gray60", size = 0.5) +
        ggplot2::geom_hline(yintercept = y_inv, color = "gray60", size = 0.5)

    # Add arcs for each SV type
    for (svtype in c("DUP", "DEL", "t2tINV", "h2hINV")) {
        sv_subset <- sv_data[sv_data$SVtype == svtype, ]

        if (nrow(sv_subset) > 0) {
            y_pos <- if (svtype %in% c("DUP", "DEL")) y_del_dup else y_inv
            curv_sign <- if (svtype %in% c("DUP", "t2tINV")) 1 else -1

            for (i in seq_len(nrow(sv_subset))) {
                p <- p + ggplot2::geom_curve(
                    data = sv_subset[i, ],
                    ggplot2::aes(x = pos1, y = y_pos, xend = pos2, yend = y_pos),
                    curvature = curv_sign * sv_subset$curvature[i],
                    color = sv_colors[svtype],
                    size = 0.4,
                    alpha = 0.7,
                    ncp = 8
                )

                # Add breakpoint markers
                p <- p +
                    ggplot2::geom_point(
                        data = sv_subset[i, ],
                        ggplot2::aes(x = pos1, y = y_pos),
                        color = sv_colors[svtype],
                        size = 0.8,
                        alpha = 0.8
                    ) +
                    ggplot2::geom_point(
                        data = sv_subset[i, ],
                        ggplot2::aes(x = pos2, y = y_pos),
                        color = sv_colors[svtype],
                        size = 0.8,
                        alpha = 0.8
                    )
            }
        }
    }

    # Add invisible points for legend
    # Create data for each SV type present in the data
    present_types <- unique(sv_data$SVtype)
    present_types <- present_types[present_types %in% names(sv_colors)]

    if (length(present_types) > 0) {
        legend_data <- data.frame(
            svtype = factor(present_types, levels = names(sv_colors)),
            x = start_pos,
            y = 1,
            stringsAsFactors = FALSE
        )

        p <- p +
            ggplot2::geom_point(
                data = legend_data,
                ggplot2::aes(x = x, y = y, color = svtype),
                size = 0, # Invisible points for legend only
                alpha = 0
            ) +
            ggplot2::scale_color_manual(
                values = sv_colors,
                name = "SV Type",
                breaks = present_types
            )
    }

    # Theme
    p <- p +
        ggplot2::scale_x_continuous(
            expand = c(0.01, 0.01),
            limits = c(start_pos, end_pos)
        ) +
        ggplot2::ylim(0, 16) +
        ggplot2::labs(x = NULL, y = NULL) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_blank(),
            panel.grid = ggplot2::element_blank(),
            legend.position = "right",
            legend.text = ggplot2::element_text(size = 9),
            plot.margin = ggplot2::margin(5, 5, 0, 5)
        )

    return(p)
}


#' Plot region copy number profile
#'
#' @keywords internal
plot_region_cn <- function(cnv_data, start_pos, end_pos) {
    if (nrow(cnv_data) == 0) {
        p <- ggplot2::ggplot() +
            ggplot2::annotate("text",
                x = (start_pos + end_pos) / 2, y = 2,
                label = "No CNV data in this region",
                size = 5, color = "gray50"
            ) +
            ggplot2::theme_void()
        return(p)
    }

    # Determine appropriate y-axis scale
    max_cn <- max(cnv_data$total_cn, na.rm = TRUE)
    min_cn <- min(cnv_data$total_cn, na.rm = TRUE)

    # gGnome-style color for CN segments
    p <- ggplot2::ggplot() +
        ggplot2::geom_segment(
            data = cnv_data,
            ggplot2::aes(x = start, xend = end, y = total_cn, yend = total_cn),
            color = "#4a4a4a",
            size = 2
        )

    # Add reference line at CN=2
    p <- p +
        ggplot2::geom_hline(
            yintercept = 2,
            linetype = "dashed",
            color = "gray40",
            size = 0.5
        )

    # Determine y-axis breaks
    if (max_cn <= 10) {
        y_breaks <- c(0, sort(unique(cnv_data$total_cn)))
        p <- p + ggplot2::scale_y_continuous(
            breaks = y_breaks,
            limits = c(max(0, min_cn - 0.5), max_cn + 0.5)
        )
    } else if (max_cn <= 100) {
        # Log2 scale for moderate CN
        cnv_data_log <- cnv_data
        cnv_data_log$total_cn[cnv_data_log$total_cn == 0] <- 0.99
        y_breaks <- c(0, 1, 2, 4, 10, 20, max_cn)
        p <- ggplot2::ggplot() +
            ggplot2::geom_segment(
                data = cnv_data_log,
                ggplot2::aes(x = start, xend = end, y = log2(total_cn), yend = log2(total_cn)),
                color = "#4a4a4a",
                size = 2
            ) +
            ggplot2::scale_y_continuous(
                breaks = log2(y_breaks),
                labels = as.character(y_breaks)
            )
    } else {
        # Log10 scale for high CN
        cnv_data_log <- cnv_data
        cnv_data_log$total_cn[cnv_data_log$total_cn == 0] <- 0.99
        y_breaks <- c(0, 1, 2, 4, 20, 50, 100, max_cn)
        p <- ggplot2::ggplot() +
            ggplot2::geom_segment(
                data = cnv_data_log,
                ggplot2::aes(x = start, xend = end, y = log10(total_cn), yend = log10(total_cn)),
                color = "#4a4a4a",
                size = 2
            ) +
            ggplot2::scale_y_continuous(
                breaks = log10(y_breaks),
                labels = as.character(y_breaks)
            )
    }

    # Theme
    p <- p +
        ggplot2::scale_x_continuous(
            expand = c(0.01, 0.01),
            limits = c(start_pos, end_pos),
            labels = function(x) {
                format_genomic_position(x)
            }
        ) +
        ggplot2::labs(x = "Position", y = "Copy Number") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text = ggplot2::element_text(size = 9),
            axis.title = ggplot2::element_text(size = 10),
            panel.grid.minor = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(0, 5, 5, 5)
        )

    return(p)
}


# Note: format_genomic_position is defined in visualization_utils.R


#' Arrange regional plots
#'
#' Combines ideogram, SV arcs, and CN profile into a single figure.
#'
#' @param plot_list List of plots from plot_chromoanagenesis_region()
#' @return Combined plot
#' @export
arrange_regional_plots <- function(plot_list) {
    if (!requireNamespace("gridExtra", quietly = TRUE)) {
        stop("Package 'gridExtra' is required.")
    }

    if (!requireNamespace("grid", quietly = TRUE)) {
        stop("Package 'grid' is required.")
    }

    title_text <- attr(plot_list, "title")
    if (is.null(title_text)) {
        title_text <- "Chromoanagenesis Region Detail"
    }

    # Determine layout based on available plots
    if (!is.null(plot_list$ideogram)) {
        gridExtra::grid.arrange(
            grobs = list(plot_list$ideogram, plot_list$sv_arcs, plot_list$cn_profile),
            heights = c(0.1, 0.5, 0.4),
            ncol = 1,
            top = grid::textGrob(
                title_text,
                gp = grid::gpar(fontsize = 14, fontface = "bold")
            )
        )
    } else {
        gridExtra::grid.arrange(
            grobs = list(plot_list$sv_arcs, plot_list$cn_profile),
            heights = c(0.6, 0.4),
            ncol = 1,
            top = grid::textGrob(
                title_text,
                gp = grid::gpar(fontsize = 14, fontface = "bold")
            )
        )
    }
}
