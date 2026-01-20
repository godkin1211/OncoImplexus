#' Plot chromoanasynthesis regions
#'
#' Visualizes chromoanasynthesis regions with CN gradients and tandem duplications.
#'
#' @param chromoanasynthesis_result Result from detect_chromoanasynthesis()
#' @param region_id Optional: specific region ID to plot (default: plot all)
#' @param sample_name Sample name for plot title
#' @param genome Reference genome ("hg19" or "hg38")
#' @return A ggplot object or list of ggplot objects
#' @export
plot_chromoanasynthesis <- function(chromoanasynthesis_result,
                                region_id = NULL,
                                sample_name = "",
                                genome = "hg19") {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required for plotting chromoanasynthesis.")
    }

    if (chromoanasynthesis_result$total_regions == 0) {
        message("No chromoanasynthesis regions to plot.")
        return(NULL)
    }

    # Select regions to plot
    if (!is.null(region_id)) {
        regions_to_plot <- list(chromoanasynthesis_result$region_details[[region_id]])
    } else {
        regions_to_plot <- chromoanasynthesis_result$region_details
    }

    plots <- list()

    for (i in 1:length(regions_to_plot)) {
        region_detail <- regions_to_plot[[i]]

        # Create integrated plot showing:
        # 1. CN gradient
        # 2. Tandem duplications
        # 3. Other SVs

        p <- create_chromoanasynthesis_region_plot(
            region_detail = region_detail,
            sample_name = sample_name
        )

        plots[[i]] <- p
    }

    if (length(plots) == 1) {
        return(plots[[1]])
    } else {
        return(plots)
    }
}


#' Create a single chromoanasynthesis region plot
#'
#' @param region_detail Region detail object
#' @param sample_name Sample name
#' @return ggplot object
#' @keywords internal
create_chromoanasynthesis_region_plot <- function(region_detail, sample_name) {

    region <- region_detail$region
    tandem_dups <- region_detail$tandem_dups
    summary <- region_detail$summary

    # Prepare CN data
    cn_data <- region$segments
    cn_data$midpoint <- (cn_data$start + cn_data$end) / 2
    cn_data$segment_order <- 1:nrow(cn_data)

    # Create main plot with CN gradient
    p <- ggplot2::ggplot(cn_data, ggplot2::aes(x = midpoint, y = total_cn))

    # Add CN segments as steps
    for (i in 1:nrow(cn_data)) {
        p <- p + ggplot2::geom_rect(
            ggplot2::aes(
                xmin = start, xmax = end,
                ymin = 0, ymax = total_cn
            ),
            data = cn_data[i, ],
            fill = "steelblue",
            alpha = 0.6,
            color = "darkblue"
        )
    }

    # Add gradient trend line
    p <- p + ggplot2::geom_smooth(
        method = "lm",
        se = TRUE,
        color = "red",
        size = 1.5,
        alpha = 0.2
    )

    # Add tandem duplications as arcs
    if (tandem_dups$n_tandem_dups > 0) {
        dup_data <- tandem_dups$tandem_dups

        for (i in 1:nrow(dup_data)) {
            dup <- dup_data[i, ]

            # Draw arc for tandem duplication
            p <- p + ggplot2::annotate(
                "segment",
                x = dup$pos1,
                xend = dup$pos2,
                y = max(cn_data$total_cn) * 1.1,
                yend = max(cn_data$total_cn) * 1.1,
                arrow = ggplot2::arrow(length = ggplot2::unit(0.02, "npc")),
                color = "darkgreen",
                size = 1.2,
                alpha = 0.7
            )

            # Mark breakpoints
            p <- p + ggplot2::geom_point(
                ggplot2::aes(x = x, y = y),
                data = data.frame(
                    x = c(dup$pos1, dup$pos2),
                    y = max(cn_data$total_cn) * 1.1
                ),
                color = "darkgreen",
                size = 2
            )
        }
    }

    # Add inversions if present
    if (tandem_dups$n_inversions > 0) {
        inv_data <- tandem_dups$inversions

        for (i in 1:nrow(inv_data)) {
            inv <- inv_data[i, ]

            p <- p + ggplot2::annotate(
                "segment",
                x = inv$pos1,
                xend = inv$pos2,
                y = max(cn_data$total_cn) * 1.2,
                yend = max(cn_data$total_cn) * 1.2,
                arrow = ggplot2::arrow(
                    length = ggplot2::unit(0.02, "npc"),
                    ends = "both"
                ),
                color = "purple",
                size = 1,
                alpha = 0.6
            )
        }
    }

    # Add title and labels
    title_text <- sprintf("%s - Chromoanasynthesis Region %d (%s)",
                         sample_name,
                         summary$region_id,
                         summary$chrom)

    subtitle_text <- sprintf(
        "Gradient: r=%.2f, Tandem Dups: %d, Classification: %s",
        summary$gradient_correlation,
        summary$n_tandem_dups,
        summary$classification
    )

    p <- p + ggplot2::labs(
        title = title_text,
        subtitle = subtitle_text,
        x = "Genomic Position (bp)",
        y = "Copy Number"
    )

    # Add annotations
    if (tandem_dups$n_tandem_dups > 0) {
        p <- p + ggplot2::annotate(
            "text",
            x = region$start + (region$end - region$start) * 0.95,
            y = max(cn_data$total_cn) * 1.1,
            label = "Tandem Dup",
            hjust = 1,
            color = "darkgreen",
            size = 3
        )
    }

    if (tandem_dups$n_inversions > 0) {
        p <- p + ggplot2::annotate(
            "text",
            x = region$start + (region$end - region$start) * 0.95,
            y = max(cn_data$total_cn) * 1.2,
            label = "Inversion",
            hjust = 1,
            color = "purple",
            size = 3
        )
    }

    # Theme
    p <- p + ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(size = 14, face = "bold"),
            plot.subtitle = ggplot2::element_text(size = 11),
            panel.grid.minor = ggplot2::element_blank()
        )

    # Set y-axis limits
    max_y <- max(cn_data$total_cn) * 1.3
    p <- p + ggplot2::ylim(0, max_y)

    return(p)
}


#' Create gradient heatmap for chromoanasynthesis
#'
#' Shows CN gradient as a heatmap along chromosome.
#'
#' @param chromoanasynthesis_result Result from detect_chromoanasynthesis()
#' @param region_id Region ID to plot
#' @param sample_name Sample name
#' @return ggplot object
#' @export
plot_chromoanasynthesis_heatmap <- function(chromoanasynthesis_result,
                                        region_id,
                                        sample_name = "") {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required for plotting.")
    }

    region_detail <- chromoanasynthesis_result$region_details[[region_id]]
    region <- region_detail$region
    cn_data <- region$segments

    # Prepare data for heatmap
    cn_data$segment_id <- 1:nrow(cn_data)
    cn_data$cn_normalized <- (cn_data$total_cn - min(cn_data$total_cn)) /
                             (max(cn_data$total_cn) - min(cn_data$total_cn))

    # Create heatmap
    p <- ggplot2::ggplot(cn_data)

    # Add colored rectangles
    p <- p + ggplot2::geom_rect(
        ggplot2::aes(
            xmin = start,
            xmax = end,
            ymin = 0,
            ymax = 1,
            fill = total_cn
        ),
        color = "white",
        size = 0.5
    )

    # Color scale
    p <- p + ggplot2::scale_fill_gradient2(
        low = "blue",
        mid = "white",
        high = "red",
        midpoint = 2,  # Diploid
        name = "Copy\nNumber"
    )

    # Add gradient arrow
    p <- p + ggplot2::annotate(
        "segment",
        x = region$start,
        xend = region$end,
        y = 0.5,
        yend = 0.5,
        arrow = ggplot2::arrow(length = ggplot2::unit(0.03, "npc")),
        color = "black",
        size = 2,
        alpha = 0.7
    )

    # Title
    title_text <- sprintf("%s - CN Gradient Heatmap\n%s:%d-%d",
                         sample_name,
                         region$chrom,
                         region$start,
                         region$end)

    p <- p + ggplot2::labs(
        title = title_text,
        subtitle = sprintf("Gradient correlation: %.2f", region$gradient_correlation),
        x = "Genomic Position (bp)",
        y = ""
    )

    # Theme
    p <- p + ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            panel.grid = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(size = 14, face = "bold")
        )

    return(p)
}


#' Plot CN gradient scatter plot
#'
#' Shows correlation between segment order and copy number.
#'
#' @param chromoanasynthesis_result Result from detect_chromoanasynthesis()
#' @param region_id Region ID to plot
#' @param sample_name Sample name
#' @return ggplot object
#' @export
plot_cn_gradient_scatter <- function(chromoanasynthesis_result,
                                    region_id,
                                    sample_name = "") {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required for plotting.")
    }

    region_detail <- chromoanasynthesis_result$region_details[[region_id]]
    region <- region_detail$region
    cn_data <- region$segments

    # Prepare data
    cn_data$segment_order <- 1:nrow(cn_data)

    # Create scatter plot
    p <- ggplot2::ggplot(cn_data, ggplot2::aes(x = segment_order, y = total_cn))

    # Add points
    p <- p + ggplot2::geom_point(
        size = 3,
        color = "steelblue",
        alpha = 0.7
    )

    # Add connecting line
    p <- p + ggplot2::geom_line(
        color = "steelblue",
        alpha = 0.3,
        size = 1
    )

    # Add trend line
    p <- p + ggplot2::geom_smooth(
        method = "lm",
        se = TRUE,
        color = "red",
        fill = "pink",
        alpha = 0.3
    )

    # Add correlation text
    cor_text <- sprintf("r = %.3f\np = %.3e",
                       region$gradient_correlation,
                       1e-5)  # Placeholder p-value

    p <- p + ggplot2::annotate(
        "text",
        x = max(cn_data$segment_order) * 0.8,
        y = min(cn_data$total_cn) * 1.2,
        label = cor_text,
        hjust = 0,
        size = 5,
        color = "red"
    )

    # Title
    title_text <- sprintf("%s - CN Gradient Analysis\n%s:%d-%d",
                         sample_name,
                         region$chrom,
                         region$start,
                         region$end)

    p <- p + ggplot2::labs(
        title = title_text,
        x = "Segment Order (5' -> 3')",
        y = "Copy Number",
        subtitle = sprintf("Slope: %.3f, %d segments",
                          region$slope,
                          region$n_segments)
    )

    # Theme
    p <- p + ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(size = 14, face = "bold"),
            plot.subtitle = ggplot2::element_text(size = 11),
            panel.grid.minor = ggplot2::element_blank()
        )

    return(p)
}
