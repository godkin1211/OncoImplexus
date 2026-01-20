#' Sort chromosomes in natural order
#'
#' @param chroms Character vector of chromosome names
#' @return Sorted character vector
#' @keywords internal
sort_chromosomes <- function(chroms) {
    # Extract numeric part
    chr_nums <- gsub("^chr", "", chroms)

    # Separate autosomes and sex chromosomes
    auto <- chr_nums[!chr_nums %in% c("X", "Y", "M", "MT")]
    sex <- chr_nums[chr_nums %in% c("X", "Y")]
    mito <- chr_nums[chr_nums %in% c("M", "MT")]

    # Sort autosomes numerically
    if (length(auto) > 0) {
        auto <- auto[order(as.numeric(auto))]
    }

    # Combine in order
    sorted_nums <- c(auto, sex, mito)

    # Add back "chr" prefix if original had it
    if (any(grepl("^chr", chroms))) {
        sorted_chroms <- paste0("chr", sorted_nums)
    } else {
        sorted_chroms <- sorted_nums
    }

    # Only return chromosomes that exist in input
    sorted_chroms[sorted_chroms %in% chroms]
}


#' Plot mechanism landscape across chromosomes
#'
#' Creates a genome-wide view showing distribution of chromoanagenesis mechanisms
#' across all chromosomes, highlighting mixed mechanism regions.
#'
#' @param mixed_mechanisms_result Result from classify_mixed_mechanisms()
#' @param sample_name Sample name for plot title
#' @param show_confidence Include confidence scores in visualization (default: TRUE)
#' @return A ggplot object
#' @export
plot_mechanism_landscape <- function(mixed_mechanisms_result,
                                    sample_name = "",
                                    show_confidence = TRUE) {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required for plotting.")
    }

    if (!inherits(mixed_mechanisms_result, "mixed_mechanisms")) {
        stop("Input must be a mixed_mechanisms result object")
    }

    locations <- mixed_mechanisms_result$mechanism_locations

    if (nrow(locations) == 0) {
        message("No mechanisms to plot.")
        return(NULL)
    }

    # Prepare chromosome order using natural sorting
    chr_order <- sort_chromosomes(unique(locations$chrom))
    locations$chrom <- factor(locations$chrom, levels = chr_order)

    # Define mechanism order to match y-axis labels
    mechanism_order <- c("chromothripsis", "chromoplexy", "chromoanasynthesis")
    locations$mechanism <- factor(locations$mechanism, levels = mechanism_order)

    # Color scheme for mechanisms
    mechanism_colors <- c(
        "chromothripsis" = "#E41A1C",
        "chromoplexy" = "#377EB8",
        "chromoanasynthesis" = "#4DAF4A"
    )

    # Create plot
    p <- ggplot2::ggplot(locations)

    # Add chromosome backgrounds (alternating)
    chr_levels <- levels(locations$chrom)
    for (i in seq_along(chr_levels)) {
        if (i %% 2 == 0) {
            p <- p + ggplot2::annotate(
                "rect",
                xmin = i - 0.4, xmax = i + 0.4,
                ymin = -Inf, ymax = Inf,
                fill = "gray95",
                alpha = 0.5
            )
        }
    }

    # Plot mechanisms as points
    if (show_confidence) {
        p <- p + ggplot2::geom_point(
            ggplot2::aes(
                x = chrom,
                y = as.numeric(mechanism),  # Use pre-factored mechanism
                color = mechanism,
                size = confidence
            ),
            position = ggplot2::position_jitter(width = 0.2, height = 0.1),
            alpha = 0.7
        )
    } else {
        p <- p + ggplot2::geom_point(
            ggplot2::aes(
                x = chrom,
                y = as.numeric(mechanism),  # Use pre-factored mechanism
                color = mechanism
            ),
            position = ggplot2::position_jitter(width = 0.2, height = 0.1),
            size = 3,
            alpha = 0.7
        )
    }

    # Add colors
    p <- p + ggplot2::scale_color_manual(
        values = mechanism_colors,
        name = "Mechanism"
    )

    # Add size scale if showing confidence
    if (show_confidence) {
        p <- p + ggplot2::scale_size_continuous(
            range = c(2, 6),
            name = "Confidence"
        )
    }

    # Customize y-axis
    p <- p + ggplot2::scale_y_continuous(
        breaks = 1:3,
        labels = c("Chromothripsis", "Chromoplexy", "Chromoanasynthesis"),
        limits = c(0.5, 3.5)
    )

    # Labels and theme
    title_text <- sprintf("Chromoanagenesis Landscape%s",
                         ifelse(sample_name != "", paste0(" - ", sample_name), ""))

    subtitle_text <- sprintf(
        "Classification: %s | Complexity: %s (%.2f)",
        mixed_mechanisms_result$sample_classification$classification,
        mixed_mechanisms_result$complexity$complexity_level,
        mixed_mechanisms_result$complexity$complexity_score
    )

    p <- p + ggplot2::labs(
        title = title_text,
        subtitle = subtitle_text,
        x = "Chromosome",
        y = "Mechanism Type"
    )

    p <- p + ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(size = 14, face = "bold"),
            plot.subtitle = ggplot2::element_text(size = 11),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            panel.grid.major.x = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            legend.position = "right"
        )

    return(p)
}


#' Plot chromosome-level mechanism classification
#'
#' Creates a heatmap showing which mechanisms are present on each chromosome
#' and highlights mixed mechanism regions.
#'
#' @param mixed_mechanisms_result Result from classify_mixed_mechanisms()
#' @param sample_name Sample name for plot title
#' @return A ggplot object
#' @export
plot_chromosome_mechanisms <- function(mixed_mechanisms_result,
                                      sample_name = "") {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required for plotting.")
    }

    chr_class <- mixed_mechanisms_result$chromosome_classification

    if (nrow(chr_class) == 0) {
        message("No chromosome classification to plot.")
        return(NULL)
    }

    locations <- mixed_mechanisms_result$mechanism_locations

    # Create presence/absence matrix using natural chromosome sorting
    chr_order <- sort_chromosomes(chr_class$chrom)

    mechanisms <- c("chromothripsis", "chromoplexy", "chromoanasynthesis")

    # Build matrix
    presence_data <- list()
    for (chr in chr_order) {
        for (mech in mechanisms) {
            is_present <- any(locations$chrom == chr & locations$mechanism == mech)

            presence_data[[length(presence_data) + 1]] <- list(
                chrom = chr,
                mechanism = mech,
                present = is_present,
                is_mixed = chr_class$is_mixed[chr_class$chrom == chr][1]
            )
        }
    }

    presence_df <- do.call(rbind, lapply(presence_data, function(x) {
        data.frame(
            chrom = as.character(x$chrom),
            mechanism = as.character(x$mechanism),
            present = as.logical(x$present),
            is_mixed = as.logical(x$is_mixed),
            stringsAsFactors = FALSE
        )
    }))

    # Set factor levels
    presence_df$chrom <- factor(presence_df$chrom, levels = chr_order)
    presence_df$mechanism <- factor(presence_df$mechanism,
                                   levels = rev(mechanisms),
                                   labels = rev(c("Chromothripsis", "Chromoplexy", "Chromoanasynthesis")))

    # Create plot
    p <- ggplot2::ggplot(presence_df, ggplot2::aes(x = chrom, y = mechanism))

    # Add tiles
    p <- p + ggplot2::geom_tile(
        ggplot2::aes(fill = present),
        color = "white",
        size = 1
    )

    # Highlight mixed chromosomes with border
    mixed_chrs <- unique(presence_df$chrom[presence_df$is_mixed])
    if (length(mixed_chrs) > 0) {
        for (chr in mixed_chrs) {
            for (i in 1:length(mechanisms)) {
                p <- p + ggplot2::annotate(
                    "rect",
                    xmin = as.numeric(factor(chr, levels = levels(presence_df$chrom))) - 0.45,
                    xmax = as.numeric(factor(chr, levels = levels(presence_df$chrom))) + 0.45,
                    ymin = i - 0.45,
                    ymax = i + 0.45,
                    fill = NA,
                    color = "red",
                    size = 1.5,
                    alpha = 0.8
                )
            }
        }
    }

    # Color scale
    p <- p + ggplot2::scale_fill_manual(
        values = c("TRUE" = "#2166AC", "FALSE" = "gray90"),
        labels = c("TRUE" = "Present", "FALSE" = "Absent"),
        name = "Status"
    )

    # Labels
    title_text <- sprintf("Chromosome-Level Mechanism Distribution%s",
                         ifelse(sample_name != "", paste0(" - ", sample_name), ""))

    subtitle_text <- sprintf(
        "Mixed chromosomes: %d/%d (highlighted with red border)",
        mixed_mechanisms_result$sample_classification$n_mixed_chromosomes,
        nrow(chr_class)
    )

    p <- p + ggplot2::labs(
        title = title_text,
        subtitle = subtitle_text,
        x = "Chromosome",
        y = ""
    )

    # Theme
    p <- p + ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(size = 14, face = "bold"),
            plot.subtitle = ggplot2::element_text(size = 10),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10),
            axis.text.y = ggplot2::element_text(size = 11),
            panel.grid = ggplot2::element_blank(),
            legend.position = "right"
        )

    return(p)
}


#' Plot complexity components breakdown
#'
#' Creates a bar plot showing the contribution of different components
#' to the overall complexity score.
#'
#' @param mixed_mechanisms_result Result from classify_mixed_mechanisms()
#' @param sample_name Sample name for plot title
#' @return A ggplot object
#' @export
plot_complexity_breakdown <- function(mixed_mechanisms_result,
                                     sample_name = "") {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required for plotting.")
    }

    components <- mixed_mechanisms_result$complexity$components

    # Create plot
    p <- ggplot2::ggplot(components, ggplot2::aes(x = component, y = weighted_score))

    # Add bars
    p <- p + ggplot2::geom_col(
        fill = "#377EB8",
        alpha = 0.8,
        width = 0.7
    )

    # Add raw score as points
    p <- p + ggplot2::geom_point(
        ggplot2::aes(y = raw_score),
        color = "#E41A1C",
        size = 4,
        shape = 18
    )

    # Add weight labels
    p <- p + ggplot2::geom_text(
        ggplot2::aes(
            label = sprintf("Weight: %.0f%%", weight * 100),
            y = 0.05
        ),
        hjust = 0,
        size = 3,
        color = "white",
        angle = 90
    )

    # Add horizontal line for total score
    total_score <- mixed_mechanisms_result$complexity$complexity_score
    p <- p + ggplot2::geom_hline(
        yintercept = total_score,
        linetype = "dashed",
        color = "darkgreen",
        size = 1
    )

    p <- p + ggplot2::annotate(
        "text",
        x = 4.3,
        y = total_score,
        label = sprintf("Total: %.3f", total_score),
        color = "darkgreen",
        hjust = 0,
        size = 4,
        fontface = "bold"
    )

    # Labels
    title_text <- sprintf("Complexity Score Breakdown%s",
                         ifelse(sample_name != "", paste0(" - ", sample_name), ""))

    subtitle_text <- sprintf(
        "Overall Complexity: %s",
        mixed_mechanisms_result$complexity$complexity_level
    )

    p <- p + ggplot2::labs(
        title = title_text,
        subtitle = subtitle_text,
        x = "",
        y = "Score",
        caption = "Blue bars = weighted score, Red diamonds = raw score"
    )

    # Theme
    p <- p + ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(size = 14, face = "bold"),
            plot.subtitle = ggplot2::element_text(size = 11),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10),
            panel.grid.major.x = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank()
        )

    # Set y-axis limits
    p <- p + ggplot2::ylim(0, max(1, total_score * 1.1))

    return(p)
}


#' Plot mechanism dominance pie chart
#'
#' Creates a pie chart showing the relative proportions of different
#' chromoanagenesis mechanisms.
#'
#' @param mixed_mechanisms_result Result from classify_mixed_mechanisms()
#' @param sample_name Sample name for plot title
#' @return A ggplot object
#' @export
plot_mechanism_dominance <- function(mixed_mechanisms_result,
                                    sample_name = "") {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required for plotting.")
    }

    proportions <- mixed_mechanisms_result$dominance$mechanism_proportions

    # Filter out mechanisms with 0 events
    proportions <- proportions[proportions$n_events > 0, ]

    if (nrow(proportions) == 0) {
        message("No events to plot.")
        return(NULL)
    }

    # Calculate percentages and labels
    proportions$percentage <- proportions$proportion * 100
    proportions$label <- sprintf("%s\n%d events\n(%.1f%%)",
                                tools::toTitleCase(proportions$mechanism),
                                proportions$n_events,
                                proportions$percentage)

    # Calculate position for labels
    proportions$pos <- cumsum(proportions$proportion) - 0.5 * proportions$proportion

    # Color scheme
    mechanism_colors <- c(
        "chromothripsis" = "#E41A1C",
        "chromoplexy" = "#377EB8",
        "chromoanasynthesis" = "#4DAF4A"
    )

    # Create pie chart (stacked bar converted to polar)
    p <- ggplot2::ggplot(proportions, ggplot2::aes(x = "", y = proportion, fill = mechanism))

    p <- p + ggplot2::geom_col(width = 1, color = "white", size = 2)

    # Add labels
    p <- p + ggplot2::geom_text(
        ggplot2::aes(y = pos, label = label),
        size = 4,
        color = "white",
        fontface = "bold"
    )

    # Convert to pie chart
    p <- p + ggplot2::coord_polar(theta = "y")

    # Colors
    p <- p + ggplot2::scale_fill_manual(
        values = mechanism_colors,
        name = "Mechanism"
    )

    # Labels
    title_text <- sprintf("Mechanism Dominance%s",
                         ifelse(sample_name != "", paste0(" - ", sample_name), ""))

    subtitle_text <- sprintf(
        "Dominant: %s | Total events: %d",
        tools::toTitleCase(mixed_mechanisms_result$dominance$dominant_mechanism),
        mixed_mechanisms_result$dominance$total_events
    )

    p <- p + ggplot2::labs(
        title = title_text,
        subtitle = subtitle_text,
        x = NULL,
        y = NULL
    )

    # Theme
    p <- p + ggplot2::theme_void() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
            plot.subtitle = ggplot2::element_text(size = 11, hjust = 0.5),
            legend.position = "none"
        )

    return(p)
}


#' Create comprehensive mechanism report plot
#'
#' Combines multiple visualizations into a single comprehensive report.
#'
#' @param mixed_mechanisms_result Result from classify_mixed_mechanisms()
#' @param sample_name Sample name for plot title
#' @return A combined plot object (requires gridExtra)
#' @export
plot_mechanism_report <- function(mixed_mechanisms_result,
                                 sample_name = "") {

    if (!requireNamespace("gridExtra", quietly = TRUE)) {
        stop("Package 'gridExtra' is required for creating reports.")
    }

    if (!requireNamespace("grid", quietly = TRUE)) {
        stop("Package 'grid' is required for creating reports.")
    }

    # Create individual plots
    p1 <- plot_mechanism_landscape(mixed_mechanisms_result, sample_name, show_confidence = TRUE)
    p2 <- plot_chromosome_mechanisms(mixed_mechanisms_result, sample_name)
    p3 <- plot_complexity_breakdown(mixed_mechanisms_result, sample_name)
    p4 <- plot_mechanism_dominance(mixed_mechanisms_result, sample_name)

    # Combine plots
    if (!is.null(p1) && !is.null(p2) && !is.null(p3) && !is.null(p4)) {
        # Create layout
        combined <- gridExtra::grid.arrange(
            p1, p2,
            p3, p4,
            ncol = 2,
            nrow = 2,
            top = grid::textGrob(
                sprintf("Integrated Chromoanagenesis Report%s",
                       ifelse(sample_name != "", paste0(" - ", sample_name), "")),
                gp = grid::gpar(fontsize = 16, fontface = "bold")
            )
        )

        return(combined)
    } else {
        message("Not all plots could be generated. Check if there are events to plot.")
        return(NULL)
    }
}
