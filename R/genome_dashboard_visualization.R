#' Genome-Wide Chromoanagenesis Dashboard Visualization
#'
#' Functions for creating comprehensive genome-wide overview of
#' chromoanagenesis events, inspired by gGnome-style integrated views.
#'
#' @name genome_dashboard
#' @keywords internal


#' Plot genome-wide chromoanagenesis dashboard
#'
#' Creates an integrated dashboard with multiple panels showing
#' genome-wide patterns of chromoanagenesis events.
#'
#' @param chromoanagenesis_result Result from detect_chromoanagenesis()
#' @param SV.sample Original SV data (SVs object or data frame)
#' @param CNV.sample Original CNV data (CNVsegs object or data frame)
#' @param sample_name Sample identifier for plot title
#' @param genome Reference genome ("hg19" or "hg38")
#' @param include_panels Vector of panel names to include:
#'   "ideogram", "sv_density", "cn_profile", "mechanism_dist"
#' @return A grid arrangement of plots (uses grid/gridExtra)
#' @export
#'
#' @examples
#' \dontrun{
#' results <- detect_chromoanagenesis(sv_data, cnv_data)
#' dashboard <- plot_genome_dashboard(results, sv_data, cnv_data,
#'                                    sample_name = "Patient_001",
#'                                    genome = "hg38")
#' print(dashboard)
#' }
plot_genome_dashboard <- function(chromoanagenesis_result,
                                  SV.sample = NULL,
                                  CNV.sample = NULL,
                                  sample_name = "",
                                  genome = "hg19",
                                  include_panels = c("ideogram", "sv_density",
                                                    "cn_profile", "mechanism_dist")) {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required for dashboard visualization.")
    }

    if (!requireNamespace("gridExtra", quietly = TRUE)) {
        stop("Package 'gridExtra' is required for dashboard layout.")
    }

    if (!inherits(chromoanagenesis_result, "chromoanagenesis")) {
        stop("Input must be a chromoanagenesis result object from detect_chromoanagenesis()")
    }

    # Convert S4 objects to data frames if needed
    if (!is.null(SV.sample) && is(SV.sample, "SVs")) {
        SV.sample <- as(SV.sample, "data.frame")
    }
    if (!is.null(CNV.sample) && is(CNV.sample, "CNVsegs")) {
        CNV.sample <- as(CNV.sample, "data.frame")
    }

    plots <- list()

    # Panel 1: Genome ideogram with mechanism overlay
    if ("ideogram" %in% include_panels) {
        plots$ideogram <- plot_genome_ideogram(chromoanagenesis_result, sample_name, genome)
    }

    # Panel 2: SV density heatmap
    if ("sv_density" %in% include_panels) {
        plots$sv_density <- plot_sv_density_heatmap(SV.sample)
    }

    # Panel 3: Genome-wide CN profile
    if ("cn_profile" %in% include_panels) {
        plots$cn_profile <- plot_genome_wide_cn(CNV.sample, chromoanagenesis_result, genome)
    }

    # Panel 4: Mechanism distribution
    if ("mechanism_dist" %in% include_panels) {
        plots$mechanism_dist <- plot_mechanism_distribution(chromoanagenesis_result)
    }

    # Arrange panels in grid layout
    if (length(plots) == 4) {
        # Default 2x2 layout with ideogram spanning top
        gridExtra::grid.arrange(
            grobs = list(plots$ideogram, plots$sv_density,
                        plots$cn_profile, plots$mechanism_dist),
            layout_matrix = rbind(c(1, 1, 2),
                                 c(3, 3, 2),
                                 c(4, 4, 4)),
            top = grid::textGrob(
                sprintf("Chromoanagenesis Dashboard%s",
                       ifelse(sample_name != "", paste0(" - ", sample_name), "")),
                gp = grid::gpar(fontsize = 16, fontface = "bold")
            )
        )
    } else {
        # Simple vertical arrangement for subset of panels
        do.call(gridExtra::grid.arrange, c(plots, ncol = 1))
    }
}


#' Plot genome ideogram with chromoanagenesis events
#'
#' Creates a karyogram-style visualization showing all chromosomes
#' with chromoanagenesis events highlighted.
#'
#' @param chromoanagenesis_result Chromoanagenesis results
#' @param sample_name Sample name for subtitle
#' @param genome Reference genome ("hg19" or "hg38")
#' @return ggplot2 object
#' @export
plot_genome_ideogram <- function(chromoanagenesis_result, sample_name = "", genome = "hg19") {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    # Extract mechanism information per chromosome
    mech_summary <- get_chromosome_mechanism_summary(chromoanagenesis_result)

    # Get all chromosome info for ideogram
    centro_data <- get_centromere_positions(genome)

    # Merge mechanism data with chromosome data
    chr_data <- merge(centro_data, mech_summary, by = "chrom", all.x = TRUE)

    # Fill in missing values for chromosomes without events
    chr_data$dominant_mechanism[is.na(chr_data$dominant_mechanism)] <- "none"
    chr_data$n_events[is.na(chr_data$n_events)] <- 0

    # Normalize sizes for visualization (scale to max = 1)
    max_size <- max(chr_data$size)
    chr_data$norm_size <- chr_data$size / max_size
    chr_data$norm_centro <- chr_data$centromere / chr_data$size

    # Create chromosome ordering (1-22, X)
    chr_order <- sort_chromosomes(chr_data$chrom)
    chr_data$chrom <- factor(chr_data$chrom, levels = chr_order)

    # Create karyogram layout (2 rows)
    chr_data$row <- ifelse(as.numeric(chr_data$chrom) <= 12 | chr_data$chrom == "X", 1, 2)
    chr_data$col <- ifelse(chr_data$row == 1,
                           as.numeric(chr_data$chrom),
                           as.numeric(chr_data$chrom) - 12)
    chr_data$col[chr_data$chrom == "X"] <- 13

    # Calculate coordinates for chromosome shapes
    chr_shapes <- list()
    for (i in seq_len(nrow(chr_data))) {
        chr <- chr_data$chrom[i]
        col <- chr_data$col[i]
        row_y <- 3 - chr_data$row[i]  # Invert so row 1 is at top

        # Create p-arm (short arm) rectangle
        p_arm <- data.frame(
            chrom = chr,
            xmin = col - 0.3,
            xmax = col + 0.3,
            ymin = row_y,
            ymax = row_y + chr_data$norm_size[i] * chr_data$norm_centro[i] * 0.8,
            arm = "p",
            mechanism = chr_data$dominant_mechanism[i],
            n_events = chr_data$n_events[i]
        )

        # Create centromere (constriction)
        centro <- data.frame(
            chrom = chr,
            xmin = col - 0.15,
            xmax = col + 0.15,
            ymin = p_arm$ymax,
            ymax = p_arm$ymax + 0.05,
            arm = "centro",
            mechanism = chr_data$dominant_mechanism[i],
            n_events = chr_data$n_events[i]
        )

        # Create q-arm (long arm) rectangle
        q_arm <- data.frame(
            chrom = chr,
            xmin = col - 0.3,
            xmax = col + 0.3,
            ymin = centro$ymax,
            ymax = centro$ymax + chr_data$norm_size[i] * (1 - chr_data$norm_centro[i]) * 0.8,
            arm = "q",
            mechanism = chr_data$dominant_mechanism[i],
            n_events = chr_data$n_events[i]
        )

        chr_shapes[[i]] <- rbind(p_arm, centro, q_arm)
    }

    chr_shapes_df <- do.call(rbind, chr_shapes)

    # Color scheme
    mechanism_colors <- c(
        "chromothripsis" = "#E41A1C",
        "chromoplexy" = "#377EB8",
        "chromoanasynthesis" = "#4DAF4A",
        "mixed" = "#984EA3",
        "none" = "gray90"
    )

    # Create ideogram plot
    p <- ggplot2::ggplot(chr_shapes_df)

    # Draw chromosome shapes
    p <- p + ggplot2::geom_rect(
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = mechanism),
        color = "gray40", size = 0.3
    )

    # Add event count labels for chromosomes with events
    label_data <- chr_data[chr_data$n_events > 0, ]
    if (nrow(label_data) > 0) {
        label_data$y_pos <- 3 - label_data$row + 0.4
        p <- p + ggplot2::geom_text(
            data = label_data,
            ggplot2::aes(x = col, y = y_pos, label = n_events),
            color = "black", fontface = "bold", size = 3.5
        )
    }

    # Add chromosome labels
    chr_data$y_label <- 3 - chr_data$row - 0.1
    p <- p + ggplot2::geom_text(
        data = chr_data,
        ggplot2::aes(x = col, y = y_label, label = chrom),
        size = 3, color = "gray30"
    )

    p <- p + ggplot2::scale_fill_manual(
        values = mechanism_colors,
        name = "Mechanism",
        na.value = "gray90"
    )

    # Theme and labels
    p <- p + ggplot2::labs(
        title = "Genome-Wide Chromoanagenesis Events (Karyogram)",
        subtitle = ifelse(sample_name != "", sample_name, "")
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
        legend.position = "right",
        plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 11, hjust = 0.5),
        plot.margin = ggplot2::margin(10, 10, 10, 10)
    ) +
    ggplot2::coord_fixed(ratio = 0.5)

    return(p)
}


#' Plot SV density heatmap by chromosome and type
#'
#' Creates a heatmap showing the density of different SV types
#' across chromosomes.
#'
#' @param SV.sample SV data (SVs object or data frame)
#' @return ggplot2 object
#' @export
plot_sv_density_heatmap <- function(SV.sample) {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    # Convert S4 object if needed
    if (!is.null(SV.sample) && is(SV.sample, "SVs")) {
        sv_data <- as(SV.sample, "data.frame")
    } else {
        sv_data <- SV.sample
    }

    if (is.null(sv_data) || nrow(sv_data) == 0) {
        p <- ggplot2::ggplot() +
            ggplot2::annotate("text", x = 0.5, y = 0.5,
                            label = "No SV data available",
                            size = 5, color = "gray50") +
            ggplot2::theme_void()
        return(p)
    }

    # Count SVs by chromosome and type
    sv_counts <- as.data.frame(table(
        chrom = sv_data$chrom1,
        svtype = sv_data$SVtype
    ))
    names(sv_counts) <- c("chrom", "svtype", "count")

    # Sort chromosomes
    chr_order <- sort_chromosomes(unique(as.character(sv_counts$chrom)))
    sv_counts$chrom <- factor(sv_counts$chrom, levels = chr_order)

    # Create heatmap
    p <- ggplot2::ggplot(sv_counts,
                        ggplot2::aes(x = chrom, y = svtype, fill = count))

    p <- p + ggplot2::geom_tile(color = "white", size = 0.5)

    p <- p + ggplot2::geom_text(
        ggplot2::aes(label = ifelse(count > 0, as.character(count), "")),
        color = "white", size = 3, fontface = "bold"
    )

    # Color scale
    p <- p + ggplot2::scale_fill_gradient(
        low = "gray90",
        high = "#d73027",
        name = "SV Count"
    )

    # Theme and labels
    p <- p + ggplot2::labs(
        title = "SV Density by Type and Chromosome",
        x = "Chromosome",
        y = "SV Type"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        plot.title = ggplot2::element_text(face = "bold", size = 12)
    )

    return(p)
}


#' Plot genome-wide copy number profile
#'
#' Creates a linear plot showing copy number across all chromosomes.
#'
#' @param CNV.sample CNV data (CNVsegs object or data frame)
#' @param chromoanagenesis_result Chromoanagenesis results (for annotation)
#' @param genome Reference genome ("hg19" or "hg38")
#' @return ggplot2 object
#' @export
plot_genome_wide_cn <- function(CNV.sample, chromoanagenesis_result = NULL, genome = "hg19") {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    # Convert S4 object if needed
    if (!is.null(CNV.sample) && is(CNV.sample, "CNVsegs")) {
        cnv_data <- as(CNV.sample, "data.frame")
    } else {
        cnv_data <- CNV.sample
    }

    if (is.null(cnv_data) || nrow(cnv_data) == 0) {
        p <- ggplot2::ggplot() +
            ggplot2::annotate("text", x = 0.5, y = 0.5,
                            label = "No CNV data available",
                            size = 5, color = "gray50") +
            ggplot2::theme_void()
        return(p)
    }

    # Sort chromosomes
    chr_order <- sort_chromosomes(unique(as.character(cnv_data$chrom)))
    cnv_data$chrom <- factor(cnv_data$chrom, levels = chr_order)

    # Calculate cumulative genomic positions for linear view
    cnv_data_linear <- calculate_linear_positions(cnv_data, genome)

    # Identify chromoanagenesis regions (if result provided)
    if (!is.null(chromoanagenesis_result)) {
        cnv_data_linear <- annotate_chromoanagenesis_regions(
            cnv_data_linear,
            chromoanagenesis_result
        )
    } else {
        # No annotation available
        cnv_data_linear$has_chromoanagenesis <- FALSE
    }

    # Create CN profile plot
    p <- ggplot2::ggplot(cnv_data_linear)

    # Add chromosome alternating backgrounds
    chr_bounds <- get_chromosome_boundaries(cnv_data_linear)
    for (i in seq_len(nrow(chr_bounds))) {
        if (i %% 2 == 0) {
            p <- p + ggplot2::annotate(
                "rect",
                xmin = chr_bounds$start[i],
                xmax = chr_bounds$end[i],
                ymin = -Inf, ymax = Inf,
                fill = "gray95", alpha = 0.5
            )
        }
    }

    # Add CN segments
    p <- p + ggplot2::geom_segment(
        ggplot2::aes(x = linear_start, xend = linear_end,
                    y = total_cn, yend = total_cn,
                    color = has_chromoanagenesis),
        size = 1.5
    )

    # Color scale
    p <- p + ggplot2::scale_color_manual(
        values = c("FALSE" = "gray50", "TRUE" = "#E41A1C"),
        name = "Chromoanagenesis",
        labels = c("FALSE" = "Normal", "TRUE" = "Event")
    )

    # Add chromosome labels
    p <- p + ggplot2::scale_x_continuous(
        breaks = chr_bounds$mid,
        labels = chr_bounds$chrom,
        expand = c(0.01, 0.01)
    )

    # Add reference line at CN=2
    p <- p + ggplot2::geom_hline(
        yintercept = 2,
        linetype = "dashed",
        color = "gray40",
        size = 0.5
    )

    # Theme and labels
    p <- p + ggplot2::labs(
        title = "Genome-Wide Copy Number Profile",
        x = "Chromosome",
        y = "Copy Number"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
        plot.title = ggplot2::element_text(face = "bold", size = 12),
        legend.position = "right"
    ) +
    ggplot2::coord_cartesian(ylim = c(0, max(cnv_data_linear$total_cn, 6)))

    return(p)
}


#' Plot mechanism distribution summary
#'
#' Creates a summary visualization of mechanism proportions.
#'
#' @param chromoanagenesis_result Chromoanagenesis results
#' @return ggplot2 object
#' @export
plot_mechanism_distribution <- function(chromoanagenesis_result) {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    # Get mechanism counts
    mech_counts <- data.frame(
        mechanism = character(),
        count = numeric(),
        classification = character(),
        stringsAsFactors = FALSE
    )

    # Chromothripsis
    if (!is.null(chromoanagenesis_result$chromothripsis)) {
        ct_class <- chromoanagenesis_result$chromothripsis$classification
        if (!is.null(ct_class) && nrow(ct_class) > 0) {
            ct_summary <- table(ct_class$classification)
            for (cls in names(ct_summary)) {
                if (ct_summary[cls] > 0 && cls %in% c("High confidence", "Low confidence")) {
                    mech_counts <- rbind(mech_counts, data.frame(
                        mechanism = "Chromothripsis",
                        count = as.numeric(ct_summary[cls]),
                        classification = cls,
                        stringsAsFactors = FALSE
                    ))
                }
            }
        }
    }

    # Chromoplexy
    if (!is.null(chromoanagenesis_result$chromoplexy)) {
        cp_summary <- chromoanagenesis_result$chromoplexy$summary
        if (!is.null(cp_summary) && nrow(cp_summary) > 0) {
            cp_counts <- table(cp_summary$classification)
            for (cls in names(cp_counts)) {
                if (cp_counts[cls] > 0 && cls %in% c("Likely chromoplexy", "Possible chromoplexy")) {
                    mech_counts <- rbind(mech_counts, data.frame(
                        mechanism = "Chromoplexy",
                        count = as.numeric(cp_counts[cls]),
                        classification = cls,
                        stringsAsFactors = FALSE
                    ))
                }
            }
        }
    }

    # Chromoanasynthesis
    if (!is.null(chromoanagenesis_result$chromoanasynthesis)) {
        cs_summary <- chromoanagenesis_result$chromoanasynthesis$summary
        if (!is.null(cs_summary) && nrow(cs_summary) > 0) {
            cs_counts <- table(cs_summary$classification)
            for (cls in names(cs_counts)) {
                if (cs_counts[cls] > 0 && cls %in% c("Likely chromoanasynthesis", "Possible chromoanasynthesis")) {
                    mech_counts <- rbind(mech_counts, data.frame(
                        mechanism = "Chromoanasynthesis",
                        count = as.numeric(cs_counts[cls]),
                        classification = cls,
                        stringsAsFactors = FALSE
                    ))
                }
            }
        }
    }

    if (nrow(mech_counts) == 0) {
        p <- ggplot2::ggplot() +
            ggplot2::annotate("text", x = 0.5, y = 0.5,
                            label = "No chromoanagenesis events detected",
                            size = 5, color = "gray50") +
            ggplot2::theme_void()
        return(p)
    }

    # Simplify classification labels
    mech_counts$class_simple <- ifelse(
        grepl("Likely", mech_counts$classification),
        "Likely",
        "Possible"
    )

    # Colors
    mechanism_colors <- c(
        "Chromothripsis" = "#E41A1C",
        "Chromoplexy" = "#377EB8",
        "Chromoanasynthesis" = "#4DAF4A"
    )

    # Create stacked bar chart
    p <- ggplot2::ggplot(mech_counts,
                        ggplot2::aes(x = mechanism, y = count, fill = mechanism, alpha = class_simple))

    p <- p + ggplot2::geom_col(color = "white", size = 0.5)

    p <- p + ggplot2::geom_text(
        ggplot2::aes(label = count),
        position = ggplot2::position_stack(vjust = 0.5),
        color = "white",
        fontface = "bold",
        size = 5
    )

    p <- p + ggplot2::scale_fill_manual(
        values = mechanism_colors,
        name = "Mechanism"
    )

    p <- p + ggplot2::scale_alpha_manual(
        values = c("Likely" = 1.0, "Possible" = 0.6),
        name = "Classification"
    )

    # Theme and labels
    p <- p + ggplot2::labs(
        title = "Mechanism Distribution",
        x = "",
        y = "Event Count"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold", size = 12),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        legend.position = "right"
    )

    return(p)
}
