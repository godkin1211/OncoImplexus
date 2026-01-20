#' Plot repair mechanisms by SV type
#'
#' Creates a stacked bar plot showing the distribution of DNA repair mechanisms
#' for each structural variant type.
#'
#' @param repair_data Result object from analyze_breakpoint_sequences()
#' @param SV.sample Original SV data used for analysis
#' @return A ggplot2 object
#' @export
plot_repair_by_svtype <- function(repair_data, SV.sample) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    # Extract repair mechanism and SV type
    if (isS4(SV.sample)) {
        sv_df <- as(SV.sample, "data.frame")
    } else {
        sv_df <- SV.sample
    }

    repair_df <- repair_data$repair_mechanisms

    # Merge
    plot_df <- data.frame(
        SVtype = sv_df$SVtype,
        Mechanism = repair_df$repair_mechanism,
        stringsAsFactors = FALSE
    )

    # Clean up names for plotting
    plot_df$SVtype <- factor(plot_df$SVtype, levels = c("DEL", "DUP", "h2hINV", "t2tINV", "TRA"))

    # Define mechanism colors
    mech_colors <- c(
        "NHEJ" = "#333333", # Dark gray/black for classic NHEJ
        "MMEJ" = "#E41A1C", # Red
        "MMBIR/FoSTeS" = "#4DAF4A", # Green
        "SSA" = "#377EB8", # Blue
        "Unknown" = "gray80"
    )

    # Create plot
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = SVtype, fill = Mechanism)) +
        ggplot2::geom_bar(position = "fill", color = "white", size = 0.2) +
        ggplot2::scale_y_continuous(labels = scales::percent) +
        ggplot2::scale_fill_manual(values = mech_colors) +
        ggplot2::labs(
            title = "DNA Repair Mechanisms by SV Type",
            subtitle = "Proportional distribution of inferred repair pathways",
            x = "Structural Variant Type",
            y = "Proportion of Breakpoints",
            fill = "Repair Pathway"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 14),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )

    return(p)
}

#' Plot inferred repair mechanisms distribution
#'
#' @param repair_data Result object from analyze_breakpoint_sequences()
#' @return A ggplot2 object
#' @export
plot_repair_mechanisms <- function(repair_data) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")

    mech_stats <- repair_data$summary$repair_mechanisms
    colnames(mech_stats) <- c("Mechanism", "Count")

    # Use consistent colors
    mech_colors <- c(
        "NHEJ" = "#333333", "MMEJ" = "#E41A1C", "MMBIR/FoSTeS" = "#4DAF4A",
        "SSA" = "#377EB8", "Unknown" = "gray80"
    )

    ggplot2::ggplot(mech_stats, ggplot2::aes(x = reorder(Mechanism, -Count), y = Count, fill = Mechanism)) +
        ggplot2::geom_col() +
        ggplot2::scale_fill_manual(values = mech_colors) +
        ggplot2::labs(title = "Inferred DNA Repair Mechanisms", x = "Repair Pathway", y = "Number of Breakpoints") +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none")
}

#' Plot microhomology length distribution
#'
#' @param repair_data Result object from analyze_breakpoint_sequences()
#' @return A ggplot2 object
#' @export
plot_microhomology_distribution <- function(repair_data) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")

    mh_df <- repair_data$microhomology
    mh_data <- mh_df[!is.na(mh_df$microhomology_length) & mh_df$microhomology_length > 0, ]

    if (nrow(mh_data) == 0) {
        return(ggplot2::ggplot() +
            ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No microhomology detected") +
            ggplot2::theme_void())
    }

    ggplot2::ggplot(mh_data, ggplot2::aes(x = microhomology_length)) +
        ggplot2::geom_histogram(binwidth = 1, fill = "#673ab7", color = "white") +
        ggplot2::scale_x_continuous(breaks = seq(0, max(mh_data$microhomology_length, 5), 2)) +
        ggplot2::labs(title = "Microhomology Lengths", x = "Length (bp)", y = "Frequency") +
        ggplot2::theme_minimal()
}

#' Generate a combined breakpoint sequence analysis report plot
#'
#' @param repair_data Result object from analyze_breakpoint_sequences()
#' @param SV.sample Original SV data
#' @return A grid arrangement of plots
#' @export
plot_breakpoint_report <- function(repair_data, SV.sample) {
    if (!requireNamespace("gridExtra", quietly = TRUE)) stop("gridExtra required")

    p1 <- plot_repair_mechanisms(repair_data)
    p2 <- plot_microhomology_distribution(repair_data)
    p3 <- plot_repair_by_svtype(repair_data, SV.sample)

    gridExtra::grid.arrange(
        p1, p2, p3,
        layout_matrix = rbind(c(1, 2), c(3, 3)),
        top = grid::textGrob("Breakpoint Sequence Analysis Summary", gp = grid::gpar(fontsize = 16, fontface = "bold"))
    )
}

#' Get representative sequence preview HTML
#'
#' Generates an HTML snippet showing the sequence-level evidence for a breakpoint.
#' Highlights microhomology and insertions.
#'
#' @param sv_id Integer ID of the structural variant
#' @param repair_data Result object from analyze_breakpoint_sequences()
#' @return Character string containing HTML
#' @param sv_id The ID of the structural variant
#' @param repair_data Result object from analyze_breakpoint_sequences()
#' @param chrom1 Optional chromosome name for fallback coordinate matching
#' @param pos1 Optional genomic position for fallback coordinate matching
#' @return An HTML string with highlighted microhomology/insertions
#' @export
get_fusion_seq_html <- function(sv_id, repair_data, chrom1 = NULL, pos1 = NULL) {
    if (is.null(repair_data$sequences) || is.null(repair_data$microhomology)) {
        return("")
    }

    target_id <- as.character(sv_id)
    seq_idx <- which(as.character(repair_data$sequences$sv_id) == target_id)

    # If no ID match, try coordinate-based match
    if (length(seq_idx) == 0 && !is.null(chrom1) && !is.null(pos1)) {
        # Support both chr-prefixed and non-prefixed matching
        c1_alt <- if (grepl("^chr", chrom1)) gsub("^chr", "", chrom1) else paste0("chr", chrom1)

        seq_idx <- which(
            (as.character(repair_data$sequences$chrom1) == as.character(chrom1) |
                as.character(repair_data$sequences$chrom1) == as.character(c1_alt)) &
                abs(as.numeric(repair_data$sequences$pos1) - as.numeric(pos1)) <= 2
        )
    }

    if (length(seq_idx) == 0) {
        return("<span class='text-muted'>Sequence unavailable</span>")
    }

    seq <- repair_data$sequences[seq_idx[1], ]
    # Also find corresponding MH and Insertion data by the SAME index to be consistent
    mh <- repair_data$microhomology[seq_idx[1], ]
    ins <- repair_data$insertions[seq_idx[1], ]

    if (!seq$has_sequence) {
        return("<span class='text-muted'>Sequence unavailable</span>")
    }

    # Extract pieces
    s1_kept <- seq$s1_kept
    s2_kept <- seq$s2_kept

    # Identify MH and Ins
    mh_len <- if (!is.na(mh$microhomology_length)) mh$microhomology_length else 0
    mh_seq <- if (!is.na(mh$microhomology_seq)) mh$microhomology_seq else ""

    ins_len <- if (!is.na(ins$insertion_length)) ins$insertion_length else 0
    ins_seq <- if (!is.na(ins$insertion_seq)) ins$insertion_seq else ""

    # Construct stylized string
    # We show [Kept1] [MH/Ins] [Kept2]
    # For visualization, we'll show the last 15bp of Side 1 and first 15bp of Side 2
    ctx <- 15
    k1_disp <- substr(s1_kept, nchar(s1_kept) - ctx + 1, nchar(s1_kept))
    k2_disp <- substr(s2_kept, 1, ctx)

    # HTML styling
    html <- "<code style='font-family: monospace; font-size: 0.9em;'>"
    html <- paste0(html, k1_disp)

    if (mh_len > 0) {
        html <- paste0(html, "<span style='background-color: #d1c4e9; font-weight: bold; border-bottom: 2px solid #673ab7;' title='Microhomology'>", mh_seq, "</span>")
    }

    html <- paste0(html, "<span style='color: red; font-weight: bold;'>|</span>")

    if (ins_len > 0) {
        html <- paste0(html, "<span style='background-color: #fff9c4; font-weight: bold; border-bottom: 2px solid #fbc02d;' title='Insertion'>", ins_seq, "</span>")
    }

    html <- paste0(html, k2_disp)
    html <- paste0(html, "</code>")

    return(html)
}
