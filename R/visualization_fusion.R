#' Plot Fusion Architecture (Arriba-style)
#' 
#' Creates a detailed transcript-level diagram showing how two genes are joined,
#' including 3D-rendered exons and junction connections.
#'
#' @param results A chromoanagenesis result object.
#' @param fusion_name Name of the fusion (e.g., "FGFR3::TACC3").
#' @param txdb A TxDb object for exon coordinates.
#' @param output_file Optional PDF output filename.
#' @export
plot_fusion_arriba <- function(results, fusion_name, txdb, output_file = NULL) {
    # [Implementation from R/fusion_graphics_engine.R with latest fixes]
    if (is.null(results$fusions)) stop("No fusion data found.")
    fus <- results$fusions[results$fusions$fusion_name == fusion_name, ]
    if (nrow(fus) == 0) stop("Fusion not found.")
    fus <- fus[1, ]

    # Internal helper: Change color brightness for 3D effect
    .change_brightness <- function(color, delta) {
        rgb_val <- grDevices::col2rgb(color)
        grDevices::rgb(
            min(255, max(0, rgb_val["red", ] + delta)),
            min(255, max(0, rgb_val["green", ] + delta)),
            min(255, max(0, rgb_val["blue", ] + delta)),
            maxColorValue = 255
        )
    }

    .draw_exon <- function(left, right, y, color, label = "", type = "exon") {
        exon_height <- 0.03
        dark_color <- .change_brightness(color, -100)
        if (type == "CDS") {
            graphics::rect(left, y + exon_height, right, y + exon_height/2, col = color, border = NA)
            graphics::rect(left, y - exon_height, right, y - exon_height/2, col = color, border = NA)
            graphics::lines(c(left, left, right, right), c(y + exon_height/2, y + exon_height, y + exon_height, y + exon_height/2), col = dark_color)
            graphics::lines(c(left, left, right, right), c(y - exon_height/2, y - exon_height, y - exon_height, y - exon_height/2), col = dark_color)
        } else {
            graphics::rect(left, y + exon_height/2, right, y - exon_height/2, col = color, border = dark_color)
        }
        graphics::text((left + right)/2, y, label, cex = 0.6)
    }

    .draw_strand <- function(left, right, y, color, strand) {
        graphics::lines(c(left, right), c(y, y), col = color, lwd = 2)
        if (abs(right - left) > 0.01) {
            points <- seq(left + 0.02, right - 0.02, length.out = max(2, as.integer(abs(right-left)/0.05)))
            for (p in points) graphics::arrows(p, y, p + 0.005 * ifelse(strand == "+", 1, -1), y, col = color, length = 0.05, lwd = 2, angle = 60)
        }
    }

    # Fetch structures
    get_tx_data <- function(gene_name) {
        if (is.null(gene_name) || is.na(gene_name) || gene_name == "NA" || gene_name == "") {
            return(NULL)
        }
        
        id <- tryCatch({
            AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, 
                                 keys = as.character(gene_name), 
                                 column = "ENTREZID", 
                                 keytype = "SYMBOL", 
                                 multiVals = "first")
        }, error = function(e) {
            return(NA)
        })
        
        if (is.na(id)) return(NULL)
        txs <- GenomicFeatures::transcriptsBy(txdb, by = "gene")[[id]]
        if (is.null(txs) || length(txs) == 0) return(NULL)
        best_tx <- as.character(txs$tx_id[which.max(GenomicRanges::width(txs))])
        return(GenomicFeatures::exonsBy(txdb, by = "tx")[[best_tx]])
    }

    exs1 <- get_tx_data(fus$gene_5p); exs2 <- get_tx_data(fus$gene_3p)
    if (is.null(exs1) || is.null(exs2)) stop("Exon models not found.")

    transform_to_plot <- function(exs, bp, x_start, x_end) {
        g_start <- min(GenomicRanges::start(exs)); g_end <- max(GenomicRanges::end(exs))
        scale <- (x_end - x_start) / max(g_end - g_start + 1, 1)
        strand <- as.character(GenomicRanges::strand(exs)[1])
        df <- as.data.frame(exs)
        if (strand == "-") {
            df$p_start <- x_start + (g_end - df$end) * scale
            df$p_end <- x_start + (g_end - df$start) * scale
            p_bp <- x_start + (g_end - bp) * scale
        } else {
            df$p_start <- x_start + (df$start - g_start) * scale
            df$p_end <- x_start + (df$end - g_start) * scale
            p_bp <- x_start + (bp - g_start) * scale
        }
        return(list(exs = df, bp = p_bp, strand = strand))
    }

    bp1_gr <- GenomicRanges::GRanges(fus$chrom1, IRanges::IRanges(fus$pos1, width=1))
    is_bp1_gene1 <- length(GenomicRanges::findOverlaps(bp1_gr, exs1, maxgap = 50000)) > 0
    p_bp1 <- if (is_bp1_gene1) fus$pos1 else fus$pos2
    p_bp2 <- if (is_bp1_gene1) fus$pos2 else fus$pos1

    g1 <- transform_to_plot(exs1, p_bp1, 0, 0.45); g2 <- transform_to_plot(exs2, p_bp2, 0.55, 1.0)

    if (!is.null(output_file)) grDevices::pdf(output_file, width = 11, height = 6)
    graphics::par(mar = c(1, 1, 4, 1))
    graphics::plot(0, 0, type = "n", xlim = c(-0.05, 1.05), ylim = c(0.1, 1.2), bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    graphics::title(main = paste("Fusion Architecture:", fusion_name), font.main = 2, cex.main = 1.5)
    graphics::mtext(paste("Status:", fus$status, "| Genome:", GenomeInfoDb::genome(txdb)[1]), line = 0.5)

    y_top <- 0.9; y_bot <- 0.4
    .draw_strand(0, 0.45, y_top, "#e5a5a5", g1$strand)
    for(j in 1:nrow(g1$exs)) .draw_exon(g1$exs$p_start[j], g1$exs$p_end[j], y_top, "#e5a5a5", j, "exon")
    .draw_strand(0.55, 1.0, y_top, "#a7c4e5", g2$strand)
    for(j in 1:nrow(g2$exs)) .draw_exon(g2$exs$p_start[j], g2$exs$p_end[j], y_top, "#a7c4e5", j, "exon")
    
    graphics::lines(c(g1$bp, g1$bp, 0.45, 0.55, g2$bp, g2$bp), c(y_top, 0.65, 0.65, 0.65, 0.65, y_top), col = "#FFC107", lwd = 1.5, lty = 3)
    graphics::text(0.5, 0.7, "Breakpoint Join", col = "#FF9800", font = 4, cex = 0.9)
    .draw_strand(0.2, 0.5, y_bot, "#e5a5a5", "+"); .draw_strand(0.5, 0.8, y_bot, "#a7c4e5", "+")
    graphics::text(0.5, y_bot - 0.15, "Resulting Chimeric Transcript", font = 3, cex = 1)
    graphics::text(0.22, y_top + 0.12, fus$gene_5p, font = 2, col = "#d32f2f", cex = 1.2)
    graphics::text(0.77, y_top + 0.12, fus$gene_3p, font = 2, col = "#1976d2", cex = 1.2)

    if (!is.null(output_file)) grDevices::dev.off()
}

#' Plot Fusion Repair Context (ggplot2 version)
#' 
#' @param results Analysis results.
#' @param fusion_name Fusion name.
#' @param gene_granges GRanges object with symbols.
#' @export
plot_fusion_repair_context <- function(results, fusion_name, gene_granges) {
    # [Implementation from R/fusion_repair_visualization.R]
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
    if (!requireNamespace("patchwork", quietly = TRUE)) stop("patchwork is required.")
    
    fus <- results$fusions[results$fusions$fusion_name == fusion_name, ]
    if (nrow(fus) == 0) stop("Fusion ", fusion_name, " not found.")
    fus <- fus[1, ]
    sv_id <- fus$sv_id
    
    repair_info <- results$repair_mechanisms$repair_mechanisms[results$repair_mechanisms$repair_mechanisms$sv_id == sv_id, ]
    seq_info <- results$repair_mechanisms$sequences[results$repair_mechanisms$sequences$sv_id == sv_id, ]
    mh_info <- results$repair_mechanisms$microhomology[results$repair_mechanisms$microhomology$sv_id == sv_id, ]
    
    g1_obj <- gene_granges[mcols(gene_granges)$symbol == fus$gene_5p | mcols(gene_granges)$symbol == fus$gene_3p]
    
    p1 <- ggplot2::ggplot() +
        ggplot2::geom_rect(data = as.data.frame(g1_obj), ggplot2::aes(xmin = start, xmax = end, ymin = -0.2, ymax = 0.2, fill = symbol), alpha = 0.7) +
        ggplot2::geom_vline(xintercept = c(fus$pos1, fus$pos2), linetype = "dashed", color = "red") +
        ggplot2::facet_wrap(~seqnames, scales = "free_x") +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = paste("Fusion Structure:", fusion_name), subtitle = paste("Mechanism:", repair_info$repair_mechanism), x = "Genomic Position", fill = "Gene") +
        ggplot2::theme(legend.position = "bottom")

    mh_len <- ifelse(is.na(mh_info$microhomology_length), 0, mh_info$microhomology_length)
    mh_seq <- ifelse(is.na(mh_info$microhomology_seq), "", mh_info$microhomology_seq)
    left <- substr(seq_info$s1_kept, nchar(seq_info$s1_kept)-14, nchar(seq_info$s1_kept))
    right <- substr(seq_info$s2_kept, nchar(seq_info$s2_kept)-14, nchar(seq_info$s2_kept))
    
    df_seq <- data.frame(x = 1, y = c(3, 2, 1), text = c(paste0("5' Partner (", fus$gene_5p, "): ", left, " [", mh_seq, "]"), 
                                                       paste0(paste(rep(" ", 30), collapse=""), "  ", paste(rep("|", mh_len), collapse="")),
                                                       paste0("3' Partner (", fus$gene_3p, "): ", " [", mh_seq, "] ", right)))
    
    p2 <- ggplot2::ggplot(df_seq, ggplot2::aes(x = x, y = y, label = text)) +
        ggplot2::geom_text(family = "mono", size = 4, hjust = 0) +
        ggplot2::xlim(1, 2) + ggplot2::ylim(0, 4) +
        ggplot2::theme_void() +
        ggplot2::labs(title = "Breakpoint Junction Evidence", subtitle = paste0("Evidence: ", repair_info$evidence))

    return(p1 / p2 + patchwork::plot_layout(heights = c(2, 1)))
}

