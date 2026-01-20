#' Circos Plot Visualization for Chromoanagenesis
#'
#' Functions for creating circos plots to visualize genome-wide
#' chromoanagenesis events using the circlize package.
#'
#' @name circos_viz
#' @keywords internal


#' Plot genome-wide chromoanagenesis as circos
#'
#' Creates a circos plot showing:
#' - Chromosome ideogram (outer track)
#' - Copy number variations (middle track)
#' - Chromoanagenesis regions (highlighted)
#' - Structural variants (inner links)
#'
#' @param chromoanagenesis_result Result from detect_chromoanagenesis()
#' @param SV.sample Original SV data
#' @param CNV.sample Original CNV data
#' @param genome Reference genome ("hg19" or "hg38")
#' @param highlight_mechanisms Which mechanisms to highlight (default: all)
#' @param sample_name Sample name for title
#' @return NULL (plot is created directly)
#' @export
plot_chromoanagenesis_circos <- function(chromoanagenesis_result,
                                        SV.sample,
                                        CNV.sample,
                                        genome = "hg19",
                                        highlight_mechanisms = c("chromothripsis", "chromoplexy", "chromoanasynthesis"),
                                        sample_name = "",
                                        verbose = TRUE) {

    if (!requireNamespace("circlize", quietly = TRUE)) {
        stop("Package 'circlize' is required for circos plots. Install with: install.packages('circlize')")
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

    # Defensive: Remove rows with NA coordinates or chromosomes
    if (!is.null(sv_data)) {
        sv_data <- sv_data[!is.na(sv_data$chrom1) & !is.na(sv_data$chrom2) & 
                          !is.na(sv_data$pos1) & !is.na(sv_data$pos2), ]
    }
    if (!is.null(cnv_data)) {
        cnv_data <- cnv_data[!is.na(cnv_data$chrom) & !is.na(cnv_data$start) & !is.na(cnv_data$end), ]
    }

    # Get chromosome sizes (base set)
    chr_sizes <- get_chromosome_sizes()

    # Check which chromosomes are actually in the data and add missing ones
    all_chroms_sv <- if (!is.null(sv_data)) unique(c(sv_data$chrom1, sv_data$chrom2)) else c()
    all_chroms_cnv <- if (!is.null(cnv_data)) unique(cnv_data$chrom) else c()
    all_chroms <- unique(na.omit(c(all_chroms_sv, all_chroms_cnv)))

    # Add Y chromosome if present in data but not in chr_sizes
    if ("Y" %in% all_chroms && !"Y" %in% names(chr_sizes)) {
        chr_sizes$Y <- 59373566  # hg19 Y chromosome size
    }

    # Filter to chromosomes present in chr_sizes
    chroms_to_use <- names(chr_sizes)[names(chr_sizes) %in% all_chroms]
    
    if (length(chroms_to_use) == 0) {
        warning("No valid chromosomes found for Circos plot.")
        return(NULL)
    }

    chr_df <- data.frame(
        chr = paste0("chr", chroms_to_use),
        start = 1,
        end = unlist(chr_sizes[chroms_to_use]),
        stringsAsFactors = FALSE
    )

    # Extract chromoanagenesis regions
    ca_regions <- extract_chromoanagenesis_regions(chromoanagenesis_result, highlight_mechanisms)

    # Initialize circos
    circlize::circos.clear()
    circlize::circos.par(
        start.degree = 90,
        gap.degree = 2,
        track.margin = c(0.01, 0.01),
        cell.padding = c(0, 0, 0, 0)
    )

    circlize::circos.initialize(
        factors = chr_df$chr,
        xlim = as.matrix(chr_df[, c("start", "end")])
    )

    # Track 1: Chromosome ideogram
    circlize::circos.track(
        ylim = c(0, 1),
        panel.fun = function(x, y) {
            chr <- circlize::get.cell.meta.data("sector.index")
            xlim <- circlize::get.cell.meta.data("xlim")
            ylim <- circlize::get.cell.meta.data("ylim")

            # Draw chromosome background
            circlize::circos.rect(
                xlim[1], 0, xlim[2], 1,
                col = "gray90",
                border = "gray60"
            )

            # Add chromosome label
            circlize::circos.text(
                mean(xlim), mean(ylim),
                labels = gsub("chr", "", chr),
                cex = 0.7,
                facing = "bending.inside",
                niceFacing = TRUE
            )
        },
        bg.border = NA,
        track.height = 0.05
    )

    # Track 2: Copy number
    if (!is.null(cnv_data) && nrow(cnv_data) > 0) {
        # Prepare CN data for circos
        cnv_circos <- cnv_data
        cnv_circos$chr <- paste0("chr", cnv_circos$chrom)

        # Normalize CN for visualization (cap at 6 for display)
        cnv_circos$cn_display <- pmin(cnv_circos$total_cn, 6)

        circlize::circos.track(
            ylim = c(0, 6),
            panel.fun = function(x, y) {
                chr <- circlize::get.cell.meta.data("sector.index")
                chr_cn <- cnv_circos[cnv_circos$chr == chr, ]
                chr_cn <- chr_cn[!is.na(chr_cn$cn_display), ]

                if (nrow(chr_cn) > 0) {
                    # Draw CN segments
                    for (i in seq_len(nrow(chr_cn))) {
                        # Color by CN state
                        cn_val <- chr_cn$cn_display[i]
                        
                        # Handle potential NAs although we filtered above
                        if (is.na(cn_val)) next
                        
                        if (cn_val < 1.5) {
                            col <- "#2166ac"  # Blue for loss
                        } else if (cn_val > 2.5) {
                            col <- "#b2182b"  # Red for gain
                        } else {
                            col <- "gray70"    # Gray for neutral
                        }

                        circlize::circos.rect(
                            chr_cn$start[i], 0,
                            chr_cn$end[i], cn_val,
                            col = col,
                            border = NA
                        )
                    }

                    # Add reference line at CN=2
                    circlize::circos.lines(
                        circlize::get.cell.meta.data("xlim"),
                        c(2, 2),
                        col = "gray40",
                        lty = 2,
                        lwd = 0.5
                    )
                }
            },
            bg.border = "gray60",
            track.height = 0.15
        )
    }

    # Track 3: Chromoanagenesis regions
    if (!is.null(ca_regions) && nrow(ca_regions) > 0) {
        mechanism_colors <- c(
            "chromothripsis" = "#E41A1C",
            "chromoplexy" = "#377EB8",
            "chromoanasynthesis" = "#4DAF4A"
        )

        circlize::circos.track(
            ylim = c(0, 1),
            panel.fun = function(x, y) {
                chr <- circlize::get.cell.meta.data("sector.index")
                chr_regions <- ca_regions[ca_regions$chr == chr, ]

                if (nrow(chr_regions) > 0) {
                    for (i in seq_len(nrow(chr_regions))) {
                        circlize::circos.rect(
                            chr_regions$start[i], 0,
                            chr_regions$end[i], 1,
                            col = mechanism_colors[chr_regions$mechanism[i]],
                            border = NA
                        )
                    }
                }
            },
            bg.col = "white",
            bg.border = "gray60",
            track.height = 0.05
        )
    }

    # Add SV links (intrachromosomal and interchromosomal)
    if (!is.null(sv_data) && nrow(sv_data) > 0) {
        # Prepare SV data
        sv_circos <- sv_data
        sv_circos$chr1 <- paste0("chr", sv_circos$chrom1)
        sv_circos$chr2 <- paste0("chr", sv_circos$chrom2)

        # Color by SV type
        sv_colors <- c(
            "DEL" = "#f4a582",
            "DUP" = "#92c5de",
            "t2tINV" = "#66c2a5",
            "h2hINV" = "#fc8d62",
            "TRA" = "#984ea3"
        )

        if (verbose) {
            cat(sprintf("\nCircos plot: Total SVs = %d\n", nrow(sv_circos)))
            cat("SV types:\n")
            print(table(sv_circos$SVtype))
            n_inter <- sum(sv_circos$chrom1 != sv_circos$chrom2)
            cat(sprintf("Interchromosomal SVs: %d\n", n_inter))
        }

        # Determine which SVs to show (limit for visibility)
        # Prioritize SVs in chromoanagenesis regions
        sv_in_regions <- merge_sv_with_regions(sv_circos, ca_regions)
        
        # Defensive: Handle NAs in in_ca_region
        sv_in_regions$in_ca_region[is.na(sv_in_regions$in_ca_region)] <- FALSE

        # Show all SVs in chromoanagenesis regions + sample of others
        sv_priority <- sv_in_regions[sv_in_regions$in_ca_region, ]
        sv_other <- sv_in_regions[!sv_in_regions$in_ca_region, ]

        if (verbose) {
            cat(sprintf("SVs in chromoanagenesis regions: %d\n", nrow(sv_priority)))
            cat(sprintf("Other SVs: %d\n", nrow(sv_other)))
        }

        if (nrow(sv_other) > 100) {
            # Sample 100 random other SVs for clarity
            sv_other <- sv_other[sample(nrow(sv_other), 100), ]
            if (verbose) {
                cat(sprintf("Sampled %d other SVs for display\n", nrow(sv_other)))
            }
        }

        sv_to_plot <- rbind(sv_priority, sv_other)

        if (verbose) {
            cat(sprintf("Total SVs to plot: %d\n", nrow(sv_to_plot)))
            n_inter_plot <- sum(sv_to_plot$chr1 != sv_to_plot$chr2)
            cat(sprintf("  - Interchromosomal: %d\n", n_inter_plot))
            cat(sprintf("  - Intrachromosomal: %d\n\n", nrow(sv_to_plot) - n_inter_plot))
        }

        # Draw links
        for (i in seq_len(nrow(sv_to_plot))) {
            sv <- sv_to_plot[i, ]
            col <- sv_colors[sv$SVtype]
            if (is.na(col)) col <- "gray50"

            # Make chromoanagenesis SVs more prominent
            if (sv$in_ca_region) {
                lwd <- 1.5
                alpha <- 0.6
            } else {
                lwd <- 0.5
                alpha <- 0.2
            }

            col_alpha <- adjustcolor(col, alpha.f = alpha)

            circlize::circos.link(
                sv$chr1, sv$pos1,
                sv$chr2, sv$pos2,
                col = col_alpha,
                lwd = lwd
            )
        }
    }

    # Add title
    title_text <- "Genome-Wide Chromoanagenesis"
    if (sample_name != "") {
        title_text <- paste0(title_text, " - ", sample_name)
    }
    title(main = title_text, cex.main = 1.5, font.main = 2)

    # Add legend
    legend(
        "bottomleft",
        legend = c("Chromothripsis", "Chromoplexy", "Chromoanasynthesis"),
        fill = c("#E41A1C", "#377EB8", "#4DAF4A"),
        border = NA,
        bty = "n",
        cex = 0.8,
        title = "Mechanism"
    )

    # Dynamic SV type legend based on actual data
    if (!is.null(sv_data) && nrow(sv_data) > 0) {
        # Get unique SV types present in the data
        present_types <- unique(sv_data$SVtype)

        # Filter sv_colors to only include present types
        legend_types <- present_types[present_types %in% names(sv_colors)]

        if (length(legend_types) > 0) {
            legend(
                "bottomright",
                legend = legend_types,
                fill = sv_colors[legend_types],
                border = NA,
                bty = "n",
                cex = 0.8,
                title = "SV Type"
            )
        }
    }

    # Clean up
    circlize::circos.clear()
}


#' Extract chromoanagenesis regions from result
#'
#' @keywords internal
extract_chromoanagenesis_regions <- function(chromoanagenesis_result, mechanisms) {

    regions_list <- list()

    # Chromothripsis
    if ("chromothripsis" %in% mechanisms && !is.null(chromoanagenesis_result$chromothripsis)) {
        ct_class <- chromoanagenesis_result$chromothripsis$classification
        if (!is.null(ct_class) && nrow(ct_class) > 0) {
            ct_regions <- ct_class[ct_class$classification %in%
                                  c("High confidence", "Low confidence"), ]
            if (nrow(ct_regions) > 0) {
                # Try to get coordinates from detection_output if available
                if (!is.null(chromoanagenesis_result$chromothripsis$detection_output) &&
                    isS4(chromoanagenesis_result$chromothripsis$detection_output) &&
                    .hasSlot(chromoanagenesis_result$chromothripsis$detection_output, "chromSummary")) {
                    
                    ct_summary <- chromoanagenesis_result$chromothripsis$detection_output@chromSummary
                    
                    # Merge to get coordinates if not already present or to update them
                    ct_with_coords <- merge(ct_regions, ct_summary[, c("chrom", "start", "end")],
                                           by = "chrom", all.x = TRUE)
                } else {
                    # Fallback: assume start/end are in classification table (e.g. from mock or modified objects)
                    if ("start" %in% colnames(ct_regions) && "end" %in% colnames(ct_regions)) {
                        ct_with_coords <- ct_regions
                    } else {
                        warning("Could not find start/end coordinates for chromothripsis regions.")
                        ct_with_coords <- NULL
                    }
                }

                if (!is.null(ct_with_coords)) {
                    regions_list$chromothripsis <- data.frame(
                        chr = paste0("chr", ct_with_coords$chrom),
                        start = ct_with_coords$start,
                        end = ct_with_coords$end,
                        mechanism = "chromothripsis",
                        stringsAsFactors = FALSE
                    )
                }
            }
        }
    }

    # Chromoplexy - extract affected regions from chains
    if ("chromoplexy" %in% mechanisms && !is.null(chromoanagenesis_result$chromoplexy)) {
        cp_summary <- chromoanagenesis_result$chromoplexy$summary
        if (!is.null(cp_summary) && nrow(cp_summary) > 0) {
            cp_likely <- cp_summary[cp_summary$classification %in%
                                   c("Likely chromoplexy", "Possible chromoplexy"), ]
            
            # Use chain_details which contains the actual SVs with coordinates
            if (nrow(cp_likely) > 0 && !is.null(chromoanagenesis_result$chromoplexy$chain_details)) {
                cp_regions_df <- extract_chromoplexy_regions(chromoanagenesis_result$chromoplexy$chain_details)
                if (!is.null(cp_regions_df)) {
                    regions_list$chromoplexy <- cp_regions_df
                }
            }
        }
    }

    # Chromoanasynthesis
    if ("chromoanasynthesis" %in% mechanisms && !is.null(chromoanagenesis_result$chromoanasynthesis)) {
        cs_summary <- chromoanagenesis_result$chromoanasynthesis$summary
        if (!is.null(cs_summary) && nrow(cs_summary) > 0) {
            cs_regions <- cs_summary[cs_summary$classification %in%
                                    c("Likely chromoanasynthesis", "Possible chromoanasynthesis"), ]
            if (nrow(cs_regions) > 0) {
                regions_list$chromoanasynthesis <- data.frame(
                    chr = paste0("chr", cs_regions$chrom),
                    start = cs_regions$start,
                    end = cs_regions$end,
                    mechanism = "chromoanasynthesis",
                    stringsAsFactors = FALSE
                )
            }
        }
    }

    if (length(regions_list) > 0) {
        return(do.call(rbind, regions_list))
    } else {
        return(NULL)
    }
}


#' Extract chromoplexy regions from chain details
#'
#' @keywords internal
extract_chromoplexy_regions <- function(chain_details) {
    regions <- list()

    for (i in seq_along(chain_details)) {
        detail <- chain_details[[i]]
        
        # Check classification
        if (!detail$summary$classification %in% c("Likely chromoplexy", "Possible chromoplexy")) {
            next
        }
        
        # Get SVs
        chain_svs <- detail$SVs
        
        if (!is.null(chain_svs) && nrow(chain_svs) > 0) {
            # Get regions for each chromosome in this chain
            for (chr in unique(c(chain_svs$chrom1, chain_svs$chrom2))) {
                # Find all breakpoints on this chromosome
                pos1 <- chain_svs$pos1[chain_svs$chrom1 == chr]
                pos2 <- chain_svs$pos2[chain_svs$chrom2 == chr]
                all_pos <- c(pos1, pos2)

                if (length(all_pos) > 0) {
                    # Add a small buffer to make the region visible if it's just a single point
                    # But typically chromoplexy involves a region between breakpoints
                    
                    regions[[length(regions) + 1]] <- data.frame(
                        chr = paste0("chr", chr),
                        start = min(all_pos),
                        end = max(all_pos),
                        mechanism = "chromoplexy",
                        stringsAsFactors = FALSE
                    )
                }
            }
        }
    }

    if (length(regions) > 0) {
        return(do.call(rbind, regions))
    } else {
        return(NULL)
    }
}


#' Merge SV data with chromoanagenesis regions
#'
#' @keywords internal
merge_sv_with_regions <- function(sv_data, ca_regions) {

    if (is.null(ca_regions) || nrow(ca_regions) == 0) {
        sv_data$in_ca_region <- FALSE
        return(sv_data)
    }

    sv_data$in_ca_region <- FALSE

    for (i in seq_len(nrow(sv_data))) {
        chr1 <- sv_data$chr1[i]
        pos1 <- sv_data$pos1[i]
        chr2 <- sv_data$chr2[i]
        pos2 <- sv_data$pos2[i]

        # Check if either breakpoint is in a chromoanagenesis region
        in_region <- any(
            (ca_regions$chr == chr1 & !is.na(ca_regions$start) & !is.na(ca_regions$end) & ca_regions$start <= pos1 & ca_regions$end >= pos1) |
            (ca_regions$chr == chr2 & !is.na(ca_regions$start) & !is.na(ca_regions$end) & ca_regions$start <= pos2 & ca_regions$end >= pos2),
            na.rm = TRUE
        )

        sv_data$in_ca_region[i] <- in_region
    }

    return(sv_data)
}
