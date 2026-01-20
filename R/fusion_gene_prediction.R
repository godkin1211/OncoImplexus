#' Predict Fusion Genes and Disruptions
#'
#' Maps structural variant breakpoints to gene models to predict potential
#' gene fusions and disruptions.
#'
#' @param SV.sample SV data.
#' @param gene_granges GRanges object with gene models.
#' @param txdb Optional TxDb object.
#' @param max_distance Maximum distance to gene (default 0).
#' @return Data frame of predicted fusions.
#' @export
predict_fusion_genes <- function(SV.sample, 
                                gene_granges, 
                                txdb = NULL,
                                max_distance = 1000) {

    if (is.null(gene_granges)) return(NULL)

    # Convert SV.sample to data frame and capture IDs
    if (inherits(SV.sample, "SVs")) {
        # Use our S4 slots directly to be safe
        sv_df <- data.frame(
            chrom1 = SV.sample@chrom1,
            pos1 = SV.sample@pos1,
            strand1 = SV.sample@strand1,
            chrom2 = SV.sample@chrom2,
            pos2 = SV.sample@pos2,
            strand2 = SV.sample@strand2,
            SVtype = SV.sample@SVtype,
            sv_id = as.character(SV.sample@sv_id),
            stringsAsFactors = FALSE
        )
    } else {
        sv_df <- as.data.frame(SV.sample)
        if (!("sv_id" %in% colnames(sv_df))) {
            sv_df$sv_id <- rownames(sv_df)
        }
        sv_df$sv_id <- as.character(sv_df$sv_id)
    }

    if (nrow(sv_df) == 0) return(NULL)
    
    # Store IDs for mapping later - FORCE CHARACTER
    actual_sv_ids <- as.character(sv_df$sv_id)

    cat("\n======================================================================\n")
    cat("     FUSION GENE & DISRUPTION PREDICTION\n")
    cat("======================================================================\n\n")

    cat(sprintf("Step 1: Mapping breakpoints to gene models (Buffer: %d bp)...\n", max_distance))
    
    target_has_chr <- any(grepl("^chr", as.character(GenomicRanges::seqnames(gene_granges))))
    norm_chroms1 <- if (target_has_chr) ifelse(grepl("^chr", sv_df$chrom1), sv_df$chrom1, paste0("chr", sv_df$chrom1)) else gsub("^chr", "", sv_df$chrom1)
    norm_chroms2 <- if (target_has_chr) ifelse(grepl("^chr", sv_df$chrom2), sv_df$chrom2, paste0("chr", sv_df$chrom2)) else gsub("^chr", "", sv_df$chrom2)

    # Clean strands
    s1_clean <- ifelse(sv_df$strand1 %in% c("+", "-"), sv_df$strand1, "*")
    s2_clean <- ifelse(sv_df$strand2 %in% c("+", "-"), sv_df$strand2, "*")

    bp1_gr <- GenomicRanges::GRanges(seqnames = norm_chroms1, ranges = IRanges::IRanges(start = as.numeric(sv_df$pos1), width = 1), strand = s1_clean)
    bp2_gr <- GenomicRanges::GRanges(seqnames = norm_chroms2, ranges = IRanges::IRanges(start = as.numeric(sv_df$pos2), width = 1), strand = s2_clean)

    symbol_col <- if ("symbol" %in% names(GenomicRanges::mcols(gene_granges))) "symbol" else "gene_name"
    hits1 <- GenomicRanges::findOverlaps(bp1_gr, gene_granges, maxgap = max_distance)
    hits2 <- GenomicRanges::findOverlaps(bp2_gr, gene_granges, maxgap = max_distance)

    # Get actual SV IDs for mapping - ALREADY PREPARED ABOVE
    # actual_sv_ids <- if ("sv_id" %in% colnames(sv_df)) as.character(sv_df$sv_id) else as.character(1:nrow(sv_df))

    map1 <- data.frame(
        sv_id = as.character(actual_sv_ids[S4Vectors::queryHits(hits1)]), 
        subject_idx = S4Vectors::subjectHits(hits1),
        gene1 = as.character(GenomicRanges::mcols(gene_granges)[[symbol_col]][S4Vectors::subjectHits(hits1)]), 
        gene1_strand = as.character(GenomicRanges::strand(gene_granges)[S4Vectors::subjectHits(hits1)]),
        gene1_start = GenomicRanges::start(gene_granges)[S4Vectors::subjectHits(hits1)],
        gene1_end = GenomicRanges::end(gene_granges)[S4Vectors::subjectHits(hits1)],
        type = if ("type" %in% names(mcols(gene_granges))) as.character(mcols(gene_granges)$type[S4Vectors::subjectHits(hits1)]) else "gene",
        stringsAsFactors = FALSE
    )
    
    # Side 2
    map2 <- data.frame(
        sv_id = as.character(actual_sv_ids[S4Vectors::queryHits(hits2)]), 
        subject_idx = S4Vectors::subjectHits(hits2),
        gene2 = as.character(GenomicRanges::mcols(gene_granges)[[symbol_col]][S4Vectors::subjectHits(hits2)]), 
        gene2_strand = as.character(GenomicRanges::strand(gene_granges)[S4Vectors::subjectHits(hits2)]), 
        gene2_start = GenomicRanges::start(gene_granges)[S4Vectors::subjectHits(hits2)],
        gene2_end = GenomicRanges::end(gene_granges)[S4Vectors::subjectHits(hits2)],
        type = if ("type" %in% names(mcols(gene_granges))) as.character(mcols(gene_granges)$type[S4Vectors::subjectHits(hits2)]) else "gene",
        stringsAsFactors = FALSE
    )

    # Helper function to get context (Intron/Exon and % position)
    get_context <- function(sv_pos, g_start, g_end, g_strand, sv_id, map_df) {
        rel_pos <- if (g_strand == "-") (g_end - sv_pos) / (g_end - g_start) else (sv_pos - g_start) / (g_end - g_start)
        rel_pos <- pmax(0, pmin(1, rel_pos))
        
        # Check if any hit for this SV_ID was an exon
        hits_for_sv <- map_df[map_df$sv_id == sv_id, ]
        is_exon <- any(hits_for_sv$type == "exon")
        
        context <- if (is_exon) "Exon" else "Intron/Unknown"
        return(list(rel_pos = rel_pos, context = context))
    }

    # 1. Detect Inter-genic Fusions (Both breakpoints hit different genes)
    # Filter map objects to genes only for merging
    map1_g <- map1[map1$type == "gene" | map1$type == "mRNA" | is.na(map1$type), ]
    map2_g <- map2[map2$type == "gene" | map2$type == "mRNA" | is.na(map2$type), ]
    
    fusions_both <- merge(map1_g, map2_g, by = "sv_id")
    fusions_both <- fusions_both[fusions_both$gene1 != fusions_both$gene2, ]
    
    # 2. Detect Disruptions
    sv_ids_hits1 <- unique(map1_g$sv_id)
    sv_ids_hits2 <- unique(map2_g$sv_id)
    
    only_side1 <- map1_g[!(map1_g$sv_id %in% sv_ids_hits2), ]
    only_side2 <- map2_g[!(map2_g$sv_id %in% sv_ids_hits1), ]

    fusions_list <- list()
    
    # Process Fusions
    if (nrow(fusions_both) > 0) {
        for (i in 1:nrow(fusions_both)) {
            row <- fusions_both[i, ]
            # Correctly find the SV by ID
            sv_match_idx <- which(sv_df$sv_id == row$sv_id)[1]
            sv <- sv_df[sv_match_idx, ]
            
            ctx1 <- get_context(sv$pos1, row$gene1_start, row$gene1_end, row$gene1_strand, row$sv_id, map1)
            ctx2 <- get_context(sv$pos2, row$gene2_start, row$gene2_end, row$gene2_strand, row$sv_id, map2)
            
            g1 <- row$gene1; g2 <- row$gene2
            g1s <- row$gene1_strand; g2s <- row$gene2_strand
            s1 <- as.character(sv$strand1); s2 <- as.character(sv$strand2)
            
            ok <- FALSE; p5 <- g1; p3 <- g2
            if (!is.na(g1s) && !is.na(g2s) && s1 != "*" && s2 != "*") {
                if ((g1s == "+" && s1 == "+") || (g1s == "-" && s1 == "-")) {
                    if ((g2s == "+" && s2 == "-") || (g2s == "-" && s2 == "+")) {
                        p5 <- g1; p3 <- g2; ok <- TRUE
                    }
                }
                if (!ok) {
                    if ((g2s == "+" && s2 == "+") || (g2s == "-" && s2 == "-")) {
                        if ((g1s == "+" && s1 == "-") || (g1s == "-" && s1 == "+")) {
                            p5 <- g2; p3 <- g1; ok <- TRUE
                        }
                    }
                }
            }
            
            # Determine detailed status
            final_status <- if (ok) "Canonical" else "Non-canonical"
            if (ctx1$context == "Exon" || ctx2$context == "Exon") {
                final_status <- paste0(final_status, " (Exon disruption)")
            } else if (ctx1$context == "Exon" && ctx2$context == "Exon") {
                final_status <- paste0(final_status, " (Exon-Exon)")
            } else {
                final_status <- paste0(final_status, " (Intronic)")
            }

            fusions_list[[length(fusions_list) + 1]] <- data.frame(
                fusion_name = paste0(p5, "::", p3),
                gene_5p = p5, gene_3p = p3,
                type = "Fusion",
                status = final_status,
                orientation_compatible = ok,
                sv_type = sv$SVtype,
                chrom1 = sv$chrom1, pos1 = sv$pos1,
                chrom2 = sv$chrom2, pos2 = sv$pos2,
                gene1_pos_pct = round(ctx1$rel_pos * 100, 1),
                gene2_pos_pct = round(ctx2$rel_pos * 100, 1),
                sv_id = row$sv_id, stringsAsFactors = FALSE
            )
        }
    }

    # Process Disruptions Side 1
    if (nrow(only_side1) > 0) {
        for (i in 1:nrow(only_side1)) {
            row <- only_side1[i, ]
            sv_match_idx <- which(sv_df$sv_id == row$sv_id)[1]
            sv <- sv_df[sv_match_idx, ]
            ctx1 <- get_context(sv$pos1, row$gene1_start, row$gene1_end, row$gene1_strand, row$sv_id, map1)
            
            fusions_list[[length(fusions_list) + 1]] <- data.frame(
                fusion_name = paste0(row$gene1, "::(intergenic)"),
                gene_5p = row$gene1, gene_3p = "intergenic",
                type = "Disruption",
                status = if(ctx1$context == "Exon") "Exon disruption" else "Intron disruption",
                orientation_compatible = NA,
                sv_type = sv$SVtype,
                chrom1 = sv$chrom1, pos1 = sv$pos1,
                chrom2 = sv$chrom2, pos2 = sv$pos2,
                gene1_pos_pct = round(ctx1$rel_pos * 100, 1),
                gene2_pos_pct = NA,
                sv_id = row$sv_id, stringsAsFactors = FALSE
            )
        }
    }

    # Process Disruptions Side 2
    if (nrow(only_side2) > 0) {
        for (i in 1:nrow(only_side2)) {
            row <- only_side2[i, ]
            sv_match_idx <- which(sv_df$sv_id == row$sv_id)[1]
            sv <- sv_df[sv_match_idx, ]
            ctx2 <- get_context(sv$pos2, row$gene2_start, row$gene2_end, row$gene2_strand, row$sv_id, map2)
            
            fusions_list[[length(fusions_list) + 1]] <- data.frame(
                fusion_name = paste0("(intergenic)::", row$gene2),
                gene_5p = "intergenic", gene_3p = row$gene2,
                type = "Disruption",
                status = if(ctx2$context == "Exon") "Exon disruption" else "Intron disruption",
                orientation_compatible = NA,
                sv_type = sv$SVtype,
                chrom1 = sv$chrom1, pos1 = sv$pos1,
                chrom2 = sv$chrom2, pos2 = sv$pos2,
                gene1_pos_pct = NA,
                gene2_pos_pct = round(ctx2$rel_pos * 100, 1),
                sv_id = row$sv_id, stringsAsFactors = FALSE
            )
        }
    }

    if (length(fusions_list) == 0) return(NULL)
    fusions_df <- do.call(rbind, fusions_list)
    
    driver_genes <- get_default_drivers()
    fusions_df$is_driver_involved <- fusions_df$gene_5p %in% driver_genes | fusions_df$gene_3p %in% driver_genes
    
    cat(sprintf("Detected %d potential events (Fusions/Disruptions).\n", nrow(fusions_df)))
    return(fusions_df)
}
