#' @include chromoanagenesis_integrated.R
NULL

load_cancer_genes <- function() {
    # 1. Check if the data is already available in the global option/cache to avoid repeated IO
    if (exists(".cancerGeneCache", envir = .GlobalEnv)) {
        return(get(".cancerGeneCache", envir = .GlobalEnv))
    }
    
    # 2. Try multiple possible paths for the TSV file
    possible_paths <- c(
        system.file("extdata", "cancerGeneList.tsv", package = "OncoImplexus"),
        system.file("extdata", "cancerGeneList.tsv", package = "OncoImplexus", mustWork = FALSE),
        "inst/extdata/cancerGeneList.tsv",
        "cancerGeneList.tsv"
    )
    
    # Filter out empty paths
    possible_paths <- possible_paths[possible_paths != ""]
    
    df <- NULL
    for (p in possible_paths) {
        if (file.exists(p)) {
            tryCatch({
                df <- read.table(p, sep = "\t", header = TRUE, 
                                stringsAsFactors = FALSE, check.names = FALSE, 
                                quote = "", comment.char = "")
                colnames(df) <- gsub("[[:space:].]+", "_", colnames(df))
                if ("Hugo_Symbol" %in% colnames(df)) {
                    message(sprintf("Successfully loaded %d genes from %s", nrow(df), p))
                    # Cache it
                    assign(".cancerGeneCache", df, envir = .GlobalEnv)
                    return(df)
                }
            }, error = function(e) {
                return(NULL)
            })
        }
    }
    
    # 3. Fallback to binary
    data_path <- system.file("data", "cancerGeneList.rda", package = "OncoImplexus")
    if (file.exists(data_path)) {
        env <- new.env()
        load(data_path, envir = env)
        if (exists("cancerGeneList", envir = env)) {
            df <- get("cancerGeneList", envir = env)
            colnames(df) <- gsub("[[:space:].]+", "_", colnames(df))
            assign(".cancerGeneCache", df, envir = .GlobalEnv)
            return(df)
        }
    }

    return(NULL)
}
#' Get a default list of cancer driver genes
#'
#' Returns a curated list of high-confidence cancer driver genes.
#' Tries to load from the comprehensive OncoKB list included in the package.
#' Falls back to a smaller hardcoded list if the file is unavailable.
#'
#' @return A character vector of gene symbols
#' @export
get_default_drivers <- function() {
    # Try to load from file first
    df <- load_cancer_genes()
    
    if (!is.null(df) && "Hugo_Symbol" %in% colnames(df)) {
        return(unique(df$Hugo_Symbol))
    }
    
    # Fallback: A curated subset of COSMIC Cancer Gene Census Tier 1 & OncoKB
    warning("Using fallback hardcoded driver list (cancerGeneList.tsv not found or invalid format).")
    c("TP53", "MYC", "PTEN", "BRCA1", "BRCA2", "EGFR", "KRAS", "NRAS", "HRAS",
      "BRAF", "PIK3CA", "APC", "VHL", "RB1", "CDKN2A", "CCND1", "MDM2", "ERBB2",
      "ATM", "ATR", "MLH1", "MSH2", "MSH6", "PMS2", "ARID1A", "SMARCA4", "PBRM1",
      "KMT2A", "KMT2D", "CREBBP", "EP300", "NOTCH1", "FBXW7", "SF3B1", "TERT",
      "ALK", "ROS1", "RET", "MET", "NTRK1", "NTRK2", "NTRK3", "FGFR1", "FGFR2",
      "FGFR3", "KIT", "PDGFRA", "FLT3", "IDH1", "IDH2", "WT1", "NPM1", "RUNX1",
      "DNMT3A", "TET2", "JAK2", "BCL2", "BCL6", "MYCN", "CDK4", "CDK6", "CTNNB1",
      "NF1", "NF2", "SMAD4", "GATA3", "FOXA1", "MED12", "SETD2", "BAP1", "CDH1")
}

#' Annotate chromoanagenesis events with gene information
#'
#' Maps chromoanagenesis events (Chromothripsis, Chromoplexy, Chromoanasynthesis)
#' to genes, identifying potentially impacted cancer driver genes.
#' 
#' Incorporates logic to check for "Mechanism Consistency":
#' - Oncogenes are expected to be Amplified (Chromoanasynthesis)
#' - Tumor Suppressors (TSGs) are expected to be Disrupted/Lost (Chromothripsis/Chromoplexy)
#'
#' @param result A chromoanagenesis result object from detect_chromoanagenesis()
#' @param gene_granges A GRanges object containing gene coordinates. 
#'   Must have a 'symbol' or 'gene_name' metadata column.
#' @param driver_genes Optional character vector of cancer driver gene symbols. 
#'   If NULL, uses get_default_drivers() (OncoKB list).
#' @return A list containing:
#'   \item{driver_hits}{Detailed table of impacted driver genes with mechanism consistency analysis}
#'   \item{all_impacted_genes}{Table of all impacted genes}
#' @export
annotate_chromoanagenesis <- function(result, gene_granges, driver_genes = NULL) {
    
    if (!inherits(result, "chromoanagenesis")) {
        stop("Input must be a chromoanagenesis result object")
    }
    
    if (is.null(gene_granges) || !inherits(gene_granges, "GRanges")) {
        stop("gene_granges must be a valid GRanges object containing gene annotations")
    }
    
    # Load cancer gene metadata for type (Oncogene vs TSG)
    cancer_genes_df <- load_cancer_genes()
    
    if (is.null(driver_genes)) {
        if (!is.null(cancer_genes_df)) {
            driver_genes <- unique(cancer_genes_df$Hugo_Symbol)
        } else {
            driver_genes <- get_default_drivers()
        }
    }
    
    # Ensure gene symbol column exists
    gene_cols <- colnames(GenomicRanges::mcols(gene_granges))
    symbol_col <- NULL
    if ("symbol" %in% gene_cols) symbol_col <- "symbol"
    else if ("gene_name" %in% gene_cols) symbol_col <- "gene_name"
    else if ("gene_id" %in% gene_cols) symbol_col <- "gene_id"
    else stop("gene_granges must have a 'symbol', 'gene_name', or 'gene_id' column")
    
    # Helper to normalize chromosome names to match gene_granges style
    # This prevents mismatches like "1" vs "chr1"
    target_style_has_chr <- any(grepl("^chr", GenomeInfoDb::seqlevels(gene_granges)))
    
    normalize_chrom <- function(chroms) {
        if (target_style_has_chr) {
            # Target has chr, ensure input has chr
            return(ifelse(grepl("^chr", chroms), chroms, paste0("chr", chroms)))
        } else {
            # Target has NO chr, remove chr from input
            return(gsub("^chr", "", chroms))
        }
    }
    
    # Helper: Get Gene Type (Oncogene/TSG)
    get_gene_type <- function(gene) {
        if (is.null(cancer_genes_df)) return("Unknown")
        idx <- which(cancer_genes_df$Hugo_Symbol == gene)
        if (length(idx) > 0) {
            return(cancer_genes_df$Gene_Type[idx[1]]) # Take first match
        }
        return("Unknown")
    }
    
    # Helper: Check Mechanism Consistency
    check_consistency <- function(gene_type, impact_type) {
        if (gene_type == "Unknown") return("Unknown")
        
        is_oncogene <- grepl("ONCOGENE", gene_type)
        is_tsg <- grepl("TSG", gene_type)
        
        if (impact_type == "Amplification_Target") {
            # Amplification fits Oncogenes
            if (is_oncogene) return("Consistent (Oncogene Amp)")
            if (is_tsg) return("Inconsistent (TSG Amp)")
        } else if (impact_type %in% c("Breakpoint_Disruption", "Region_Overlap")) {
            # Disruption fits TSGs (Loss of function)
            if (is_tsg) return("Consistent (TSG Loss)")
            if (is_oncogene) return("Potential Fusion/Disruption") # Could be fusion (activating) or loss (inconsistent)
        }
        
        return("Indeterminate")
    }
    
    # Helper to format output
    format_hit <- function(gene, mechanism, type, region) {
        g_type <- get_gene_type(gene)
        consistency <- check_consistency(g_type, type)
        
        data.frame(
            Gene = gene,
            Gene_Type = g_type,
            Mechanism = mechanism,
            Impact_Type = type,
            Mechanism_Consistency = consistency,
            Region = region,
            Is_Driver = gene %in% driver_genes,
            stringsAsFactors = FALSE
        )
    }
    
    hits_list <- list()
    
    # 1. Annotate Chromothripsis
    chromoth_hits <- data.frame()
    if (!is.null(result$chromothripsis) && result$chromothripsis$n_high_confidence > 0) {
        cls <- result$chromothripsis$classification
        high_conf <- cls[cls$classification == "High confidence", ]
        
        if (nrow(high_conf) > 0) {
            for (i in 1:nrow(high_conf)) {
                chrom <- normalize_chrom(high_conf$chrom[i])
                
                # Robust extraction of coordinates
                start_val <- if ("start" %in% colnames(high_conf)) high_conf$start[i] else NA
                end_val <- if ("end" %in% colnames(high_conf)) high_conf$end[i] else NA
                
                # Handle zero-length or NA
                start <- if (length(start_val) == 0 || is.na(start_val)) 1 else start_val
                end <- if (length(end_val) == 0 || is.na(end_val)) 500e6 else end_val
                
                region_gr <- GenomicRanges::GRanges(seqnames = chrom, ranges = IRanges::IRanges(start = as.numeric(start), end = as.numeric(end)))
                ov <- GenomicRanges::findOverlaps(region_gr, gene_granges)
                impacted_indices <- S4Vectors::subjectHits(ov)
                
                if (length(impacted_indices) > 0) {
                    genes <- GenomicRanges::mcols(gene_granges)[[symbol_col]][impacted_indices]
                    for (g in unique(genes)) {
                        hits_list[[length(hits_list) + 1]] <- format_hit(
                            g, "Chromothripsis", "Region_Overlap", 
                            sprintf("%s:%d-%d", chrom, start, end)
                        )
                    }
                }
            }
        }
    }
    
    # 2. Annotate Chromoplexy
    if (!is.null(result$chromoplexy) && result$chromoplexy$likely_chromoplexy > 0) {
        chain_details <- result$chromoplexy$chain_details
        
        for (i in seq_along(chain_details)) {
            detail <- chain_details[[i]]
            
            if (detail$summary$classification == "Likely chromoplexy") {
                svs <- detail$SVs
                if (is.null(svs) || nrow(svs) == 0) next
                
                # Check both breakpoints for each SV
                for (j in 1:nrow(svs)) {
                    # Breakpoint 1 with 1kb buffer
                    bp1_gr <- GenomicRanges::GRanges(seqnames = normalize_chrom(svs$chrom1[j]), 
                                      ranges = IRanges::IRanges(start = svs$pos1[j] - 500, width = 1001))
                    
                    ov1 <- GenomicRanges::findOverlaps(bp1_gr, gene_granges)
                    
                    if (length(ov1) > 0) {
                        genes <- unique(GenomicRanges::mcols(gene_granges)[[symbol_col]][S4Vectors::subjectHits(ov1)])
                        
                        for (g in genes) {
                            hits_list[[length(hits_list) + 1]] <- format_hit(
                                g, "Chromoplexy", "Breakpoint_Disruption",
                                sprintf("%s:%d", svs$chrom1[j], svs$pos1[j])
                            )
                        }
                    }
                    
                    # Breakpoint 2 with 1kb buffer
                    bp2_gr <- GenomicRanges::GRanges(seqnames = normalize_chrom(svs$chrom2[j]), 
                                      ranges = IRanges::IRanges(start = svs$pos2[j] - 500, width = 1001))
                    
                    ov2 <- GenomicRanges::findOverlaps(bp2_gr, gene_granges)
                    if (length(ov2) > 0) {
                        genes <- unique(GenomicRanges::mcols(gene_granges)[[symbol_col]][S4Vectors::subjectHits(ov2)])
                        for (g in genes) {
                            hits_list[[length(hits_list) + 1]] <- format_hit(
                                g, "Chromoplexy", "Breakpoint_Disruption",
                                sprintf("%s:%d", svs$chrom2[j], svs$pos2[j])
                            )
                        }
                    }
                }
            }
        }
    }
    
    # 3. Annotate Chromoanasynthesis
    if (!is.null(result$chromoanasynthesis) && result$chromoanasynthesis$likely_chromoanasynthesis > 0) {
        summary <- result$chromoanasynthesis$summary
        likely <- summary[summary$classification == "Likely chromoanasynthesis", ]
        
        if (nrow(likely) > 0) {
            for (i in 1:nrow(likely)) {
                chrom <- normalize_chrom(likely$chrom[i])
                start <- likely$start[i]
                end <- likely$end[i]
                
                region_gr <- GenomicRanges::GRanges(seqnames = chrom, ranges = IRanges::IRanges(start = start, end = end))
                ov <- GenomicRanges::findOverlaps(region_gr, gene_granges)
                genes <- unique(GenomicRanges::mcols(gene_granges)[[symbol_col]][S4Vectors::subjectHits(ov)])
                
                for (g in genes) {
                    hits_list[[length(hits_list) + 1]] <- format_hit(
                        g, "Chromoanasynthesis", "Amplification_Target",
                        sprintf("%s:%d-%d", chrom, start, end)
                    )
                }
            }
        }
    }
    
    # Combine results
    all_hits <- do.call(rbind, hits_list)
    if (is.null(all_hits) || nrow(all_hits) == 0) { # Check for empty all_hits as well
        return(list(
            driver_hits = data.frame(
                Gene = character(), Gene_Type = character(), 
                Mechanism = character(), Impact_Type = character(),
                Mechanism_Consistency = character(), Region = character(),
                Is_Driver = logical()
            ),
            all_impacted_genes = data.frame()
        ))
    }
    
    # Extract drivers
    driver_hits <- all_hits[all_hits$Is_Driver, ]
    
    return(list(
        driver_hits = driver_hits,
        all_impacted_genes = all_hits
    ))
}