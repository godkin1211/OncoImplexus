' Analyze breakpoint sequence features
#'
#' Extracts and analyzes sequence features at structural variant breakpoints,
#' including microhomology, insertions, and inferred DNA repair mechanisms.
#'
#' @param SV.sample An instance of class SVs or data frame with SV data
#' @param genome Reference genome object (BSgenome) or path to FASTA file
#' @param flank_size Size of flanking sequence to extract (default: 50 bp)
#' @param min_microhomology Minimum microhomology length to report (default: 2 bp)
#' @param max_microhomology Maximum microhomology length to search (default: 25 bp)
#' @return A list containing breakpoint sequence analysis results
#' @details
#' This function analyzes breakpoint junction sequences to identify:
#'
#' 1. **Microhomology**: Short homologous sequences (1-25 bp) at breakpoint junctions,
#'    indicative of microhomology-mediated end joining (MMEJ) or alt-EJ
#'
#' 2. **Insertions**: Non-templated or templated insertions at junctions
#'
#' 3. **DNA Repair Mechanisms**:
#'    - NHEJ (Non-Homologous End Joining): No microhomology, blunt ends
#'    - MMEJ (Microhomology-Mediated End Joining): 2-25 bp microhomology
#'    - HR (Homologous Recombination): Long homology (>25 bp)
#'    - FoSTeS/MMBIR: Serial replication with microhomology
#'
#' The analysis requires a reference genome to extract sequences around
#' breakpoint positions. Supports BSgenome objects or FASTA files.
#'
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg19)
#'
#' # Analyze breakpoint sequences
#' bp_analysis <- analyze_breakpoint_sequences(
#'   SV_data,
#'   genome = BSgenome.Hsapiens.UCSC.hg19
#' )
#'
#' # View results
#' print(bp_analysis)
#' summary(bp_analysis)
#'
#' # Plot repair mechanism distribution
#' plot_repair_mechanisms(bp_analysis)
#' }
#'
#' @export
analyze_breakpoint_sequences <- function(SV.sample,
                                        genome,
                                        flank_size = 50,
                                        min_microhomology = 2,
                                        max_microhomology = 25) {

    # Input validation
    if (missing(genome)) {
        stop("Reference genome is required for sequence analysis.\n",
             "Provide a BSgenome object or path to FASTA file.")
    }

    # Convert SV.sample to data frame if needed
    if (inherits(SV.sample, "SVs")) {
        sv_df <- data.frame(
            chrom1 = SV.sample@chrom1,
            pos1 = SV.sample@pos1,
            strand1 = SV.sample@strand1,
            chrom2 = SV.sample@chrom2,
            pos2 = SV.sample@pos2,
            strand2 = SV.sample@strand2,
            SVtype = SV.sample@SVtype,
            stringsAsFactors = FALSE
        )
        # Try to include sv_id if it exists in the object (PRIORITIZE STABLE ID)
        if (length(SV.sample@sv_id) > 0) {
            sv_df$sv_id <- as.character(SV.sample@sv_id)
        } else {
            # Fallback to row names or sequential IDs ONLY IF missing
            sv_df$sv_id <- rownames(sv_df)
        }
    } else {
        sv_df <- as.data.frame(SV.sample)
        if (!("sv_id" %in% colnames(sv_df))) {
            # If no ID column, use row names which are more stable than 1:nrow()
            sv_df$sv_id <- rownames(sv_df)
        } else {
            sv_df$sv_id <- as.character(sv_df$sv_id)
        }
    }

    if (nrow(sv_df) == 0) {
        stop("No structural variants provided")
    }

    cat("\n")
    cat(rep("=", 70), "\n", sep = "")
    cat("     BREAKPOINT SEQUENCE ANALYSIS\n")
    cat(rep("=", 70), "\n\n", sep = "")

    cat(sprintf("Analyzing %d structural variants...\n", nrow(sv_df)))
    cat(sprintf("Flank size: %d bp\n", flank_size))
    cat(sprintf("Microhomology range: %d-%d bp\n\n", min_microhomology, max_microhomology))

    results <- list()

    # 1. Extract breakpoint sequences
    cat("Step 1: Extracting breakpoint sequences...\n")
    bp_sequences <- extract_breakpoint_sequences(
        sv_df,
        genome,
        flank_size
    )
    results$sequences <- bp_sequences

    # 2. Detect microhomology
    cat("Step 2: Detecting microhomology...\n")
    microhomology <- detect_microhomology(
        bp_sequences,
        sv_df,
        min_length = min_microhomology,
        max_length = max_microhomology
    )
    results$microhomology <- microhomology

    # 3. Detect insertions
    cat("Step 3: Detecting insertions...\n")
    insertions <- detect_insertions(bp_sequences, sv_df)
    results$insertions <- insertions

    # 4. Classify repair mechanisms
    cat("Step 4: Classifying DNA repair mechanisms...\n")
    repair_classification <- classify_repair_mechanisms(
        microhomology,
        insertions,
        sv_df
    )
    results$repair_mechanisms <- repair_classification

    # 5. Summary statistics
    cat("Step 5: Generating summary statistics...\n")
    summary_stats <- create_sequence_summary(
        microhomology,
        insertions,
        repair_classification
    )
    results$summary <- summary_stats

    cat("\n")
    cat(rep("=", 70), "\n", sep = "")
    cat("ANALYSIS COMPLETE\n")
    cat(rep("=", 70), "\n\n", sep = "")

    class(results) <- c("breakpoint_sequences", "list")
    return(results)
}


#' Extract sequences around breakpoints
#'
#' @param sv_df Data frame with SV information
#' @param genome Reference genome
#' @param flank_size Flanking sequence size
#' @return Data frame with extracted sequences
#' @keywords internal
extract_breakpoint_sequences <- function(sv_df, genome, flank_size) {

    # Check if genome is BSgenome object
    is_bsgenome <- inherits(genome, "BSgenome")
    if (!is_bsgenome) {
        warning("Reference genome must be a BSgenome object for sequence analysis.")
        return(NULL)
    }

    # Normalize chromosome names in genome
    gen_chroms <- GenomeInfoDb::seqnames(genome)
    has_chr_prefix <- any(grepl("^chr", gen_chroms))

    sequences <- list()

    for (i in 1:nrow(sv_df)) {
        sv <- sv_df[i, ]
        
        # Use the ID column we prepared in the parent function
        current_id <- as.character(sv$sv_id)

        # Normalize sample chromosome names to match genome
        c1 <- as.character(sv$chrom1)
        c2 <- as.character(sv$chrom2)
        
        if (has_chr_prefix) {
            if (!grepl("^chr", c1)) c1 <- paste0("chr", c1)
            if (!grepl("^chr", c2)) c2 <- paste0("chr", c2)
        } else {
            c1 <- gsub("^chr", "", c1)
            c2 <- gsub("^chr", "", c2)
        }

        # Validate chromosomes exist
        if (!(c1 %in% gen_chroms) || !(c2 %in% gen_chroms)) {
            sequences[[i]] <- list(sv_id = current_id, has_sequence = FALSE)
            next
        }

        # Side 1 (P1, S1):
        if (sv$strand1 == "+") {
            s1_kept <- extract_sequence_region(genome, c1, sv$pos1 - flank_size + 1, sv$pos1, TRUE)
            s1_disc <- extract_sequence_region(genome, c1, sv$pos1 + 1, sv$pos1 + flank_size, TRUE)
        } else {
            s1_kept <- extract_sequence_region(genome, c1, sv$pos1, sv$pos1 + flank_size - 1, TRUE)
            s1_kept <- reverse_complement(s1_kept)
            s1_disc <- extract_sequence_region(genome, c1, sv$pos1 - flank_size, sv$pos1 - 1, TRUE)
            s1_disc <- reverse_complement(s1_disc)
        }

        # Side 2 (P2, S2):
        if (sv$strand2 == "-") {
            s2_kept <- extract_sequence_region(genome, c2, sv$pos2, sv$pos2 + flank_size - 1, TRUE)
            s2_disc <- extract_sequence_region(genome, c2, sv$pos2 - flank_size, sv$pos2 - 1, TRUE)
        } else {
            s2_kept <- extract_sequence_region(genome, c2, sv$pos2 - flank_size + 1, sv$pos2, TRUE)
            s2_kept <- reverse_complement(s2_kept)
            s2_disc <- extract_sequence_region(genome, c2, sv$pos2 + 1, sv$pos2 + flank_size, TRUE)
            s2_disc <- reverse_complement(s2_disc)
        }

        sequences[[i]] <- list(
            sv_id = current_id,
            chrom1 = as.character(sv$chrom1),
            pos1 = as.numeric(sv$pos1),
            chrom2 = as.character(sv$chrom2),
            pos2 = as.numeric(sv$pos2),
            has_sequence = TRUE,
            s1_kept = s1_kept,
            s1_disc = s1_disc,
            s2_kept = s2_kept,
            s2_disc = s2_disc
        )
    }

    # Convert to data frame
    seq_df <- do.call(rbind, lapply(sequences, function(x) {
        if (!x$has_sequence) {
            return(data.frame(sv_id = as.character(x$sv_id), has_sequence = FALSE, 
                             chrom1=NA, pos1=NA, chrom2=NA, pos2=NA,
                             s1_kept="", s1_disc="", s2_kept="", s2_disc="", stringsAsFactors=FALSE))
        }
        data.frame(
            sv_id = as.character(x$sv_id),
            chrom1 = as.character(x$chrom1),
            pos1 = as.numeric(x$pos1),
            chrom2 = as.character(x$chrom2),
            pos2 = as.numeric(x$pos2),
            has_sequence = TRUE,
            s1_kept = x$s1_kept,
            s1_disc = x$s1_disc,
            s2_kept = x$s2_kept,
            s2_disc = x$s2_disc,
            stringsAsFactors = FALSE
        )
    }))

    return(seq_df)
}




#' Extract sequence from genome region
#'
#' @param genome Genome object or file path
#' @param chrom Chromosome name
#' @param start Start position
#' @param end End position
#' @param is_bsgenome Whether genome is BSgenome object
#' @return Sequence string
#' @keywords internal
extract_sequence_region <- function(genome, chrom, start, end, is_bsgenome) {

    if (is_bsgenome) {
        # Use BSgenome
        tryCatch({
            seq <- as.character(genome[[chrom]][start:end])
            return(seq)
        }, error = function(e) {
            return(NA)
        })
    } else {
        return(NA)
    }
}


#' Reverse complement of DNA sequence
#'
#' @param seq DNA sequence string
#' @return Reverse complement
#' @keywords internal
reverse_complement <- function(seq) {
    if (is.na(seq)) return(NA)

    # Complement mapping
    comp_map <- c("A" = "T", "T" = "A", "G" = "C", "C" = "G",
                 "a" = "t", "t" = "a", "g" = "c", "c" = "g",
                 "N" = "N", "n" = "n")

    # Split sequence
    bases <- strsplit(seq, "")[[1]]

    # Complement
    comp_bases <- sapply(bases, function(b) {
        if (b %in% names(comp_map)) comp_map[b] else "N"
    })

    # Reverse
    rev_comp <- paste(rev(comp_bases), collapse = "")

    return(rev_comp)
}


#' Detect microhomology at breakpoint junctions
#'
#' @param bp_sequences Data frame with breakpoint sequences.
#' @param sv_df Data frame with original SV information.
#' @param min_length Minimum microhomology length to report.
#' @param max_length Maximum microhomology length to search.
#' @return Data frame with microhomology information.
#' @keywords internal
detect_microhomology <- function(bp_sequences, sv_df, min_length = 2, max_length = 25) {

    microhomology <- list()

    for (i in 1:nrow(bp_sequences)) {
        bp <- bp_sequences[i, ]
        sv <- sv_df[i, ]

        if (!bp$has_sequence) {
            microhomology[[i]] <- list(sv_id = bp$sv_id, has_microhomology = NA, microhomology_length = NA, microhomology_seq = NA)
            next
        }

        # Use the fixed logic based on strands
        mh_result <- find_microhomology(
            bp$s1_kept, bp$s2_kept, bp$s1_disc, bp$s2_disc,
            strand1 = sv$strand1, strand2 = sv$strand2,
            min_length = min_length, max_length = max_length
        )

        microhomology[[i]] <- list(
            sv_id = bp$sv_id,
            has_microhomology = mh_result$found,
            microhomology_length = mh_result$length,
            microhomology_seq = mh_result$sequence,
            position = "junction"
        )
    }

    # Convert to data frame
    mh_df <- do.call(rbind, lapply(microhomology, function(x) {
        data.frame(
            sv_id = as.character(x$sv_id),
            has_microhomology = as.logical(x$has_microhomology),
            microhomology_length = as.numeric(x$microhomology_length),
            microhomology_seq = as.character(x$microhomology_seq),
            stringsAsFactors = FALSE
        )
    }))

    return(mh_df)
}

#' Internal: Correctly compare microhomology
#' @keywords internal
find_microhomology <- function(s1_kept, s2_kept, s1_disc, s2_disc, 
                               strand1, strand2, min_length, max_length) {
    
    best_len <- 0
    best_seq <- NA
    
    n1 <- nchar(s1_disc)
    n2 <- nchar(s2_disc)
    limit <- min(n1, n2, max_length)
    
    if (limit < 1) return(list(found=FALSE, length=0, sequence=NA))

    # Single, unified loop for all strand combinations
    for (k in 1:limit) {
        f1 <- toupper(substr(s1_disc, 1, k))
        f2 <- toupper(substr(s2_disc, n2 - k + 1, n2))
        
        if (f1 == f2) {
            best_len <- k
            best_seq <- f1
        } else {
            break # Must be continuous from the junction
        }
    }
    
    return(list(
        found = (best_len >= min_length),
        length = best_len,
        sequence = best_seq
    ))
}


#' Detect insertions at breakpoint junctions
#'
#' @param bp_sequences Data frame with breakpoint sequences
#' @param sv_df Data frame with original SV information
#' @return Data frame with insertion information
#' @keywords internal
detect_insertions <- function(bp_sequences, sv_df) {

    insertions <- list()

    for (i in 1:nrow(bp_sequences)) {
        bp <- bp_sequences[i, ]
        sv <- sv_df[i, ]

        has_ins <- FALSE
        ins_len <- 0
        ins_seq <- NA
        ins_type <- "none"

        # Check for INSSEQ (Standard in many callers like Manta, GRIDSS)
        if ("INSSEQ" %in% names(sv) && !is.na(sv$INSSEQ) && nchar(as.character(sv$INSSEQ)) > 0) {
            has_ins <- TRUE
            ins_seq <- as.character(sv$INSSEQ)
            ins_len <- nchar(ins_seq)
            ins_type <- "non-templated"
        } else if ("ALT" %in% names(sv) && !is.na(sv$ALT) && grepl("^[ACGTN]+$", as.character(sv$ALT))) {
            alt_seq <- as.character(sv$ALT)
            ref_seq <- if ("REF" %in% names(sv) && !is.na(sv$REF)) as.character(sv$REF) else ""
            
            if (nchar(ref_seq) > 0 && nchar(alt_seq) > nchar(ref_seq)) {
                has_ins <- TRUE
                ins_seq <- alt_seq
                ins_len <- nchar(alt_seq) - nchar(ref_seq)
                ins_type <- "vcf_alt"
            }
        }

        insertions[[i]] <- list(
            sv_id = bp$sv_id,
            has_insertion = has_ins,
            insertion_length = ins_len,
            insertion_seq = ins_seq,
            insertion_type = ins_type
        )
    }

    # Convert to data frame
    ins_df <- do.call(rbind, lapply(insertions, function(x) {
        data.frame(
            sv_id = as.character(x$sv_id),
            has_insertion = as.logical(x$has_insertion),
            insertion_length = as.numeric(x$insertion_length),
            insertion_seq = as.character(x$insertion_seq),
            insertion_type = as.character(x$insertion_type),
            stringsAsFactors = FALSE
        )
    }))

    return(ins_df)
}


#' Classify DNA repair mechanisms based on sequence features
#'
#' @param microhomology Microhomology data frame
#' @param insertions Insertions data frame
#' @param sv_df Original SV data frame
#' @return Data frame with repair mechanism classification
#' @keywords internal
classify_repair_mechanisms <- function(microhomology, insertions, sv_df) {

    mechanisms <- list()

    for (i in 1:nrow(microhomology)) {
        mh <- microhomology[i, ]
        ins <- insertions[i, ]
        sv <- sv_df[i, ]
        
        # Robust ID retrieval
        current_id <- as.character(mh$sv_id)

        # Advanced Classification Logic
        mechanism <- "Unknown"
        confidence <- "Low"
        
        has_mh <- !is.na(mh$has_microhomology) && mh$has_microhomology
        mh_len <- ifelse(is.na(mh$microhomology_length), 0, mh$microhomology_length)
        has_ins <- !is.na(ins$has_insertion) && ins$has_insertion
        ins_len <- ifelse(is.na(ins$insertion_length), 0, ins$insertion_length)

        if (has_mh) {
            if (mh_len >= 2 && mh_len <= 20) {
                mechanism <- "MMEJ"
                confidence <- "High"
                
                if (sv$SVtype == "DUP" && mh_len >= 2 && mh_len <= 5) {
                    mechanism <- "MMBIR/FoSTeS"
                    confidence <- "Moderate"
                }
            } else if (mh_len > 20) {
                mechanism <- "SSA"
                confidence <- "High"
            } else if (mh_len == 1) {
                mechanism <- "NHEJ"
                confidence <- "Moderate"
            }
        } else {
            if (ins_len > 0) {
                mechanism <- "NHEJ"
                confidence <- "High"
            } else {
                mechanism <- "NHEJ"
                confidence <- "Moderate"
            }
        }

        mechanisms[[i]] <- list(
            sv_id = current_id,
            repair_mechanism = mechanism,
            confidence = confidence,
            evidence = sprintf("MH: %d bp, Ins: %d bp", mh_len, ins_len)
        )
    }

    # Convert to data frame
    mech_df <- do.call(rbind, lapply(mechanisms, function(x) {
        data.frame(
            sv_id = as.character(x$sv_id),
            repair_mechanism = as.character(x$repair_mechanism),
            confidence = as.character(x$confidence),
            evidence = as.character(x$evidence),
            stringsAsFactors = FALSE
        )
    }))

    return(mech_df)
}


#' Summarize repair mechanisms for a specific genomic region
#'
#' @param repair_data Result object from analyze_breakpoint_sequences()
#' @param sv_ids IDs of SVs belonging to the region
#' @return A concise string summarizing the dominant repair signatures
#' @export
summarize_repair_per_region <- function(repair_data, sv_ids) {
    if (is.null(repair_data) || length(sv_ids) == 0) return("Unknown")
    
    # Ensure IDs are characters for robust matching
    target_ids <- as.character(sv_ids)
    avail_ids <- as.character(repair_data$repair_mechanisms$sv_id)
    
    sub_repair <- repair_data$repair_mechanisms[avail_ids %in% target_ids, ]
    if (nrow(sub_repair) == 0) return("Unknown")
    
    counts <- table(sub_repair$repair_mechanism)
    props <- counts / sum(counts)
    
    sorted_props <- sort(props, decreasing = TRUE)
    top_mechs <- names(sorted_props)
    
    summary_parts <- c()
    for (i in 1:min(2, length(top_mechs))) {
        summary_parts <- c(summary_parts, sprintf("%s (%.0f%%)", top_mechs[i], sorted_props[i] * 100))
    }
    
    return(paste(summary_parts, collapse = ", "))
}


#' Create summary statistics for sequence analysis
#'
#' @param microhomology Microhomology data frame
#' @param insertions Insertions data frame
#' @param repair_classification Repair mechanism data frame
#' @return List with summary statistics
#' @keywords internal
create_sequence_summary <- function(microhomology, insertions, repair_classification) {

    mh_present <- sum(microhomology$has_microhomology, na.rm = TRUE)
    mh_total <- sum(!is.na(microhomology$has_microhomology))
    mh_proportion <- if (mh_total > 0) mh_present / mh_total else 0

    mh_lengths <- microhomology$microhomology_length[!is.na(microhomology$microhomology_length)]
    mh_mean_length <- if (length(mh_lengths) > 0) mean(mh_lengths) else NA
    mh_median_length <- if (length(mh_lengths) > 0) median(mh_lengths) else NA

    mech_table <- table(repair_classification$repair_mechanism)

    summary <- list(
        n_breakpoints = nrow(microhomology),
        microhomology = list(
            n_with_mh = mh_present,
            proportion = mh_proportion,
            mean_length = mh_mean_length,
            median_length = mh_median_length,
            length_range = if (length(mh_lengths) > 0) range(mh_lengths) else c(NA, NA)
        ),
        repair_mechanisms = as.data.frame(mech_table, stringsAsFactors = FALSE),
        dominant_mechanism = if (length(mech_table) > 0) names(mech_table)[which.max(mech_table)] else "Unknown"
    )

    return(summary)
}


#' Print method for breakpoint sequence analysis
#'
#' @param x Breakpoint sequence analysis result
#' @param ... Additional arguments
#' @export
print.breakpoint_sequences <- function(x, ...) {
    cat("\n")
    cat(rep("=", 70), "\n", sep = "")
    cat("         BREAKPOINT SEQUENCE ANALYSIS SUMMARY\n")
    cat(rep("=", 70), "\n\n")

    cat(sprintf("Total breakpoints analyzed: %d\n\n", x$summary$n_breakpoints))

    cat("MICROHOMOLOGY:\n")
    cat(sprintf("  Breakpoints with microhomology: %d (%.1f%%)\n",
               x$summary$microhomology$n_with_mh,
               x$summary$microhomology$proportion * 100))
    if (!is.na(x$summary$microhomology$mean_length)) {
        cat(sprintf("  Mean length: %.1f bp\n", x$summary$microhomology$mean_length))
        cat(sprintf("  Median length: %.0f bp\n", x$summary$microhomology$median_length))
        cat(sprintf("  Range: %d-%d bp\n",
                   x$summary$microhomology$length_range[1],
                   x$summary$microhomology$length_range[2]))
    }
    cat("\n")

    cat("DNA REPAIR MECHANISMS:\n")
    cat(sprintf("  Dominant mechanism: %s\n\n", x$summary$dominant_mechanism))
    cat("  Distribution:\n")
    if (nrow(x$summary$repair_mechanisms) > 0) {
        for (i in 1:nrow(x$summary$repair_mechanisms)) {
            mech <- x$summary$repair_mechanisms$Var1[i]
            freq <- x$summary$repair_mechanisms$Freq[i]
            prop <- freq / x$summary$n_breakpoints * 100
            cat(sprintf("    - %s: %d (%.1f%%)\n", mech, freq, prop))
        }
    } else {
        cat("    - No data available\n")
    }

    cat("\n")
    cat(rep("=", 70), "\n")
    cat("\n")

    invisible(x)
}