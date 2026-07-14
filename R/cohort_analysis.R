#' @include chromoanagenesis_integrated.R
NULL

#' Class OncoImplexusCohort
#'
#' A class to represent a cohort of samples with chromoanagenesis analysis results
#' and associated clinical data.
#'
#' @slot results A list of chromoanagenesis result objects (one per sample)
#' @slot clinical A data frame containing clinical metadata (rownames must match sample IDs)
#' @slot sample_ids A character vector of sample IDs
#' @export
setClass("OncoImplexusCohort",
         slots = c(
             results = "list",
             clinical = "data.frame",
             sample_ids = "character"
         )
)

#' Create a cohort object for comparative analysis
#'
#' Integrates chromoanagenesis results from multiple samples with clinical metadata
#' to enable cohort-level comparative analysis.
#'
#' @param results_list A named list of chromoanagenesis result objects. Names must correspond to sample IDs.
#' @param clinical_data A data frame containing clinical variables (e.g., stage, subtype, survival).
#'   Must have a column named "sample_id" or rownames matching the names in results_list.
#' @return An OncoImplexusCohort object
#' @examples
#' \dontrun{
#' # Assuming results_list is a list of results from detect_chromoanagenesis
#' # and clinical_df contains metadata
#' cohort <- create_cohort(results_list, clinical_df)
#'
#' # Compare mechanisms by tumor stage
#' stats <- compare_mechanisms(cohort, group_by = "Stage")
#' }
#' @export
create_cohort <- function(results_list, clinical_data) {

    # Validation
    if (!is.list(results_list) || length(results_list) == 0) {
        stop("results_list must be a non-empty list of chromoanagenesis results")
    }

    if (is.null(names(results_list))) {
        stop("results_list must be a named list (names = sample IDs)")
    }

    # Handle clinical data sample IDs
    if (!"sample_id" %in% colnames(clinical_data)) {
        if (!is.null(rownames(clinical_data))) {
            clinical_data$sample_id <- rownames(clinical_data)
        } else {
            stop("clinical_data must have a 'sample_id' column or valid rownames")
        }
    }

    # Intersect samples
    common_samples <- intersect(names(results_list), clinical_data$sample_id)

    if (length(common_samples) == 0) {
        stop("No matching sample IDs found between results_list and clinical_data")
    }

    if (length(common_samples) < length(results_list)) {
        warning(sprintf("Only %d samples matched between results and clinical data (out of %d total results)",
                        length(common_samples), length(results_list)))
    }

    # Subset and order data
    filtered_results <- results_list[common_samples]
    
    # Use drop=FALSE to prevent conversion to vector when subsetting
    filtered_clinical <- clinical_data[clinical_data$sample_id %in% common_samples, , drop = FALSE]
    rownames(filtered_clinical) <- filtered_clinical$sample_id
    filtered_clinical <- filtered_clinical[common_samples, , drop = FALSE]

    # Verify result object types
    for (i in seq_along(filtered_results)) {
        # Relax inheritance check to support list-wrapped results as well
        if (!inherits(filtered_results[[i]], "chromoanagenesis") && !is.list(filtered_results[[i]])) {
            warning(sprintf("Result for sample %s is not a valid result object. Skipping.",
                            names(filtered_results)[i]))
        }
    }

    new("OncoImplexusCohort",
        results = filtered_results,
        clinical = filtered_clinical,
        sample_ids = common_samples
    )
}

#' Compare chromoanagenesis mechanisms between groups
#'
#' Performs statistical tests to determine if the frequency of chromoanagenesis mechanisms
#' differs significantly between clinical groups (e.g., subtypes, stages).
#'
#' @param cohort An OncoImplexusCohort object
#' @param group_by Column name in clinical data to group samples by
#' @param stringency Analysis stringency. "strict" (default) considers only High confidence/Likely events.
#'   "inclusive" considers High+Low confidence and Likely+Possible events.
#' @return A list containing:
#'   \item{summary_table}{Frequency of each mechanism per group}
#'   \item{stats_table}{P-values from Fisher's exact test or Chi-squared test}
#' @export
compare_mechanisms <- function(cohort, group_by, stringency = "strict") {

    if (!group_by %in% colnames(cohort@clinical)) {
        stop(sprintf("Column '%s' not found in clinical data", group_by))
    }
    
    if (!stringency %in% c("strict", "inclusive")) {
        stop("stringency must be either 'strict' or 'inclusive'")
    }

    groups <- cohort@clinical[[group_by]]
    unique_groups <- unique(groups)
    unique_groups <- unique_groups[!is.na(unique_groups)]

    if (length(unique_groups) < 2) {
        stop("group_by variable must have at least 2 levels for comparison")
    }

    # Extract mechanism status for all samples
    mechanism_matrix <- do.call(rbind, lapply(cohort@sample_ids, function(sid) {
        res <- cohort@results[[sid]]
        
        # Check Chromothripsis
        has_chromothripsis <- 0
        if (!is.null(res$chromothripsis)) {
            if (stringency == "strict") {
                if (res$chromothripsis$n_high_confidence > 0) has_chromothripsis <- 1
            } else {
                if ((res$chromothripsis$n_high_confidence + res$chromothripsis$n_low_confidence) > 0) has_chromothripsis <- 1
            }
        }

        # Check Chromoplexy
        has_chromoplexy <- 0
        if (!is.null(res$chromoplexy)) {
            if (stringency == "strict") {
                if (res$chromoplexy$likely_chromoplexy > 0) has_chromoplexy <- 1
            } else {
                if ((res$chromoplexy$likely_chromoplexy + res$chromoplexy$possible_chromoplexy) > 0) has_chromoplexy <- 1
            }
        }

        # Check Chromoanasynthesis
        has_chromoanasynthesis <- 0
        if (!is.null(res$chromoanasynthesis)) {
            if (stringency == "strict") {
                if (res$chromoanasynthesis$likely_chromoanasynthesis > 0) has_chromoanasynthesis <- 1
            } else {
                if ((res$chromoanasynthesis$likely_chromoanasynthesis + res$chromoanasynthesis$possible_chromoanasynthesis) > 0) has_chromoanasynthesis <- 1
            }
        }

        return(c(Chromothripsis = has_chromothripsis,
                 Chromoplexy = has_chromoplexy,
                 Chromoanasynthesis = has_chromoanasynthesis))
    }))

    mechanism_df <- as.data.frame(mechanism_matrix)
    mechanism_df$Group <- groups

    # Remove samples with NA group
    mechanism_df <- mechanism_df[!is.na(mechanism_df$Group), ]

    # Initialize results
    stats_results <- data.frame(
        Mechanism = character(0),
        Test = character(0),
        P_value = numeric(0),
        Significant = character(0),
        stringsAsFactors = FALSE
    )

    mechanisms <- c("Chromothripsis", "Chromoplexy", "Chromoanasynthesis")

    cat(sprintf("Comparing mechanisms by '%s' (%d groups: %s)\n", 
                group_by, length(unique_groups), paste(unique_groups, collapse=", ")))
    cat(sprintf("Stringency: %s\n\n", stringency))

    for (mech in mechanisms) {
        # Create contingency table
        tbl <- table(mechanism_df$Group, mechanism_df[[mech]])
        
        # Ensure table has 2 columns (0 and 1)
        if (ncol(tbl) == 1) {
            # Handle case where all are 0 or all are 1
            cat(sprintf("  - %s: Cannot test (no variation in outcome)\n", mech))
            next
        }

        # Perform test
        # Use Fisher's exact test for 2x2, Chi-squared for larger
        test_name <- ""
        p_val <- NA

        if (nrow(tbl) == 2) {
            test_res <- fisher.test(tbl)
            test_name <- "Fisher's Exact"
            p_val <- test_res$p.value
        } else {
            # Check expected counts for Chi-squared validity
            test_res <- chisq.test(tbl, simulate.p.value = TRUE, B = 2000) # Use simulation for robustness
            test_name <- "Chi-squared"
            p_val <- test_res$p.value
        }

        sig_star <- ""
        if (p_val < 0.001) sig_star <- "***"
        else if (p_val < 0.01) sig_star <- "**"
        else if (p_val < 0.05) sig_star <- "*"

        stats_results[nrow(stats_results) + 1, ] <- list(
            Mechanism = mech,
            Test = test_name,
            P_value = p_val,
            Significant = sig_star
        )
        
        # Print summary
        cat(sprintf("  - %s: p = %.4g %s (%s)\n", mech, p_val, sig_star, test_name))
        
        # Print rates per group
        rates <- prop.table(tbl, margin = 1)[, "1"] * 100
        for (g in names(rates)) {
            cat(sprintf("    %s: %.1f%%\n", g, rates[g]))
        }
        cat("\n")
    }

    return(list(
        summary_table = table(mechanism_df$Group, mechanism_df$Chromothripsis), # Simplified summary for now
        stats_table = stats_results
    ))
}

#' Process a cohort of samples from a sample sheet
#'
#' Batch processes multiple samples defined in a CSV/TSV file or data frame.
#' Reads SV and CNV files for each sample and runs comprehensive chromoanagenesis detection.
#' Supports parallel processing for faster execution.
#'
#' @param sample_sheet Path to a CSV/TSV file or a data frame containing sample information.
#'   Must contain column: 'sv_file'. Optional column: 'cnv_file' for full SV+CNV analysis.
#'   Optional columns: 'sample_id' (if missing, row numbers or filenames are used).
#' @param output_dir Optional directory to save individual result objects (.rds files).
#'   Useful for caching results of large cohorts.
#' @param num_cores Number of cores to use for parallel processing (default: 1).
#'   Requires 'parallel' package.
#' @param ... Additional arguments passed to detect_chromoanagenesis (e.g., genome="hg38").
#' @return A named list of chromoanagenesis results, ready for create_cohort().
#' @examples
#' \dontrun{
#' # Create a simple sample sheet
#' sheet <- data.frame(
#'   sample_id = c("Sample1", "Sample2"),
#'   sv_file = c("path/to/s1.sv.vcf", "path/to/s2.sv.vcf"),
#'   cnv_file = c("path/to/s1.cnv.vcf", "path/to/s2.cnv.vcf")
#' )
#' 
#' # Run batch analysis
#' results <- process_cohort_files(sheet, genome = "hg19", num_cores = 2)
#' 
#' # Create cohort object with clinical data
#' cohort <- create_cohort(results, clinical_data)
#' }
#' @export
process_cohort_files <- function(sample_sheet, 
                                output_dir = NULL, 
                                num_cores = 1, 
                                ...) {
    
    # 1. Parse sample sheet
    if (is.character(sample_sheet)) {
        if (!file.exists(sample_sheet)) {
            stop(sprintf("Sample sheet file not found: %s", sample_sheet))
        }
        # Auto-detect separator based on extension
        if (grepl("\\.csv$", sample_sheet, ignore.case = TRUE)) {
            df <- read.csv(sample_sheet, stringsAsFactors = FALSE)
        } else {
            df <- read.table(sample_sheet, header = TRUE, stringsAsFactors = FALSE)
        }
    } else if (is.data.frame(sample_sheet)) {
        df <- sample_sheet
    } else {
        stop("sample_sheet must be a file path or data frame")
    }
    
    # 2. Validate columns
    required_cols <- c("sv_file")
    missing_cols <- setdiff(required_cols, colnames(df))
    if (length(missing_cols) > 0) {
        stop(sprintf("Sample sheet missing required columns: %s", paste(missing_cols, collapse = ", ")))
    }
    if (!"cnv_file" %in% colnames(df)) {
        df$cnv_file <- NA_character_
    }
    
    # Handle sample IDs
    if (!"sample_id" %in% colnames(df)) {
        warning("No 'sample_id' column found. Using SV filenames as IDs.")
        df$sample_id <- basename(df$sv_file)
    }
    
    # Ensure IDs are unique
    if (any(duplicated(df$sample_id))) {
        stop("Duplicate sample IDs found in sample sheet.")
    }
    
    # 3. Setup Output Directory
    if (!is.null(output_dir)) {
        if (!dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE)
        }
        cat(sprintf("Results will be cached in: %s\n", output_dir))
    }
    
    # 4. Define Processing Function (per sample)
    process_single_sample <- function(i) {
        sid <- df$sample_id[i]
        sv_path <- df$sv_file[i]
        cnv_path <- df$cnv_file[i]
        
        # Check cache first
        if (!is.null(output_dir)) {
            cache_file <- file.path(output_dir, paste0(sid, "_chromoanagenesis.rds"))
            if (file.exists(cache_file)) {
                message(sprintf("[%d/%d] Loading cached result for %s...", i, nrow(df), sid))
                return(readRDS(cache_file))
            }
        }
        
        message(sprintf("[%d/%d] Processing sample: %s...", i, nrow(df), sid))
        
        # Check input files
        if (!file.exists(sv_path)) {
            warning(sprintf("SV file not found for sample %s: %s", sid, sv_path))
            return(NULL)
        }
        has_cnv <- !is.na(cnv_path) && nzchar(cnv_path) && file.exists(cnv_path)
        if (!has_cnv) {
            message(sprintf("  CNV file not available for %s. Running SV-only chromoplexy analysis.", sid))
        }
        
        tryCatch({
            # Read Data
            # Note: We assume the package functions are available in the worker environment
            # When using parallel::mclapply (forking), this is automatic.
            # For snow/cluster, libraries need to be loaded on workers.
            
            # Using min_sv_size=0 by default for robustness as per our testing
            sv_data <- read_sv_vcf(sv_path, min_sv_size = 0)
            cnv_data <- if (has_cnv) read_cnv_vcf(cnv_path) else NULL
            
            # Run Detection
            # Pass ellipses (...) arguments to detect_chromoanagenesis
            res <- detect_chromoanagenesis(
                SV.sample = sv_data,
                CNV.sample = cnv_data,
                verbose = FALSE,
                ...
            )
            
            # Save Cache
            if (!is.null(output_dir)) {
                cache_file <- file.path(output_dir, paste0(sid, "_chromoanagenesis.rds"))
                saveRDS(res, cache_file)
            }
            
            return(res)
            
        }, error = function(e) {
            warning(sprintf("Error processing sample %s: %s", sid, e$message))
            return(NULL)
        })
    }
    
    # 5. Run Processing (Serial or Parallel)
    results_list <- list()
    
    if (num_cores > 1 && requireNamespace("parallel", quietly = TRUE)) {
        message(sprintf("Running in parallel with %d cores...", num_cores))
        
        # Use mclapply for forking (Unix-like) or parLapply for Windows
        if (.Platform$OS.type == "unix") {
            results_list <- parallel::mclapply(1:nrow(df), process_single_sample, mc.cores = num_cores)
        } else {
            # Windows setup
            cl <- parallel::makeCluster(num_cores)
            # Export necessary libraries/functions if needed
            parallel::clusterEvalQ(cl, {
                library(OncoImplexus)
            })
            results_list <- parallel::parLapply(cl, 1:nrow(df), process_single_sample)
            parallel::stopCluster(cl)
        }
    } else {
        # Serial processing
        results_list <- lapply(1:nrow(df), process_single_sample)
    }
    
    # 6. Format Output
    names(results_list) <- df$sample_id
    
    # Remove NULLs (failed samples)
    failed_count <- sum(sapply(results_list, is.null))
    if (failed_count > 0) {
        warning(sprintf("%d samples failed processing.", failed_count))
        results_list <- results_list[!sapply(results_list, is.null)]
    }
    
    message(sprintf("Successfully processed %d samples.", length(results_list)))
    return(results_list)
}

#' Summarize cohort results into a CSV file
#'
#' Scans a directory of chromoanagenesis result objects (.rds) and creates
#' a summary CSV table with key metrics for each sample.
#'
#' @param results_dir Directory containing .rds files.
#' @param output_csv Path to save the summary CSV.
#' @return A data frame containing the summary.
#' @export
summarize_cohort_results <- function(results_dir, output_csv = "cohort_summary.csv") {
    files <- list.files(results_dir, pattern = "\\.rds$", full.names = TRUE)
    if (length(files) == 0) stop("No RDS files found.")
    
    summary_list <- list()
    for (f in files) {
        sid <- gsub("\\.rds$", "", basename(f))
        tryCatch({
            obj <- readRDS(f)
            # Handle different wrapper structures
            res <- if ("results" %in% names(obj)) obj$results else obj
            summary_info <- res$integrated_summary
            
            # Robust extraction of counts
            n_fusions <- 0
            if (!is.null(res$fusions) && is.data.frame(res$fusions)) {
                n_fusions <- nrow(res$fusions)
            }
            
            row <- data.frame(
                sample_id = sid,
                classification = if(!is.null(summary_info$overall_classification)) as.character(summary_info$overall_classification) else "Unknown",
                ct_hc = if(!is.null(summary_info$chromothripsis_high_confidence)) as.numeric(summary_info$chromothripsis_high_confidence) else 0,
                cp_likely = if(!is.null(summary_info$chromoplexy_likely)) as.numeric(summary_info$chromoplexy_likely) else 0,
                cs_likely = if(!is.null(summary_info$chromoanasynthesis_likely)) as.numeric(summary_info$chromoanasynthesis_likely) else 0,
                n_fusions = n_fusions,
                stringsAsFactors = FALSE
            )
            summary_list[[sid]] <- row
        }, error = function(e) {
            warning("Failed to summarize ", sid, ": ", e$message)
        })
    }
    
    final_df <- do.call(rbind, summary_list)
    write.csv(final_df, output_csv, row.names = FALSE)
    message("Summary saved to ", output_csv)
    return(final_df)
}

#' Summarize collapsed chromoplexy events across a cohort
#'
#' Aggregates event-level chromoplexy calls produced by
#' \code{collapse_chromoplexy_chains()} across many samples. The output is
#' designed for recurrent gene, recurrent breakpoint region, and per-sample
#' burden analysis.
#'
#' @param results A named list of chromoanagenesis results, an
#'   \code{OncoImplexusCohort} object, or a directory containing result
#'   \code{.rds} files.
#' @param output_dir Optional directory. If supplied, summary tables are written
#'   as TSV files.
#' @param breakpoint_window Window size in bp for recurrent breakpoint region
#'   aggregation.
#' @return A list containing sample, event, gene, chromosome, and breakpoint
#'   region summary tables.
#' @export
summarize_chromoplexy_cohort_events <- function(results,
                                                output_dir = NULL,
                                                breakpoint_window = 1e6) {
    results_list <- coerce_chromoplexy_results_list(results)
    if (length(results_list) == 0) {
        stop("No chromoanagenesis results found")
    }

    event_rows <- list()
    gene_rows <- list()
    breakpoint_rows <- list()
    chromosome_rows <- list()
    sample_rows <- list()

    for (sid in names(results_list)) {
        res <- unwrap_chromoanagenesis_result(results_list[[sid]])
        cp <- if (!is.null(res$chromoplexy)) res$chromoplexy else res
        ce <- cp$collapsed_events
        if (is.null(ce) || is.null(ce$event_summary)) {
            ce <- if (!is.null(cp$summary)) collapse_chromoplexy_chains(cp) else empty_collapsed_chromoplexy_events()
        }

        events <- ce$event_summary
        if (!is.null(events) && nrow(events) > 0) {
            events$sample_id <- sid
            events$cohort_event_id <- paste(sid, events$collapsed_event_id, sep = ":")
            event_rows[[sid]] <- events

            event_chroms <- unique(unlist(strsplit(as.character(events$chromosomes_involved), ","), use.names = FALSE))
            event_chroms <- event_chroms[nzchar(event_chroms)]
            if (length(event_chroms) > 0) {
                chromosome_rows[[sid]] <- data.frame(
                    sample_id = sid,
                    chrom = event_chroms,
                    stringsAsFactors = FALSE
                )
            }
        }

        if (!is.null(ce$gene_detail) && nrow(ce$gene_detail) > 0) {
            genes <- unique(ce$gene_detail[, intersect(c(
                "collapsed_event_id", "gene_id", "symbol", "gene_type",
                "is_driver", "chrom", "breakpoint_id"
            ), colnames(ce$gene_detail)), drop = FALSE])
            if (!"gene_id" %in% colnames(genes)) genes$gene_id <- genes$symbol
            if (!"gene_type" %in% colnames(genes)) genes$gene_type <- "Unknown"
            if (!"is_driver" %in% colnames(genes)) genes$is_driver <- FALSE
            if (!"chrom" %in% colnames(genes)) genes$chrom <- ""
            genes$sample_id <- sid
            genes$cohort_event_id <- paste(sid, genes$collapsed_event_id, sep = ":")
            gene_rows[[sid]] <- genes
        } else if (!is.null(events) && nrow(events) > 0 && "genes" %in% colnames(events)) {
            gene_rows[[sid]] <- event_summary_gene_rows(events, sid)
        }

        if (!is.null(ce$event_breakpoints) && nrow(ce$event_breakpoints) > 0) {
            bp <- ce$event_breakpoints
            bp$sample_id <- sid
            bp$cohort_event_id <- paste(sid, bp$collapsed_event_id, sep = ":")
            bp$window_start <- floor((as.numeric(bp$pos) - 1) / breakpoint_window) * breakpoint_window + 1
            bp$window_end <- bp$window_start + breakpoint_window - 1
            breakpoint_rows[[sid]] <- bp
        }

        sample_rows[[sid]] <- summarize_single_chromoplexy_sample(sid, res, ce)
    }

    sample_summary <- do.call(rbind, sample_rows)
    event_summary <- bind_rows_or_empty(event_rows)
    gene_detail <- bind_rows_or_empty(gene_rows)
    breakpoint_detail <- bind_rows_or_empty(breakpoint_rows)
    chromosome_detail <- bind_rows_or_empty(chromosome_rows)

    gene_summary <- summarize_recurrent_chromoplexy_genes(gene_detail, event_summary, length(results_list))
    event_summary <- add_cohort_recurrence_to_events(event_summary, gene_detail, gene_summary)
    chromosome_summary <- summarize_recurrent_chromosomes(chromosome_detail)
    breakpoint_region_summary <- summarize_recurrent_breakpoint_regions(breakpoint_detail)

    out <- list(
        sample_summary = sample_summary,
        event_summary = event_summary,
        gene_summary = gene_summary,
        chromosome_summary = chromosome_summary,
        breakpoint_region_summary = breakpoint_region_summary,
        gene_detail = gene_detail,
        breakpoint_detail = breakpoint_detail,
        parameters = list(
            n_samples = length(results_list),
            breakpoint_window = breakpoint_window
        )
    )

    if (!is.null(output_dir)) {
        if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
        write_summary_tsv(out$sample_summary, file.path(output_dir, "chromoplexy_sample_summary.tsv"))
        write_summary_tsv(out$event_summary, file.path(output_dir, "chromoplexy_collapsed_event_summary.tsv"))
        write_summary_tsv(out$gene_summary, file.path(output_dir, "chromoplexy_recurrent_gene_summary.tsv"))
        write_summary_tsv(out$chromosome_summary, file.path(output_dir, "chromoplexy_recurrent_chromosome_summary.tsv"))
        write_summary_tsv(out$breakpoint_region_summary, file.path(output_dir, "chromoplexy_recurrent_breakpoint_regions.tsv"))
        write_summary_tsv(out$gene_detail, file.path(output_dir, "chromoplexy_gene_detail.tsv"))
        write_summary_tsv(out$breakpoint_detail, file.path(output_dir, "chromoplexy_breakpoint_detail.tsv"))
    }

    out
}

#' Flag chromoplexy events that recur at near-identical breakpoints across
#' unrelated cohort samples
#'
#' Chromoplexy is, by definition, a one-off catastrophic event private to a
#' single tumor clone. If the *same* breakpoint (within \code{artifact_window}
#' bp) is called in multiple, unrelated samples, that is evidence of a
#' systematic artifact (e.g. mis-mapping in a segmental-duplication/repeat
#' region, or a common structural polymorphism) rather than independent
#' patient-specific rearrangements. Somatic cohorts usually filter this class
#' of noise via matched-normal subtraction or a panel-of-normals; germline
#' cohorts have no such step, so this check is especially important there.
#'
#' This reuses \code{summarize_chromoplexy_cohort_events()}'s breakpoint
#' aggregation with a much tighter window than its cancer-hotspot default
#' (1e6 bp), and, unlike that function's \code{recurrence_score} (which
#' treats gene-level recurrence as supporting evidence, appropriate for
#' recurrent cancer driver events), treats breakpoint-level cross-sample
#' recurrence as a red flag: matching events are marked
#' \code{artifact_suspected = TRUE} and their confidence is downgraded rather
#' than boosted.
#'
#' @param results A named list of chromoanagenesis results, an
#'   \code{OncoImplexusCohort} object, or a directory of result \code{.rds}
#'   files.
#' @param artifact_window Window size in bp for considering two breakpoints
#'   from different samples "the same". Default 1000.
#' @param min_samples Minimum number of unrelated samples sharing a
#'   breakpoint region before events there are flagged. Default 2.
#' @return A list with \code{event_summary} (per-sample collapsed chromoplexy
#'   events, with \code{artifact_suspected}, \code{n_samples_sharing_breakpoint}
#'   and a downgraded \code{event_confidence} for flagged rows) and
#'   \code{artifact_regions} (the recurrent breakpoint regions driving the
#'   flags, for manual review e.g. against a segmental-duplication track).
#' @export
flag_recurrent_chromoplexy_artifacts <- function(results, artifact_window = 1000, min_samples = 2) {
    cohort <- summarize_chromoplexy_cohort_events(results, breakpoint_window = artifact_window)
    events <- cohort$event_summary
    regions <- cohort$breakpoint_region_summary

    if (nrow(events) == 0) {
        events$artifact_suspected <- logical(0)
        events$n_samples_sharing_breakpoint <- integer(0)
        return(list(event_summary = events, artifact_regions = regions))
    }

    events$artifact_suspected <- FALSE
    events$n_samples_sharing_breakpoint <- 1L

    if (nrow(regions) > 0) {
        artifact_regions <- regions[regions$n_samples >= min_samples, , drop = FALSE]
        if (nrow(artifact_regions) > 0) {
            for (i in seq_len(nrow(artifact_regions))) {
                flagged_ids <- unlist(strsplit(artifact_regions$cohort_event_ids[i], ",", fixed = TRUE))
                match_idx <- events$cohort_event_id %in% flagged_ids
                events$artifact_suspected[match_idx] <- TRUE
                events$n_samples_sharing_breakpoint[match_idx] <- pmax(
                    events$n_samples_sharing_breakpoint[match_idx],
                    artifact_regions$n_samples[i]
                )
            }
        }
    } else {
        artifact_regions <- regions
    }

    if ("event_confidence" %in% colnames(events)) {
        events$event_confidence_original <- events$event_confidence
        events$event_confidence[events$artifact_suspected] <- "Low (recurrent cross-sample artifact)"
    }

    events <- events[order(-events$artifact_suspected, events$sample_id), , drop = FALSE]

    list(event_summary = events, artifact_regions = artifact_regions)
}

coerce_chromoplexy_results_list <- function(results) {
    if (inherits(results, "OncoImplexusCohort")) {
        out <- results@results
    } else if (is.character(results) && length(results) == 1 && dir.exists(results)) {
        files <- list.files(results, pattern = "\\.rds$", full.names = TRUE)
        out <- lapply(files, readRDS)
        names(out) <- clean_result_sample_ids(files)
    } else if (is.list(results)) {
        out <- results
    } else {
        stop("results must be a named list, OncoImplexusCohort object, or RDS directory")
    }

    if (is.null(names(out)) || any(!nzchar(names(out)))) {
        names(out) <- paste0("Sample_", seq_along(out))
    }
    out
}

clean_result_sample_ids <- function(files) {
    ids <- sub("\\.rds$", "", basename(files))
    ids <- sub("_chromoanagenesis$", "", ids)
    ids <- sub("_result$", "", ids)
    ids
}

unwrap_chromoanagenesis_result <- function(obj) {
    if (is.list(obj) && "results" %in% names(obj)) obj$results else obj
}

summarize_single_chromoplexy_sample <- function(sample_id, res, collapsed_events) {
    cp <- if (!is.null(res$chromoplexy)) res$chromoplexy else res
    events <- collapsed_events$event_summary
    gene_summary <- collapsed_events$gene_event_summary

    data.frame(
        sample_id = sample_id,
        analysis_mode = if (!is.null(res$analysis_mode)) res$analysis_mode else if (!is.null(cp$analysis_mode)) cp$analysis_mode else NA_character_,
        total_chains = if (!is.null(cp$total_chains)) cp$total_chains else 0,
        likely_chains = if (!is.null(cp$likely_chromoplexy)) cp$likely_chromoplexy else 0,
        possible_chains = if (!is.null(cp$possible_chromoplexy)) cp$possible_chromoplexy else 0,
        n_collapsed_events = if (!is.null(events)) nrow(events) else 0,
        n_high_confidence_events = if (!is.null(events) && nrow(events) > 0) sum(events$event_confidence == "High", na.rm = TRUE) else 0,
        mean_event_qc_score = if (!is.null(events) && nrow(events) > 0) mean(events$event_qc_score, na.rm = TRUE) else NA_real_,
        max_event_qc_score = if (!is.null(events) && nrow(events) > 0) max(events$event_qc_score, na.rm = TRUE) else NA_real_,
        n_genes = if (!is.null(gene_summary) && nrow(gene_summary) > 0) length(unique(gene_summary$symbol)) else 0,
        n_driver_genes = if (!is.null(gene_summary) && nrow(gene_summary) > 0 && "is_driver" %in% colnames(gene_summary)) {
            length(unique(gene_summary$symbol[gene_summary$is_driver]))
        } else 0,
        driver_genes = if (!is.null(gene_summary) && nrow(gene_summary) > 0 && "is_driver" %in% colnames(gene_summary)) {
            paste(sort(unique(gene_summary$symbol[gene_summary$is_driver])), collapse = ",")
        } else "",
        stringsAsFactors = FALSE
    )
}

event_summary_gene_rows <- function(events, sample_id) {
    rows <- list()
    for (i in seq_len(nrow(events))) {
        genes <- unlist(strsplit(as.character(events$genes[i]), ","), use.names = FALSE)
        genes <- genes[nzchar(genes)]
        if (length(genes) == 0) next
        rows[[length(rows) + 1]] <- data.frame(
            collapsed_event_id = events$collapsed_event_id[i],
            gene_id = genes,
            symbol = genes,
            gene_type = "Unknown",
            is_driver = genes %in% get_default_drivers(),
            sample_id = sample_id,
            cohort_event_id = paste(sample_id, events$collapsed_event_id[i], sep = ":"),
            stringsAsFactors = FALSE
        )
    }
    bind_rows_or_empty(rows)
}

summarize_recurrent_chromoplexy_genes <- function(gene_detail, event_summary, n_samples_total) {
    if (nrow(gene_detail) == 0) return(data.frame())
    if (nrow(event_summary) > 0) {
        event_scores <- event_summary[, intersect(c("cohort_event_id", "event_qc_score", "event_priority_score"), colnames(event_summary)), drop = FALSE]
        gene_detail <- merge(gene_detail, event_scores, by = "cohort_event_id", all.x = TRUE)
    }
    symbols <- sort(unique(gene_detail$symbol[nzchar(gene_detail$symbol)]))
    rows <- lapply(symbols, function(sym) {
        sub <- gene_detail[gene_detail$symbol == sym, , drop = FALSE]
        samples <- sort(unique(sub$sample_id))
        events <- sort(unique(sub$cohort_event_id))
        gene_type_value <- if ("gene_type" %in% colnames(sub) && any(!is.na(sub$gene_type))) {
            sub$gene_type[which(!is.na(sub$gene_type))[1]]
        } else {
            "Unknown"
        }
        is_driver_value <- if ("is_driver" %in% colnames(sub)) any(sub$is_driver, na.rm = TRUE) else FALSE
        data.frame(
            symbol = sym,
            gene_id = paste(sort(unique(sub$gene_id)), collapse = ","),
            gene_type = gene_type_value,
            is_driver = is_driver_value,
            n_samples = length(samples),
            recurrence_fraction = length(samples) / n_samples_total,
            recurrence_score = min(length(samples) / 3, 1),
            n_collapsed_events = length(events),
            samples = paste(samples, collapse = ","),
            cohort_event_ids = paste(events, collapse = ","),
            chromosomes = if ("chrom" %in% colnames(sub)) paste(sort(unique(sub$chrom)), collapse = ",") else "",
            mean_event_qc_score = if ("event_qc_score" %in% colnames(sub)) mean_or_na(sub$event_qc_score) else NA_real_,
            max_event_qc_score = if ("event_qc_score" %in% colnames(sub)) max_or_na(sub$event_qc_score) else NA_real_,
            stringsAsFactors = FALSE
        )
    })
    out <- do.call(rbind, rows)
    out[order(-out$n_samples, -out$n_collapsed_events, -out$is_driver, out$symbol), , drop = FALSE]
}

add_cohort_recurrence_to_events <- function(event_summary, gene_detail, gene_summary) {
    if (nrow(event_summary) == 0) return(event_summary)
    event_summary$recurrence_score <- 0
    if (nrow(gene_detail) > 0 && nrow(gene_summary) > 0) {
        gs <- gene_summary[, c("symbol", "recurrence_score"), drop = FALSE]
        gd <- merge(unique(gene_detail[, c("cohort_event_id", "symbol"), drop = FALSE]), gs, by = "symbol", all.x = TRUE)
        recurrence_by_event <- tapply(gd$recurrence_score, gd$cohort_event_id, max_or_na)
        recurrence_by_event[is.na(recurrence_by_event)] <- 0
        hit <- match(event_summary$cohort_event_id, names(recurrence_by_event))
        event_summary$recurrence_score[!is.na(hit)] <- recurrence_by_event[hit[!is.na(hit)]]
    }
    if (!"driver_impact_score" %in% colnames(event_summary)) {
        event_summary$driver_impact_score <- 0
    }
    event_summary$event_cohort_priority_score <- (
        event_summary$event_qc_score * 0.70 +
        event_summary$driver_impact_score * 0.20 +
        event_summary$recurrence_score * 0.10
    )
    event_summary[order(-event_summary$event_cohort_priority_score, event_summary$sample_id), , drop = FALSE]
}

summarize_recurrent_chromosomes <- function(chromosome_detail) {
    if (nrow(chromosome_detail) == 0) return(data.frame())
    chroms <- sort(unique(chromosome_detail$chrom))
    rows <- lapply(chroms, function(chr) {
        sub <- chromosome_detail[chromosome_detail$chrom == chr, , drop = FALSE]
        data.frame(
            chrom = chr,
            n_samples = length(unique(sub$sample_id)),
            samples = paste(sort(unique(sub$sample_id)), collapse = ","),
            stringsAsFactors = FALSE
        )
    })
    out <- do.call(rbind, rows)
    out[order(-out$n_samples, out$chrom), , drop = FALSE]
}

summarize_recurrent_breakpoint_regions <- function(breakpoint_detail) {
    if (nrow(breakpoint_detail) == 0) return(data.frame())
    keys <- unique(breakpoint_detail[, c("chrom", "window_start", "window_end"), drop = FALSE])
    rows <- lapply(seq_len(nrow(keys)), function(i) {
        key <- keys[i, ]
        sub <- breakpoint_detail[
            breakpoint_detail$chrom == key$chrom &
                breakpoint_detail$window_start == key$window_start &
                breakpoint_detail$window_end == key$window_end,
            ,
            drop = FALSE
        ]
        data.frame(
            chrom = key$chrom,
            window_start = key$window_start,
            window_end = key$window_end,
            n_samples = length(unique(sub$sample_id)),
            n_collapsed_events = length(unique(sub$cohort_event_id)),
            n_breakpoints = length(unique(paste(sub$sample_id, sub$breakpoint_id, sep = ":"))),
            samples = paste(sort(unique(sub$sample_id)), collapse = ","),
            cohort_event_ids = paste(sort(unique(sub$cohort_event_id)), collapse = ","),
            stringsAsFactors = FALSE
        )
    })
    out <- do.call(rbind, rows)
    out[order(-out$n_samples, -out$n_collapsed_events, out$chrom, out$window_start), , drop = FALSE]
}

bind_rows_or_empty <- function(rows) {
    rows <- Filter(function(x) !is.null(x) && nrow(x) > 0, rows)
    if (length(rows) == 0) return(data.frame())
    # Per-sample tables can have different optional columns (e.g.
    # n_gene_overlapping_breakpoints is only present when a sample has gene
    # overlaps), so a plain rbind() would error; fill missing columns with NA.
    as.data.frame(data.table::rbindlist(rows, fill = TRUE, use.names = TRUE))
}

write_summary_tsv <- function(x, file) {
    utils::write.table(x, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
}
