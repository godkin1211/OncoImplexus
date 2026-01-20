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
#'   Must contain columns: 'sv_file' and 'cnv_file'. 
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
    required_cols <- c("sv_file", "cnv_file")
    missing_cols <- setdiff(required_cols, colnames(df))
    if (length(missing_cols) > 0) {
        stop(sprintf("Sample sheet missing required columns: %s", paste(missing_cols, collapse = ", ")))
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
        if (!file.exists(cnv_path)) {
            warning(sprintf("CNV file not found for sample %s: %s", sid, cnv_path))
            return(NULL)
        }
        
        tryCatch({
            # Read Data
            # Note: We assume the package functions are available in the worker environment
            # When using parallel::mclapply (forking), this is automatic.
            # For snow/cluster, libraries need to be loaded on workers.
            
            # Using min_sv_size=0 by default for robustness as per our testing
            sv_data <- read_sv_vcf(sv_path, min_sv_size = 0)
            cnv_data <- read_cnv_vcf(cnv_path)
            
            # Run Detection
            # Pass ellipses (...) arguments to detect_chromoanagenesis
            res <- detect_chromoanagenesis(sv_data, cnv_data, verbose = FALSE, ...)
            
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
