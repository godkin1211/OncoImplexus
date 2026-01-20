#' Comprehensive chromoanagenesis detection
#'
#' Performs integrated analysis to detect chromothripsis, chromoplexy, and
#' chromoanasynthesis from the same dataset, providing a complete view of complex
#' chromosomal rearrangements.
#'
#' @param SV.sample An instance of class SVs or data frame with SV data
#' @param CNV.sample An instance of class CNVsegs or data frame with CNV data
#' @param genome Reference genome ("hg19" or "hg38", default: "hg19")
#' @param detect_chromothripsis Detect chromothripsis events (default: TRUE)
#' @param detect_chromoplexy Detect chromoplexy events (default: TRUE)
#' @param detect_chromoanasynthesis Detect chromoanasynthesis events (default: TRUE)
#' @param min_chromothripsis_size Minimum SV cluster size for chromothripsis (default: 1)
#' @param min_chromoplexy_chromosomes Minimum number of chromosomes involved in a chromoplexy chain
#' @param min_chromoanasynthesis_tandem_dups Minimum number of tandem duplications in a chromoanasynthesis region
#' @param max_path_search Maximum number of paths to search during chromoplexy detection
#' @param max_neighbors Maximum number of neighbors to consider for each node in the translocation graph
#' @param max_region_size Maximum genomic size for a chromoanasynthesis candidate region
#' @param gene_granges Optional GRanges object containing gene models for driver impact annotation
#' @param verbose Logical, whether to print detailed progress messages
#' @return A chromoanagenesis object containing detection results for all mechanisms
#' @details
#' This function provides a comprehensive analysis of chromoanagenesis events
#' by detecting chromothripsis, chromoplexy, and chromoanasynthesis in the same sample.
#'
#' Chromothripsis characteristics:
#' - Localized to one or few chromosomes
#' - Many clustered breakpoints
#' - Oscillating copy number patterns
#' - Random fragment joins
#'
#' Chromoplexy characteristics:
#' - Involves multiple chromosomes (3-8)
#' - Chained translocations
#' - Minimal copy number changes
#' - May form cycles
#'
#' Chromoanasynthesis characteristics:
#' - Replication-based mechanism (FoSTeS/MMBIR)
#' - Gradual copy number increases
#' - Tandem duplications and inversions
#' - Localized to specific regions
#'
#' The function returns integrated results allowing comparison and
#' classification of all complex rearrangement patterns.
#'
#' @examples
#' \dontrun{
#' library(ShatterSeek)
#' data(DO17373)
#'
#' # Prepare data
#' SV_data <- SVs(chrom1=SV_DO17373$chrom1, ...)
#' CN_data <- CNVsegs(chrom=SCNA_DO17373$chromosome, ...)
#'
#' # Run comprehensive analysis
#' results <- detect_chromoanagenesis(SV_data, CN_data)
#'
#' # View results
#' print(results)
#' summary(results)
#'
#' # Access specific results
#' results$chromothripsis
#' results$chromoplexy
#' }
#'
#' @export
detect_chromoanagenesis <- function(SV.sample,
                                   CNV.sample,
                                   genome = "hg19",
                                   detect_chromothripsis = TRUE,
                                   detect_chromoplexy = TRUE,
                                   detect_chromoanasynthesis = TRUE,
                                   gene_granges = NULL, # New parameter
                                   min_chromothripsis_size = 1,
                                   min_chromoplexy_chromosomes = 3,
                                   min_chromoanasynthesis_tandem_dups = 3,
                                   max_path_search = 50,
                                   max_region_size = 10e6,
                                   max_neighbors = 3,
                                   verbose = TRUE) {

    if (verbose) {
        cat("\n")
        cat(rep("=", 70), "\n", sep = "")
        cat("     COMPREHENSIVE CHROMOANAGENESIS ANALYSIS\n")
        cat(rep("=", 70), "\n\n", sep = "")
    }

    results <- list()

    # 1. Data quality check
    if (verbose) cat("Step 1: Checking data quality...\n")
    quality <- check_data_quality(SV.sample, CNV.sample, verbose = FALSE)

    if (quality$has_issues && verbose) {
        cat("  WARNING: Data quality issues detected.\n")
        cat(sprintf("  Number of warnings: %d\n", length(quality$warnings)))
    } else if (verbose) {
        cat("  Data quality: OK\n")
    }

    results$quality_check <- quality

    # 2. Detect chromothripsis
    if (detect_chromothripsis) {
        if (verbose) {
            cat("\nStep 2: Detecting chromothripsis...\n")
        }

        chromothripsis_result <- detect_chromothripsis(
            SV.sample = SV.sample,
            seg.sample = CNV.sample,
            min.Size = min_chromothripsis_size,
            genome = genome
        )

        # Calculate scores and classifications
        if (verbose) cat("  Calculating confidence scores...\n")
        chromothripsis_classification <- classify_chromothripsis(chromothripsis_result)

        results$chromothripsis <- list(
            detection_output = chromothripsis_result,
            classification = chromothripsis_classification,
            n_high_confidence = sum(chromothripsis_classification$classification == "High confidence"),
            n_low_confidence = sum(chromothripsis_classification$classification == "Low confidence")
        )

        if (verbose) {
            cat(sprintf("  Found %d high confidence and %d low confidence chromothripsis events.\n",
                       results$chromothripsis$n_high_confidence,
                       results$chromothripsis$n_low_confidence))
        }
    }

    # 3. Detect chromoplexy
    if (detect_chromoplexy) {
        if (verbose) {
            cat("\nStep 3: Detecting chromoplexy...\n")
        }

        chromoplexy_result <- detect_chromoplexy(
            SV.sample = SV.sample,
            CNV.sample = CNV.sample,
            min_chromosomes = min_chromoplexy_chromosomes,
            max_path_search = max_path_search, # Pass new parameter
            max_neighbors = max_neighbors # Pass new parameter
        )

        results$chromoplexy <- chromoplexy_result

        if (verbose) {
            cat(sprintf("  Found %d likely and %d possible chromoplexy events.\n",
                       chromoplexy_result$likely_chromoplexy,
                       chromoplexy_result$possible_chromoplexy))
        }
    }

    # 4. Detect chromoanasynthesis
    if (detect_chromoanasynthesis) {
        if (verbose) {
            cat("\nStep 4: Detecting chromoanasynthesis...\n")
        }

        chromoanasynthesis_result <- detect_chromoanasynthesis(
            SV.sample = SV.sample,
            CNV.sample = CNV.sample,
            min_tandem_dups = min_chromoanasynthesis_tandem_dups,
            max_region_size = max_region_size # Pass new parameter
        )

        results$chromoanasynthesis <- chromoanasynthesis_result

        if (verbose) {
            cat(sprintf("  Found %d likely and %d possible chromoanasynthesis events.\n",
                       chromoanasynthesis_result$likely_chromoanasynthesis,
                       chromoanasynthesis_result$possible_chromoanasynthesis))
        }
    }

    # 4.5 Detect fusion genes
    if (!is.null(gene_granges)) {
        if (verbose) {
            cat("\nStep 4.5: Predicting fusion genes...\n")
        }
        
        # Robust call with existence check and try-catch
        if (exists("predict_fusion_genes")) {
            results$fusions <- tryCatch({
                predict_fusion_genes(
                    SV.sample = SV.sample,
                    gene_granges = gene_granges
                )
            }, error = function(e) {
                if (verbose) message("  Warning: Fusion prediction failed: ", e$message)
                NULL
            })
        } else {
            if (verbose) message("  Warning: predict_fusion_genes() not found. Skipping.")
            results$fusions <- NULL
        }
    }

    # 5. Integrated classification
    if (verbose) cat("\nStep 5: Integrated classification...\n")

    integrated_summary <- create_integrated_summary(results)
    results$integrated_summary <- integrated_summary

    if (verbose) {
        cat("\n")
        cat(rep("=", 70), "\n", sep = "")
        cat("ANALYSIS COMPLETE\n")
        cat(rep("=", 70), "\n\n", sep = "")
    }

    class(results) <- c("chromoanagenesis", "list")
    return(results)
}


#' Create integrated summary of chromoanagenesis results
#'
#' @param results Results list from detect_chromoanagenesis
#' @return Integrated summary data frame
#' @keywords internal
create_integrated_summary <- function(results) {

    summary <- list()

    # Chromothripsis summary
    if (!is.null(results$chromothripsis)) {
        summary$chromothripsis_high_confidence <- results$chromothripsis$n_high_confidence
        summary$chromothripsis_low_confidence <- results$chromothripsis$n_low_confidence

        ct_class <- results$chromothripsis$classification
        if (!is.null(ct_class) && all(c("classification", "chrom") %in% colnames(ct_class))) {
            high_conf_chroms <- ct_class[ct_class$classification == "High confidence", "chrom"]
            low_conf_chroms <- ct_class[ct_class$classification == "Low confidence", "chrom"]

            if (length(high_conf_chroms) > 0) {
                summary$chromothripsis_chromosomes <- paste(unique(high_conf_chroms), collapse = ", ")
                if (length(low_conf_chroms) > 0) {
                    summary$chromothripsis_chromosomes <- paste0(
                        summary$chromothripsis_chromosomes, " (high); ",
                        paste(unique(low_conf_chroms), collapse = ", "), " (low)"
                    )
                }
            } else if (length(low_conf_chroms) > 0) {
                summary$chromothripsis_chromosomes <- paste0(
                    paste(unique(low_conf_chroms), collapse = ", "), " (low confidence)"
                )
            } else {
                summary$chromothripsis_chromosomes <- "None"
            }
        } else {
            summary$chromothripsis_chromosomes <- "Analysis incomplete"
        }
    } else {
        summary$chromothripsis_high_confidence <- NA
        summary$chromothripsis_low_confidence <- NA
        summary$chromothripsis_chromosomes <- "Not analyzed"
    }

    # Chromoplexy summary
    if (!is.null(results$chromoplexy)) {
        summary$chromoplexy_likely <- results$chromoplexy$likely_chromoplexy
        summary$chromoplexy_possible <- results$chromoplexy$possible_chromoplexy
        summary$chromoplexy_chains <- results$chromoplexy$total_chains

        cp_sum <- results$chromoplexy$summary
        if (!is.null(cp_sum) && "classification" %in% colnames(cp_sum) && results$chromoplexy$likely_chromoplexy > 0) {
            likely_chains <- cp_sum[cp_sum$classification == "Likely chromoplexy", ]
            if (nrow(likely_chains) > 0) {
                chr_list <- unique(unlist(strsplit(as.character(likely_chains$chromosomes_involved), ",")))
                summary$chromoplexy_chromosomes <- paste(chr_list, collapse = ", ")
            } else {
                summary$chromoplexy_chromosomes <- "None"
            }
        } else {
            summary$chromoplexy_chromosomes <- "None"
        }
    } else {
        summary$chromoplexy_likely <- NA
        summary$chromoplexy_possible <- NA
        summary$chromoplexy_chains <- NA
        summary$chromoplexy_chromosomes <- "Not analyzed"
    }

    # Chromoanasynthesis summary
    if (!is.null(results$chromoanasynthesis)) {
        summary$chromoanasynthesis_likely <- results$chromoanasynthesis$likely_chromoanasynthesis
        summary$chromoanasynthesis_possible <- results$chromoanasynthesis$possible_chromoanasynthesis
        summary$chromoanasynthesis_regions <- results$chromoanasynthesis$total_regions

        if (results$chromoanasynthesis$likely_chromoanasynthesis > 0) {
            likely_regions <- results$chromoanasynthesis$summary[
                results$chromoanasynthesis$summary$classification == "Likely chromoanasynthesis",
            ]
            chr_list <- unique(likely_regions$chrom)
            summary$chromoanasynthesis_chromosomes <- paste(chr_list, collapse = ", ")
        } else {
            summary$chromoanasynthesis_chromosomes <- "None"
        }
    } else {
        summary$chromoanasynthesis_likely <- NA
        summary$chromoanasynthesis_possible <- NA
        summary$chromoanasynthesis_regions <- NA
        summary$chromoanasynthesis_chromosomes <- "Not analyzed"
    }

    # Overall assessment
    has_chromothripsis <- !is.null(results$chromothripsis) &&
                         (results$chromothripsis$n_high_confidence > 0 ||
                          results$chromothripsis$n_low_confidence > 0)
    has_chromoplexy <- !is.null(results$chromoplexy) &&
                      results$chromoplexy$likely_chromoplexy > 0
    has_chromoanasynthesis <- !is.null(results$chromoanasynthesis) &&
                          results$chromoanasynthesis$likely_chromoanasynthesis > 0

    # Create detailed classification
    mechanisms <- c()
    if (has_chromothripsis) mechanisms <- c(mechanisms, "chromothripsis")
    if (has_chromoplexy) mechanisms <- c(mechanisms, "chromoplexy")
    if (has_chromoanasynthesis) mechanisms <- c(mechanisms, "chromoanasynthesis")

    if (length(mechanisms) == 0) {
        summary$overall_classification <- "No chromoanagenesis detected"
    } else if (length(mechanisms) == 1) {
        summary$overall_classification <- paste(tools::toTitleCase(mechanisms[1]), "only")
    } else {
        summary$overall_classification <- paste("Multiple mechanisms:",
                                               paste(mechanisms, collapse = ", "))
    }

    return(as.data.frame(summary, stringsAsFactors = FALSE))
}


#' Print method for chromoanagenesis results
#'
#' @param x Chromoanagenesis result object
#' @param ... Additional arguments
#' @export
print.chromoanagenesis <- function(x, ...) {
    cat("\n")
    cat(rep("=", 70), "\n", sep = "")
    cat("         CHROMOANAGENESIS ANALYSIS SUMMARY\n")
    cat(rep("=", 70), "\n\n", sep = "")

    cat("Overall Classification:", x$integrated_summary$overall_classification, "\n\n")

    # Chromothripsis results
    if (!is.null(x$chromothripsis)) {
        cat("CHROMOTHRIPSIS:\n")
        cat(sprintf("  - High confidence: %d\n", x$chromothripsis$n_high_confidence))
        cat(sprintf("  - Low confidence:  %d\n", x$chromothripsis$n_low_confidence))
        cat(sprintf("  - Chromosomes:     %s\n",
                   x$integrated_summary$chromothripsis_chromosomes))
        cat("\n")
    }

    # Chromoplexy results
    if (!is.null(x$chromoplexy)) {
        cat("CHROMOPLEXY:\n")
        cat(sprintf("  - Likely events:   %d\n", x$chromoplexy$likely_chromoplexy))
        cat(sprintf("  - Possible events: %d\n", x$chromoplexy$possible_chromoplexy))
        cat(sprintf("  - Total chains:    %d\n", x$chromoplexy$total_chains))
        cat(sprintf("  - Chromosomes:     %s\n",
                   x$integrated_summary$chromoplexy_chromosomes))
        cat("\n")
    }

    # Chromoanasynthesis results
    if (!is.null(x$chromoanasynthesis)) {
        cat("CHROMOANASYNTHESIS:\n")
        cat(sprintf("  - Likely events:   %d\n", x$chromoanasynthesis$likely_chromoanasynthesis))
        cat(sprintf("  - Possible events: %d\n", x$chromoanasynthesis$possible_chromoanasynthesis))
        cat(sprintf("  - Total regions:   %d\n", x$chromoanasynthesis$total_regions))
        cat(sprintf("  - Chromosomes:     %s\n",
                   x$integrated_summary$chromoanasynthesis_chromosomes))
        cat("\n")
    }

    # Data quality
    if (!is.null(x$quality_check) && x$quality_check$has_issues) {
        cat("DATA QUALITY WARNINGS:\n")
        for (w in x$quality_check$warnings) {
            cat(sprintf("  - %s\n", w))
        }
        cat("\n")
    }

    cat(rep("=", 70), "\n", sep = "")
    cat("\n")

    invisible(x)
}


#' Summary method for chromoanagenesis results
#'
#' @param object Chromoanagenesis result object
#' @param ... Additional arguments
#' @export
summary.chromoanagenesis <- function(object, ...) {

    cat("\nIntegrated Summary Table:\n")
    print(object$integrated_summary)

    cat("\n")

    if (!is.null(object$chromothripsis)) {
        cat("Chromothripsis Details:\n")
        print(object$chromothripsis$classification[, c("chrom", "classification",
                                                       "confidence_score",
                                                       "clusterSize")])
        cat("\n")
    }

    if (!is.null(object$chromoplexy) && object$chromoplexy$total_chains > 0) {
        cat("Chromoplexy Details:\n")
        print(object$chromoplexy$summary[, c("chain_id", "n_chromosomes",
                                            "n_translocations", "classification")])
        cat("\n")
    }

    if (!is.null(object$chromoanasynthesis) && object$chromoanasynthesis$total_regions > 0) {
        cat("Chromoanasynthesis Details:\n")
        # Identify available columns to show
        sum_cols <- colnames(object$chromoanasynthesis$summary)
        show_cols <- intersect(c("region_id", "chrom", "start", "end", "classification", "combined_score"), sum_cols)
        print(object$chromoanasynthesis$summary[, show_cols])
        cat("\n")
    }

    invisible(object)
}
