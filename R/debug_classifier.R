#' Debug function to inspect chromoanagenesis results
#'
#' @param chromoanagenesis_result Result from detect_chromoanagenesis()
#' @export
debug_chromoanagenesis_result <- function(chromoanagenesis_result) {
    cat("\n=== DEBUGGING CHROMOANAGENESIS RESULT ===\n\n")

    # Check chromothripsis
    cat("1. CHROMOTHRIPSIS:\n")
    if (!is.null(chromoanagenesis_result$chromothripsis)) {
        cat("  - Not NULL\n")
        if (!is.null(chromoanagenesis_result$chromothripsis$classification)) {
            cat(sprintf("  - Classification rows: %d\n",
                       nrow(chromoanagenesis_result$chromothripsis$classification)))
            if (nrow(chromoanagenesis_result$chromothripsis$classification) > 0) {
                cat("  - Column names:", paste(colnames(chromoanagenesis_result$chromothripsis$classification), collapse = ", "), "\n")
                cat("  - First row:\n")
                print(chromoanagenesis_result$chromothripsis$classification[1, ])
            }
        } else {
            cat("  - Classification is NULL\n")
        }
    } else {
        cat("  - NULL\n")
    }

    # Check chromoplexy
    cat("\n2. CHROMOPLEXY:\n")
    if (!is.null(chromoanagenesis_result$chromoplexy)) {
        cat("  - Not NULL\n")
        if (!is.null(chromoanagenesis_result$chromoplexy$summary)) {
            cat(sprintf("  - Summary rows: %d\n",
                       nrow(chromoanagenesis_result$chromoplexy$summary)))
            if (nrow(chromoanagenesis_result$chromoplexy$summary) > 0) {
                cat("  - Column names:", paste(colnames(chromoanagenesis_result$chromoplexy$summary), collapse = ", "), "\n")
                cat("  - First row:\n")
                print(chromoanagenesis_result$chromoplexy$summary[1, ])
                cat("  - Classification values:\n")
                print(table(chromoanagenesis_result$chromoplexy$summary$classification))
            }
        } else {
            cat("  - Summary is NULL\n")
        }
    } else {
        cat("  - NULL\n")
    }

    # Check chromoanasynthesis
    cat("\n3. CHROMOANASYNTHESIS:\n")
    if (!is.null(chromoanagenesis_result$chromoanasynthesis)) {
        cat("  - Not NULL\n")
        if (!is.null(chromoanagenesis_result$chromoanasynthesis$summary)) {
            cat(sprintf("  - Summary rows: %d\n",
                       nrow(chromoanagenesis_result$chromoanasynthesis$summary)))
            if (nrow(chromoanagenesis_result$chromoanasynthesis$summary) > 0) {
                cat("  - Column names:", paste(colnames(chromoanagenesis_result$chromoanasynthesis$summary), collapse = ", "), "\n")
                cat("  - First row:\n")
                print(chromoanagenesis_result$chromoanasynthesis$summary[1, ])
                cat("  - Classification values:\n")
                print(table(chromoanagenesis_result$chromoanasynthesis$summary$classification))
            }
        } else {
            cat("  - Summary is NULL\n")
        }
    } else {
        cat("  - NULL\n")
    }

    cat("\n=== END DEBUG ===\n\n")
    invisible(NULL)
}
