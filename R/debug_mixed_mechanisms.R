#' Debug function to inspect mixed mechanisms classification
#'
#' @param mixed_class Result from classify_mixed_mechanisms()
#' @export
debug_mixed_mechanisms <- function(mixed_class) {
    cat("\n=== DEBUGGING MIXED MECHANISMS CLASSIFICATION ===\n\n")

    # Check mechanism locations
    cat("1. MECHANISM LOCATIONS:\n")
    cat(sprintf("  - Rows: %d\n", nrow(mixed_class$mechanism_locations)))
    if (nrow(mixed_class$mechanism_locations) > 0) {
        cat("  - Column names:", paste(colnames(mixed_class$mechanism_locations), collapse = ", "), "\n")
        cat("  - Mechanisms found:", paste(unique(mixed_class$mechanism_locations$mechanism), collapse = ", "), "\n")
        cat("  - Sample:\n")
        print(head(mixed_class$mechanism_locations))
    } else {
        cat("  - EMPTY! No mechanism locations found\n")
    }

    # Check overlaps
    cat("\n2. OVERLAPS:\n")
    cat(sprintf("  - Number of overlaps: %d\n", mixed_class$overlaps$n_overlaps))
    if (mixed_class$overlaps$n_overlaps > 0) {
        cat("  - Sample:\n")
        print(head(mixed_class$overlaps$overlaps))
    }

    # Check chromosome classification
    cat("\n3. CHROMOSOME CLASSIFICATION:\n")
    cat(sprintf("  - Rows: %d\n", nrow(mixed_class$chromosome_classification)))
    if (nrow(mixed_class$chromosome_classification) > 0) {
        cat("  - Sample:\n")
        print(head(mixed_class$chromosome_classification))
    } else {
        cat("  - EMPTY!\n")
    }

    # Check sample classification
    cat("\n4. SAMPLE CLASSIFICATION:\n")
    cat(sprintf("  - Classification: %s\n", mixed_class$sample_classification$classification))
    cat(sprintf("  - Category: %s\n", mixed_class$sample_classification$category))
    cat(sprintf("  - N mechanisms: %d\n", mixed_class$sample_classification$n_mechanisms))

    # Check complexity
    cat("\n5. COMPLEXITY:\n")
    cat(sprintf("  - Level: %s\n", mixed_class$complexity$complexity_level))
    cat(sprintf("  - Score: %.3f\n", mixed_class$complexity$complexity_score))
    cat(sprintf("  - Components rows: %d\n", nrow(mixed_class$complexity$components)))

    # Check dominance
    cat("\n6. DOMINANCE:\n")
    cat(sprintf("  - Dominant mechanism: %s\n", mixed_class$dominance$dominant_mechanism))
    cat("  - Proportions:\n")
    print(mixed_class$dominance$mechanism_proportions)

    cat("\n=== END DEBUG ===\n\n")
    invisible(NULL)
}
