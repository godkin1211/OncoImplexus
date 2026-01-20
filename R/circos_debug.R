#' Debug circos visualization
#'
#' Helper function to debug circos plot issues
#'
#' @param sv_data SV data (data frame)
#' @export
debug_circos_svs <- function(sv_data) {

    if (is(sv_data, "SVs")) {
        sv_data <- as(sv_data, "data.frame")
    }

    cat("\n=== Circos SV Debug Information ===\n\n")

    # Check SV types
    cat("SV Types:\n")
    print(table(sv_data$SVtype))
    cat("\n")

    # Check chromosomes involved
    cat("Chromosomes in chrom1:\n")
    print(sort(unique(sv_data$chrom1)))
    cat("\n")

    cat("Chromosomes in chrom2:\n")
    print(sort(unique(sv_data$chrom2)))
    cat("\n")

    # Check for interchromosomal SVs
    inter_chr <- sv_data[sv_data$chrom1 != sv_data$chrom2, ]
    cat(sprintf("Total SVs: %d\n", nrow(sv_data)))
    cat(sprintf("Intrachromosomal SVs: %d\n", sum(sv_data$chrom1 == sv_data$chrom2)))
    cat(sprintf("Interchromosomal SVs: %d\n", nrow(inter_chr)))
    cat("\n")

    if (nrow(inter_chr) > 0) {
        cat("Interchromosomal SV types:\n")
        print(table(inter_chr$SVtype))
        cat("\n")

        cat("Interchromosomal SV chromosome pairs (first 20):\n")
        inter_pairs <- data.frame(
            chrom1 = inter_chr$chrom1,
            chrom2 = inter_chr$chrom2,
            type = inter_chr$SVtype
        )
        print(head(inter_pairs, 20))
    }

    invisible(NULL)
}
