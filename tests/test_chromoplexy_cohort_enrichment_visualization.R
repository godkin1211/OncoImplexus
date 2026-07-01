library(OncoImplexus)
library(GenomicRanges)
library(IRanges)

make_result <- function(offset = 0) {
    sv <- SVs(
        chrom1 = c("1", "2", "3"),
        pos1 = c(1000000, 1100000, 1100000) + offset,
        chrom2 = c("2", "3", "1"),
        pos2 = c(1000000, 1000000, 1100000) + offset,
        strand1 = c("+", "+", "+"),
        strand2 = c("-", "-", "-"),
        SVtype = c("TRA", "TRA", "TRA"),
        sv_id = paste0("sv", seq_len(3), "_", offset)
    )

    genes <- GRanges(
        seqnames = c("1", "2", "3"),
        ranges = IRanges(
            start = c(999500, 1099500, 999500) + offset,
            end = c(1000500, 1100500, 1000500) + offset
        )
    )
    mcols(genes)$symbol <- c("TP53", "EGFR", "KRAS")
    mcols(genes)$gene_id <- c("ENSGTP53", "ENSGEGFR", "ENSGKRAS")

    detect_chromoplexy(
        SV.sample = sv,
        CNV.sample = NULL,
        use_statistical_testing = FALSE,
        gene_granges = genes,
        verbose = FALSE
    )
}

res1 <- make_result(0)
res2 <- make_result(1000)

cohort <- summarize_chromoplexy_cohort_events(list(S1 = res1, S2 = res2))

stopifnot(nrow(cohort$sample_summary) == 2)
stopifnot(nrow(cohort$event_summary) >= 2)
stopifnot(all(c(
    "recurrence_score",
    "event_cohort_priority_score"
) %in% colnames(cohort$event_summary)))
stopifnot("TP53" %in% cohort$gene_summary$symbol)
stopifnot(max(cohort$gene_summary$n_samples) >= 2)
stopifnot(nrow(cohort$breakpoint_region_summary) > 0)

genes <- extract_chromoplexy_genes(cohort)
stopifnot(all(c("TP53", "EGFR", "KRAS") %in% genes))

p <- plot_collapsed_chromoplexy_event(res1, event_id = "CE001", sample_name = "S1")
stopifnot(inherits(p, "ggplot"))

if (requireNamespace("clusterProfiler", quietly = TRUE) &&
    requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    enrich <- run_chromoplexy_enrichment(
        genes,
        run_go = FALSE,
        run_kegg = FALSE,
        min_genes = 1
    )
    stopifnot(nrow(enrich$id_mapping) >= 1)
    stopifnot(all(c("SYMBOL", "ENTREZID") %in% colnames(enrich$id_mapping)))
}
