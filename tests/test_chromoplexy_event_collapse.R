library(OncoImplexus)
library(GenomicRanges)
library(IRanges)

sv <- SVs(
    chrom1 = c("1", "2", "3"),
    pos1 = c(1000000, 1100000, 1100000),
    chrom2 = c("2", "3", "1"),
    pos2 = c(1000000, 1000000, 1100000),
    strand1 = c("+", "+", "+"),
    strand2 = c("-", "-", "-"),
    SVtype = c("TRA", "TRA", "TRA"),
    sv_id = c("sv1", "sv2", "sv3")
)

genes <- GRanges(
    seqnames = c("1", "2", "3"),
    ranges = IRanges(
        start = c(999500, 1099500, 999500),
        end = c(1000500, 1100500, 1000500)
    )
)
mcols(genes)$symbol <- c("TP53", "GENE2", "GENE3")
mcols(genes)$gene_id <- c("ENSGTP53", "ENSG000002", "ENSG000003")

res <- detect_chromoplexy(
    SV.sample = sv,
    CNV.sample = NULL,
    use_statistical_testing = FALSE,
    gene_granges = genes,
    verbose = FALSE
)

stopifnot(!is.null(res$collapsed_events))
stopifnot(nrow(res$collapsed_events$event_summary) >= 1)
stopifnot(nrow(res$collapsed_events$chain_to_event) >= 1)
stopifnot(nrow(res$collapsed_events$event_breakpoints) >= 1)
stopifnot(all(c(
    "collapsed_event_id",
    "event_qc_score",
    "event_confidence",
    "n_unique_svs",
    "n_breakpoints"
) %in% colnames(res$collapsed_events$event_summary)))

score <- res$collapsed_events$event_summary$event_qc_score[1]
stopifnot(!is.na(score), score >= 0, score <= 1)

stopifnot(nrow(res$collapsed_events$gene_detail) >= 1)
stopifnot("TP53" %in% res$collapsed_events$gene_detail$symbol)
stopifnot(any(res$collapsed_events$gene_detail$is_driver))

annotated <- annotate_chromoplexy_events(
    chromoplexy_result = res,
    gene_granges = genes,
    breakpoint_padding = 1000
)
stopifnot(nrow(annotated$event_summary) == nrow(res$collapsed_events$event_summary))
stopifnot(nrow(annotated$gene_event_summary) >= 1)
