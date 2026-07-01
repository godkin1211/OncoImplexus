library(OncoImplexus)

sv <- SVs(
    chrom1 = c("1", "2", "3"),
    pos1 = c(1000000, 1200000, 1200000),
    chrom2 = c("2", "3", "1"),
    pos2 = c(1000000, 1000000, 1200000),
    strand1 = c("+", "+", "+"),
    strand2 = c("-", "-", "-"),
    SVtype = c("TRA", "TRA", "TRA"),
    sv_id = c("sv1", "sv2", "sv3")
)

sv_only <- detect_chromoanagenesis(
    SV.sample = sv,
    CNV.sample = NULL,
    genome = "hg19",
    verbose = FALSE
)

stopifnot(identical(sv_only$analysis_mode, "SV-only chromoplexy"))
stopifnot(is.null(sv_only$chromothripsis))
stopifnot(is.null(sv_only$chromoanasynthesis))
stopifnot(!is.null(sv_only$chromoplexy))
stopifnot(sv_only$chromoplexy$total_chains >= 1)
stopifnot(all(sv_only$chromoplexy$summary$evidence_mode == "SV-only"))

cnv <- CNVsegs(
    chrom = rep(c("1", "2", "3"), each = 4),
    start = rep(c(1, 1000001, 2000001, 3000001), 3),
    end = rep(c(1000000, 2000000, 3000000, 4000000), 3),
    total_cn = rep(2, 12)
)

full <- detect_chromoanagenesis(
    SV.sample = sv,
    CNV.sample = cnv,
    genome = "hg19",
    verbose = FALSE
)

stopifnot(identical(full$analysis_mode, "SV+CNV integrated"))
stopifnot(!is.null(full$chromothripsis))
stopifnot(!is.null(full$chromoanasynthesis))
stopifnot(!is.null(full$chromoplexy))
stopifnot(all(full$chromoplexy$summary$evidence_mode == "SV+CNV"))
