library(OncoImplexus)

vcf_file <- tempfile(fileext = ".vcf")
writeLines(c(
    "##fileformat=VCFv4.2",
    "##source=Sniffles2_2.6.3",
    "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"SV type\">",
    "##contig=<ID=chr1,length=248956422>",
    "##contig=<ID=chr2,length=242193529>",
    "##contig=<ID=chr3,length=198295559>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t1000000\tbnd1\tN\tN[CHR2:1000000[\t60\tPASS\tSVTYPE=BND",
    "chr2\t1000000\tbnd2\tN\tN[CHR3:1000000[\t60\tPASS\tSVTYPE=BND",
    "chr3\t1000000\tbnd3\tN\tN[CHR1:1000000[\t60\tPASS\tSVTYPE=BND"
), vcf_file)

sv <- suppressWarnings(read_sv_vcf(vcf_file, genome = "hg38", min_sv_size = 0))
sv_df <- as(sv, "data.frame")

stopifnot(nrow(sv_df) == 3)
stopifnot(all(sv_df$SVtype == "TRA"))
stopifnot(sum(as.character(sv_df$chrom1) != as.character(sv_df$chrom2)) == 3)

res <- suppressWarnings(detect_chromoanagenesis(
    SV.sample = sv,
    CNV.sample = NULL,
    genome = "hg38",
    verbose = FALSE
))

stopifnot(identical(res$analysis_mode, "SV-only chromoplexy"))
stopifnot(!is.null(res$chromoplexy))
if (res$chromoplexy$total_chains > 0) {
    stopifnot(all(res$chromoplexy$summary$evidence_mode == "SV-only"))
}
