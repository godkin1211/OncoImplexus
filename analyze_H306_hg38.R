library(OncoImplexus)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

sample_id <- commandArgs(TRUE)[1]
data_dir <- "WGS_results/my_data"
sv_file <- file.path(data_dir, paste0(sample_id, ".sv.vcf.gz"))
cnv_file <- file.path(data_dir, paste0(sample_id, ".cnv.vcf.gz"))

if (!file.exists(sv_file)) {
    stop("SV VCF file not found for ", sample_id)
}

# 1. Load Genome Info & Genes
message(">>> Loading hg38 gene models...")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes_gr <- readRDS(system.file("extdata", "hg38_genes.rds", package = "OncoImplexus"))

# 2. Read Data with min_sv_size=1000
message(">>> Reading VCF data (min_sv_size=1000)...")
sv_data <- read_sv_vcf(sv_file, min_sv_size = 1000)
cnv_data <- if (file.exists(cnv_file)) {
    read_cnv_vcf(cnv_file)
} else {
    message(">>> CNV VCF not found. Running SV-only chromoplexy analysis.")
    NULL
}

# 3. Run Analysis
message(">>> Running detect_chromoanagenesis...")
results <- detect_chromoanagenesis(
    SV.sample = sv_data, 
    CNV.sample = cnv_data, 
    gene_granges = genes_gr,
    genome = "hg38",
    verbose = TRUE
)

# 4. Generate Report
message(">>> Generating report...")
output_file <- paste0(sample_id, "_hg38_Report.html")
output_dir <- "reports/sample_reports"

report_res <- generate_interactive_report(
    result = results,
    SV.sample = sv_data,
    CNV.sample = cnv_data,
    output_file = output_file,
    output_dir = output_dir,
    sample_name = sample_id,
    genome = "hg38",
    gene_granges = genes_gr,
    txdb = txdb
)

message("Report successfully generated at: ", report_res$report_path)
