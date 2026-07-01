sampleId <- commandArgs(TRUE)[1]
inputDir <- commandArgs(TRUE)[2]

library(OncoImplexus)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

# Set paths
SV_VCF <- paste0(inputDir,"/",sampleId,".sv.vcf.gz")
CNV_VCF <- paste0(inputDir,"/",sampleId,".cnv.vcf.gz")
GENE_RDS <- system.file("extdata", "hg38_genes.rds", package = "OncoImplexus")
if (GENE_RDS == "" && file.exists("inst/extdata/hg38_genes.rds")) {
  GENE_RDS <- "inst/extdata/hg38_genes.rds"
}

cat("=======================================================\n")
cat("Starting Analysis for",sampleId,"\n")
cat(sprintf("SV VCF: %s\n", SV_VCF))
cat(sprintf("CNV VCF: %s\n", ifelse(file.exists(CNV_VCF), CNV_VCF, "not provided")))
cat("=======================================================\n")

# 1. Load Data
if (!file.exists(SV_VCF)) {
  stop("SV VCF file not found: ", SV_VCF)
}

cat("[Step 1] Loading SVs")
if (file.exists(CNV_VCF)) {
  cat(" and CNVs")
}
cat("...\n")
sv_data <- read_sv_vcf(SV_VCF, genome = "hg38")
cnv_data <- if (file.exists(CNV_VCF)) {
  read_cnv_vcf(CNV_VCF)
} else {
  cat("  CNV file not found. Running SV-only chromoplexy analysis.\n")
  NULL
}
if (GENE_RDS == "" || !file.exists(GENE_RDS)) {
  stop("Gene annotation file not found: hg38_genes.rds")
}
gene_data <- readRDS(GENE_RDS)

# 2. Run Detection
cat("\n[Step 2] Running Detection...\n")
res <- detect_chromoanagenesis(
  SV.sample = sv_data,
  CNV.sample = cnv_data,
  genome = "hg38",
  gene_granges = gene_data
)

# 3. Print Summary
cat("\n[Step 3] Summary of Results...\n")
print(summary(res))

# 4. Generate Report
cat("\n[Step 4] Generating Report...\n")
output_file <- paste0(sampleId,"_Test_Report.html")
generate_interactive_report(
  result = res,
  SV.sample = sv_data,
  CNV.sample = cnv_data,
  output_file = output_file,
  output_dir = "reports",
  sample_name = sampleId,
  genome = "hg38",
  gene_granges = gene_data
)

cat(sprintf("\nSUCCESS: Report generated at reports/%s\n", output_file))
