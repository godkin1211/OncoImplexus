# OncoImplexus: Integrated Detection of Complex Chromosomal Rearrangements in Cancer

**OncoImplexus** is a comprehensive R package designed for the integrated detection, classification, and biological interpretation of **chromoanagenesis** events (Chromothripsis, Chromoplexy, and Chromoanasynthesis) in cancer genomes.

By unifying structural variation (SV) topology with copy number (CN) oscillations and DNA repair signatures, OncoImplexus provides a publication-quality landscape of genomic catastrophes across diverse cancer types.

---

## 🔬 Core Scientific Value

### 1. Multi-Mechanism Integration
Unlike traditional tools that often focus on a single phenomenon, OncoImplexus simultaneously identifies and differentiates:
- **Chromothripsis (CT)**: Localized shattering and random rejoining (NHEJ-driven).
- **Chromoplexy (CP)**: Inter-chromosomal translocation chains with copy number stability.
- **Chromoanasynthesis (CS)**: Replication-based rearrangements (MMBIR-driven) with characteristic CN gradients.

### 2. Superior Sensitivity & Robustness
- **Corrected Statistical Models**: Addresses limitations in random-join tests and coordinate biases found in legacy frameworks.
- **Data Self-Healing**: Adaptive parsing of SV/Strand conflicts ensures compatibility with diverse clinical VCFs (e.g., DRAGEN, Manta).
- **Optimized Signal-to-Noise**: employs a standardized **1kb minimum SV size filter** to remove background noise that confounds statistical tests.
- **Benchmark Performance**: Demonstrated significantly improved detection prevalence in benchmarks (e.g., PCAWG Bladder-TCC) compared to legacy tools.

### 3. Biological & Clinical Depth
- **DNA Repair Fingerprinting**: Quantifies NHEJ, MMBIR, and MMEJ proportions at breakpoints.
- **Driver Impact Ranking**: Automatically maps catastrophes to core oncogenes and TSGs.
- **Collapsed Chromoplexy Events**: Collapses redundant chain enumerations into event-level components with QC scores, breakpoint tables, and optional gene impact annotation.
- **SV-only Long-read Mode**: Supports chromoplexy analysis from SV VCF alone when CNV VCF is unavailable, including ONT/Sniffles2 BND ALT parsing.
- **Clinical Correlation**: Supports integrated survival analysis and Whole Genome Duplication (WGD) association tests.

---

## � Installation

You can install the development version of OncoImplexus from GitHub:

```r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("godkin1211/OncoImplexus")
```

**Prerequisites**: R >= 4.0.0 and standard Bioconductor packages (`GenomicRanges`, `VariantAnnotation`, etc.).

---

## 🚀 Quick Start Tutorial

OncoImplexus includes a bundled example dataset (`DO17373`) to help you get started quickly. This sample exhibits chromothripsis on chromosomes 2, 21, and X.

### 1. Load Example Data & Run Analysis

```r
library(OncoImplexus)

# === Load Bundled Example Data ===
load(system.file("extdata", "DO17373.RData", package = "OncoImplexus"))

# === Create SVs Object ===
sv_obj <- SVs(
  chrom1 = SV_DO17373$chrom1,
  pos1 = SV_DO17373$start1,
  chrom2 = SV_DO17373$chrom2,
  pos2 = SV_DO17373$start2,
  strand1 = SV_DO17373$strand1,
  strand2 = SV_DO17373$strand2,
  SVtype = SV_DO17373$svclass,
  sv_id = SV_DO17373$sv_id
)

# === Create CNVsegs Object ===
cnv_obj <- CNVsegs(
  chrom = as.character(SCNA_DO17373$chromosome),
  start = SCNA_DO17373$start,
  end = SCNA_DO17373$end,
  total_cn = SCNA_DO17373$total_cn
)

# === Run Integrated Detection ===
results <- detect_chromoanagenesis(sv_obj, cnv_obj, genome = "hg19")
```

### 2. Generate Interactive Report

Visualization is a key component of OncoImplexus. You can generate a comprehensive HTML report as follows:

```r
# === Load Gene Annotations ===
gene_granges <- readRDS(
  system.file("extdata", "hg19_genes.rds", package = "OncoImplexus")
)

# === Generate HTML Report ===
generate_interactive_report(
  result = results,
  SV.sample = sv_obj,
  CNV.sample = cnv_obj,
  output_file = "DO17373_Report.html",
  output_dir = "reports",
  sample_name = "DO17373",
  genome = "hg19",
  gene_granges = gene_granges
)
```

### 3. Working with VCF Files

OncoImplexus supports direct import from VCF files (Manta, Delly, DRAGEN, etc.):

```r
# Load SVs and CNVs directly from VCF
sv_data <- read_sv_vcf("tumor.sv.vcf.gz")
cnv_data <- read_cnv_vcf("tumor.cnv.vcf.gz")

# Run integrated detection
results <- detect_chromoanagenesis(sv_data, cnv_data, genome = "hg38")
```

If only an SV VCF is available, OncoImplexus runs SV-only chromoplexy detection
and skips CNV-dependent chromothripsis/chromoanasynthesis modules:

```r
sv_data <- read_sv_vcf("tumor.sv.vcf.gz", genome = "hg38")
gene_granges <- readRDS(
  system.file("extdata", "hg38_genes.rds", package = "OncoImplexus")
)

results <- detect_chromoanagenesis(
  SV.sample = sv_data,
  CNV.sample = NULL,
  genome = "hg38",
  gene_granges = gene_granges
)

# Chain-level calls
results$chromoplexy$summary

# Event-level collapsed calls
results$chromoplexy$collapsed_events$event_summary
results$chromoplexy$collapsed_events$gene_event_summary
```

For an existing chromoplexy result, event collapse and gene annotation can also
be run directly:

```r
events <- collapse_chromoplexy_chains(
  chromoplexy_result = results$chromoplexy,
  gene_granges = gene_granges,
  breakpoint_padding = 1000
)
```

For a more detailed walkthrough, please refer to the tutorial located at `inst/tutorial/USER_GUIDE_END_TO_END.Rmd`.

---

## 📂 Repository Structure

The package is organized as follows:

*   **`R/`**: Core R functions and algorithms for detection and visualization.
*   **`inst/`**:
    *   `tutorial/`: Detailed user guides and tutorials (e.g., `USER_GUIDE_END_TO_END.Rmd`).
    *   `extdata/`: Bundled example data and gene annotations.
*   **`ForPaper/`**: Contains scripts and analyses related to the manuscript.
    *   `develop/scripts/`: specialized analysis scripts such as `analyze_pcawg_cohorts.R`, `analyze_utuc_hg38.R`, and visualization tools used for the paper's figures.
*   **`tests/`**: Unit tests and validation scripts.

---

## 📧 Author & Contact

**Author**: Chia-Chun Chiu  
**Maintainer Email**: n28111021@gs.ncku.edu.tw
