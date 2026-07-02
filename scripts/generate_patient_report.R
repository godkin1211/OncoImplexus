#!/usr/bin/env Rscript

usage <- function(exit_status = 0) {
    cat("
Generate a complete single-patient OncoImplexus HTML report.

Required:
  --sample-id ID          Sample/patient ID used in output filenames
  --sv-vcf FILE          Structural variant VCF (.vcf or .vcf.gz)

Optional:
  --cnv-vcf FILE         Copy-number VCF (.vcf or .vcf.gz). If omitted, SV-only chromoplexy is run
  --genome hg19|hg38    Reference genome [default: hg38]
  --out-dir DIR         Output directory [default: reports/<sample-id>]
  --gene-rds FILE       Gene annotation GRanges RDS. Defaults to bundled <genome>_genes.rds
  --min-sv-size N       Minimum intrachromosomal SV size for VCF import [default: 0]
  --max-path-search N   Chromoplexy max path search [default: 50]
  --max-neighbors N     Chromoplexy graph neighbors per node [default: 3]
  --prefix TEXT         Output filename prefix [default: sample ID]
  --enrichment          Run GO BP and KEGG enrichment if clusterProfiler is installed
  --no-tables           Do not export TSV summary tables
  --quiet               Suppress verbose detection messages
  --help                Show this help

Outputs:
  <prefix>_Report.html
  <prefix>_result.rds
  <prefix>_integrated_summary.tsv
  <prefix>_chromoplexy_chain_summary.tsv
  <prefix>_collapsed_chromoplexy_events.tsv
  <prefix>_chromoplexy_gene_summary.tsv
  <prefix>_chromoplexy_breakpoints.tsv
  <prefix>_GO_BP_enrichment.tsv and <prefix>_KEGG_enrichment.tsv if --enrichment succeeds

Examples:
  Rscript scripts/generate_patient_report.R \\
    --sample-id Patient01 \\
    --sv-vcf input/Patient01.sv.vcf.gz \\
    --cnv-vcf input/Patient01.cnv.vcf.gz \\
    --genome hg38 \\
    --out-dir reports/Patient01

  Rscript scripts/generate_patient_report.R \\
    --sample-id ONT01 \\
    --sv-vcf bam_pass.wf_sv.vcf.gz \\
    --genome hg38 \\
    --enrichment
")
    quit(status = exit_status, save = "no")
}

parse_args <- function(args) {
    if (length(args) == 0 || any(args %in% c("--help", "-h"))) {
        usage(0)
    }

    opts <- list(
        cnv_vcf = NULL,
        genome = "hg38",
        out_dir = NULL,
        gene_rds = NULL,
        min_sv_size = 0,
        max_path_search = 50,
        max_neighbors = 3,
        prefix = NULL,
        enrichment = FALSE,
        export_tables = TRUE,
        verbose = TRUE
    )

    i <- 1
    while (i <= length(args)) {
        key <- args[[i]]
        needs_value <- function() {
            if (i == length(args) || startsWith(args[[i + 1]], "--")) {
                stop("Missing value for ", key, call. = FALSE)
            }
            args[[i + 1]]
        }

        if (key == "--sample-id") {
            opts$sample_id <- needs_value()
            i <- i + 2
        } else if (key == "--sv-vcf") {
            opts$sv_vcf <- needs_value()
            i <- i + 2
        } else if (key == "--cnv-vcf") {
            opts$cnv_vcf <- needs_value()
            i <- i + 2
        } else if (key == "--genome") {
            opts$genome <- needs_value()
            i <- i + 2
        } else if (key == "--out-dir") {
            opts$out_dir <- needs_value()
            i <- i + 2
        } else if (key == "--gene-rds") {
            opts$gene_rds <- needs_value()
            i <- i + 2
        } else if (key == "--min-sv-size") {
            opts$min_sv_size <- as.numeric(needs_value())
            i <- i + 2
        } else if (key == "--max-path-search") {
            opts$max_path_search <- as.integer(needs_value())
            i <- i + 2
        } else if (key == "--max-neighbors") {
            opts$max_neighbors <- as.integer(needs_value())
            i <- i + 2
        } else if (key == "--prefix") {
            opts$prefix <- needs_value()
            i <- i + 2
        } else if (key == "--enrichment") {
            opts$enrichment <- TRUE
            i <- i + 1
        } else if (key == "--no-tables") {
            opts$export_tables <- FALSE
            i <- i + 1
        } else if (key == "--quiet") {
            opts$verbose <- FALSE
            i <- i + 1
        } else {
            stop("Unknown argument: ", key, call. = FALSE)
        }
    }

    if (is.null(opts$sample_id) || !nzchar(opts$sample_id)) {
        stop("--sample-id is required", call. = FALSE)
    }
    if (is.null(opts$sv_vcf) || !nzchar(opts$sv_vcf)) {
        stop("--sv-vcf is required", call. = FALSE)
    }
    if (!opts$genome %in% c("hg19", "hg38")) {
        stop("--genome must be hg19 or hg38", call. = FALSE)
    }
    if (is.null(opts$out_dir)) {
        opts$out_dir <- file.path("reports", opts$sample_id)
    }
    if (is.null(opts$prefix)) {
        opts$prefix <- opts$sample_id
    }

    opts
}

configure_library_path <- function() {
    cmd <- commandArgs(trailingOnly = FALSE)
    file_arg <- cmd[startsWith(cmd, "--file=")]
    if (length(file_arg) == 0) {
        return(invisible(NULL))
    }

    script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)
    repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
    local_lib <- file.path(repo_root, ".r-lib")
    if (dir.exists(local_lib)) {
        .libPaths(unique(c(local_lib, .libPaths())))
    }

    invisible(NULL)
}

check_package_version <- function(min_version = "1.3.0") {
    pkg_version <- as.character(utils::packageVersion("OncoImplexus"))
    if (utils::compareVersion(pkg_version, min_version) < 0) {
        stop(
            "OncoImplexus >= ", min_version, " is required for the complete report workflow. ",
            "Installed version is ", pkg_version, " at ", find.package("OncoImplexus"), ". ",
            "Install the latest package or run from a checkout with a repo-local .r-lib.",
            call. = FALSE
        )
    }
    invisible(pkg_version)
}

load_gene_granges <- function(genome, gene_rds = NULL) {
    if (!is.null(gene_rds)) {
        if (!file.exists(gene_rds)) {
            stop("Gene annotation RDS not found: ", gene_rds, call. = FALSE)
        }
        return(readRDS(gene_rds))
    }

    gene_file <- system.file("extdata", paste0(genome, "_genes.rds"), package = "OncoImplexus")
    if (gene_file == "" && file.exists(file.path("inst", "extdata", paste0(genome, "_genes.rds")))) {
        gene_file <- file.path("inst", "extdata", paste0(genome, "_genes.rds"))
    }
    if (gene_file == "" || !file.exists(gene_file)) {
        stop("Bundled gene annotation not found for ", genome,
             ". Use --gene-rds to provide a GRanges RDS.", call. = FALSE)
    }

    readRDS(gene_file)
}

write_tsv <- function(x, path) {
    if (is.null(x) || !is.data.frame(x)) {
        return(invisible(FALSE))
    }
    utils::write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE)
    invisible(TRUE)
}

safe_export_tables <- function(result, output_dir, prefix) {
    write_tsv(result$integrated_summary,
              file.path(output_dir, paste0(prefix, "_integrated_summary.tsv")))

    if (!is.null(result$chromoplexy)) {
        write_tsv(result$chromoplexy$summary,
                  file.path(output_dir, paste0(prefix, "_chromoplexy_chain_summary.tsv")))

        ce <- result$chromoplexy$collapsed_events
        if (!is.null(ce)) {
            write_tsv(ce$event_summary,
                      file.path(output_dir, paste0(prefix, "_collapsed_chromoplexy_events.tsv")))
            write_tsv(ce$chain_to_event,
                      file.path(output_dir, paste0(prefix, "_chromoplexy_chain_to_event.tsv")))
            write_tsv(ce$event_breakpoints,
                      file.path(output_dir, paste0(prefix, "_chromoplexy_breakpoints.tsv")))
            write_tsv(ce$gene_detail,
                      file.path(output_dir, paste0(prefix, "_chromoplexy_gene_detail.tsv")))
            write_tsv(ce$gene_event_summary,
                      file.path(output_dir, paste0(prefix, "_chromoplexy_gene_summary.tsv")))
        }
    }

    if (!is.null(result$fusions) && is.data.frame(result$fusions)) {
        write_tsv(result$fusions,
                  file.path(output_dir, paste0(prefix, "_fusion_candidates.tsv")))
    }
}

run_enrichment_if_requested <- function(result, output_dir, prefix) {
    if (!requireNamespace("clusterProfiler", quietly = TRUE) ||
        !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
        message("  Enrichment skipped: clusterProfiler and/or org.Hs.eg.db not installed.")
        return(invisible(NULL))
    }

    genes <- tryCatch(
        extract_chromoplexy_genes(result),
        error = function(e) character(0)
    )
    if (length(genes) < 3) {
        message("  Enrichment skipped: fewer than 3 chromoplexy-affected genes.")
        return(invisible(NULL))
    }

    enrich_dir <- file.path(output_dir, "enrichment")
    message("  Running GO BP / KEGG enrichment on ", length(genes), " genes...")
    tryCatch({
        enrich <- run_chromoplexy_enrichment(
            result,
            organism_db = "org.Hs.eg.db",
            kegg_organism = "hsa",
            output_dir = enrich_dir,
            prefix = prefix
        )
        write_tsv(summarize_chromoplexy_pathways(enrich, top_n = 20),
                  file.path(enrich_dir, paste0(prefix, "_top_pathways.tsv")))
    }, error = function(e) {
        message("  Enrichment failed: ", conditionMessage(e))
    })
}

main <- function() {
    opts <- parse_args(commandArgs(trailingOnly = TRUE))

    configure_library_path()

    suppressPackageStartupMessages({
        library(OncoImplexus)
    })
    check_package_version()

    if (!file.exists(opts$sv_vcf)) {
        stop("SV VCF not found: ", opts$sv_vcf, call. = FALSE)
    }
    has_cnv <- !is.null(opts$cnv_vcf) && nzchar(opts$cnv_vcf) && file.exists(opts$cnv_vcf)
    if (!is.null(opts$cnv_vcf) && nzchar(opts$cnv_vcf) && !has_cnv) {
        message("CNV VCF not found; continuing in SV-only chromoplexy mode: ", opts$cnv_vcf)
    }

    if (!dir.exists(opts$out_dir)) {
        dir.create(opts$out_dir, recursive = TRUE)
    }

    cat("=======================================================\n")
    cat("OncoImplexus single-patient report\n")
    cat("Sample ID: ", opts$sample_id, "\n", sep = "")
    cat("Genome:    ", opts$genome, "\n", sep = "")
    cat("SV VCF:    ", opts$sv_vcf, "\n", sep = "")
    cat("CNV VCF:   ", if (has_cnv) opts$cnv_vcf else "not provided", "\n", sep = "")
    cat("Output:    ", opts$out_dir, "\n", sep = "")
    cat("=======================================================\n")

    cat("[1/6] Loading gene annotation...\n")
    gene_granges <- load_gene_granges(opts$genome, opts$gene_rds)

    cat("[2/6] Reading SV VCF...\n")
    sv_data <- read_sv_vcf(
        opts$sv_vcf,
        genome = opts$genome,
        min_sv_size = opts$min_sv_size
    )

    cnv_data <- NULL
    if (has_cnv) {
        cat("[3/6] Reading CNV VCF...\n")
        cnv_data <- read_cnv_vcf(opts$cnv_vcf)
    } else {
        cat("[3/6] CNV VCF unavailable; SV-only chromoplexy mode.\n")
    }

    cat("[4/6] Running chromoanagenesis analysis...\n")
    result <- detect_chromoanagenesis(
        SV.sample = sv_data,
        CNV.sample = cnv_data,
        genome = opts$genome,
        gene_granges = gene_granges,
        max_path_search = opts$max_path_search,
        max_neighbors = opts$max_neighbors,
        verbose = opts$verbose
    )

    rds_file <- file.path(opts$out_dir, paste0(opts$prefix, "_result.rds"))
    saveRDS(result, rds_file)

    if (opts$export_tables) {
        cat("[5/6] Exporting summary tables...\n")
        safe_export_tables(result, opts$out_dir, opts$prefix)
    } else {
        cat("[5/6] TSV export disabled.\n")
    }

    if (opts$enrichment) {
        run_enrichment_if_requested(result, opts$out_dir, opts$prefix)
    }

    cat("[6/6] Rendering HTML report...\n")
    report_file <- paste0(opts$prefix, "_Report.html")
    generate_interactive_report(
        result = result,
        SV.sample = sv_data,
        CNV.sample = cnv_data,
        output_file = report_file,
        output_dir = opts$out_dir,
        sample_name = opts$sample_id,
        genome = opts$genome,
        gene_granges = gene_granges
    )

    cat("\nSUCCESS\n")
    cat("HTML report: ", file.path(opts$out_dir, report_file), "\n", sep = "")
    cat("Result RDS:  ", rds_file, "\n", sep = "")
    if (opts$export_tables) {
        cat("Tables:      ", opts$out_dir, "\n", sep = "")
    }
}

tryCatch(main(), error = function(e) {
    cat("\nERROR: ", conditionMessage(e), "\n", sep = "", file = stderr())
    quit(status = 1, save = "no")
})
