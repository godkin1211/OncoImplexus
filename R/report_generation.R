#' @include chromoanagenesis_integrated.R circos_visualization.R gene_annotation_utils.R quality_control.R
NULL

#' Generate an interactive HTML report for Chromoanagenesis analysis
#' 
#' @param result A chromoanagenesis result object from detect_chromoanagenesis()
#' @param SV.sample The original SV data used for analysis
#' @param CNV.sample The original CNV data used for analysis
#' @param output_file Output HTML filename
#' @param output_dir Output directory
#' @param sample_name Sample name
#' @param genome Reference genome
#' @param gene_granges GRanges for annotation
#' @param txdb TxDb for fusion plotting
#' @export
generate_interactive_report <- function(result,
                                      SV.sample,
                                      CNV.sample,
                                      output_file = "report.html",
                                      output_dir = "reports/sample_reports",
                                      sample_name = "Sample",
                                      genome = "hg19",
                                      gene_granges = NULL,
                                      txdb = NULL) {

    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    # Locate template file
    template_file <- system.file("report_templates", "main_report.Rmd", package = "OncoImplexus")
    if (template_file == "" || !file.exists(template_file)) template_file <- "inst/report_templates/main_report.Rmd"

    # Prepare data with ULTRA-SHORT relative path to avoid connection issues
    tmp_data_file <- "tmp_rep.rds"
    
    # 1. Annotation
    # If gene_granges is provided, re-annotate. 
    # Otherwise, try to reuse existing annotations from the result object.
    annotation_data <- NULL
    if (!is.null(gene_granges)) {
        annotation_data <- annotate_chromoanagenesis(result, gene_granges)
        
        # Also predict fusions if missing or requested
        if (is.null(result$fusions) || nrow(result$fusions) == 0) {
            result$fusions <- tryCatch({
                predict_fusion_genes(SV.sample, gene_granges)
            }, error = function(e) {
                message("  Warning: Fusion prediction failed during report generation: ", e$message)
                return(NULL)
            })
        }
    } else {
        # Failover to existing annotations in the result object
        if (!is.null(result$gene_annotations)) {
            annotation_data <- result$gene_annotations
        } else if (!is.null(result$annotation)) {
            annotation_data <- result$annotation
        }
    }

    # 2. Repair Mechanism Analysis
    repair_data <- NULL
    if (requireNamespace("BSgenome", quietly = TRUE)) {
        try({
            g_obj <- BSgenome::getBSgenome(genome)
            repair_data <- analyze_breakpoint_sequences(SV.sample, g_obj)
            if (!is.null(repair_data) && !is.null(result$fusions)) {
                # ENSURE CHARACTER MATCHING
                fus_ids <- as.character(result$fusions$sv_id)
                rep_ids <- as.character(repair_data$sequences$sv_id)
                
                # Check for overlap and provide feedback
                match_count <- sum(fus_ids %in% rep_ids)
                message(sprintf("  -> Sequence matching: %d/%d fusions linked to sequence data.", 
                                match_count, length(fus_ids)))
                
                if (match_count == 0 && length(fus_ids) > 0) {
                    message("  [CRITICAL] ID Mismatch Alert!")
                    message(sprintf("  [CRITICAL] First fusion ID: '%s'", fus_ids[1]))
                    message(sprintf("  [CRITICAL] First repair ID: '%s'", rep_ids[1]))
                }
                
                result$fusions$sequence_preview <- sapply(1:nrow(result$fusions), function(i) {
                    row <- result$fusions[i, ]
                    get_fusion_seq_html(row$sv_id, repair_data, chrom1 = row$chrom1, pos1 = row$pos1)
                })
            }
        })
    }

    # 3. Diagnostics
    qc_metrics <- assess_data_quality(SV.sample, CNV.sample, genome)

    # 4. TxDb Info
    txdb_pkg <- if (!is.null(txdb)) (if (is.character(txdb)) txdb else "provided_txdb") else NULL

    # 5. Create report list (S4-safe)
    sv_report_df <- if(isS4(SV.sample)) as(SV.sample, "data.frame") else as.data.frame(SV.sample)
    if (isS4(SV.sample) && !("sv_id" %in% colnames(sv_report_df))) {
        sv_report_df$sv_id <- as.character(SV.sample@sv_id)
    }

    report_list <- list(
        result = result,
        SV.sample = sv_report_df,
        CNV.sample = if(isS4(CNV.sample)) as(CNV.sample, "data.frame") else CNV.sample,
        sample_name = sample_name,
        genome = genome,
        annotation = annotation_data,
        repair_mechanisms = repair_data,
        qc_metrics = qc_metrics,
        txdb_name = txdb_pkg
    )
    
    # Prepare rendering environment to inject data directly
    render_env <- new.env(parent = .GlobalEnv)
    render_env$report_data <- report_list
    
    # Render
    message("Rendering HTML report...")
    render_status <- tryCatch({
        rmarkdown::render(
            input = template_file,
            output_file = output_file,
            output_dir = output_dir,
            envir = render_env, # Direct injection!
            quiet = TRUE
        )
        TRUE
    }, error = function(e) {
        message("Error during report rendering: ", e$message)
        FALSE
    })

    if (render_status) {
        message(sprintf("Report successfully generated: %s", file.path(output_dir, output_file)))
    }

    return(list(report_path = file.path(output_dir, output_file), result_object = result))
}