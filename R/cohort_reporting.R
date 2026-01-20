#' Generate a comprehensive cohort-level analysis report
#'
#' Aggregates results from multiple chromoanagenesis analysis files and 
#' generates an interactive HTML summary report.
#'
#' @param results_dir Directory containing .rds result files from detect_chromoanagenesis()
#' @param output_file Output filename (default: "cohort_analysis_report.html")
#' @param clinical_data Optional data frame with clinical metadata
#' @param title Report title
#' @return Path to the generated report
#' @export
generate_cohort_report <- function(results_dir,
                                  output_file = "cohort_analysis_report.html",
                                  clinical_data = NULL,
                                  title = "Cohort Chromoanagenesis Summary") {

    output_dir <- "reports/cohort_reports"
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    if (!dir.exists(results_dir)) stop("Results directory not found.")

    # 1. Load all RDS files in the directory
    cat("Loading result files from", results_dir, "...\n")
    files <- list.files(results_dir, pattern = "\\.rds$", full.names = TRUE)
    
    if (length(files) == 0) stop("No .rds files found in the specified directory.")

    results_list <- list()
    for (f in files) {
        sample_id <- gsub("\\.rds$", "", basename(f))
        tryCatch({
            obj <- readRDS(f)
            # Support both full objects and list-wrapped objects
            if ("results" %in% names(obj)) {
                results_list[[sample_id]] <- obj$results
                # Sync annotation if it exists in the outer list
                if (!is.null(obj$annotation)) results_list[[sample_id]]$annotation <- obj$annotation
            } else {
                results_list[[sample_id]] <- obj
            }
        }, error = function(e) {
            warning(sprintf("Failed to load %s: %s", f, e$message))
        })
    }

    if (length(results_list) == 0) stop("Could not load any valid result objects.")

    # 2. Create OncoImplexusCohort object
    cat("Creating cohort object...\n")
    # If no clinical data, create dummy
    if (is.null(clinical_data)) {
        clinical_data <- data.frame(
            sample_id = names(results_list),
            row.names = names(results_list),
            stringsAsFactors = FALSE
        )
    }
    
    cohort <- create_cohort(results_list, clinical_data)

    # 3. Render Report
    template_file <- system.file("report_templates", "cohort_report.Rmd", package = "OncoImplexus")
    if (template_file == "" || !file.exists(template_file)) {
        template_file <- "inst/report_templates/cohort_report.Rmd"
    }
    
    # Save temporary data for Rmd
    tmp_file <- tempfile(fileext = ".rds")
    saveRDS(list(cohort = cohort, title = title), tmp_file)

    cat("Rendering cohort HTML report...\n")
    tryCatch({
        rmarkdown::render(
            input = template_file,
            output_file = output_file,
            output_dir = output_dir,
            params = list(data_file = tmp_file),
            quiet = TRUE
        )
        cat(sprintf("Cohort report successfully generated: %s\n", file.path(output_dir, output_file)))
    }, finally = {
        if (file.exists(tmp_file)) unlink(tmp_file)
    })

    return(file.path(output_dir, output_file))
}
