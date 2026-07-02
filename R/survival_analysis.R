#' Build sample-level survival features from OncoImplexus cohort results
#'
#' Converts chromoanagenesis and collapsed chromoplexy calls into one row per
#' sample. The output is designed for Kaplan-Meier and Cox regression analyses
#' and can be merged with clinical metadata.
#'
#' @param results A named list of chromoanagenesis results, an
#'   \code{OncoImplexusCohort} object, or a directory containing result
#'   \code{.rds} files.
#' @param clinical_data Optional clinical data frame. If \code{results} is an
#'   \code{OncoImplexusCohort} object and this is NULL, cohort clinical data are
#'   used automatically.
#' @param sample_id_col Column containing sample IDs in clinical data.
#' @param pathway_gene_sets Optional named list of gene vectors. For each gene
#'   set, \code{has_pathway_*} and \code{n_pathway_*} features are added.
#' @param output_file Optional TSV path for writing the feature table.
#' @return A data frame with one row per sample.
#' @export
build_chromoanagenesis_survival_features <- function(results,
                                                     clinical_data = NULL,
                                                     sample_id_col = "sample_id",
                                                     pathway_gene_sets = NULL,
                                                     output_file = NULL) {
    if (inherits(results, "OncoImplexusCohort") && is.null(clinical_data)) {
        clinical_data <- results@clinical
    }

    results_list <- coerce_chromoplexy_results_list(results)
    if (length(results_list) == 0) {
        stop("No chromoanagenesis results found")
    }

    feature_rows <- lapply(names(results_list), function(sample_id) {
        res <- unwrap_chromoanagenesis_result(results_list[[sample_id]])
        build_single_survival_feature_row(sample_id, res, pathway_gene_sets)
    })
    features <- do.call(rbind, feature_rows)
    feature_cols <- setdiff(colnames(features), c(
        sample_id_col, "analysis_mode", "affected_chromoplexy_genes",
        "chromoplexy_driver_genes"
    ))

    if (!is.null(clinical_data)) {
        features <- merge_clinical_data(features, clinical_data, sample_id_col)
    }
    attr(features, "onco_feature_cols") <- intersect(feature_cols, colnames(features))

    if (!is.null(output_file)) {
        output_dir <- dirname(output_file)
        if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
        utils::write.table(features, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    }

    features
}

#' Run chromoanagenesis survival association analysis
#'
#' Runs log-rank tests from dichotomized features and Cox regression from the
#' original numeric features. If covariates are supplied, both univariable and
#' multivariable Cox models are returned.
#'
#' @param features Sample-level feature table from
#'   \code{build_chromoanagenesis_survival_features()}.
#' @param clinical_data Optional clinical data frame to merge before analysis.
#' @param time_col Survival time column.
#' @param event_col Survival event column. Events should be coded as 1/0,
#'   TRUE/FALSE, or common text labels such as "event"/"censored".
#' @param sample_id_col Sample ID column.
#' @param feature_cols Optional feature columns to test. If NULL, OncoImplexus
#'   burden, complexity, driver, and pathway features are selected.
#' @param covariates Optional covariates for multivariable Cox regression.
#' @param dichotomize Whether to run KM/log-rank tests for continuous features.
#' @param cutpoint Dichotomization cutpoint. Currently "median" or a numeric
#'   scalar.
#' @param min_group_n Minimum samples per KM group.
#' @param min_events Minimum events required for a survival test.
#' @param p_adjust_method P-value adjustment method.
#' @param output_dir Optional directory for TSV outputs.
#' @param prefix Output prefix when \code{output_dir} is supplied.
#' @return A list with analysis data, Cox table, KM table, and KM group data.
#' @export
run_chromoanagenesis_survival <- function(features,
                                          clinical_data = NULL,
                                          time_col,
                                          event_col,
                                          sample_id_col = "sample_id",
                                          feature_cols = NULL,
                                          covariates = NULL,
                                          dichotomize = TRUE,
                                          cutpoint = "median",
                                          min_group_n = 5,
                                          min_events = 3,
                                          p_adjust_method = "BH",
                                          output_dir = NULL,
                                          prefix = "chromoanagenesis_survival") {
    require_survival_package()

    if (!is.data.frame(features)) {
        stop("features must be a data frame")
    }
    analysis_data <- if (!is.null(clinical_data)) {
        merge_clinical_data(features, clinical_data, sample_id_col)
    } else {
        features
    }

    required_cols <- c(sample_id_col, time_col, event_col)
    missing_cols <- setdiff(required_cols, colnames(analysis_data))
    if (length(missing_cols) > 0) {
        stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
    }

    analysis_data$.survival_time <- suppressWarnings(as.numeric(analysis_data[[time_col]]))
    analysis_data$.survival_event <- coerce_survival_event(analysis_data[[event_col]])
    keep <- !is.na(analysis_data$.survival_time) &
        analysis_data$.survival_time > 0 &
        !is.na(analysis_data$.survival_event)
    analysis_data <- analysis_data[keep, , drop = FALSE]

    if (nrow(analysis_data) < 2) {
        stop("Fewer than 2 samples have usable survival time/event data")
    }
    if (sum(analysis_data$.survival_event == 1, na.rm = TRUE) < min_events) {
        warning("Total survival events are fewer than min_events; most tests may be skipped.")
    }

    if (is.null(feature_cols)) {
        feature_cols <- select_default_survival_features(
            analysis_data,
            sample_id_col = sample_id_col,
            time_col = time_col,
            event_col = event_col,
            covariates = covariates,
            attr_features = attr(features, "onco_feature_cols")
        )
    }
    feature_cols <- intersect(feature_cols, colnames(analysis_data))
    if (length(feature_cols) == 0) {
        stop("No usable feature columns found")
    }
    missing_covars <- setdiff(covariates, colnames(analysis_data))
    if (length(missing_covars) > 0) {
        stop("Missing covariate columns: ", paste(missing_covars, collapse = ", "))
    }

    cox_rows <- list()
    km_rows <- list()
    km_group_data <- list()

    for (feature in feature_cols) {
        x <- normalize_survival_feature(analysis_data[[feature]])
        if (is.null(x) || length(unique(x[!is.na(x)])) < 2) {
            next
        }

        univ <- fit_survival_cox(
            analysis_data = analysis_data,
            feature_values = x,
            feature = feature,
            model = "univariable",
            covariates = NULL,
            min_events = min_events
        )
        if (!is.null(univ)) cox_rows[[length(cox_rows) + 1]] <- univ

        if (length(covariates) > 0) {
            multiv <- fit_survival_cox(
                analysis_data = analysis_data,
                feature_values = x,
                feature = feature,
                model = "multivariable",
                covariates = covariates,
                min_events = min_events
            )
            if (!is.null(multiv)) cox_rows[[length(cox_rows) + 1]] <- multiv
        }

        if (isTRUE(dichotomize)) {
            km <- fit_survival_km(
                analysis_data = analysis_data,
                feature_values = x,
                feature = feature,
                cutpoint = cutpoint,
                min_group_n = min_group_n,
                min_events = min_events,
                sample_id_col = sample_id_col
            )
            if (!is.null(km)) {
                km_rows[[length(km_rows) + 1]] <- km$table
                km_group_data[[feature]] <- km$group_data
            }
        }
    }

    cox_table <- bind_rows_or_empty(cox_rows)
    km_table <- bind_rows_or_empty(km_rows)
    if (nrow(cox_table) > 0) {
        cox_table$fdr <- adjust_p_by_group(cox_table$p_value, cox_table$model, p_adjust_method)
        cox_table <- cox_table[order(cox_table$model, cox_table$p_value), , drop = FALSE]
    }
    if (nrow(km_table) > 0) {
        km_table$fdr <- stats::p.adjust(km_table$p_value, method = p_adjust_method)
        km_table <- km_table[order(km_table$p_value), , drop = FALSE]
    }

    out <- list(
        features = features,
        analysis_data = analysis_data,
        feature_cols = feature_cols,
        cox_table = cox_table,
        km_table = km_table,
        km_group_data = km_group_data,
        parameters = list(
            time_col = time_col,
            event_col = event_col,
            sample_id_col = sample_id_col,
            covariates = covariates,
            dichotomize = dichotomize,
            cutpoint = cutpoint,
            min_group_n = min_group_n,
            min_events = min_events,
            p_adjust_method = p_adjust_method
        )
    )
    class(out) <- c("onco_survival", "list")

    if (!is.null(output_dir)) {
        write_survival_outputs(out, output_dir, prefix)
    }

    out
}

#' Generate an HTML cohort survival report
#'
#' @param survival_result Result from \code{run_chromoanagenesis_survival()}.
#' @param output_file Output HTML filename.
#' @param output_dir Output directory.
#' @param title Report title.
#' @return Path to the generated HTML report.
#' @export
generate_survival_report <- function(survival_result,
                                     output_file = "cohort_survival_report.html",
                                     output_dir = "reports/cohort_reports",
                                     title = "Cohort Chromoanagenesis Survival Analysis") {
    if (!inherits(survival_result, "onco_survival")) {
        stop("survival_result must come from run_chromoanagenesis_survival()")
    }
    if (!requireNamespace("rmarkdown", quietly = TRUE)) {
        stop("rmarkdown is required to generate the survival report")
    }
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    template_file <- system.file("report_templates", "survival_report.Rmd", package = "OncoImplexus")
    if (template_file == "" || !file.exists(template_file)) {
        template_file <- "inst/report_templates/survival_report.Rmd"
    }
    if (!file.exists(template_file)) {
        stop("Survival report template not found")
    }

    tmp_file <- tempfile(fileext = ".rds")
    saveRDS(list(survival_result = survival_result, title = title), tmp_file)
    on.exit(if (file.exists(tmp_file)) unlink(tmp_file), add = TRUE)

    rmarkdown::render(
        input = template_file,
        output_file = output_file,
        output_dir = output_dir,
        params = list(data_file = tmp_file),
        quiet = TRUE
    )

    file.path(output_dir, output_file)
}

#' Plot a Cox regression forest plot
#'
#' @param survival_result Result from \code{run_chromoanagenesis_survival()}.
#' @param model "multivariable" or "univariable".
#' @param top_n Number of features to show.
#' @return A ggplot object.
#' @export
plot_survival_forest <- function(survival_result,
                                 model = c("multivariable", "univariable"),
                                 top_n = 20) {
    model <- match.arg(model)
    cox_table <- survival_result$cox_table
    if (is.null(cox_table) || nrow(cox_table) == 0) {
        return(empty_survival_plot("No Cox regression results available"))
    }
    df <- cox_table[cox_table$model == model, , drop = FALSE]
    if (nrow(df) == 0 && model == "multivariable") {
        df <- cox_table[cox_table$model == "univariable", , drop = FALSE]
        model <- "univariable"
    }
    df <- df[is.finite(df$hazard_ratio) &
                 is.finite(df$conf_low) &
                 is.finite(df$conf_high) &
                 df$hazard_ratio > 0 &
                 df$conf_low > 0 &
                 df$conf_high > 0, , drop = FALSE]
    if (nrow(df) == 0) {
        return(empty_survival_plot("No finite Cox estimates available"))
    }
    df <- df[order(df$p_value), , drop = FALSE]
    df <- head(df, top_n)
    df$feature_label <- factor(df$feature, levels = rev(df$feature))

    ggplot2::ggplot(df, ggplot2::aes(x = hazard_ratio, y = feature_label)) +
        ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
        ggplot2::geom_errorbarh(
            ggplot2::aes(xmin = conf_low, xmax = conf_high),
            height = 0.18,
            color = "grey35"
        ) +
        ggplot2::geom_point(size = 2.4, color = "#2b6cb0") +
        ggplot2::scale_x_log10() +
        ggplot2::theme_minimal() +
        ggplot2::labs(
            x = "Hazard ratio (log scale)",
            y = NULL,
            title = paste("Top Cox Associations -", model)
        )
}

#' Plot a Kaplan-Meier curve for one feature
#'
#' @param survival_result Result from \code{run_chromoanagenesis_survival()}.
#' @param feature Feature name present in \code{survival_result$km_group_data}.
#' @return A ggplot object.
#' @export
plot_survival_km <- function(survival_result, feature) {
    require_survival_package()
    group_data <- survival_result$km_group_data[[feature]]
    if (is.null(group_data) || nrow(group_data) == 0) {
        return(empty_survival_plot(paste("No KM data available for", feature)))
    }
    fit <- survival::survfit(survival::Surv(time, event) ~ group, data = group_data)
    s <- summary(fit)
    if (length(s$time) == 0) {
        return(empty_survival_plot(paste("No survival curve available for", feature)))
    }
    plot_df <- data.frame(
        time = s$time,
        survival = s$surv,
        group = sub("^group=", "", as.character(s$strata)),
        stringsAsFactors = FALSE
    )
    starts <- unique(plot_df["group"])
    starts$time <- 0
    starts$survival <- 1
    plot_df <- rbind(starts[, c("time", "survival", "group")], plot_df)

    km_row <- survival_result$km_table[survival_result$km_table$feature == feature, , drop = FALSE]
    subtitle <- if (nrow(km_row) > 0) {
        sprintf("log-rank p = %.3g, FDR = %.3g", km_row$p_value[1], km_row$fdr[1])
    } else {
        NULL
    }

    ggplot2::ggplot(plot_df, ggplot2::aes(x = time, y = survival, color = group)) +
        ggplot2::geom_step(linewidth = 0.9) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
            x = survival_result$parameters$time_col,
            y = "Survival probability",
            color = feature,
            title = paste("Kaplan-Meier:", feature),
            subtitle = subtitle
        ) +
        ggplot2::coord_cartesian(ylim = c(0, 1))
}

#' @export
print.onco_survival <- function(x, ...) {
    cat("OncoImplexus survival analysis\n")
    cat("  Samples: ", nrow(x$analysis_data), "\n", sep = "")
    cat("  Features tested: ", length(x$feature_cols), "\n", sep = "")
    cat("  Cox rows: ", if (!is.null(x$cox_table)) nrow(x$cox_table) else 0, "\n", sep = "")
    cat("  KM/log-rank rows: ", if (!is.null(x$km_table)) nrow(x$km_table) else 0, "\n", sep = "")
    invisible(x)
}

build_single_survival_feature_row <- function(sample_id, res, pathway_gene_sets = NULL) {
    integrated <- if (!is.null(res$integrated_summary)) res$integrated_summary else list()
    cp <- if (!is.null(res$chromoplexy)) res$chromoplexy else NULL
    ce <- get_survival_collapsed_events(cp)
    events <- if (!is.null(ce$event_summary)) ce$event_summary else data.frame()
    gene_table <- get_survival_gene_table(ce)
    genes <- extract_symbols_from_survival_gene_table(gene_table, driver_only = FALSE)
    driver_genes <- extract_symbols_from_survival_gene_table(gene_table, driver_only = TRUE)

    ct_high <- num_or_default(integrated$chromothripsis_high_confidence,
                              num_or_default(res$chromothripsis$n_high_confidence, 0))
    ct_low <- num_or_default(integrated$chromothripsis_low_confidence,
                             num_or_default(res$chromothripsis$n_low_confidence, 0))
    cp_likely <- num_or_default(integrated$chromoplexy_likely,
                                num_or_default(cp$likely_chromoplexy, 0))
    cp_possible <- num_or_default(integrated$chromoplexy_possible,
                                  num_or_default(cp$possible_chromoplexy, 0))
    cs_likely <- num_or_default(integrated$chromoanasynthesis_likely,
                                num_or_default(res$chromoanasynthesis$likely_chromoanasynthesis, 0))
    cs_possible <- num_or_default(integrated$chromoanasynthesis_possible,
                                  num_or_default(res$chromoanasynthesis$possible_chromoanasynthesis, 0))

    row <- data.frame(
        sample_id = sample_id,
        analysis_mode = char_or_default(res$analysis_mode, NA_character_),
        is_sv_only = as.integer(grepl("SV-only", char_or_default(res$analysis_mode, ""), fixed = TRUE)),
        ct_high_confidence = ct_high,
        ct_low_confidence = ct_low,
        has_chromothripsis = as.integer((ct_high + ct_low) > 0),
        cp_likely_chains = cp_likely,
        cp_possible_chains = cp_possible,
        cp_total_chains = num_or_default(cp$total_chains, cp_likely + cp_possible),
        has_chromoplexy = as.integer((cp_likely + cp_possible) > 0 || nrow(events) > 0),
        cs_likely_events = cs_likely,
        cs_possible_events = cs_possible,
        has_chromoanasynthesis = as.integer((cs_likely + cs_possible) > 0),
        has_chromoanagenesis = as.integer((ct_high + ct_low + cp_likely + cp_possible + cs_likely + cs_possible) > 0 || nrow(events) > 0),
        n_chromoplexy_events = nrow(events),
        n_high_confidence_chromoplexy_events = count_matching(events, "event_confidence", "High"),
        n_moderate_confidence_chromoplexy_events = count_matching(events, "event_confidence", "Moderate"),
        n_low_confidence_chromoplexy_events = count_matching(events, "event_confidence", "Low"),
        has_high_confidence_chromoplexy_event = as.integer(count_matching(events, "event_confidence", "High") > 0),
        mean_event_qc_score = mean_col_or_na(events, "event_qc_score"),
        max_event_qc_score = max_col_or_na(events, "event_qc_score"),
        total_event_svs = sum_col_or_zero(events, "n_unique_svs"),
        total_event_breakpoints = sum_col_or_zero(events, "n_breakpoints"),
        max_n_svs = max_col_or_zero(events, "n_unique_svs"),
        max_n_breakpoints = max_col_or_zero(events, "n_breakpoints"),
        max_n_chromosomes = max_col_or_zero(events, "n_chromosomes"),
        has_multichromosomal_chromoplexy_event = as.integer(max_col_or_zero(events, "n_chromosomes") >= 2),
        has_cycle_chromoplexy_event = as.integer(any_col_true(events, "has_cycle")),
        n_chromoplexy_genes = length(genes),
        n_chromoplexy_driver_genes = length(driver_genes),
        has_chromoplexy_driver_gene = as.integer(length(driver_genes) > 0),
        n_fusions = if (!is.null(res$fusions) && is.data.frame(res$fusions)) nrow(res$fusions) else 0,
        affected_chromoplexy_genes = paste(genes, collapse = ","),
        chromoplexy_driver_genes = paste(driver_genes, collapse = ","),
        stringsAsFactors = FALSE
    )

    if (!is.null(pathway_gene_sets)) {
        pathway_features <- build_pathway_survival_features(genes, pathway_gene_sets)
        row <- cbind(row, pathway_features)
    }

    row
}

get_survival_collapsed_events <- function(cp) {
    if (is.null(cp)) return(empty_collapsed_chromoplexy_events())
    if (!is.null(cp$collapsed_events)) return(cp$collapsed_events)
    tryCatch(collapse_chromoplexy_chains(cp), error = function(e) empty_collapsed_chromoplexy_events())
}

get_survival_gene_table <- function(ce) {
    if (!is.null(ce$gene_event_summary) && nrow(ce$gene_event_summary) > 0) return(ce$gene_event_summary)
    if (!is.null(ce$gene_detail) && nrow(ce$gene_detail) > 0) return(ce$gene_detail)
    data.frame()
}

extract_symbols_from_survival_gene_table <- function(gene_table, driver_only = FALSE) {
    if (is.null(gene_table) || nrow(gene_table) == 0 || !"symbol" %in% colnames(gene_table)) {
        return(character(0))
    }
    x <- gene_table
    if (isTRUE(driver_only) && "is_driver" %in% colnames(x)) {
        x <- x[x$is_driver %in% TRUE, , drop = FALSE]
    } else if (isTRUE(driver_only)) {
        return(character(0))
    }
    genes <- unique(as.character(x$symbol))
    sort(genes[!is.na(genes) & nzchar(genes)])
}

build_pathway_survival_features <- function(genes, pathway_gene_sets) {
    if (!is.list(pathway_gene_sets) || is.null(names(pathway_gene_sets))) {
        stop("pathway_gene_sets must be a named list of gene vectors")
    }
    out <- data.frame(row.names = 1)
    genes_upper <- toupper(genes)
    for (nm in names(pathway_gene_sets)) {
        key <- sanitize_feature_name(nm)
        pathway_genes <- unique(toupper(as.character(pathway_gene_sets[[nm]])))
        n_hit <- sum(genes_upper %in% pathway_genes)
        out[[paste0("has_pathway_", key)]] <- as.integer(n_hit > 0)
        out[[paste0("n_pathway_", key, "_genes")]] <- n_hit
    }
    out
}

merge_clinical_data <- function(features, clinical_data, sample_id_col) {
    if (!is.data.frame(clinical_data)) stop("clinical_data must be a data frame")
    clinical <- clinical_data
    if (!sample_id_col %in% colnames(clinical)) {
        if (!is.null(rownames(clinical))) {
            clinical[[sample_id_col]] <- rownames(clinical)
        } else {
            stop("clinical_data must contain ", sample_id_col, " or have rownames")
        }
    }
    if (any(duplicated(clinical[[sample_id_col]]))) {
        stop("clinical_data contains duplicated sample IDs")
    }
    if (!sample_id_col %in% colnames(features)) {
        stop("features must contain ", sample_id_col)
    }
    idx <- seq_len(nrow(features))
    features$.feature_order <- idx
    merged <- merge(features, clinical, by = sample_id_col, all.x = TRUE, sort = FALSE)
    merged <- merged[order(merged$.feature_order), , drop = FALSE]
    merged$.feature_order <- NULL
    rownames(merged) <- NULL
    merged
}

select_default_survival_features <- function(analysis_data,
                                             sample_id_col,
                                             time_col,
                                             event_col,
                                             covariates = NULL,
                                             attr_features = NULL) {
    base_features <- c(
        "is_sv_only",
        "ct_high_confidence", "ct_low_confidence", "has_chromothripsis",
        "cp_likely_chains", "cp_possible_chains", "cp_total_chains", "has_chromoplexy",
        "cs_likely_events", "cs_possible_events", "has_chromoanasynthesis",
        "has_chromoanagenesis",
        "n_chromoplexy_events", "n_high_confidence_chromoplexy_events",
        "n_moderate_confidence_chromoplexy_events", "n_low_confidence_chromoplexy_events",
        "has_high_confidence_chromoplexy_event",
        "mean_event_qc_score", "max_event_qc_score",
        "total_event_svs", "total_event_breakpoints",
        "max_n_svs", "max_n_breakpoints", "max_n_chromosomes",
        "has_multichromosomal_chromoplexy_event", "has_cycle_chromoplexy_event",
        "n_chromoplexy_genes", "n_chromoplexy_driver_genes",
        "has_chromoplexy_driver_gene", "n_fusions"
    )
    pathway_features <- grep("^(has_pathway_|n_pathway_)", colnames(analysis_data), value = TRUE)
    candidates <- unique(c(attr_features, base_features, pathway_features))
    excluded <- c(sample_id_col, time_col, event_col, covariates,
                  ".survival_time", ".survival_event")
    candidates <- setdiff(intersect(candidates, colnames(analysis_data)), excluded)
    candidates[vapply(analysis_data[candidates], is_survival_testable_feature, logical(1))]
}

is_survival_testable_feature <- function(x) {
    x <- normalize_survival_feature(x)
    !is.null(x) && length(unique(x[!is.na(x)])) >= 2
}

normalize_survival_feature <- function(x) {
    if (is.logical(x)) return(as.numeric(x))
    if (is.integer(x) || is.numeric(x)) return(as.numeric(x))
    if (is.factor(x)) x <- as.character(x)
    if (is.character(x)) {
        vals <- unique(x[!is.na(x) & nzchar(x)])
        if (length(vals) == 2) return(as.numeric(factor(x, levels = sort(vals))) - 1)
    }
    NULL
}

fit_survival_cox <- function(analysis_data,
                             feature_values,
                             feature,
                             model,
                             covariates = NULL,
                             min_events = 3) {
    dat <- data.frame(
        .time = analysis_data$.survival_time,
        .event = analysis_data$.survival_event,
        .feature = feature_values
    )
    if (length(covariates) > 0) {
        for (i in seq_along(covariates)) {
            dat[[paste0(".cov", i)]] <- analysis_data[[covariates[[i]]]]
        }
    }
    dat <- dat[stats::complete.cases(dat), , drop = FALSE]
    if (nrow(dat) < 2 || sum(dat$.event == 1, na.rm = TRUE) < min_events ||
        length(unique(dat$.feature)) < 2) {
        return(NULL)
    }

    rhs <- ".feature"
    if (length(covariates) > 0) {
        rhs <- paste(c(rhs, paste0(".cov", seq_along(covariates))), collapse = " + ")
    }
    f <- stats::as.formula(paste("survival::Surv(.time, .event) ~", rhs))
    fit <- tryCatch(
        suppressWarnings(survival::coxph(f, data = dat)),
        error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)
    s <- summary(fit)
    coefs <- s$coefficients
    conf <- s$conf.int
    if (is.null(coefs) || !".feature" %in% rownames(coefs)) return(NULL)

    data.frame(
        feature = feature,
        model = model,
        n = nrow(dat),
        n_events = sum(dat$.event == 1, na.rm = TRUE),
        hazard_ratio = unname(conf[".feature", "exp(coef)"]),
        conf_low = unname(conf[".feature", "lower .95"]),
        conf_high = unname(conf[".feature", "upper .95"]),
        p_value = unname(coefs[".feature", "Pr(>|z|)"]),
        covariates = if (length(covariates) > 0) paste(covariates, collapse = ",") else "",
        stringsAsFactors = FALSE
    )
}

fit_survival_km <- function(analysis_data,
                            feature_values,
                            feature,
                            cutpoint,
                            min_group_n,
                            min_events,
                            sample_id_col) {
    group_info <- make_survival_groups(feature_values, cutpoint)
    if (is.null(group_info)) return(NULL)
    dat <- data.frame(
        sample_id = analysis_data[[sample_id_col]],
        time = analysis_data$.survival_time,
        event = analysis_data$.survival_event,
        feature_value = feature_values,
        group = group_info$group,
        stringsAsFactors = FALSE
    )
    dat <- dat[stats::complete.cases(dat[, c("time", "event", "group")]), , drop = FALSE]
    if (nrow(dat) < 2 || length(unique(dat$group)) < 2) return(NULL)
    group_counts <- table(dat$group)
    if (any(group_counts < min_group_n)) return(NULL)
    if (sum(dat$event == 1, na.rm = TRUE) < min_events) return(NULL)

    diff <- tryCatch(
        survival::survdiff(survival::Surv(time, event) ~ group, data = dat),
        error = function(e) NULL
    )
    fit <- tryCatch(
        survival::survfit(survival::Surv(time, event) ~ group, data = dat),
        error = function(e) NULL
    )
    if (is.null(diff) || is.null(fit)) return(NULL)
    p_value <- stats::pchisq(diff$chisq, df = length(diff$n) - 1, lower.tail = FALSE)

    group_events <- stats::aggregate(event ~ group, dat, sum)
    names(group_events)[2] <- "events"
    group_ns <- data.frame(group = names(group_counts), n = as.integer(group_counts))
    group_stats <- merge(group_ns, group_events, by = "group", all.x = TRUE)
    table_row <- data.frame(
        feature = feature,
        n = nrow(dat),
        n_events = sum(dat$event == 1, na.rm = TRUE),
        cutpoint = as.character(group_info$cutpoint),
        groups = paste(sprintf("%s n=%d events=%d",
                               group_stats$group, group_stats$n, group_stats$events),
                       collapse = "; "),
        chisq = unname(diff$chisq),
        p_value = p_value,
        stringsAsFactors = FALSE
    )
    list(table = table_row, group_data = dat)
}

make_survival_groups <- function(x, cutpoint = "median") {
    vals <- x[!is.na(x)]
    uniq <- sort(unique(vals))
    if (length(uniq) < 2) return(NULL)
    if (length(uniq) == 2 && all(uniq %in% c(0, 1))) {
        group <- ifelse(x == 1, "Present", "Absent")
        return(list(group = factor(group, levels = c("Absent", "Present")),
                    cutpoint = "binary"))
    }
    if (identical(cutpoint, "median")) {
        cp <- stats::median(vals, na.rm = TRUE)
    } else if (is.numeric(cutpoint) && length(cutpoint) == 1) {
        cp <- cutpoint
    } else {
        stop("cutpoint must be 'median' or a numeric scalar")
    }
    if (!is.finite(cp) || all(vals <= cp) || all(vals > cp)) return(NULL)
    group <- ifelse(x > cp, "High", "Low")
    list(group = factor(group, levels = c("Low", "High")), cutpoint = cp)
}

coerce_survival_event <- function(x) {
    if (is.logical(x)) return(as.integer(x))
    if (is.numeric(x) || is.integer(x)) {
        out <- as.numeric(x)
        out[!(out %in% c(0, 1))] <- NA_real_
        return(out)
    }
    vals <- tolower(trimws(as.character(x)))
    out <- rep(NA_real_, length(vals))
    out[vals %in% c("1", "true", "t", "yes", "y", "event", "dead", "deceased",
                    "progressed", "progression", "relapsed")] <- 1
    out[vals %in% c("0", "false", "f", "no", "n", "censored", "alive",
                    "none", "no event")] <- 0
    numeric_vals <- suppressWarnings(as.numeric(vals))
    fill <- is.na(out) & numeric_vals %in% c(0, 1)
    out[fill] <- numeric_vals[fill]
    out
}

write_survival_outputs <- function(x, output_dir, prefix) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    write_survival_tsv(x$features, file.path(output_dir, paste0(prefix, "_features.tsv")))
    write_survival_tsv(x$analysis_data, file.path(output_dir, paste0(prefix, "_analysis_data.tsv")))
    write_survival_tsv(x$cox_table, file.path(output_dir, paste0(prefix, "_cox.tsv")))
    write_survival_tsv(x$km_table, file.path(output_dir, paste0(prefix, "_km_logrank.tsv")))
    saveRDS(x, file.path(output_dir, paste0(prefix, "_result.rds")))
}

write_survival_tsv <- function(x, file) {
    if (is.null(x) || !is.data.frame(x)) return(invisible(FALSE))
    utils::write.table(x, file, sep = "\t", quote = FALSE, row.names = FALSE)
    invisible(TRUE)
}

adjust_p_by_group <- function(p, group, method) {
    out <- rep(NA_real_, length(p))
    for (g in unique(group)) {
        idx <- which(group == g & !is.na(p))
        if (length(idx) > 0) out[idx] <- stats::p.adjust(p[idx], method = method)
    }
    out
}

require_survival_package <- function() {
    if (!requireNamespace("survival", quietly = TRUE)) {
        stop("The survival package is required. Install it with install.packages('survival').")
    }
}

empty_survival_plot <- function(label) {
    ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0, y = 0, label = label) +
        ggplot2::theme_void()
}

num_or_default <- function(x, default = 0) {
    if (is.null(x) || length(x) == 0) return(default)
    out <- suppressWarnings(as.numeric(x[1]))
    if (is.na(out)) default else out
}

char_or_default <- function(x, default = "") {
    if (is.null(x) || length(x) == 0 || is.na(x[1])) return(default)
    as.character(x[1])
}

count_matching <- function(df, col, value) {
    if (is.null(df) || nrow(df) == 0 || !col %in% colnames(df)) return(0)
    sum(df[[col]] == value, na.rm = TRUE)
}

mean_col_or_na <- function(df, col) {
    if (is.null(df) || nrow(df) == 0 || !col %in% colnames(df)) return(NA_real_)
    vals <- suppressWarnings(as.numeric(df[[col]]))
    if (all(is.na(vals))) NA_real_ else mean(vals, na.rm = TRUE)
}

max_col_or_na <- function(df, col) {
    if (is.null(df) || nrow(df) == 0 || !col %in% colnames(df)) return(NA_real_)
    vals <- suppressWarnings(as.numeric(df[[col]]))
    if (all(is.na(vals))) NA_real_ else max(vals, na.rm = TRUE)
}

max_col_or_zero <- function(df, col) {
    out <- max_col_or_na(df, col)
    if (is.na(out)) 0 else out
}

sum_col_or_zero <- function(df, col) {
    if (is.null(df) || nrow(df) == 0 || !col %in% colnames(df)) return(0)
    vals <- suppressWarnings(as.numeric(df[[col]]))
    sum(vals, na.rm = TRUE)
}

any_col_true <- function(df, col) {
    if (is.null(df) || nrow(df) == 0 || !col %in% colnames(df)) return(FALSE)
    any(df[[col]] %in% TRUE, na.rm = TRUE)
}

sanitize_feature_name <- function(x) {
    x <- tolower(gsub("[^A-Za-z0-9]+", "_", x))
    x <- gsub("^_+|_+$", "", x)
    if (!nzchar(x)) "geneset" else x
}
