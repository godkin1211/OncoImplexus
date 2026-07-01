#' Perform driver gene enrichment analysis for chromoanagenesis events
#'
#' @param cohort_summary A list containing summarized data for each sample, 
#'                       including detected mechanisms and affected genes.
#' @param target_mechanism The mechanism to analyze (e.g., "CT_High", "CP_Likely").
#' @return A data frame containing Fisher's exact test results for each gene.
#' @export
calculate_gene_enrichment <- function(cohort_summary, target_mechanism = "CT_High") {
  
  # 1. Flatten results to sample-level data frame
  all_samples <- unlist(cohort_summary, recursive = FALSE)
  sample_ids <- names(all_samples)
  
  if (length(all_samples) == 0) stop("Cohort summary is empty")
  
  # Determine which samples have the mechanism
  has_mech <- sapply(all_samples, function(x) {
    val <- x$mechs[[target_mechanism]]
    if (is.null(val)) return(FALSE)
    return(val > 0)
  })
  
  mech_positive_samples <- sample_ids[has_mech]
  mech_negative_samples <- sample_ids[!has_mech]
  
  n_pos <- length(mech_positive_samples)
  n_neg <- length(mech_negative_samples)
  
  if (n_pos == 0) {
    warning(paste("No samples found with mechanism:", target_mechanism))
    return(NULL)
  }
  
  # 2. Get the list of all disrupted driver genes
  all_affected_genes <- unique(unlist(lapply(all_samples, function(x) x$genes)))
  
  if (length(all_affected_genes) == 0) return(NULL)
  
  # 3. Perform Fisher's Exact Test for each gene
  results <- lapply(all_affected_genes, function(gene) {
    # Count occurrences
    a <- sum(sapply(all_samples[mech_positive_samples], function(x) gene %in% x$genes))
    b <- sum(sapply(all_samples[mech_negative_samples], function(x) gene %in% x$genes))
    
    # Contingency Table:
    #            Mech+ | Mech-
    # Gene+ |    a    |   b
    # Gene- | n_pos-a | n_neg-b
    
    mat <- matrix(c(a, b, n_pos - a, n_neg - b), nrow = 2)
    ft <- fisher.test(mat, alternative = "greater")
    
    data.frame(
      Gene = gene,
      Mech_Pos_Hits = a,
      Mech_Neg_Hits = b,
      Mech_Pos_Freq = a / n_pos,
      Mech_Neg_Freq = b / n_neg,
      OddsRatio = ft$estimate,
      PValue = ft$p.value,
      stringsAsFactors = FALSE
    )
  })
  
  enrich_df <- rbindlist(results)
  if (nrow(enrich_df) > 0) {
    enrich_df$FDR <- p.adjust(enrich_df$PValue, method = "BH")
    return(enrich_df[order(enrich_df$PValue), ])
  }
  
  return(enrich_df)
}

#' Extract genes affected by collapsed chromoplexy events
#'
#' @param x A chromoanagenesis result, chromoplexy result, collapsed event list,
#'   cohort summary from \code{summarize_chromoplexy_cohort_events()}, data frame,
#'   or character vector.
#' @param driver_only If TRUE, return only driver genes where driver annotation is
#'   available.
#' @return A unique character vector of gene symbols.
#' @export
extract_chromoplexy_genes <- function(x, driver_only = FALSE) {
  if (is.character(x)) {
    genes <- x
  } else if (is.data.frame(x)) {
    genes <- extract_gene_symbols_from_table(x, driver_only = driver_only)
  } else if (is.list(x) && !is.null(x$gene_summary)) {
    genes <- extract_gene_symbols_from_table(x$gene_summary, driver_only = driver_only)
  } else if (is.list(x) && !is.null(x$gene_event_summary)) {
    genes <- extract_gene_symbols_from_table(x$gene_event_summary, driver_only = driver_only)
  } else if (is.list(x) && !is.null(x$gene_detail)) {
    genes <- extract_gene_symbols_from_table(x$gene_detail, driver_only = driver_only)
  } else if (is.list(x) && !is.null(x$collapsed_events)) {
    genes <- extract_chromoplexy_genes(x$collapsed_events, driver_only = driver_only)
  } else if (is.list(x) && !is.null(x$chromoplexy)) {
    genes <- extract_chromoplexy_genes(x$chromoplexy, driver_only = driver_only)
  } else {
    stop("Unsupported input for chromoplexy gene extraction")
  }

  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  sort(genes)
}

#' Run GO BP and KEGG enrichment for chromoplexy-affected genes
#'
#' This is a thin, optional wrapper around clusterProfiler. clusterProfiler and
#' an organism annotation package are required only when this function is used.
#'
#' @param x Gene vector, chromoplexy result, chromoanagenesis result, collapsed
#'   event list, or cohort summary from \code{summarize_chromoplexy_cohort_events()}.
#' @param organism_db Annotation package name for SYMBOL-to-ENTREZ mapping.
#' @param kegg_organism KEGG organism code.
#' @param ontology GO ontology. Defaults to biological process ("BP").
#' @param pvalue_cutoff P-value cutoff passed to clusterProfiler.
#' @param qvalue_cutoff Q-value cutoff passed to clusterProfiler.
#' @param universe Optional background gene symbols.
#' @param min_genes Minimum mapped genes required before enrichment is attempted.
#' @param driver_only If TRUE, enrich only driver genes.
#' @param run_go Whether to run GO enrichment.
#' @param run_kegg Whether to run KEGG enrichment.
#' @param output_dir Optional directory for TSV outputs.
#' @param prefix Prefix for TSV output files.
#' @return A list with input genes, ID mapping, enrichment objects, and tables.
#' @export
run_chromoplexy_enrichment <- function(x,
                                       organism_db = "org.Hs.eg.db",
                                       kegg_organism = "hsa",
                                       ontology = "BP",
                                       pvalue_cutoff = 0.05,
                                       qvalue_cutoff = 0.20,
                                       universe = NULL,
                                       min_genes = 3,
                                       driver_only = FALSE,
                                       run_go = TRUE,
                                       run_kegg = TRUE,
                                       output_dir = NULL,
                                       prefix = "chromoplexy") {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("clusterProfiler is required. Install it with BiocManager::install('clusterProfiler').")
  }
  if (!requireNamespace(organism_db, quietly = TRUE)) {
    stop(sprintf("%s is required for gene ID mapping.", organism_db))
  }

  genes <- extract_chromoplexy_genes(x, driver_only = driver_only)
  genes <- unique(genes[!is.na(genes) & nzchar(genes)])
  orgdb <- getExportedValue(organism_db, organism_db)
  mapping <- map_symbols_to_entrez(genes, orgdb)

  universe_mapping <- NULL
  if (!is.null(universe)) {
    universe_mapping <- map_symbols_to_entrez(unique(as.character(universe)), orgdb)
  }

  entrez <- unique(mapping$ENTREZID)
  if (length(entrez) < min_genes) {
    warning(sprintf("Only %d genes mapped to ENTREZ IDs; enrichment was not run.", length(entrez)))
    out <- empty_chromoplexy_enrichment(genes, mapping, universe_mapping)
    if (!is.null(output_dir)) write_chromoplexy_enrichment_outputs(out, output_dir, prefix)
    return(out)
  }

  go_result <- NULL
  kegg_result <- NULL

  if (run_go) {
    go_result <- clusterProfiler::enrichGO(
      gene = entrez,
      universe = if (!is.null(universe_mapping)) unique(universe_mapping$ENTREZID) else NULL,
      OrgDb = orgdb,
      keyType = "ENTREZID",
      ont = ontology,
      pvalueCutoff = pvalue_cutoff,
      qvalueCutoff = qvalue_cutoff,
      readable = TRUE
    )
  }

  if (run_kegg) {
    kegg_result <- clusterProfiler::enrichKEGG(
      gene = entrez,
      universe = if (!is.null(universe_mapping)) unique(universe_mapping$ENTREZID) else NULL,
      organism = kegg_organism,
      pvalueCutoff = pvalue_cutoff,
      qvalueCutoff = qvalue_cutoff
    )
  }

  out <- list(
    input_genes = sort(genes),
    id_mapping = mapping,
    universe_mapping = universe_mapping,
    GO_BP = go_result,
    KEGG = kegg_result,
    GO_BP_table = enrichment_to_table(go_result),
    KEGG_table = enrichment_to_table(kegg_result),
    parameters = list(
      organism_db = organism_db,
      kegg_organism = kegg_organism,
      ontology = ontology,
      pvalue_cutoff = pvalue_cutoff,
      qvalue_cutoff = qvalue_cutoff,
      driver_only = driver_only,
      run_go = run_go,
      run_kegg = run_kegg
    )
  )

  if (!is.null(output_dir)) write_chromoplexy_enrichment_outputs(out, output_dir, prefix)
  out
}

#' Summarize chromoplexy pathway enrichment results
#'
#' @param enrichment_result Result from \code{run_chromoplexy_enrichment()}.
#' @param top_n Number of rows to retain per database.
#' @return A combined GO/KEGG table ordered by adjusted p-value.
#' @export
summarize_chromoplexy_pathways <- function(enrichment_result, top_n = 10) {
  tables <- list()
  if (!is.null(enrichment_result$GO_BP_table) && nrow(enrichment_result$GO_BP_table) > 0) {
    go <- enrichment_result$GO_BP_table
    go$source <- "GO_BP"
    tables$GO_BP <- go
  }
  if (!is.null(enrichment_result$KEGG_table) && nrow(enrichment_result$KEGG_table) > 0) {
    kegg <- enrichment_result$KEGG_table
    kegg$source <- "KEGG"
    tables$KEGG <- kegg
  }
  if (length(tables) == 0) return(data.frame())

  out <- do.call(rbind, lapply(tables, function(tbl) {
    tbl <- tbl[order(tbl$p.adjust, tbl$pvalue), , drop = FALSE]
    tbl[seq_len(min(nrow(tbl), top_n)), , drop = FALSE]
  }))
  show_cols <- intersect(c("source", "ID", "Description", "GeneRatio", "BgRatio",
                           "pvalue", "p.adjust", "qvalue", "geneID", "Count"),
                         colnames(out))
  out[, show_cols, drop = FALSE]
}

extract_gene_symbols_from_table <- function(tbl, driver_only = FALSE) {
  if (driver_only && "is_driver" %in% colnames(tbl)) {
    tbl <- tbl[tbl$is_driver %in% TRUE, , drop = FALSE]
  }
  if ("symbol" %in% colnames(tbl)) {
    return(tbl$symbol)
  }
  if ("Gene" %in% colnames(tbl)) {
    return(tbl$Gene)
  }
  if ("genes" %in% colnames(tbl)) {
    return(unlist(strsplit(paste(tbl$genes, collapse = ","), ","), use.names = FALSE))
  }
  if ("driver_genes" %in% colnames(tbl)) {
    return(unlist(strsplit(paste(tbl$driver_genes, collapse = ","), ","), use.names = FALSE))
  }
  character(0)
}

map_symbols_to_entrez <- function(symbols, orgdb) {
  symbols <- unique(symbols[!is.na(symbols) & nzchar(symbols)])
  if (length(symbols) == 0) {
    return(data.frame(SYMBOL = character(), ENTREZID = character(), stringsAsFactors = FALSE))
  }
  suppressMessages(clusterProfiler::bitr(
    symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = orgdb
  ))
}

enrichment_to_table <- function(enrichment) {
  if (is.null(enrichment)) return(data.frame())
  as.data.frame(enrichment)
}

empty_chromoplexy_enrichment <- function(genes, mapping, universe_mapping = NULL) {
  list(
    input_genes = sort(genes),
    id_mapping = mapping,
    universe_mapping = universe_mapping,
    GO_BP = NULL,
    KEGG = NULL,
    GO_BP_table = data.frame(),
    KEGG_table = data.frame(),
    parameters = list(run_go = FALSE, run_kegg = FALSE)
  )
}

write_chromoplexy_enrichment_outputs <- function(enrichment_result, output_dir, prefix) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  write_summary_tsv(enrichment_result$id_mapping, file.path(output_dir, paste0(prefix, "_input_gene_mapping.tsv")))
  write_summary_tsv(enrichment_result$GO_BP_table, file.path(output_dir, paste0(prefix, "_GO_BP_enrichment.tsv")))
  write_summary_tsv(enrichment_result$KEGG_table, file.path(output_dir, paste0(prefix, "_KEGG_enrichment.tsv")))
  summary_tbl <- summarize_chromoplexy_pathways(enrichment_result)
  write_summary_tsv(summary_tbl, file.path(output_dir, paste0(prefix, "_pathway_summary.tsv")))
}
