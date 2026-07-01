#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(S4Vectors)
})

args <- commandArgs(trailingOnly = TRUE)
result_file <- if (length(args) >= 1) args[[1]] else "analysis_output/ont_bam_pass/bam_pass_sv_only_chromoplexy_result.rds"
out_dir <- if (length(args) >= 2) args[[2]] else "analysis_output/ont_bam_pass"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

res <- readRDS(result_file)
genes <- readRDS(system.file("extdata", "hg38_genes.rds", package = "OncoImplexus"))

summary_df <- res$chromoplexy$summary
likely_ids <- summary_df$chain_id[summary_df$classification == "Likely chromoplexy"]
details <- res$chromoplexy$chain_details[
    vapply(res$chromoplexy$chain_details, function(x) x$summary$chain_id %in% likely_ids, logical(1))
]
chain_ids <- vapply(details, function(x) as.integer(x$summary$chain_id), integer(1))

all_nodes <- sort(unique(unlist(lapply(details, function(x) x$chain$nodes))))
node_index <- setNames(seq_along(all_nodes), all_nodes)
parent <- seq_along(all_nodes)

find_root <- function(x) {
    while (parent[x] != x) {
        parent[x] <<- parent[parent[x]]
        x <- parent[x]
    }
    x
}

union_nodes <- function(a, b) {
    ra <- find_root(a)
    rb <- find_root(b)
    if (ra != rb) parent[rb] <<- ra
}

for (detail in details) {
    nodes <- detail$chain$nodes
    if (length(nodes) >= 2) {
        for (i in seq_len(length(nodes) - 1L)) {
            union_nodes(node_index[[nodes[i]]], node_index[[nodes[i + 1L]]])
        }
    }
}

node_root <- vapply(seq_along(all_nodes), find_root, integer(1))
root_key <- setNames(node_root, all_nodes)
chain_root <- vapply(details, function(x) {
    roots <- unique(root_key[x$chain$nodes])
    if (length(roots) != 1) {
        stop("Chain spans multiple collapsed components: ", x$summary$chain_id)
    }
    roots[[1]]
}, integer(1))

raw_roots <- sort(unique(chain_root))
root_rows <- lapply(raw_roots, function(root) {
    cids <- chain_ids[chain_root == root]
    s <- summary_df[match(cids, summary_df$chain_id), , drop = FALSE]
    best_i <- which.max(s$combined_score)
    data.frame(
        raw_root = root,
        representative_chain_id = s$chain_id[best_i],
        representative_score = s$combined_score[best_i],
        n_likely_chains = length(cids),
        stringsAsFactors = FALSE
    )
})
root_df <- do.call(rbind, root_rows)
root_df <- root_df[order(-root_df$representative_score, -root_df$n_likely_chains, root_df$representative_chain_id), ]
root_df$collapsed_event_id <- sprintf("CE%03d", seq_len(nrow(root_df)))
root_to_event <- setNames(root_df$collapsed_event_id, root_df$raw_root)

chain_map <- data.frame(
    chain_id = chain_ids,
    collapsed_event_id = unname(root_to_event[as.character(chain_root)]),
    stringsAsFactors = FALSE
)
chain_map <- merge(
    chain_map,
    summary_df[, c("chain_id", "chromosomes_involved", "n_translocations", "combined_score", "pvalue", "fdr", "classification")],
    by = "chain_id",
    all.x = TRUE
)
chain_map <- chain_map[order(chain_map$collapsed_event_id, -chain_map$combined_score, chain_map$chain_id), ]

event_rows <- list()
breakpoint_rows <- list()
bp_idx <- 1L

for (event_id in root_df$collapsed_event_id) {
    cids <- chain_map$chain_id[chain_map$collapsed_event_id == event_id]
    edetails <- details[chain_ids %in% cids]
    esum <- summary_df[match(cids, summary_df$chain_id), , drop = FALSE]
    sv_df <- do.call(rbind, lapply(edetails, function(x) x$SVs))
    sv_df <- sv_df[!duplicated(sv_df$sv_id), , drop = FALSE]
    event_nodes <- sort(unique(unlist(lapply(edetails, function(x) x$chain$nodes))))
    chroms <- sort(unique(c(as.character(sv_df$chrom1), as.character(sv_df$chrom2))))
    rep_chain <- root_df$representative_chain_id[root_df$collapsed_event_id == event_id]

    event_rows[[event_id]] <- data.frame(
        collapsed_event_id = event_id,
        n_likely_chains = length(cids),
        representative_chain_id = rep_chain,
        representative_combined_score = max(esum$combined_score, na.rm = TRUE),
        min_pvalue = min(esum$pvalue, na.rm = TRUE),
        min_fdr = min(esum$fdr, na.rm = TRUE),
        n_unique_svs = nrow(sv_df),
        n_unique_breakpoint_nodes = length(event_nodes),
        n_chromosomes = length(chroms),
        chromosomes_involved = paste(chroms, collapse = ","),
        chain_ids = paste(sort(cids), collapse = ","),
        sv_ids = paste(sort(unique(sv_df$sv_id)), collapse = ","),
        stringsAsFactors = FALSE
    )

    for (i in seq_len(nrow(sv_df))) {
        breakpoint_rows[[bp_idx]] <- data.frame(
            collapsed_event_id = event_id,
            sv_id = sv_df$sv_id[i],
            breakpoint_side = "left",
            chrom = paste0("chr", sv_df$chrom1[i]),
            pos = as.integer(sv_df$pos1[i]),
            partner_chrom = paste0("chr", sv_df$chrom2[i]),
            partner_pos = as.integer(sv_df$pos2[i]),
            SVtype = sv_df$SVtype[i],
            stringsAsFactors = FALSE
        )
        bp_idx <- bp_idx + 1L
        breakpoint_rows[[bp_idx]] <- data.frame(
            collapsed_event_id = event_id,
            sv_id = sv_df$sv_id[i],
            breakpoint_side = "right",
            chrom = paste0("chr", sv_df$chrom2[i]),
            pos = as.integer(sv_df$pos2[i]),
            partner_chrom = paste0("chr", sv_df$chrom1[i]),
            partner_pos = as.integer(sv_df$pos1[i]),
            SVtype = sv_df$SVtype[i],
            stringsAsFactors = FALSE
        )
        bp_idx <- bp_idx + 1L
    }
}

event_summary <- do.call(rbind, event_rows)
event_bp <- do.call(rbind, breakpoint_rows)
event_bp$breakpoint_id <- paste(event_bp$sv_id, event_bp$breakpoint_side, sep = "|")

bp_gr <- GRanges(seqnames = event_bp$chrom, ranges = IRanges(start = event_bp$pos, width = 1), strand = "*")
hits <- findOverlaps(bp_gr, genes, ignore.strand = TRUE)
if (length(hits) > 0) {
    bph <- event_bp[queryHits(hits), , drop = FALSE]
    gh <- subjectHits(hits)
    gm <- mcols(genes)
    event_gene_detail <- cbind(
        bph,
        gene_id = as.character(gm$gene_id[gh]),
        symbol = as.character(gm$symbol[gh]),
        gene_chrom = as.character(seqnames(genes)[gh]),
        gene_start = as.integer(start(genes)[gh]),
        gene_end = as.integer(end(genes)[gh]),
        gene_strand = as.character(strand(genes)[gh]),
        stringsAsFactors = FALSE
    )
} else {
    event_gene_detail <- data.frame()
}

if (nrow(event_gene_detail) > 0) {
    event_gene_summary <- do.call(rbind, lapply(sort(unique(event_gene_detail$collapsed_event_id)), function(event_id) {
        sub <- event_gene_detail[event_gene_detail$collapsed_event_id == event_id, , drop = FALSE]
        data.frame(
            collapsed_event_id = event_id,
            n_genes = length(unique(sub$gene_id)),
            genes = paste(sort(unique(sub$symbol)), collapse = ","),
            gene_ids = paste(sort(unique(sub$gene_id)), collapse = ","),
            n_gene_overlapping_svs = length(unique(sub$sv_id)),
            n_gene_overlapping_breakpoints = length(unique(sub$breakpoint_id)),
            stringsAsFactors = FALSE
        )
    }))

    gene_event_summary <- do.call(rbind, lapply(sort(unique(event_gene_detail$gene_id)), function(gid) {
        sub <- event_gene_detail[event_gene_detail$gene_id == gid, , drop = FALSE]
        data.frame(
            gene_id = gid,
            symbol = unique(sub$symbol)[1],
            n_collapsed_events = length(unique(sub$collapsed_event_id)),
            collapsed_event_ids = paste(sort(unique(sub$collapsed_event_id)), collapse = ","),
            n_unique_svs = length(unique(sub$sv_id)),
            n_unique_breakpoints = length(unique(sub$breakpoint_id)),
            chromosomes = paste(sort(unique(sub$chrom)), collapse = ","),
            sv_ids = paste(sort(unique(sub$sv_id)), collapse = ","),
            stringsAsFactors = FALSE
        )
    }))
    gene_event_summary <- gene_event_summary[order(-gene_event_summary$n_collapsed_events, gene_event_summary$symbol), ]

    event_summary <- merge(event_summary, event_gene_summary, by = "collapsed_event_id", all.x = TRUE)
    event_summary$n_genes[is.na(event_summary$n_genes)] <- 0
    event_summary$genes[is.na(event_summary$genes)] <- ""
    event_summary$gene_ids[is.na(event_summary$gene_ids)] <- ""
    event_summary$n_gene_overlapping_svs[is.na(event_summary$n_gene_overlapping_svs)] <- 0
    event_summary$n_gene_overlapping_breakpoints[is.na(event_summary$n_gene_overlapping_breakpoints)] <- 0
} else {
    gene_event_summary <- data.frame()
    event_summary$n_genes <- 0
    event_summary$genes <- ""
    event_summary$gene_ids <- ""
    event_summary$n_gene_overlapping_svs <- 0
    event_summary$n_gene_overlapping_breakpoints <- 0
}

event_summary <- event_summary[order(event_summary$collapsed_event_id), ]

write.table(chain_map, file.path(out_dir, "bam_pass_likely_chromoplexy_chain_to_collapsed_event.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(event_summary, file.path(out_dir, "bam_pass_likely_chromoplexy_collapsed_events.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(event_gene_detail, file.path(out_dir, "bam_pass_likely_chromoplexy_collapsed_event_gene_detail.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gene_event_summary, file.path(out_dir, "bam_pass_likely_chromoplexy_collapsed_gene_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("Likely chains:", length(likely_ids), "\n")
cat("Collapsed events:", nrow(event_summary), "\n")
cat("Collapsed event sizes:\n")
print(event_summary[, c(
    "collapsed_event_id", "n_likely_chains", "representative_chain_id",
    "representative_combined_score", "n_unique_svs", "n_chromosomes",
    "chromosomes_involved", "n_genes", "genes"
)], row.names = FALSE)
cat("\nCollapsed gene summary:\n")
print(gene_event_summary[, c(
    "symbol", "gene_id", "n_collapsed_events", "collapsed_event_ids",
    "n_unique_svs", "n_unique_breakpoints", "chromosomes"
)], row.names = FALSE)
