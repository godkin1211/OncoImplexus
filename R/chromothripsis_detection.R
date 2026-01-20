#' Class to store SV data
#' @param chrom1 (character): chromosome for the first breakpoint
#' @param pos1 (character): position for the first breakpoint
#' @param chrom2 (character): chromosome for the second breakpoint
#' @param pos2 (character): position for the second breakpoint
#' @param SVtype (character): type of SV, encoded as: DEL (deletion-like; +/-), DUP (duplication-like; -/+), h2hINV (head-to-head inversion; +/+), and t2tINV (tail-to-tail inversion; -/-).
#' @param strand1 (e.g. + for DEL)
#' @param strand2 (e.g. - for DEL)
#' @return an instance of the class 'SVs' that contains SV data. Required format by the function detect_chromothripsis
#' @export

SVs <- setClass(
    "SVs",
    representation(
        chrom1 = "character",
        pos1 = "numeric",
        chrom2 = "character",
        pos2 = "numeric",
        SVtype = "character",
        strand1 = "character",
        strand2 = "character",
        sv_id = "character"
    )
)

setMethod("initialize", "SVs", function(.Object, ...) {
    .Object <- callNextMethod()
    if (length(.Object@chrom1) != length(.Object@pos1) || length(.Object@chrom1) != length(.Object@pos2)) {
        stop("slots lengths are not all equal")
    }
    if (length(.Object@chrom2) == 0) .Object@chrom2 <- .Object@chrom1
    if (length(.Object@SVtype) == 0) .Object@SVtype <- rep("", length(.Object@chrom1))

    # CRITICAL: Preserve stable IDs. Only generate sequential if absolutely empty
    if (length(.Object@sv_id) == 0) {
        .Object@sv_id <- as.character(1:length(.Object@chrom1))
    } else if (length(.Object@sv_id) != length(.Object@chrom1)) {
        # This handles subsetting issues where sv_id might get out of sync
        warning("sv_id length mismatch. Resetting to sequential.")
        .Object@sv_id <- as.character(1:length(.Object@chrom1))
    }

    if (length(.Object@chrom1) != length(.Object@SVtype) || length(.Object@chrom1) != length(.Object@chrom2)) stop("slots lengths are not all equal")

    # Only swap positions for intra-chromosomal SVs
    ind <- which(.Object@chrom1 == .Object@chrom2 & .Object@pos1 > .Object@pos2)
    if (length(ind) > 0) {
        # Swap chromosomes (same anyway, but for consistency)
        tmp_chr <- .Object@chrom2[ind]
        .Object@chrom2[ind] <- .Object@chrom1[ind]
        .Object@chrom1[ind] <- tmp_chr

        # Swap positions
        tmp_pos <- .Object@pos2[ind]
        .Object@pos2[ind] <- .Object@pos1[ind]
        .Object@pos1[ind] <- tmp_pos

        # Swap strands (CRITICAL: strands must follow their breakpoints)
        if (length(.Object@strand1) > 0 && length(.Object@strand2) > 0) {
            tmp_s <- .Object@strand2[ind]
            .Object@strand2[ind] <- .Object@strand1[ind]
            .Object@strand1[ind] <- tmp_s
        }
    }
    # .Object@numSV = length(.Object@chrom1)
    ind.match.sv1 <- match(.Object@chrom1, chromNames)
    ind.match.sv2 <- match(.Object@chrom2, chromNames)

    if (sum(is.na(ind.match.sv1)) != 0 | sum(is.na(ind.match.sv2)) != 0) {
        stop(paste("chromosome name must be in \"", paste(chromNames, collapse = " ")), "\"")
    }

    .Object
})

# setMethod("show","SVs",function(object){
# 			  print(paste("SVs with", object@numSV, "structural variations"))
# 				})


setAs("SVs", "data.frame", function(from, to) {
    to <- data.frame(
        chrom1 = from@chrom1,
        pos1 = from@pos1,
        chrom2 = from@chrom2,
        pos2 = from@pos2,
        strand1 = from@strand1,
        strand2 = from@strand2,
        SVtype = from@SVtype,
        sv_id = from@sv_id,
        stringsAsFactors = FALSE
    )
})

#' Class to store CNV data
#'
#' @param chrom (character): chromosome (also in Ensembl notation)
#' @param start (numeric): start position for the CN segment
#' @param end (numeric): end position for the CN segment
#' @param CN (numeric): integer total copy number (e.g. 2 for unaltered chromosomal regions)
#' @return an instance of the class 'CNVsegs' that contains CNV data. Required format by the function detect_chromothripsis
#' @export
CNVsegs <- setClass(
    "CNVsegs",
    representation(
        chrom = "character",
        start = "numeric",
        end = "numeric",
        total_cn = "numeric",
        numSegs = "numeric"
    )
)

setMethod("initialize", "CNVsegs", function(.Object, ...) {
    .Object <- callNextMethod()
    if (length(.Object@chrom) != length(.Object@start) || length(.Object@chrom) != length(.Object@total_cn) || length(.Object@chrom) != length(.Object@end)) {
        stop("slots lengths are not all equal")
    }
    .Object@numSegs <- length(.Object@chrom)
    ind.match.cnv <- match(.Object@chrom, chromNames)
    if (sum(is.na(ind.match.cnv)) != 0) {
        stop(paste("chromosome name must be in \"", paste(chromNames, collapse = " ")), "\"")
    }

    .Object
})

setMethod("show", "CNVsegs", function(object) {
    print(paste("CNVsegs with", object@numSegs, "segments"))
})

setAs("CNVsegs", "data.frame", function(from, to) {
    to <- data.frame(chrom = from@chrom, start = from@start, end = from@end, total_cn = from@total_cn, stringsAsFactors = FALSE)
})


#' Class to store chromothripsis detection results
#'
#' @slot chromSummary A data.frame summarizing results per chromosome
#' @slot detail A list containing detailed statistical results
#' @export
chromoth <- setClass(
    "chromoth",
    representation(
        chromSummary = "data.frame",
        detail = "list"
    )
)

setMethod("show", "chromoth", function(object) {
    print(object@chromSummary)
})

setAs("chromoth", "data.frame", function(from, to) {
    to <- from@chromSummary
})

## if chromNames is not NULL, return the maximum number of cluster size in each chromsome specified in chromNames
cluster.SV <- function(SV.sample, min.Size = 1, chromNames) {
    SV.df <- as.data.frame(SV.sample)
    SV.df$chromothEvent <- rep(0, nrow(SV.df))

    # Ensure required packages are available
    if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
        stop("Package 'GenomicRanges' needed for this function to work. Please install it.", call. = FALSE)
    }
    if (!requireNamespace("S4Vectors", quietly = TRUE)) {
        stop("Package 'S4Vectors' needed for this function to work. Please install it.", call. = FALSE)
    }
    if (!requireNamespace("IRanges", quietly = TRUE)) {
        stop("Package 'IRanges' needed for this function to work. Please install it.", call. = FALSE)
    }

    chromNames <- as.character(chromNames)
    numSV.byChrom <- data.frame(chrom = chromNames, numbSV = rep(0, length(chromNames)), stringsAsFactors = FALSE)
    maxCluster.byChrom <- data.frame(chrom = chromNames, clusterSize = rep(0, length(chromNames)), stringsAsFactors = FALSE)

    ind.chr <- which(SV.df$chrom1 == SV.df$chrom2)
    if (length(ind.chr) < 1) {
        rt.value <- list(
            SV = SV.df, graph = NULL, connComp = list(), num.chromth = 0,
            # chromothripsis=list(),
            # chromothripsis.chr=c(),
            maxSVs = 0, degree = list()
        )
        rt.value$numSVByChrom <- numSV.byChrom
        rt.value$maxClusterSize <- maxCluster.byChrom
        return(rt.value)
    }
    tmp.numSV.byChrom <- table(SV.df$chrom1[ind.chr])
    ind.numSV.byChrom <- match(names(tmp.numSV.byChrom), numSV.byChrom$chrom)
    numSV.byChrom$numbSV[ind.numSV.byChrom[!is.na(ind.numSV.byChrom)]] <- tmp.numSV.byChrom[!is.na(ind.numSV.byChrom)]


    gr.chr <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(SV.df$chrom1[ind.chr]), ranges = IRanges::IRanges(SV.df$pos1[ind.chr] - 1, SV.df$pos2[ind.chr] + 1), strand = S4Vectors::Rle("+", length(ind.chr)))
    ovlp.chr <- as.data.frame(GenomicRanges::findOverlaps(gr.chr, gr.chr))
    # ovlp.chr = data.frame(ovlp.chr[[1]],ovlp.chr[[2]])
    ind.ovlp.chr <- which(ovlp.chr[, 1] < ovlp.chr[, 2])
    # to remove redundancies: e.g. 1 -2 ; 2 -1 both indicate that the SVs 1 and 2 overlap
    ovlp.chr <- ovlp.chr[ind.ovlp.chr, ]


    ovlp.chr.witin <- as.data.frame(GenomicRanges::findOverlaps(gr.chr, gr.chr, type = "within"))
    ind.rm <- which(ovlp.chr.witin[, 1] == ovlp.chr.witin[, 2])
    if (length(ind.rm) > 0) ovlp.chr.witin <- ovlp.chr.witin[-ind.rm, ]
    # ovlp.chr.witin = ovlp.chr.witin[-ind.rm,]
    if (nrow(ovlp.chr.witin) > 0 & nrow(ovlp.chr) > 0) {
        ovlp.chr.witin1 <- ovlp.chr.witin
        ind.tmp.within <- which(ovlp.chr.witin[, 2] < ovlp.chr.witin[, 1])
        if (length(ind.tmp.within) > 0) {
            ovlp.chr.witin[ind.tmp.within, 1] <- ovlp.chr.witin1[ind.tmp.within, 2]
            ovlp.chr.witin[ind.tmp.within, 2] <- ovlp.chr.witin1[ind.tmp.within, 1]
        }

        # Robust removal using string keys instead of IRanges on indices
        keys_all <- paste(ovlp.chr[, 1], ovlp.chr[, 2], sep = "_")
        keys_within <- paste(ovlp.chr.witin[, 1], ovlp.chr.witin[, 2], sep = "_")

        keep_idx <- which(!keys_all %in% keys_within)
        ovlp.chr <- ovlp.chr[keep_idx, , drop = FALSE]
    }
    # Use sequential local indices for graph nodes to ensure total consistency
    n_intra <- length(ind.chr)
    node_indices <- 1:n_intra

    adjMatrix <- matrix(0, nrow = n_intra, ncol = n_intra)
    if (nrow(ovlp.chr) > 0) {
        adjMatrix[ovlp.chr[, 1] + (ovlp.chr[, 2] - 1) * n_intra] <- 1
        adjMatrix[ovlp.chr[, 2] + (ovlp.chr[, 1] - 1) * n_intra] <- 1
    }
    rownames(adjMatrix) <- node_indices
    colnames(adjMatrix) <- node_indices
    g1 <- graphAM(adjMat = adjMatrix)
    gn <- as(g1, "graphNEL")
    cmpnt <- connComp(gn) ## finds the clusters

    ind.LargComp <- which(sapply(cmpnt, length) >= min.Size)
    chromothripsis.Reg <- cmpnt[ind.LargComp]

    degree.chromtheripsis <- list()
    maxSVs <- 0
    if (length(ind.LargComp) > 0) {
        for (i in 1:length(chromothripsis.Reg)) {
            gn.sub <- subGraph(chromothripsis.Reg[[i]], gn)
            degree.chromtheripsis[[i]] <- degree(gn.sub)
        }
        # Correctly map cluster membership back to the SV dataframe using local indices
        tmp <- as.numeric(unlist(chromothripsis.Reg))
        SV.df$chromothEvent[tmp] <- 1
        maxSVs <- max(sapply(chromothripsis.Reg, length))
    }

    rt.value <- list(
        SV = SV.df, graph = gn, connComp = cmpnt, num.chromth = length(ind.LargComp),
        # chromothripsis = chromothripsis.Reg,
        # chromothripsis.chr=chromothripsis.Reg.chr,
        maxSVs = maxSVs, degree = degree.chromtheripsis
    )
    rt.value$numSVByChrom <- numSV.byChrom

    tmp.chr <- as.character(sapply(cmpnt, FUN = function(v) {
        SV.df$chrom1[as.numeric(v[1])]
    }))
    tmp.Size <- sapply(cmpnt, length)
    tmp.maxSize <- aggregate(tmp.Size, by = list(tmp.chr), max)
    ind.match <- match(chromNames, tmp.maxSize[, 1])
    maxCluster.byChrom[, 2][!is.na(ind.match)] <- tmp.maxSize[ind.match[!is.na(ind.match)], 2]
    maxCluster.byChrom[maxCluster.byChrom[, 2] < min.Size, 2] <- 0
    rt.value$maxClusterSize <- maxCluster.byChrom
    # print(rt.value)
    return(rt.value)
}

#' Detect chromothripsis events from structural variation data
#' Identifies clusters of interleaved SVs and calculates statistical metrics for each chromosome (chromosomes 1-22 and X)
#' @param SV.sample an instance of class SVs
#' @param seg.sample an instance of class CNVsegs
#' @param min.Size minimum number of interleaved SVs required to report a cluster. Default is 1
#' @param genome reference genome (hg19 or hg38)
#' @export
detect_chromothripsis <- function(SV.sample, seg.sample, min.Size = 1, genome = "hg19") {
    cat("Running..\n\n\n")
    if (!is(SV.sample, "SVs")) {
        stop("SV.sample must be a SVs object")
    }
    if (!missing(seg.sample)) {
        if (!is(seg.sample, "CNVsegs")) {
            stop("seg.sample must be a CNVsegs object")
        }
    }
    if (!(as.character(genome) %in% c("hg19", "hg38"))) {
        stop("Reference genome assembly is not supported (Use hg19 or hg38)")
    }

    SV.sample <- as(SV.sample, "data.frame")
    # check that the strand info is correct
    if (sum(!(SV.sample$strand1 %in% c("+", "-")) > 0)) {
        stop("Error in the strand1 column. The strand values can only be + or -")
    }
    if (sum(!(SV.sample$strand2 %in% c("+", "-")) > 0)) {
        stop("Error in the strand2 column. The strand values can only be + or -")
    }

    # check SVtype-strand consistency
    if (nrow(SV.sample) > 0 & "SVtype" %in% names(SV.sample)) {
        # For intrachromosomal SVs, check consistency
        intra_idx <- which(SV.sample$chrom1 == SV.sample$chrom2)
        if (length(intra_idx) > 0) {
            for (i in intra_idx) {
                svtype <- SV.sample$SVtype[i]
                strand1 <- SV.sample$strand1[i]
                strand2 <- SV.sample$strand2[i]
                strand_combo <- paste0(strand1, strand2)

                # Skip empty SVtype or TRA (translocations can have any strand combination)
                if (svtype == "" | svtype == "TRA") next

                # Check consistency
                is_consistent <- FALSE
                if (svtype == "DEL" & strand_combo == "+-") is_consistent <- TRUE
                if (svtype == "DUP" & strand_combo == "-+") is_consistent <- TRUE
                if (svtype == "h2hINV" & strand_combo == "++") is_consistent <- TRUE
                if (svtype == "t2tINV" & strand_combo == "--") is_consistent <- TRUE

                if (!is_consistent) {
                    # Auto-correct based on strands (the physical reality)
                    inferred <- .infer_sv_type_from_strands(
                        SV.sample$chrom1[i], SV.sample$chrom2[i],
                        strand1, strand2
                    )
                    SV.sample$SVtype[i] <- inferred
                    # No longer stopping, just proceed with the corrected type
                }
            }
        }
    }

    # Check CNV data quality and completeness
    if (!missing(seg.sample)) {
        seg.sample <- as(seg.sample, "data.frame")

        # Warn if CNV data is missing or insufficient
        if (nrow(seg.sample) == 0) {
            warning(
                "CNV data is empty. Copy number oscillation analysis will be skipped.\n",
                "Chromothripsis detection accuracy may be reduced."
            )
        } else if (nrow(seg.sample) < 10) {
            warning(
                sprintf("Very few CNV segments (%d). CN oscillation detection may be unreliable.\n", nrow(seg.sample)),
                "Consider using a CNV caller with higher resolution."
            )
        }

        # Check for truly adjacent segments with same CN (which should be merged)
        # Only check segments that are genomically adjacent (end[i] == start[i+1])
        if (nrow(seg.sample) > 1) {
            adjacent_same_cn <- 0
            chrom_list <- unique(seg.sample$chrom)

            for (chr in chrom_list) {
                chr_cnv <- seg.sample[seg.sample$chrom == chr, ]
                if (nrow(chr_cnv) > 1) {
                    chr_cnv <- chr_cnv[order(chr_cnv$start), ]
                    for (i in 1:(nrow(chr_cnv) - 1)) {
                        # Check if segments are truly adjacent (no gap) AND have same CN
                        if (chr_cnv$end[i] >= chr_cnv$start[i + 1] &&
                            chr_cnv$total_cn[i] == chr_cnv$total_cn[i + 1]) {
                            adjacent_same_cn <- adjacent_same_cn + 1
                        }
                    }
                }
            }

            if (adjacent_same_cn > 0) {
                warning(
                    sprintf("Found %d pairs of truly adjacent CNV segments with identical copy numbers.\n", adjacent_same_cn),
                    "These should be merged. See README for merge code example.\n",
                    "Unmerged segments may affect CN oscillation detection."
                )
            }
        }
    } else {
        warning(
            "CNV data not provided. Copy number oscillation analysis will be skipped.\n",
            "Chromothripsis detection will rely solely on SV clustering patterns."
        )
        seg.sample <- data.frame(
            chrom = character(0), start = numeric(0),
            end = numeric(0), total_cn = numeric(0)
        )
    }

    chromothSample <- cluster.SV(SV.sample[SV.sample$chrom1 == SV.sample$chrom2, ], min.Size = min.Size, chromNames = chromNames) ## pass only intra
    # chromothSample$SV is already correctly set by cluster.SV with markings
    inter_df <- SV.sample[SV.sample$chrom1 != SV.sample$chrom2, ]
    # CRITICAL: Add chromothEvent column to match $SV structure for rbind compatibility
    if (nrow(inter_df) > 0) {
        inter_df$chromothEvent <- 0
    }
    chromothSample$SVinter <- inter_df

    chromSummary <- data.frame(chromothSample$maxClusterSize) # ,SVpvalue=chromothSample$SVpvalue[,2])
    chromothSample$CNV <- seg.sample

    out <- chromoth(chromSummary = chromSummary, detail = chromothSample)
    cat("Evaluating the statistical criteria\n")
    out@chromSummary <- suppressWarnings(statistical_criteria(out, genome))
    cat("Successfully finished!\n")
    return(out)
}
