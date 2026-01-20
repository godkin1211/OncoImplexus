#' Internal Monte Carlo KS test for Exponential Distribution
#' Corrects for parameter estimation bias (Lilliefors-style)
#' @keywords internal
.ks_test_exp_mc <- function(x, n_sim = 1000) {
    n <- length(x)
    if (n < 5) {
        return(NA)
    }

    # 1. Observed D statistic with estimated rate
    rate_hat <- 1 / mean(x)
    d_obs <- ks.test(x, "pexp", rate = rate_hat)$statistic

    # 2 & 3. Simulations (using rate=1 due to scale invariance)
    d_sims <- replicate(n_sim, {
        x_sim <- rexp(n, rate = 1)
        rate_hat_sim <- 1 / mean(x_sim)
        ks.test(x_sim, "pexp", rate = rate_hat_sim)$statistic
    })

    # 4. P-value calculation with +1 correction
    p_val <- (sum(d_sims >= d_obs) + 1) / (n_sim + 1)
    return(as.numeric(p_val))
}

statistical_criteria <- function(input, genome) {
    if (!requireNamespace("MASS", quietly = TRUE)) {
        stop("Package 'MASS' needed for this function to work. Please install it.", call. = FALSE)
    }

    summary <- data.frame(
        chrom = input@chromSummary$chrom,
        clusterSize = input@chromSummary$clusterSize
    )
    candidate_chrs <- input@chromSummary$chrom
    l_candidate_chrs <- length(candidate_chrs)
    # contains the IDs for the SVs that are in a cluster.
    cluster_sizes <- sapply(input@detail$connComp, length)
    summary$start <- rep(NA, l_candidate_chrs)
    summary$end <- rep(NA, l_candidate_chrs)
    summary$number_DEL <- rep(0, l_candidate_chrs)
    summary$number_DUP <- rep(0, l_candidate_chrs)
    summary$number_h2hINV <- rep(0, l_candidate_chrs)
    summary$number_t2tINV <- rep(0, l_candidate_chrs)
    summary$number_TRA <- rep(0, l_candidate_chrs)
    summary$clusterSize_including_TRA <- rep(0, l_candidate_chrs)
    summary$number_SVs_sample <- rep(0, l_candidate_chrs)
    summary$number_CNV_segments <- rep(NA, l_candidate_chrs)
    summary$pval_fragment_joins <- rep(NA, l_candidate_chrs)
    summary$chr_breakpoint_enrichment <- rep(NA, l_candidate_chrs)
    summary$pval_exp_chr <- rep(NA, l_candidate_chrs)
    summary$pval_exp_cluster <- rep(NA, l_candidate_chrs)
    # oscillating
    summary$max_number_oscillating_CN_segments_2_states <- rep(NA, l_candidate_chrs)
    summary$max_number_oscillating_CN_segments_3_states <- rep(NA, l_candidate_chrs)
    summary$number_CN_segments_chr <- rep(NA, l_candidate_chrs)
    summary$max_number_oscillating_CN_segments_2_states_chr <- rep(NA, l_candidate_chrs)
    summary$max_number_oscillating_CN_segments_3_states_chr <- rep(NA, l_candidate_chrs)
    #---------------------------------------------------------------------------
    #---------------------------------------------------------------------------
    l_candidate_chrs <- length(candidate_chrs)
    summary_inter <- data.frame( # main_chrom=rep(0,l_candidate_chrs),
        number_DEL = rep(0, l_candidate_chrs),
        number_h2hINV = rep(0, l_candidate_chrs),
        number_t2tINV = rep(0, l_candidate_chrs),
        # number_TRA=rep(0,l_candidate_chrs),
        number_DUP = rep(0, l_candidate_chrs),
        pval_fragment_joins = rep(NA, l_candidate_chrs)
    )

    other_chroms <- rep("", length(candidate_chrs))
    other_chroms_coords_all <- rep(" ", length(candidate_chrs))

    # iterate over candidate chromosomes
    for (cand in candidate_chrs) {
        separator_index <- 1
        cand_clust_size <- input@chromSummary$clusterSize[input@chromSummary$chrom == cand]
        idx <- which(cluster_sizes == cand_clust_size)
        SVsnow <- input@detail$SV[as.numeric(unlist(input@detail$connComp[idx])), ]
        index_chromosome <- which(summary$chrom == cand)
        SVsnow <- SVsnow[SVsnow$chrom1 == cand, ] # remove if there are more
        # get the copy number data
        CNVsnow <- input@detail$CNV
        # get CNV data for the current chromosome
        CNVsnow <- CNVsnow[which(CNVsnow$chrom == cand), ]

        #-------------------------------------------
        # CN oscillations at the chromosome level
        #-------------------------------------------
        if (nrow(CNVsnow) >= 4) {
            k <- nrow(CNVsnow)
            # Round CN to handle floating point noise (e.g., 2.01 -> 2)
            cn_vals <- round(CNVsnow$total_cn)
            v <- (cn_vals[1:(k - 2)] - cn_vals[3:k])
            v_3states <- 1 * (abs(v) %in% c(0, 1))
            v <- 1 * (v == 0)
            # i_i2 = fmax(1,v)
            number_CN_segments <- length(v)
            max_number_oscillating_CN_segments_2_states <- fmaxmax(1, v)
            max_number_oscillating_CN_segments_3_states <- fmaxmax(1, v_3states)
            summary$number_CN_segments_chr[index_chromosome] <- number_CN_segments
            summary$max_number_oscillating_CN_segments_2_states_chr[index_chromosome] <- max_number_oscillating_CN_segments_2_states + 2
            summary$max_number_oscillating_CN_segments_3_states_chr[index_chromosome] <- max_number_oscillating_CN_segments_3_states + 2
        }

        if (nrow(SVsnow) != 0) {
            # get coordinates for the cluster of interleaved SVs
            min_now <- as.numeric(min(c(SVsnow$pos1, SVsnow$pos2)))
            max_now <- max(c(SVsnow$pos1, SVsnow$pos2))
            summary$start[index_chromosome] <- min_now
            summary$end[index_chromosome] <- max_now

            if (nrow(CNVsnow) > 0) {
                idxa <- which(CNVsnow$start <= summary$start[index_chromosome])
                idxb <- which(CNVsnow$end >= summary$end[index_chromosome])
                if (length(idxa) == 0) {
                    tmp_idx <- which(CNVsnow$start >= summary$start[index_chromosome])
                    if (length(tmp_idx) > 0) idxa <- min(tmp_idx) else idxa <- numeric(0)
                }
                if (length(idxb) == 0) {
                    tmp_idx <- which(CNVsnow$end <= summary$end[index_chromosome])
                    if (length(tmp_idx) > 0) idxb <- max(tmp_idx) else idxb <- numeric(0)
                }

                if (length(idxa) != 0 && length(idxb) != 0 && max(idxa) <= min(idxb)) {
                    CNVsnow <- CNVsnow[seq(max(idxa), min(idxb), 1), ]
                    summary$number_CNV_segments[index_chromosome] <- nrow(CNVsnow)
                    i_i2_sequential <- 0
                    # CN oscillations in the cluster region
                    if (nrow(CNVsnow) >= 4) {
                        k <- nrow(CNVsnow)
                        # Round CN to handle floating point noise
                        cn_vals <- round(CNVsnow$total_cn)
                        v <- (cn_vals[1:(k - 2)] - cn_vals[3:k])
                        v_3states <- 1 * (abs(v) %in% c(0, 1))
                        v <- 1 * (v == 0)
                        max_number_oscillating_CN_segments_3_states <- fmaxmax(1, v_3states)
                        number_CN_segments <- length(v)
                        max_number_oscillating_CN_segments_2_states <- fmaxmax(1, v)

                        summary$max_number_oscillating_CN_segments_3_states[index_chromosome] <- max_number_oscillating_CN_segments_3_states + 2
                        summary$max_number_oscillating_CN_segments_2_states[index_chromosome] <- max_number_oscillating_CN_segments_2_states + 2
                    }
                }
            }
        }

        #-----------------------------------------------------
        # Consider interchromosomal SVs
        #-----------------------------------------------------
        inter <- input@detail$SVinter
        SVs_all <- rbind(input@detail$SV, inter)

        if (nrow(SVsnow) != 0 & nrow(inter) > 0) {
            # get whether any tranlocation involves the cand chromo and is within the SV cluster previously identified
            min_now <- as.numeric(min(c(SVsnow$pos1, SVsnow$pos2)))
            max_now <- max(c(SVsnow$pos1, SVsnow$pos2))
            idx_inter1 <- which(inter$chrom1 == cand & inter$pos1 >= (min_now - 10000) & inter$pos1 <= (max_now + 10000))
            idx_inter2 <- which(inter$chrom2 == cand & inter$pos2 >= (min_now - 10000) & inter$pos2 <= (max_now + 10000))
            idx_inter <- c(idx_inter1, idx_inter2)
            window <- c() # save the SVs in the window comprising multiple chrs
            CNV_window <- c()
            selection_chrs <- c()
            nb_TRA_tot <- c()
            selection_chr_coords <- c()

            if (length(idx_inter) >= 1) { # i.e., there are interchr SVs mapped to the SV cluster
                inter <- inter[idx_inter, ]
                cand_inter <- unique(c(inter$chrom1, inter$chrom2))
                cand_inter <- cand_inter[which(cand_inter != cand)]

                for (chr_inter in sort(cand_inter)) {
                    # for each of the cand_inter, look whether the translocation connects to SV clusters
                    cand_clust_size_window <- input@chromSummary$clusterSize[input@chromSummary$chrom == chr_inter]
                    idx <- which(cluster_sizes == cand_clust_size_window)
                    SVsnow_window <- input@detail$SV[as.numeric(unlist(input@detail$connComp[idx])), ]
                    coords_window <- paste(SVsnow_window$chrom1, SVsnow_window$pos1, SVsnow_window$chrom2, SVsnow_window$pos2)
                    coords_window_from_SVs_all <- paste(SVs_all$chrom1, SVs_all$pos1, SVs_all$chrom2, SVs_all$pos2)

                    idx <- which(coords_window_from_SVs_all %in% coords_window)
                    SVsnow_window <- SVs_all[idx, ]
                    min_now_window <- as.numeric(min(c(SVsnow_window$pos1[which(SVsnow_window$chrom1 == chr_inter)], SVsnow_window$pos2[which(SVsnow_window$chrom2 == chr_inter)])))
                    max_now_window <- max(c(SVsnow_window$pos1[which(SVsnow_window$chrom1 == chr_inter)], SVsnow_window$pos2[which(SVsnow_window$chrom2 == chr_inter)]))
                    # check if the tranlocations are mapped to the SV cluster
                    inter_chr_inter1 <- inter[which(inter$chrom1 == chr_inter), 1:2]
                    inter_chr_inter2 <- inter[which(inter$chrom2 == chr_inter), 3:4]
                    names(inter_chr_inter1) <- c("chr", "pos")
                    names(inter_chr_inter2) <- c("chr", "pos")
                    inter_chr_inter <- rbind(inter_chr_inter1, inter_chr_inter2)
                    min_trans <- min(c(inter_chr_inter$pos))
                    max_trans <- max(c(inter_chr_inter$pos))

                    if ((min_trans >= (min_now_window - 10000) | max_trans <= (max_now_window + 10000)) & nrow(inter_chr_inter) >= 2) {
                        CNV_window_now <- input@detail$CNV
                        idxa <- which(CNV_window_now$chrom == chr_inter & CNV_window_now$start <= (min_now_window - 10000))
                        if (length(idxa) == 0) {
                            idxa <- which(CNV_window_now$chrom == chr_inter & CNV_window_now$start <= (min_now_window))
                        }
                        idxb <- which(CNV_window_now$chrom == chr_inter & CNV_window_now$end >= (max_now_window + 10000))
                        if (length(idxb) == 0) {
                            idxb <- which(CNV_window_now$chrom == chr_inter & CNV_window_now$end >= (max_now_window + 10000))
                        }
                        if (length(idxb) == 0) {
                            idxb <- max(which(CNV_window_now$chrom == chr_inter & CNV_window_now$end <= (max_now_window)))
                        }

                        if (length(idxb) != 0 & length(idxa) != 0) {
                            nb_oscill <- fmax2(CNV_window_now[seq(max(idxa), min(idxb), 1), ], cutoff = 1)
                            if (nb_oscill >= 0) {
                                CNV_window <- rbind(CNV_window, CNV_window_now[seq(max(idxa), min(idxb), 1), ])
                                selection_chrs <- c(selection_chrs, chr_inter)
                                if (separator_index %% 2 == 0) {
                                    separator <- "\n"
                                } else {
                                    separator <- ";"
                                }
                                separator_index <- separator_index + 1
                                selection_chr_coords <- c(selection_chr_coords, paste(chr_inter, ":", min_now_window, "-", max_now_window, separator, sep = ""))
                                window <- rbind(window, SVsnow_window)
                            }
                        }
                    }
                }
                if (!is.null(window)) {
                    names(inter)[1:4] <- c("chrom1", "pos1", "chrom2", "pos2")
                    window <- rbind(window, inter[which(inter$chrom1 %in% selection_chrs | inter$chrom2 %in% selection_chrs), ])
                    #---------------------------------------------
                    # Randomness of joins
                    #---------------------------------------------
                    # check multinomial distribution for SV types
                    # add SV types for translocations
                    window$SVtype2 <- window$SVtype
                    nb_TRA <- length(which(window$SVtype == "TRA"))
                    window$SVtype2[which(window$SVtype == "TRA" & window$strand1 == "+" & window$strand2 == "-")] <- "DEL"
                    window$SVtype2[which(window$SVtype == "TRA" & window$strand1 == "+" & window$strand2 == "+")] <- "h2hINV"
                    window$SVtype2[which(window$SVtype == "TRA" & window$strand1 == "-" & window$strand2 == "-")] <- "t2tINV"
                    window$SVtype2[which(window$SVtype == "TRA" & window$strand1 == "-" & window$strand2 == "+")] <- "DUP"
                    obs <- c(0, 0, 0, 0)
                    names(obs) <- c("del", "inv", "inv2", "dup")
                    obs[1] <- sum(window$SVtype2 == "DEL")
                    obs[2] <- sum(window$SVtype2 %in% c("h2hINV"))
                    obs[3] <- sum(window$SVtype2 %in% c("t2tINV"))
                    obs[4] <- sum(window$SVtype2 == "DUP")
                    idxx <- which(candidate_chrs == cand)
                    summary_inter$number_DEL[idxx] <- obs[1]
                    summary_inter$number_h2hINV[idxx] <- obs[2]
                    summary_inter$number_t2tINV[idxx] <- obs[3]
                    summary_inter$number_DUP[idxx] <- obs[4]
                    # summary_inter$number_TRA[idxx] = nb_TRA
                    other_chroms[idxx] <- as.character(paste(as.vector(selection_chrs), collapse = "_"))
                    other_chroms_coords_all[idxx] <- as.character(paste(as.vector(selection_chr_coords), collapse = ""))
                    # Use simulated p-value for robustness, especially with low counts
                    signif <- chisq.test(obs, p = rep(1 / 4, 4), simulate.p.value = TRUE, B = 2000)$p.value
                    summary_inter$pval_fragment_joins[idxx] <- signif
                }
            }
        }

        summary_inter$other_chroms <- other_chroms
        summary_inter$other_chroms_coords_all <- other_chroms_coords_all

        #-----------------------------------------------------
        # Exponential distribution (CLUSTER SVS)
        #-----------------------------------------------------
        # The inter-breakpoint distances should follow an exponential distribution
        # under the null hypothesis of random breakage (Poisson process).
        # We test observed distances against the theoretical CDF.
        SVsnow_exp <- inter
        minnow <- summary$start[index_chromosome]
        maxnow <- summary$end[index_chromosome]

        idx_inter1 <- which(SVsnow_exp$chrom1 == cand & SVsnow_exp$pos1 >= (minnow) & SVsnow_exp$pos1 <= (maxnow))
        breaks_inter1 <- SVsnow_exp$pos1[idx_inter1]
        idx_inter2 <- which(SVsnow_exp$chrom2 == cand & SVsnow_exp$pos2 >= (minnow) & SVsnow_exp$pos2 <= (maxnow))
        breaks_inter2 <- SVsnow_exp$pos2[idx_inter2]

        SVsnow_exp <- SVsnow_exp[which(SVsnow_exp$chrom1 == SVsnow_exp$chrom2), ]
        cand_clust_size <- input@chromSummary$clusterSize[input@chromSummary$chrom == cand]
        idx <- which(cluster_sizes == cand_clust_size)
        SVsnow_exp <- input@detail$SV[as.numeric(unlist(input@detail$connComp[idx])), ]
        SVsnow_exp <- SVsnow_exp[SVsnow_exp$chrom1 == cand, ]

        breaks <- sort(unique(c(SVsnow_exp$pos1, SVsnow_exp$pos2, breaks_inter1, breaks_inter2)))

        if (length(breaks) >= 6) {
            # 1. Calculate observed distances
            distances <- as.numeric(diff(breaks))
            # Filter out zero distances from microhomology if any
            distances <- distances[distances > 0]

            if (length(distances) >= 5) {
                # Use robust Monte Carlo KS test to handle parameter estimation bias
                pval_exp <- .ks_test_exp_mc(distances, n_sim = 1000)
                summary$pval_exp_cluster[index_chromosome] <- pval_exp
            } else {
                summary$pval_exp_cluster[index_chromosome] <- NA
            }
        } else {
            summary$pval_exp_cluster[index_chromosome] <- NA
        }

        #-----------------------------------------------------
        # Exponential distribution to test clustering of breakpoints in a given chromosome
        #-----------------------------------------------------
        SVsnow_exp <- inter

        idx_inter1 <- which(SVsnow_exp$chrom1 == cand)
        breaks1 <- SVsnow_exp$pos1[idx_inter1]
        idx_inter2 <- which(SVsnow_exp$chrom2 == cand)
        breaks2 <- SVsnow_exp$pos2[idx_inter2]
        breaks <- sort(unique(c(breaks1, breaks2)))

        if (length(breaks) >= 6) {
            distances <- as.numeric(diff(breaks))
            distances <- distances[distances > 0]

            if (length(distances) >= 5) {
                pval_exp <- .ks_test_exp_mc(distances, n_sim = 1000)
                summary$pval_exp_chr[index_chromosome] <- pval_exp
            } else {
                summary$pval_exp_chr[index_chromosome] <- NA
            }
        } else {
            summary$pval_exp_chr[index_chromosome] <- NA
        }

        #-----------------------------------------------------
        # Are there more SVs in a chrs than expected by chance?
        #-----------------------------------------------------
        SVsnow_exp <- inter
        idx_inter1 <- which(SVsnow_exp$chrom1 == cand)
        idx_inter2 <- which(SVsnow_exp$chrom2 == cand)
        nb_SVs_cand <- nrow(SVsnow_exp[unique(c(idx_inter1, idx_inter2)), ])

        if (nb_SVs_cand != 0) {
            nb_SVs_all_sample <- nrow(SVsnow_exp)

            # Robust mapping probability calculation
            target_mappa <- if (genome == "hg38") info_mappa_hg38 else info_mappa

            # Handle naming mismatch (chr1 vs 1)
            cand_clean <- gsub("^chr", "", cand)
            idx_mappa <- which(target_mappa$V1 == cand | target_mappa$V1 == cand_clean | target_mappa$V1 == paste0("chr", cand_clean))

            if (length(idx_mappa) > 0) {
                prob_cand <- target_mappa$tot[idx_mappa[1]] / sum(as.numeric(target_mappa$tot))
                chr_enrich <- binom.test(nb_SVs_cand, nb_SVs_all_sample, p = prob_cand)$p.value
                summary$chr_breakpoint_enrichment[index_chromosome] <- chr_enrich
            } else {
                summary$chr_breakpoint_enrichment[index_chromosome] <- NA
            }
        } else {
            summary$chr_breakpoint_enrichment[index_chromosome] <- NA
        }
        #--------------------------------------------------
        # Randomness of DNA fragment joins
        #--------------------------------------------------
        # check multinomial distribution for SV types
        cand_clust_size <- input@chromSummary$clusterSize[input@chromSummary$chrom == cand]
        idx <- which(cluster_sizes == cand_clust_size) # ojo, puede haber clusts del mismo size en chrs diferentes
        SVsnow <- input@detail$SV[as.numeric(unlist(input@detail$connComp[idx])), ]
        SVsnow <- SVsnow[SVsnow$chrom1 == cand, ] # remove if there are more

        obs <- c(0, 0, 0, 0)
        names(obs) <- c("del", "inv", "inv2", "dup")
        obs[1] <- sum(SVsnow$SVtype == "DEL")
        obs[2] <- sum(SVsnow$SVtype == "h2hINV")
        obs[3] <- sum(SVsnow$SVtype == "t2tINV")
        obs[4] <- sum(SVsnow$SVtype == "DUP")
        summary$number_SVs_sample[index_chromosome] <- nrow(SVs_all)
        obs2 <- obs
        if (nrow(inter) > 0) {
            # get whether any tranlocation involves the cand chromo and is within the cluster of SVs previously identified
            min_now <- as.numeric(min(c(SVsnow$pos1, SVsnow$pos2)))
            max_now <- max(c(SVsnow$pos1, SVsnow$pos2))
            idx_inter1 <- which(inter$chrom1 == cand & inter$pos1 >= (min_now) & inter$pos1 <= (max_now))
            idx_inter2 <- which(inter$chrom2 == cand & inter$pos2 >= (min_now) & inter$pos2 <= (max_now))
            idx_inter <- unique(c(idx_inter1, idx_inter2))
            inter <- inter[idx_inter, ]
            obs2[1] <- obs2[1] + length(which(inter$strand1 == "+" & inter$strand2 == "-"))
            obs2[2] <- obs2[2] + length(which(inter$strand1 == "+" & inter$strand2 == "+"))
            obs2[3] <- obs2[3] + length(which(inter$strand1 == "-" & inter$strand2 == "-"))
            obs2[4] <- obs2[4] + length(which(inter$strand1 == "-" & inter$strand2 == "+"))
            summary$number_TRA[index_chromosome] <- length(idx_inter)
        }

        summary$number_DEL[index_chromosome] <- obs[1]
        summary$number_h2hINV[index_chromosome] <- obs[2]
        summary$number_t2tINV[index_chromosome] <- obs[3]
        summary$number_DUP[index_chromosome] <- obs[4]
        summary$clusterSize_including_TRA[index_chromosome] <- sum(obs2)

        if (nrow(SVsnow) != 0 & sum(obs2) > 0) {
            # Use simulated p-value for robustness with small samples
            signif <- chisq.test(obs2, p = rep(1 / 4, 4), simulate.p.value = TRUE, B = 2000)$p.value
            summary$pval_fragment_joins[index_chromosome] <- signif
        } else {
            summary$pval_fragment_joins[index_chromosome] <- NA
        }
    }


    names(summary_inter) <- paste0("inter_", names(summary_inter))
    return(cbind(summary, summary_inter))
}
