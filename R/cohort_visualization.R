#' @include cohort_analysis.R
NULL

#' Plot mechanism frequency comparison between groups
#'
#' Creates a bar plot comparing the prevalence of chromoanagenesis mechanisms
#' across different clinical groups.
#'
#' @param cohort An OncoImplexusCohort object
#' @param group_by Column name in clinical data to group samples by
#' @param stringency Analysis stringency ("strict" or "inclusive")
#' @param show_pvalue Logical, whether to show p-values on the plot (default: TRUE)
#' @return A ggplot object
#' @export
plot_mechanism_comparison <- function(cohort, group_by, stringency = "strict", show_pvalue = TRUE) {
    
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' needed for this function to work. Please install it.")
    }
    
    # Get statistics using existing function (suppress output)
    capture.output({
        res <- compare_mechanisms(cohort, group_by, stringency)
    })
    
    data <- res$data
    stats <- res$stats
    
    # Reshape data for plotting
    # Calculate percentages manually to ensure correct grouping
    plot_data <- data.frame(
        Group = character(0),
        Mechanism = character(0),
        Count = numeric(0),
        Total = numeric(0),
        Percentage = numeric(0),
        stringsAsFactors = FALSE
    )
    
    mechanisms <- c("Chromothripsis", "Chromoplexy", "Chromoanasynthesis")
    groups <- unique(data$Group)
    
    for (g in groups) {
        for (m in mechanisms) {
            sub_data <- data[data$Group == g, ]
            count <- sum(sub_data[[m]])
            total <- nrow(sub_data)
            pct <- (count / total) * 100
            
            plot_data[nrow(plot_data) + 1, ] <- list(
                Group = g,
                Mechanism = m,
                Count = count,
                Total = total,
                Percentage = pct
            )
        }
    }
    
    # Merge p-values
    plot_data$Significance <- ""
    for (i in 1:nrow(plot_data)) {
        m <- plot_data$Mechanism[i]
        sig <- stats$Significant[stats$Mechanism == m]
        if (length(sig) > 0) plot_data$Significance[i] <- sig
    }
    
    # Create Plot
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Group, y = Percentage, fill = Mechanism)) +
        ggplot2::geom_bar(stat = "identity", position = "dodge") +
        ggplot2::theme_minimal() +
        ggplot2::labs(
            title = paste("Prevalence of Chromoanagenesis Mechanisms by", group_by),
            subtitle = paste("Stringency:", tools::toTitleCase(stringency)),
            y = "Prevalence (%)",
            x = group_by
        ) +
        ggplot2::scale_fill_brewer(palette = "Set2") +
        ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 14),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            legend.position = "bottom"
        )
    
    # Add counts on top of bars
    p <- p + ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%d/%d", Count, Total)),
        position = ggplot2::position_dodge(width = 0.9),
        vjust = -0.5,
        size = 3
    )
    
    # Add statistical significance (if p-value requested)
    if (show_pvalue) {
        # Create a separate data frame for annotations to avoid overplotting
        # We put one annotation per mechanism per facet (if faceting)
        # Here we just add text to the plot title or subtitle
        sig_text <- paste(stats$Mechanism, ": p =", format(stats$P_value, digits = 3), stats$Significant)
        p <- p + ggplot2::labs(caption = paste("Statistical Significance (Fisher/Chi-sq):\n", paste(sig_text, collapse = "\n")))
    }
    
    # Option: Facet by Mechanism for clearer comparison
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Group, y = Percentage, fill = Group)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::facet_wrap(~Mechanism) +
        ggplot2::theme_bw() +
        ggplot2::labs(
            title = paste("Prevalence of Chromoanagenesis by", group_by),
            subtitle = paste("Stringency:", stringency),
            y = "Prevalence (%)",
            x = NULL
        ) +
        ggplot2::geom_text(
            ggplot2::aes(label = sprintf("%.1f%%\n(n=%d)", Percentage, Count)),
            vjust = -0.2,
            size = 3
        ) +
        ggplot2::scale_y_continuous(limits = c(0, 115)) + # Add space for labels
        ggplot2::theme(
            legend.position = "none",
            strip.text = ggplot2::element_text(face = "bold", size = 12)
        )
        
    # Add significance label to each facet
    # We need a dataframe for geom_text with facet information
    ann_text <- data.frame(
        Mechanism = stats$Mechanism,
        Group = groups[1], # Dummy x position
        Percentage = 110,  # y position
        Label = ifelse(is.na(stats$P_value), "", 
                      paste0("p=", format(stats$P_value, digits=2), " ", stats$Significant)),
        stringsAsFactors = FALSE
    )
    
    # Filter out mechanisms not in plot_data just in case
    ann_text <- ann_text[ann_text$Mechanism %in% plot_data$Mechanism, ]
    
    p <- p + ggplot2::geom_text(
        data = ann_text,
        ggplot2::aes(x = Inf, y = Inf, label = Label, fill = NULL),
        hjust = 1.1, vjust = 1.5,
        size = 3.5, fontface = "italic", color = "red"
    )

    return(p)
}

#' Plot cohort landscape (OncoPrint style)
#'
#' Creates a heatmap visualizing the presence of chromoanagenesis mechanisms
#' across all samples in the cohort, grouped by clinical variables.
#'
#' @param cohort An OncoImplexusCohort object
#' @param group_by Column name in clinical data to group samples by
#' @return A ggplot object
#' @export
plot_cohort_landscape <- function(cohort, group_by) {
    
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' needed.")
    }
    
    # Prepare data matrix
    # We want 3 rows per sample: Chromothripsis, Chromoplexy, Chromoanasynthesis
    # And columns for status: High, Low, Likely, Possible, None
    
    plot_data <- data.frame(
        Sample = character(0),
        Group = character(0),
        Mechanism = character(0),
        Status = character(0),
        Score = numeric(0), # For ordering
        stringsAsFactors = FALSE
    )
    
    groups <- cohort@clinical[[group_by]]
    names(groups) <- cohort@sample_ids
    
    for (sid in cohort@sample_ids) {
        res <- cohort@results[[sid]]
        grp <- groups[sid]
        
        # Chromothripsis
        status_ct <- "None"
        score_ct <- 0
        if (!is.null(res$chromothripsis)) {
            if (res$chromothripsis$n_high_confidence > 0) {
                status_ct <- "High confidence"
                score_ct <- 2
            } else if (res$chromothripsis$n_low_confidence > 0) {
                status_ct <- "Low confidence"
                score_ct <- 1
            }
        }
        
        # Chromoplexy
        status_cp <- "None"
        score_cp <- 0
        if (!is.null(res$chromoplexy)) {
            if (res$chromoplexy$likely_chromoplexy > 0) {
                status_cp <- "Likely"
                score_cp <- 2
            } else if (res$chromoplexy$possible_chromoplexy > 0) {
                status_cp <- "Possible"
                score_cp <- 1
            }
        }
        
        # Chromoanasynthesis
        status_cs <- "None"
        score_cs <- 0
        if (!is.null(res$chromoanasynthesis)) {
            if (res$chromoanasynthesis$likely_chromoanasynthesis > 0) {
                status_cs <- "Likely"
                score_cs <- 2
            } else if (res$chromoanasynthesis$possible_chromoanasynthesis > 0) {
                status_cs <- "Possible"
                score_cs <- 1
            }
        }
        
        plot_data <- rbind(plot_data, data.frame(
            Sample = sid, Group = grp, Mechanism = "Chromothripsis", Status = status_ct, Score = score_ct
        ))
        plot_data <- rbind(plot_data, data.frame(
            Sample = sid, Group = grp, Mechanism = "Chromoplexy", Status = status_cp, Score = score_cp
        ))
        plot_data <- rbind(plot_data, data.frame(
            Sample = sid, Group = grp, Mechanism = "Chromoanasynthesis", Status = status_cs, Score = score_cs
        ))
    }
    
    # Ordering samples:
    # 1. By Group
    # 2. By presence of mechanisms (High > Low > None)
    
    # Calculate sample score sum
    sample_scores <- aggregate(Score ~ Sample, plot_data, sum)
    sample_order <- order(groups[sample_scores$Sample], -sample_scores$Score)
    ordered_samples <- sample_scores$Sample[sample_order]
    
    plot_data$Sample <- factor(plot_data$Sample, levels = ordered_samples)
    plot_data$Mechanism <- factor(plot_data$Mechanism, 
                                 levels = c("Chromoanasynthesis", "Chromoplexy", "Chromothripsis")) # Reverse order for plotting
    
    # Remove NA groups
    plot_data <- plot_data[!is.na(plot_data$Group), ]
    
    # Define colors
    status_colors <- c(
        "High confidence" = "#E41A1C", # Red
        "Likely" = "#E41A1C",
        "Low confidence" = "#FF7F00",  # Orange
        "Possible" = "#FF7F00",
        "None" = "#F0F0F0"             # Light Gray
    )
    
    # Plot
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Sample, y = Mechanism, fill = Status)) +
        ggplot2::geom_tile(color = "white", size = 0.2) +
        ggplot2::scale_fill_manual(values = status_colors) +
        ggplot2::facet_grid(~Group, scales = "free_x", space = "free_x") +
        ggplot2::theme_minimal() +
        ggplot2::labs(
            title = "Cohort Landscape of Chromoanagenesis",
            x = NULL,
            y = NULL
        ) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_blank(), # Hide sample names if too many
            panel.grid = ggplot2::element_blank(),
            strip.background = ggplot2::element_rect(fill = "gray90", color = NA),
            strip.text = ggplot2::element_text(face = "bold")
        )
        
    return(p)
}
