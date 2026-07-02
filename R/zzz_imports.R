#' @importFrom grDevices adjustcolor
#' @importFrom graphics legend text title
#' @importFrom methods .hasSlot as callNextMethod is new
#' @importFrom stats aggregate ave binom.test chisq.test coef cor.test end fisher.test ks.test lm median na.omit p.adjust pchisq reorder rexp sd start
#' @importFrom utils capture.output head read.csv read.table write.csv
#' @importFrom data.table rbindlist
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_hline geom_curve geom_segment geom_rect theme theme_bw element_text element_blank element_line unit xlab ylab xlim ylim scale_x_continuous scale_y_continuous scale_colour_manual annotate ggplotGrob coord_cartesian
#' @importFrom grid grid.draw unit.c
#' @importFrom gridExtra grid.arrange tableGrob
#' @importFrom GenomicRanges GRanges mcols
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors Rle
NULL

if (getRversion() >= "2.15.1") {
    utils::globalVariables(c(
        "conf_high", "conf_low", "feature_label", "group",
        "hazard_ratio", "survival", "time"
    ))
}
