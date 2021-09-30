
#' Quantitative plot of transcripts
#'
#' @param splice_df same as the input to draw_splice_picture, a data frame of
#' transcripts. Must contain columns named variant, score and Transcript_id.
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' NA
plot_quant <- function(splice_df){
    ggplot2::ggplot(splice_df, ggplot2::aes(x=variant, y=score, label=Transcript_id)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color="grey") +
    ggplot2::geom_point(size=4) +
    #ggtitle(gene_name)+
    ggrepel::geom_text_repel()
}

#' Draw a picture of splice patterns
#'
#' Creates plot showing splice variants
#'
#' @param splice_data dataframe containing columns named variant, strand, score,
#' Transcript_id and splice_pattern. Intended to be the output of parse_nexons_gtf
#' function.
#' @param title_text optional title for plot (probably gene name)
#' @param quant TRUE or FALSE whether to show a quantitative plot of transcipt
#' counts alongside the splice variants
#'
#' @return
#' @export
#'
#' @examples
#' NA
#' @import ggplot2
draw_splice_picture <- function(splice_data, title_text = "", quant = FALSE) {

  splice_plot_data <- add_exon_loci(splice_data)

  splice_segment_data <- splice_plot_data %>%
    dplyr::group_by(variant) %>%
    dplyr::summarise(start=min(start),end=max(end), strand = strand)

  if(all(splice_segment_data$strand == "+")){
    strand_colour <- "red3"
  } else if(all(splice_segment_data$strand == "-")) {
    strand_colour <- "blue3"
  } else {
    strand_colour <- c("blue3", "red3")
  }

  p <- splice_plot_data %>%
    ggplot(aes(xmin = start, xmax = end, ymin = variant+0.05,
               ymax = variant+0.95, fill = strand, colour = strand)) +
    geom_segment(data=splice_segment_data,
                 aes(x=start, xend=end, y=variant+0.5, yend=variant+0.5))+
    geom_rect()+#fill=strand,colour=shade) +
    scale_y_continuous(
      breaks=splice_plot_data$variant+0.5,
      labels = paste(splice_plot_data$Transcript_id)) +
    ylab("") +
    xlab("Genomic Position") +
    ggtitle(title_text) +
    scale_fill_manual(values = strand_colour) +
    scale_colour_manual(values = strand_colour)

  if(quant){
    quant_plot <- plot_quant(splice_data)
    return(invisible(gridExtra::grid.arrange(p, quant_plot, nrow=1)))
  } else {
    return(p)
  }
}
