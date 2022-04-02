
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
    #ggplot2::geom_abline(slope = 1, intercept = 0, color="grey") +
    ggplot2::geom_point(size=4) +
    #ggtitle(gene_name)+
    ggrepel::geom_text_repel()
}


#Plot for labeling the truncated variants in the quantification plots:

plot_quant_trunc <- function(splice_df){
  qtp<- ggplot2::ggplot(splice_df, ggplot2::aes(x=variant, y=score, label=Transcript_id))+
  ggplot2::geom_point(size=4, shape=21,fill= ifelse(splice_df$truncation_origin=="none", "black",ifelse(splice_df$truncation_origin!="none", "white","black"))) +
  ggrepel::geom_text_repel(colour= "black")
  
  return(qtp)
}

#' Add more useful ids to unknown transcripts.
#' Append a number to any Transcript_ids labelled with the id "unknown".
#'
#' @param parsed_splices dataframe containing columns named variant, strand, score,
#' Transcript_id and splice_pattern. Intended to be the output of parse_nexons_gtf
#' function.
#' @param Transcript_id_col name of the Transcript id column
#' @param variant_col name of column containing the variant id
#'
#' @return tibble
#' @export
#'
#' @examples
#' NA
add_unknown_ids <- function(parsed_splices, Transcript_id_col = "Transcript_id", variant_col = "variant") {
  n_ids <- sum(parsed_splices[[Transcript_id_col]] != "unknown")
  parsed_splices %>%
    dplyr::mutate(unknown_index = .data[[variant_col]] - n_ids) %>%
    dplyr::mutate(Transcript_id = dplyr::if_else(.data[[Transcript_id_col]] == "unknown", paste("unknown", unknown_index), .data[[Transcript_id_col]]))
}


#' Draw a picture of splice patterns
#'
#' Creates plot showing splice variants
#'
#' @param splice_data dataframe containing columns named variant, strand, score,
#' Transcript_id and splice_pattern. Intended to be the output of parse_nexons_gtf
#' function.
#' @param order_splices (one of 'score', 'name', NULL) default NULL. How to order the splices on the y axis.
#' 'score' will sort data by score (highest at the top), 'name' will sort data
#' alphabetically (with unknowns at the bottom)
#' @param gene gene name/id to filter for
#' @param gene_id_col name of the column containing the gene names/ids
#' @param title_text optional title for plot (probably gene name)
#' @param quant TRUE or FALSE whether to show a quantitative plot of transcipt
#' counts alongside the splice variants
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' file <- system.file("extdata", "nexons_sirv5_f15.gtf", package = "nexonsAnalysis")
#' nexons_output <- readr::read_delim(file)
#' parsed_splices <- parse_nexons_gtf(nexons_output, min_count = 3)
#' draw_splice_picture(parsed_splices, quant = TRUE, order_splices = "score", gene="SIRV5")

#' @import ggplot2
draw_splice_picture <- function(splice_data, order_splices = NULL, gene = "", gene_id_col = "Gene_id", title_text = "", quant = FALSE) {

  if(title_text == "") title_text <- gene

  splice_data_filt <- if(gene %in% splice_data[[gene_id_col]]) {
    dplyr::filter(splice_data, .data[[gene_id_col]] == gene)
  } else if (nchar(gene) > 0){
    warning(message = paste("couldn't find gene name", gene, "in dataset, check it matches exactly. \n Continuing without filtering"))
    splice_data
  } else splice_data

  splice_data_filt <- add_unknown_ids(splice_data_filt)
  splice_plot_data <- add_exon_loci(splice_data_filt)

  # sort data by score (highest at the top)
  if(is.null(order_splices)) {
    print("No ordering of y axis specified")
  }
  else {
    if(order_splices == "score") {
      splice_plot_data <- splice_plot_data %>%
        dplyr::arrange(score, Transcript_id) %>%
        dplyr::mutate(Transcript_id = forcats::as_factor(Transcript_id)) %>%
        dplyr::mutate(variant = as.integer(Transcript_id))
    } else if (order_splices == "name") {
      # sort data alphabetically (with unknowns at the bottom)
      splice_plot_data <- splice_plot_data %>%
        dplyr::arrange(desc(startsWith(Transcript_id, "unknown")), desc(Transcript_id)) %>%
        dplyr::mutate(Transcript_id = forcats::as_factor(Transcript_id)) %>%
        dplyr::mutate(variant = as.integer(Transcript_id))
    }
  }
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
    quant_plot <- plot_quant(splice_data_filt)
    return(invisible(gridExtra::grid.arrange(p, quant_plot, nrow=1)))
  } else {
    return(p)
  }
}



#Draw variants and quantification with the truncated ones annotated with different shades:

draw_splice_picture_trunc <- function(splice_data, order_splices = NULL, gene = "", gene_id_col = "Gene_id", title_text = "", quant = FALSE) {
  
  if(title_text == "") title_text <- gene
  
  splice_data_filt <- if(gene %in% splice_data[[gene_id_col]]) {
    dplyr::filter(splice_data, .data[[gene_id_col]] == gene)
  } else if (nchar(gene) > 0){
    warning(message = paste("couldn't find gene name", gene, "in dataset, check it matches exactly. \n Continuing without filtering"))
    splice_data
  } else splice_data
  
  splice_data_filt <- add_unknown_ids(splice_data_filt)
  splice_plot_data <- add_exon_loci(splice_data_filt)
  
  # sort data by score (highest at the top)
  if(is.null(order_splices)) {
    print("No ordering of y axis specified")
  }
  else {
    if(order_splices == "score") {
      splice_plot_data <- splice_plot_data %>%
        dplyr::arrange(score, Transcript_id) %>%
        dplyr::mutate(Transcript_id = forcats::as_factor(Transcript_id)) %>%
        dplyr::mutate(variant = as.integer(Transcript_id))
    } else if (order_splices == "name") {
      # sort data alphabetically (with unknowns at the bottom)
      splice_plot_data <- splice_plot_data %>%
        dplyr::arrange(desc(startsWith(Transcript_id, "unknown")), desc(Transcript_id)) %>%
        dplyr::mutate(Transcript_id = forcats::as_factor(Transcript_id)) %>%
        dplyr::mutate(variant = as.integer(Transcript_id))
    }
  }
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
    geom_rect(alpha=ifelse(splice_plot_data$truncation_origin=="none", 1.0,ifelse(splice_plot_data$truncation_origin!="none", 0.4,1.0)))+#fill=strand,colour=shade) +
    scale_y_continuous(
      breaks=splice_plot_data$variant+0.5,
      labels = paste(splice_plot_data$Transcript_id)) +
    ylab("") +
    xlab("Genomic Position") +
    ggtitle(title_text) +
    scale_fill_manual(values = strand_colour) +
    scale_colour_manual(values = strand_colour)
  
  if(quant){
    quant_plot <- plot_quant_trunc(splice_data_filt)
    return(invisible(gridExtra::grid.arrange(p, quant_plot, nrow=1)))
  } else {
    return(p)
  }
}
