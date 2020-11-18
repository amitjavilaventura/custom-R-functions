##################
### BAR ANNO 3 ###
##################

#' @title barAnno3
#' @author amitjavilaventura
#'
#' @usage barAnno(anno_list, names, names_order, protein = NULL, main = NULL, subtitle = NULL, ylab = "Proportion", xlab = NULL,palette = "Set1", legend_position = "right", is_anno = F, is_chip = F, anno_num = 3){
#'
#' Function for ChIP-seq and ATAC-seq.
#' It must be used after the function annotatePeak() from the R package ChIPseeker. @seealso \code{\link{annotatePeak}}
#'
#' It takes a list of annotation objects that come as output of annotatePeak() and changes the features to "Promoter", "Distal" and "Gene body" (or to "Promoter" and "Distal"). Finally it plots a bargraph with the distribution of all the proportions
#' As a ggplot2-based function, it allows to add more layers to format the plot.
#'
#' @param anno_list List of annotation objects that come from annotatePeak().
#' @param names Charachter vector of the same length as 'anno_list'. Names that will be given to each of the objects in anno_list. Not that will be the names plotted in the bargraph
#' @param names_order Character vector with the same entries as 'names' with the order wanted to plot the data.
#' @param protein Character vector of the same length as 'anno_list' with the protein chipped in each dataframe of anno_list. This names will be passed through `facet_wrap()`. Only works if is_chip is set to TRUE. Default: NULL.
#' @param main Character of lenght 1. Title of the plot. Default: NULL.
#' @param subtitle Character of lenght 1. Subtitle of the plot. Default: NULL.
#' @param ylab Character of lenght 1. Title of the Y axis. Default: NULL.
#' @param xlab Character of lenght 1. Title of the X axis. Default: NULL
#' @param palette Character of lenght 1. Color palette used to color the bars through the function `scale_fill_brewer()`. Default: "Set1".
#' @param legend_position Character of lenght 1. Position of the legend. One of c("none", "bottom", "right", "left," "top"). Default: "right"
#' @param is_anno Logical. If TRUE, takes the 'anno_list' as a list of annotation objects from annotatePeak. If FALSE, takes 'anno_list' as a list of `annotatePeak()` output that have already been formatted to data.frame (i.e. read from file)... this is done because of the behaviour of our snakemake pipeline.
#' @param is_chip Logical. If TRUE, it means that the input data comes from a ChiP experiment, which will separate the plot by proteins. It can be set to FALSE even if data comes from ChIPseq, but remember that if you want to separate proteins, the names of each protein group must be different. It can be also set to TRUE if data comes from ATACseq, then in 'protein' you can write the desired grouping variable. Default: FALSE.
#' @param anno_num Integer number. Number of annotations to plot, either 2 (Promoter/Distal) or 3 (Promoter/Gene body/Distal). Default: 3.
#' @param position_fill Logical. If TRUE (default), it the plotted bars will represent proportion of peaks in each feature. If FALSE, the bars will have the height of the total number of peaks with the correspondent feature color.

barAnno  <- function(anno_list, names, names_order,
                     protein = NULL,
                     main = NULL,
                     subtitle = NULL,
                     ylab = NULL, xlab = NULL,
                     palette = "Set1",
                     legend_position = "right",
                     is_anno = F, is_chip = F, anno_num = 3,
                     position_fill = T){

  # Load packages
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(ggpubr)
  library(magrittr)

  # Create list for annotations from the list of annotatePeak object provided as input
  anno <- list()

  # Format promoter annnotation names to take parentheses out for each annotatePeak object in the list
  if(is_anno){
    for(i in 1:length(anno_list)){
      anno[[i]] <- sub(" \\(.*\\)", "", anno_list[[i]]@anno$annotation)
    }
  }
  else{
    for(i in 1:length(anno_list)){
      anno[[i]] <- sub(" \\(.*\\)", "", anno_list[[i]]$annotation)
    }
  }

  # Set the names of each annotatePeak object in the list with the vector names
  anno_df <- set_names(x = anno, nm = names) %>%

    # Convert the annotatePeak objects to dataframe
    purrr::map(~as.data.frame(x = .x)) %>%

    # Write an extra column to each dataframe with the name of the dataframe (provided in names)
    purrr::imap(~mutate(.data = .x, sample = as.character(.y))) %>%

    # Set column names
    purrr::map(~set_colnames(x = .x, c(c("annotation", "sample")))) %>%
    purrr::map(~dplyr::mutate(.data = .x, annotation = as.character(annotation)))

  if(is_chip){
    anno_df <- set_names(x = anno_df, nm = protein) %>%

      # Write an extra column to each dataframe with the name of the dataframe (provided in names)
      purrr::imap(~mutate(.data = .x, protein = as.character(.y))) %>%

      # Set column names
      purrr::map(~set_colnames(x = .x, c(c("annotation", "sample", "protein"))))
  }

  # Bind dataframes by rows
  anno_df <- anno_df %>% bind_rows() %>%

    # Mutate the names in order to
    dplyr::mutate(sample = factor(sample, levels = names_order))

  # Rewrite annotation as distal or gene body depending on the feature
  if(anno_num == 2){
    anno_df <- anno_df %>%
      dplyr::mutate(annotation = dplyr::recode(annotation,
                                               "Distal Intergenic" = "Distal",
                                               "Downstream" = "Distal",
                                               "Intron" = "Distal",
                                               "Exon" = "Distal",
                                               "5' UTR" = "Distal",
                                               "3' UTR" = "Distal"))
  }
  else if(anno_num == 3){
    anno_df <- anno_df %>%
      dplyr::mutate(annotation = dplyr::recode(annotation,
                                               "Distal Intergenic" = "Distal",
                                               "Downstream" = "Distal",
                                               "Intron" = "Gene body",
                                               "Exon" = "Gene body",
                                               "5' UTR" = "Gene body",
                                               "3' UTR" = "Gene body"))
  }

  # Create gglpot2-based barplot
  if(position_fill){
    if(is_chip){
      g <- ggplot(anno_df, aes(sample, fill = annotation)) + geom_bar(position = "fill") + facet_wrap(~protein)
    }
    else{
      g <- ggplot(anno_df, aes(sample, fill = annotation)) + geom_bar(position = "fill")
    }
  }

  else{
    if(is_chip){
      g <- ggplot(anno_df, aes(sample, fill = annotation)) + geom_bar() + facet_wrap(~protein)
    }
    else{
      g <- ggplot(anno_df, aes(sample, fill = annotation)) + geom_bar()
    }
  }
    # Write plot title, subtitle and axis labels
  g <- g + ggtitle(main, subtitle = subtitle) + ylab(ylab) + xlab(xlab) +

    # Format colors with ggplot2 function scale_fill_brewer
    scale_fill_brewer(palette = palette) +

    # Basic formatting
    theme_pubr() +

    theme(axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1),
          legend.title = element_blank(),
          legend.position = legend_position)

  # Draw plot
  return(g)
}


