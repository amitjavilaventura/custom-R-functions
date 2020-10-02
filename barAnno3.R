##################
### BAR ANNO 3 ###
##################

#' @title barAnno3
#' @author amitjavilaventura
#'
#' @usage barAnno3(anno_list, names, main = NULL, subtitle = NULL, ylab = "Proportion", xlab = NULL, palette = "Set1", legend_position = "bottom")
#' 
#' Function for ChIP-seq and ATAC-seq.
#' It must be used after the function annotatePeak() from the R package ChIPseeker. @seealso \code{\link{annotatePeak}}
#'
#' It takes a list of annotation objects that come as output of annotatePeak() and changes the features to "Promoter", "Distal" and "Gene body". Finally it plots a bargraph with the distribution of all the proportions
#' As a ggplot2-based function, it allows to add more layers to format the plot.
#' 
#' @param anno_list List of annotation objects that come from annotatePeak().
#' @param names Charachter vector of the same length as 'anno_list'. Names that will be given to each of the objects in anno_list. Not that will be the names plotted in the bargraph
#' @param main Character of lenght 1. Title of the plot. Default: NULL.
#' @param subtitle Character of lenght 1. Subtitle of the plot. Default: NULL.
#' @param ylab Character of lenght 1. Title of the Y axis. Default: "Proportion".
#' @param xlab Character of lenght 1. Title of the X axis. Default: NULL
#' @param palette Character of lenght 1. Color palette used to color the bars through the function `scale_fill_brewer()`. Default: "Set1".
#' @param legend_position Character of lenght 1. Position of the legend. One of c("none", "bottom", "right", "left," "top"). Default: "right"
#' @param is_anno Logical. If TRUE, takes the 'anno_list' as a list of annotation objects from annotatePeak. If FALSE, takes 'anno_list' as a list of df that have already been formatted (i.e. read from file)
#' 

barAnno3 <- function(anno_list, names, 
                     main = NULL, 
                     subtitle = NULL,
                     ylab = "Proportion", xlab = NULL,
                     palette = "Set1",
                     legend_position = "right", is_anno = T){
  
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
    purrr::map(~mutate(.data = .x, annotation = as.character(annotation))) %>%
  
    # Bind dataframes by rows
    bind_rows() %>% 

    # Rewrite annotation as distal or gene body depending on the feature 
    dplyr::mutate(annotation = dplyr::recode(annotation,
                                              "Distal Intergenic" = "Distal",
                                              "Downstream" = "Distal",
                                              "Intron" = "Gene body",
                                              "Exon" = "Gene body",
                                              "5' UTR" = "Gene body",
                                              "3' UTR" = "Gene body"))
  
  # Create gglpot2-based barplot
  g <- ggplot(anno_df, aes(sample, fill = factor(annotation))) + geom_bar(position = "fill") +
    
    # Write plot title, subtitle and axis labels
    ggtitle(main, subtitle = subtitle) + ylab(ylab) + xlab(xlab) +
    
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

