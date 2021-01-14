### ============= ###
### AnnoFunctions ###
### ============= ###

### filterAnno() -----
### ============================================================================================== #

#' @title filterAnno
#' @author amitjavilaventura
#'
#' @usage filterAnno(df)
#'
#' Function for ChIP-seq and ATAC-seq.
#' It must be used after the function annotatePeak() from the R package ChIPseeker. @seealso \code{\link{annotatePeak}}
#'
#' It takes the annotation of the output of annotatePeak() and changes the features to "Promoter" (in the case of the different promoter regions),
#' "Distal" (all the features that Distal intergenic or Downstream regions) and "Gene body" regions (the other regions).
#'
#' @param df Data frame obtained as output of annotatePeak().
#' @param anno_num Integer number. Number of annotations to plot, either 2 (Promoter/Distal) or 3 (Promoter/Gene body/Distal). Default: 3.

filterAnno <- function(anno, anno_num = 2){

  # load packages
  require(dplyr)

  # remove all the different promoter categories (<=1kb, 1-2kb, 2-3kb) to have only one promoter category called "Promoter"
  if(class(anno) == "csAnno"){
    anno@anno$annotation <- sub(" \\(.*\\)", "", anno@anno$annotation)
    anno <- anno@anno
  }
  else if(class(anno) == "data.frame"){ anno$annotation <- sub(" \\(.*\\)", "", anno$annotation) }
  else{ stop("The 'anno' object must be an object of class 'csAnno' or 'data.frame'")}

  # change feature names
  if(anno_num == 3){
    anno <- anno %>% dplyr::mutate(annotation = recode(annotation,
                                                "Distal Intergenic" = "Distal", "Downstream" = "Distal",
                                                "Intron" = "Gene body", "Exon" = "Gene body",
                                                "5' UTR" = "Gene body", "3' UTR" = "Gene body"))
  }
  else if(anno_num == 2){
    anno <- anno %>% dplyr::mutate(annotation = recode(annotation,
                                                "Distal Intergenic" = "Distal", "Downstream" = "Distal",
                                                "Intron" = "Distal", "Exon" = "Distal",
                                                "5' UTR" = "Distal", "3' UTR" = "Distal"))
  }
  else{ stop("'anno_num' must be 2 or 3.") }

  return(anno)
}



### pieAnno() -----
### ============================================================================================== #

#' @title pieAnno
#' @author amitjavilaventura
#'
#' @usage pieAnno <- function(df, main = NULL, size.main = 13, sub = NULL, size.sub = 11, legend = TRUE, size.legend = 9 , pos.legend = "bottom", fill.colors = c("Steelblue", "Lightsalmon3"), labels = TRUE, decimal = 3, size.label = 9)
#'
#' Function for ChIP-seq and ATAC-seq.
#' It must be used after the function annotatePeak() from the R package ChIPseeker. @seealso \code{\link{annotatePeak}}
#'
#' It takes the annotation of the output of annotatePeak() and changes the features to "Promoter" (in the case of the different promoter regions) and "Distal" (all the features that are not promoter regions).
#' Then, it creates a dataframe with the percentage of each, distal and promoter features.
#' After that, it plots a ggplot2-based pie chart with these percentages
#' Finally, it returns the pie chart returns a pie chart.
#'
#' @param df Data frame obtained as output of annotatePeak().
#' @param main Character. Title of the pie chart. Default: NULL.
#' @param size.main Numerical. Font size of the title. It works only if main is not NULL. Default: 13.
#' @param sub Character.Subtitle of the piechart. It works only if main is not NULL. Default: NULL.
#' @param size.sub Numerical. Font size of the subtitle. It works only if main is not NULL. Default: 11.
#' @param legend Logical. If TRUE, legend is plotted. Default: TRUE.
#' @param size.legend Numerical. Font size of the legend text. It works only if legend is set to TRUE. Default: 9.
#' @param pos.legend Character. Position of the legend in the graph. One of "none", "left", "right", "bottom", "top". It works only if legend is set to TRUE. Default: "bottom".
#' @param labels Logical. If TRUE, labels of the percentages are ploted in the chart. Default: TRUE.
#' @param decimal Integer. Number of decimals to round the percentage labels inside the pie chart. It only works if labels is set to TRUE. Default: 3.
#' @param size.label Numerical. Font size of the percentage labels inside the pie chart.  It only works if labels is set to TRUE. Default: 9.
#' @param fill.palette Character. Palette passet to 'fill_palette()' to fill the areas. Colors to use on each of the regions. Default: "Set2"
#' @param anno_num Integer number. Number of annotations to plot, either 2 (Promoter/Distal) or 3 (Promoter/Gene body/Distal). Default: 3.

pieAnno <- function(anno, main = NULL, size.main = 12, sub = NULL, size.sub = 10, legend = TRUE, size.legend = 8 ,
                    pos.legend = "bottom", fill.palette = "Set2", labels = TRUE, decimal = 2, size.label = 7, anno_num = 2){

  # ----- LOAD PACKAGES ----- #
  require(ggplot2)
  require(ggpubr)
  require(dplyr)

  # ----- FORMATTING THE DATA ----- #

  # call filterAnno
  anno <- filterAnno(anno, anno_num = anno_num)

  # calculate the % of promoter or distal fields.
  anno_100 <- anno %>%
    count(annotation) %>%
    group_by(annotation) %>%
    summarise(n = n, percent = n/nrow(anno)*100) %>%
    arrange(percent) %>%
    mutate(ypos = cumsum(percent)- 0.5*percent )

  # ----- DRAW THE PIE CHART ----- #

  # generate the ggplot graph
  p <- ggplot(data = anno_100, aes(x = "", y = percent, fill = annotation)) +

    # do a stacked bargraph
    geom_bar(width = 1, stat = "identity", color = "black", show.legend = legend) +

    # turn the bargraph to a pie chart
    coord_polar("y", start = 0) +

    # color the areas
    scale_fill_brewer(palette = fill.palette) +

    #stablish the format for the graph, axis, legend, panel.
    theme_pubr(legend = pos.legend) +
    theme(axis.title = element_blank(), axis.text = element_blank(),
          panel.background = element_blank(), panel.grid=element_blank(), panel.border = element_blank(),
          plot.tag = element_blank(), legend.title = element_blank())

  # write the percentages of each piece of the pie
  if(labels==TRUE){ p <- p + geom_text(aes(y = ypos, label = c(round(percent, digits = decimal))), size = size.label) }

  # stablish the format for title and subtitle
  if(!is.null(main)){
    p <- p + ggtitle(label = main, subtitle = sub) +
      theme(plot.title = element_text(face = "bold", size = size.main, hjust = 0.5),
            plot.subtitle = element_text(face = "bold.italic", size = size.sub, hjust = 0.5))
  }

  #stablish the colors of the piechart
  if(!is.null(fill.colors)){
    p <- p + scale_fill_manual(values = fill.colors)
  }

  # RETURN THE PIE CHART
  return(p)
}

