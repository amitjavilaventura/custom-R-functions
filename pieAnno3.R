#####################
### FILTER ANNO 3 ###
#####################

#' @title filterAnno3
#' @author amitjavilaventura
#' 
#' @usage filterAnno3(df)
#' 
#' Function for ChIP-seq and ATAC-seq.
#' It must be used after the function annotatePeak() from the R package ChIPseeker. @seealso \code{\link{annotatePeak}}
#' 
#' It takes the annotation of the output of annotatePeak() and changes the features to "Promoter" (in the case of the different promoter regions), 
#' "Distal" (all the features that Distal intergenic or Downstream regions) and "Gene body" regions (the other regions). 
#' 
#' @param df Data frame obtained as output of annotatePeak().

filterAnno3 <- function(df){
  
  # ----- FORMATTING THE DATA ----- #
  
  # remove all the different promoter categories (<=1kb, 1-2kb, 2-3kb) to have only one promoter category called "Promoter"
  anno <- sub(" \\(.*\\)", "", df@anno$annotation)
  
  # change feature names
  for(i in 1:length(anno)){
    
    # "Distal Intergenic" and "Downstream" to "Distal"
    if(anno[i] == "Distal Intergenic" | anno[i] == "Downstream"){
      anno[i] <- "Distal"
    }
    
    # If the feature is not "Distal" or "Promoter", change it to GeneBody
    else if(anno[i] != "Distal" & anno[i] != "Promoter"){
      anno[i] <- "Genebody"
    }
    
  }
  
  return(anno)
}
  
  

#####################
### PLOT PIE ANNO ###
#####################

#' @title pieAnno2
#' @author amitjavilaventura
#' 
#' @usage pieAnno2 <- function(df, main = NULL, size.main = 13, sub = NULL, size.sub = 11, legend = TRUE, size.legend = 9 , pos.legend = "bottom", fill.colors = c("Steelblue", "Lightsalmon3"), labels = TRUE, decimal = 3, size.label = 9)
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
#' @param fill.colors. Character vector of length 3. Colors to use on each of the regions. Default: c("Steelblue", "Lightsalmon3", "Lightgreen")

pieAnno3 <- function(df, main = NULL, size.main = 12, sub = NULL, size.sub = 10, legend = TRUE, size.legend = 8 , 
                     pos.legend = "bottom", fill.colors = c("Steelblue", "Lightsalmon3", "Lightgreen"), 
                     labels = TRUE, decimal = 2, size.label = 7){
  
  # ----- LOAD PACKAGES ----- #
  require(ggplot2)
  
  # ----- FORMATTING THE DATA ----- #
  
  # call filterAnno3"
  anno <- filterAnno3(df)
  
  # subset the promoter and distal peaks
  prom <- subset(anno, anno=="Promoter")
  gene <- subset(anno, anno=="Genebody")
  dist <- subset(anno, anno=="Distal")
  
  #put the number of promoters and distal peaks to a data frame
  peak.num <- as.data.frame(c(length(prom), length(dist), length(gene)))
  colnames(peak.num) <- "Number of peaks"
  rownames(peak.num) <- c("Promoter", "Distal", "Genebody")
  
  # calculate the % of promoter or distal fields. 
  prom.100 <- length(prom)/length(anno)*100
  dist.100 <- length(dist)/length(anno)*100
  gene.100 <- length(gene)/length(anno)*100
  
  # put the promoter and distal percentages in a data.frame
  df.100 <- as.data.frame(matrix(c(prom.100, gene.100, dist.100), ncol = 1))
  df.100$Category <- c("Promoter", "Genebody", "Distal")
  colnames(df.100) <- c("Per100","Category")
  
  
  # ----- DRAW THE PIE CHART ----- #
  
  # generate the ggplot graph
  p <- ggplot(data = df.100, aes(x = "", y = Per100, fill = Category)) +
    
    # do a stacked bargraph
    geom_bar(width = 1, stat="identity", color = "gray15", show.legend = TRUE) +
    
    # turn the bargraph to a pie chart
    coord_polar("y", start = 0) +
    
    #stablish the format for the graph, axis, legend, panel.
    theme_minimal() + 
    theme(axis.title = element_blank(), axis.text = element_blank(),
          panel.background = element_blank(), panel.grid=element_blank(), panel.border = element_blank(),
          plot.tag = element_blank())
  
  # write the percentages of each piece of the pie
  if(labels==TRUE){
    p <- p +  geom_text(aes(y = df.100$Per100/2 + c(0, cumsum(df.100$Per100)[1], cumsum(df.100$Per100)[2]), 
                            label = c(round(df.100$Per100[1], digits = decimal), 
                                      round(df.100$Per100[2], digits = decimal), 
                                      round(df.100$Per100[3], digits = decimal))), 
                        size = size.label)
  }
 
  # stablish the format for title and subtitle
  if(!is.null(main)){
    p <- p + ggtitle(label = main, subtitle = sub) + 
      theme(plot.title = element_text(face = "bold", size = size.main, hjust = 0.5), 
            plot.subtitle = element_text(face = "bold.italic", size = size.sub, hjust = 0.5))
  }
  
  #stablish a format for the legend
  if(legend){
    p <- p + theme(legend.title = element_blank(), legend.text = element_text(size = size.legend), legend.key.size = unit(5, "mm"), legend.position = pos.legend)
  }
  
  #stablish the colors of the piechart
  if(!is.null(fill.colors)){
    p <- p + scale_fill_manual(values = fill.colors)
  }
  
  # RETURN THE PIE CHART
  return(p)
}

