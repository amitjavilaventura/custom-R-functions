##############################
###     VOLCANO PLOT 2     ###
##############################

# ggplot2-based fucntion made by Dani. I modified it so if I put all together in the same panel, it is easier to modify.
# Added:
#   -Change colors of points
#   -Change the axis label names and sizes
#   -Change the title size
#   -Change the legend position
#   -Changed color and position of the number of up/down DEGs

VolcanoPlot2 <- function(df, xlim = c(-10,10), ylim = c(0,30), main = NULL, title_size = 9, 
                         labelSize = 7, labelColor = c("green", "darkred"), labelPos = 0,
                         pval = 0.05, log2FC = 1.5, xlab = "log2(FC)", ylab = "-log10(pval)", axis_label_size = 7, axis_text_size = 7, 
                         point.color = c("darkgreen", "gray", "red"), legend_title = FALSE, legend_pos = "bottom", 
                         degsLabelNum = 5, degsLabel = F) {
  
  #load package ggplot2
  require(ggplot2)
  require(dplyr)
  
  #generate the graph with ggplot(), stablish the foldchange and pvalue to colour the DEGs and stablish the point format.
  p <-  ggplot(data = df, aes(x=log2FoldChange, y=-log10(padj), colour=DEG) ) +
    geom_point(size=1.5) +
    
    #annotate the number of up and downregulated DEGs
    annotate("text", label = sum(df$DEG == "Upregulated"), color = labelColor[2], y = labelPos, x = xlim[2], 
             vjust=0.5,hjust="inward", size = labelSize) +
    annotate("text", label = sum(df$DEG == "Downregulated"), color = labelColor[1], y = labelPos, x = xlim[1],
             vjust=0.5,hjust="inward", size = labelSize) + 
    
    #stablish a predefined theme
    theme_classic() +
    
    #write and format the graph title, can be nothing. 
    ggtitle(main) +
    theme(plot.title = element_text(face="bold", hjust = .5, size = title_size)) +
    
    #stablish the x and y axes ranges. 
    xlim(xlim) + ylim(ylim) + 
    
    #put an horizontal line in the -log10(pval) value and two vertival lines in the -logFC and logFC values. 
    geom_hline(yintercept = -log10(pval), linetype = 2) +
    geom_vline(xintercept = c(-log2FC, log2FC), linetype = 2) +
    
    #format the axis names and sizes
    xlab(xlab) + ylab(ylab) + theme(axis.title = element_text(size = axis_label_size)) +
    
    #format the color of the points
    scale_colour_manual(values=point.color) +
    
    #format the axis values
    theme(axis.text = element_text(size = axis_text_size)) +
    
    #decide the position of the legend (default: "bottom")
    theme(legend.position = legend_pos) 

  #Decide if legend title is writen or not. Default: not writen.  
  if(legend_title==FALSE){ p <- p + theme(legend.title = element_blank()) }
  
  #Put labels of the genes
  if(degsLabel==TRUE) {
    #load ggrepel
    require(ggrepel)
    
    # organaize and retrieve most upregulated and most downregulated genes
    df$abs_log2fc <- sqrt(df$log2FoldChange**2)
    df  <- df %>% dplyr::arrange(abs_log2fc)
    df2 <- rbind(tail(na.omit(df),degsLabelNum)) %>% as.data.frame()
    #put labels in the plot
    p <- p + geom_text_repel(data = df2, mapping = aes(x = log2FoldChange, y = -log10(padj), label = Geneid), size = 5, color = "Black")
    
  }
  
  #draw the graph.
  p
  
}
