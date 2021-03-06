
### volcanoPlot2()
### ===============================================================================================

#' @title volcanoPlot2()
#' @author amitjavilaventura & dfernandezperez
#'
#' @usage colcanoPlot2(df, xlim = c(-10,10), ylim = c(0,30), pval = 0.05, log2FC = 1.5, main = NULL, mainSize = 9, sub = NULL, subSize = 8, labelSize = 7, labelColor = c("darkgreen", "red"), labelPos = 0, xlab = "log2(FC)", ylab = "-log10(pval)", axisLabelSize = 7, axisTextSize = 7, pointColor = c("darkgreen", "gray", "red"), legendTitle = FALSE, legendPos = "bottom", degsLabel=F, degsLabelNum=5, ggrastr = F)
#'
#' Function that draws a volcano plot with DEGs.
#' As input, it takes the output of DESeq2 and adding a DEG information column (Downregulated, Upregulated, NotDE). Columns: Geneid, baseMean, log2FoldChange, lfcSE, pvalue, padj, DEG.
#' It also counts the number of up and downregulated genes and writes them in the plot.
#' The points that are outside of the limits of the graph are drawn as triangles at the borders.
#'
#' @param df Dataframe. Output from DESeq2 with an extra DEG column. Columns: Geneid, baseMean, log2FoldChange, lfcSE, pvalue, padj, DEG.
#' @param xlim Numerical of length 2. The limits of the x axis where the log2FC is plotted. Default: c(-10,10).
#' @param ylim Numerical of length 2. The limits of the y axis where the -log10(adj_pval) is plotted.  Default: c(0,30).
#' @param pval Numerical. p-value threshold used to call de DEGs. It is used to draw a dashed line in this value. Default: 0.05.
#' @param log2FC Numerical. log2FC threshold (in absolute value) used to call the DEGs. It is used to draw two lines in this value (positive and negative). Default: 1.5.
#' @param main Character. The title of the volcano plot. Default: NULL.
#' @param mainSize Numerical. The size of the title. Default: 9.
#' @param sub Character. The subtitle of the volcano plot. Default: NULL
#' @param subSize. Numerical. The size of the subtitle. Default: 8.
#' @param xlab Character. label of the X axis. Default: "log2(FC)".
#' @param ylab Character. Label of the Y axis. Default: "-log10(pval)".
#' @param labelSize Numerical. The size of the numbers of up and downregulated genes. Default: 7.
#' @param labelColor Character of length 2. Colors of the numbers of downregulated and upregulated genes. Default: c("darkgreen", "red").
#' @param labelPos Numerical. Position of the numbers of up and downregulated genes along the y axis. Default: 0.
#' @param axisLabelSize Numerical. Size of the axis labels. Default: 7.
#' @param axisText_Size Numerical. Size of the text in the axis. Default: 7.
#' @param pointColor Character of length 3. Colors of the downregulated, notDE and upregulated points. Default: c("darkgreen", "gray", "red")
#' @param legendTitle Logical. If FALSE, the title of the legend is not plotted. Default: FALSE.
#' @param legendPos Character. Position of the legend. One of: "bottom", "top", "right", "left", "none". Default: "bottom".
#' @param degsLabel Logical. If TRUE, the genes with highest log2FC in absolute are labelled. Default: FALSE.
#' @param degsLabelNum Numerical. Number of most expressed genes to label in the plot. The double of this number will be labelled (once for DEGs with lowest p-value and once for the DEGs with highest log2fFC in absolute value). Default: 5.
#' @param degsLabelSize Numerical. Size of the labels of the DEGs Default: 3.
#' @param ggrastr Logical. If TRUE, the points of the volcano plot are drawn with geom_points_rast instead of geom_point. Default: FALSE.
#' @param plotly Logical. If TRUE, the ggplot object will be passed to 'ggplotlify()' to transform it to an interactive plotly graph.

volcanoPlot2 <- function(df, xlim = c(-10,10), ylim = c(0,30),
                         pval = 0.05, log2FC = 1.5,
                         main = NULL, mainSize = 9, sub = NULL, subSize = 8,
                         labelSize = 7, labelColor = c("darkgreen", "red"), labelPos = 0,
                         xlab = bquote(~Log[2]~ "FC"), ylab = (bquote(~-Log[10]~italic(P))) , axisLabelSize = 7, axisTextSize = 7,
                         pointColor = c("darkgreen", "gray", "red"), legendTitle = FALSE, legendPos = "bottom",
                         degsLabel = F , degsLabelNum=5, degsLabelSize = 3,
                         ggrastr = F, plotly = F) {

  #load packages
  require(ggplot2)
  require(dplyr)

  # Mutate the dataframe to include the shape of the points.
  df <- df %>%

    # Shape of al the points
    dplyr::mutate(shape = "circle") %>%

    # Change the shape of the point to triangle if the -log10(padj) is greater than the ylim variable
    dplyr::mutate(shape = ifelse(test = -log10(padj) > ylim[2], yes = "triangle", no = shape)) %>%

    # Change the shape of the points to triangle if the log2FoldChange is greater or lower than the xlim
    dplyr::mutate(shape = ifelse(test = log2FoldChange > xlim[2], yes = "triangle", no = shape)) %>%
    dplyr::mutate(shape = ifelse(test = log2FoldChange < xlim[1], yes = "triangle", no = shape)) %>%

    # Change the padj to the value that plots the points at the top of the graph.
    dplyr::mutate(padj  = ifelse(test = -log10(padj) > ylim[2], yes = 10^-ylim[2], no = padj)) %>%

    # Change the log2FoldChange of the points with log2FoldChange than the limits, so they will be plotted in the limits of the graph
    dplyr::mutate(log2FoldChange = ifelse(test = log2FoldChange > xlim[2], yes = xlim[2], no = log2FoldChange)) %>%
    dplyr::mutate(log2FoldChange = ifelse(test = log2FoldChange < xlim[1], yes = xlim[1], no = log2FoldChange))


  # Start plot with ggrastr points. This is useful when using workflowr.
  if(ggrastr){

    # Load ggrastr
    require(ggrastr)

    p <- ggplot(data = na.omit(df), aes(x=log2FoldChange, y=-log10(padj), colour=DEG, shape=shape)) +

      # Draw points with geom_point_rast from ggrastr, because it will draw points as images, making the graph lighter.
      geom_point_rast(alpha=0.7, size=1.7, raster.height = 5.15, raster.width = 6, raster.dpi = 400)
  }

  # Start plot normally
  else{
    p <- ggplot(data = na.omit(df), aes(x=log2FoldChange, y=-log10(padj), colour=DEG, shape=shape)) +

      geom_point(alpha=0.7, size=1.7)
  }


  # Annotate the number of up and downregulated DEGs
  p <- p +
    annotate("text", label = sum(df$DEG == "Upregulated"), color = labelColor[2], y = labelPos, x = xlim[2],
             vjust=0.5,hjust="inward", size = labelSize) +
    annotate("text", label = sum(df$DEG == "Downregulated"), color = labelColor[1], y = labelPos, x = xlim[1],
             vjust=0.5,hjust="inward", size = labelSize)

  # Basic formatting
  p <- p +

    # Stablish a predefined theme
    theme_classic() +

    # Write and format the graph title, can be nothing.
    ggtitle(label = main, subtitle = sub) +
    theme(plot.title = element_text(face="bold", hjust = .5, size = mainSize),
          plot.subtitle = element_text(hjust = .5, size = subSize)) +

    # Stablish the x and y axes ranges.
    coord_cartesian(xlim = xlim, ylim = ylim) +

    # Put an horizontal line in the -log10(pval) value and two vertival lines in the -logFC and logFC values.
    geom_hline(yintercept = -log10(pval), linetype = 2) +
    geom_vline(xintercept = c(-log2FC, log2FC), linetype = 2) +

    # Format the axis names and sizes
    xlab(xlab) + ylab(ylab) + theme(axis.title = element_text(size = axisLabelSize)) +

    # Format the color of the points
    scale_colour_manual(values=c("Downregulated" = pointColor[1], "NS" = pointColor[2], "Upregulated" = pointColor[3]),
                        labels = c("Downregulated" = "Downregulated", "NS" = "NS", "Upregulated" = "Upregulated"),
                        drop = FALSE) +

    # Remove the legend for shape
    guides(shape=FALSE) +

    # Format the axis values
    theme(axis.text = element_text(size = axisTextSize)) +

    # Decide the position of the legend (default: "bottom")
    theme(legend.position = legendPos)

  # Decide if legend title is writen or not. Default: not writen.
  if(!legendTitle){
    p <- p + theme(legend.title = element_blank())
  }

  # Write names of the most DE genes in terms of lowest adjusted p-value
  if(degsLabel) {

    # Load ggrepel
    require(ggrepel)

    # Organaize and retrieve lowest p-value genes
    degs <- df %>%

      # Filter non significant genes
      dplyr::filter(DEG!="NS") %>%

      # Arrange by ascendent order of padjusted
      dplyr::arrange(padj)

    # Create a dataframe with the labels of the DEGs with highest abs(log2FC).
    degs <- head(na.omit(degs),degsLabelNum) %>% as.data.frame()

    # Put labels in the plot
    p <- p + geom_text_repel(data = degs, mapping = aes(x = log2FoldChange, y = -log10(padj), label = Geneid), size = degsLabelSize, color = "Black")
  }

  if(plotly == T){
    require(plotly)
    p <- ggplotly(p)
  }

  # Draw the graph.
  return(p)

}
