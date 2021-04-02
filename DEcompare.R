### DEcompare()
### ===============================================================================================

DEcompare <- function(df1, df2, measure = "Log2FC", threshold = 1.5,
                      xlim = c(-10, 10), ylim = c(-10, 10),
                      xlab = "Contrast1", ylab = "Contrast2",
                      title = paste("Comparison of", measure),
                      subtitle = paste(xlab, "vs", ylab),
                      color_corners = c("pink", "lightgreen", "cornflowerblue", "yellow"),
                      alpha_corners = c(.7)){

  data <- inner_join(x = df1 %>% select(Geneid, log2FoldChange, padj, DEG),
                     y = df2 %>% select(Geneid, log2FoldChange, padj, DEG),
                     by="Geneid")

  if(measure == "Log2FC"){

    data <- data %>%
      mutate(color = if_else(log2FoldChange.x >= threshold & log2FoldChange.y >= threshold, "positivecorr",
                             if_else(log2FoldChange.x >= threshold & log2FoldChange.y <= -threshold, "negativecorr",
                                     if_else(log2FoldChange.x <= -threshold & log2FoldChange.y <= -threshold, "positivecorr",
                                             if_else(log2FoldChange.x <= -threshold & log2FoldChange.y >= threshold, "negativecorr", "ns"))))) %>%
      tibble::column_to_rownames("Geneid")

    g <- ggplot(data, aes(log2FoldChange.x, log2FoldChange.y, color = color))

    if(!is.null(color_corners)){
      g <- g +
        annotate(geom = "rect",
                 xmin = c(xlim[1], xlim[1], 0, 0),
                 xmax = c(0, 0, xlim[2], xlim[2]),
                 ymin = c(0,  ylim[1],  0,  ylim[1]),
                 ymax = c(ylim[2], 0, ylim[2], 0),
                 fill = color_corners, alpha = alpha_corners)
    }

    g <- g + geom_point() + coord_cartesian(xlim = xlim, ylim = ylim, expand = F) +

      geom_hline(yintercept = threshold, linetype = 2, color = "black") +
      geom_hline(yintercept = -threshold, linetype = 2, color = "black") +
      geom_vline(xintercept = threshold, linetype = 2, color = "black") +
      geom_vline(xintercept = -threshold, linetype = 2, color = "black") +

      scale_colour_manual(values=c("positivecorr" = "darkred", "negativecorr" = "darkgreen", "ns" = "gray50"),
                          drop = T) +



      ggtitle(title, subtitle) +
      xlab(xlab) + ylab(ylab) +

      theme_pubr(border = T, legend = "none")

  } else if(measure == "pvalue"){
    # to do
  }


  # return
  return(g)
}
