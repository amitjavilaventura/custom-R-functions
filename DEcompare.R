### DEcompare()
### ===============================================================================================

DEcompare <- function(deg_list, threshold = 1.5,
                      genes = c("Apc"),
                      xlim = c(-10, 10), ylim = c(-10, 10),
                      xlab = "Contrast1", ylab = "Contrast2",
                      main = paste("Comparison of Log2FC"),
                      subtitle = paste(xlab, "vs", ylab),
                      color_corners = c("pink", "lightgreen", "cornflowerblue", "yellow"),
                      alpha_corners = c(.7)){
  # Load packages
  require(dplyr)
  require(tibble)
  require(ggplot2)
  require(ggrepel)
  require(ggpubr)

  # Check that inputs are OK
  if(!is.list(deg_list) | length(deg_list) != 2){ stop("'deg_list' must be a list of 2 data frames with the columns 'Geneid', 'log2FoldChange', 'padj' and 'DEG'") }

  data <- inner_join(x = deg_list[[1]] %>% select(Geneid, log2FoldChange, padj, DEG),
                     y = deg_list[[2]] %>% select(Geneid, log2FoldChange, padj, DEG),
                     by="Geneid")

  data <- data %>%
    mutate(color = if_else(log2FoldChange.x >= threshold & log2FoldChange.y >= threshold, "positivecorr",
                           if_else(log2FoldChange.x >= threshold & log2FoldChange.y <= -threshold, "negativecorr",
                                   if_else(log2FoldChange.x <= -threshold & log2FoldChange.y <= -threshold, "positivecorr",
                                           if_else(log2FoldChange.x <= -threshold & log2FoldChange.y >= threshold, "negativecorr", "ns")))))

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

  g <- g  +
    geom_point(alpha = 1) +
    geom_text_repel(data = data %>% filter(Geneid %in% genes),
                    mapping = aes(label = Geneid, x = log2FoldChange.x, y = log2FoldChange.y),
                    color = "black", size = 4) +

    coord_cartesian(xlim = xlim, ylim = ylim, expand = F) +
    geom_hline(yintercept = threshold, linetype = 2, color = "black") +
    geom_hline(yintercept = -threshold, linetype = 2, color = "black") +
    geom_vline(xintercept = threshold, linetype = 2, color = "black") +
    geom_vline(xintercept = -threshold, linetype = 2, color = "black") +
    scale_colour_manual(values=c("positivecorr" = "darkred", "negativecorr" = "darkgreen", "ns" = "gray50"),
                        drop = T) +

    ggtitle(main, subtitle) +
    xlab(xlab) + ylab(ylab) +

    theme_pubr(border = T, legend = "none") +
    theme(legend.title = element_blank(),
          plot.title = element_text(face="bold"),
          plot.subtitle = element_text(face="italic"),
          axis.title = element_text(face="bold"))


  # return
  return(g)
}

