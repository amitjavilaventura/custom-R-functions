######################
### VENN DIAGRAM 2 ###
######################

#' @title vennDiagram2
#' @author amitjavilaventura
#'
#' @usage vennDiagram2(granges1, granges2, names = c("Peaks1", "Peaks2"), fill = c("Blue", "Red"), col = c("Darkblue", "Darkred"),cat.col = c("Darkblue", "Darkred"), cat.cex = c(1,1), cat.pos = c(1,11),label.col = c("Black", "Black", "Black"), alpha = 0.2)
#' 
#' Function for ChIP-seq and ATAC-seq.
#'
#' It takes two GRanges objects and compute the overlaps between them using findOverlapsOfPeaks() and draw.pairwise.venn(),
#' After that it will return an object that will have to be plotted with ggdraw(), ggarrange(), plot_grid() or similar.
#' 
#' @param granges1 GRanges object.
#' @param granges2 GRanges object.
#' @param names Character vector of length 2. Names of the 2 sets of peaks.
#' @param fill Caracter vector of length 2. Colors of the two circle areas.
#' @param col Caracter vector of length 2. Colors of the two circle circumferences.
#' @param cat.col Caracter vector of length 2. Colors of the two peak labels.
#' @param cat.cex Numeric vector of length 2. Size of the two peak labels.
#' @param cat.pos Numeric vector of length 2. Position of the two peak labels respect to their circles (see draw.pairwise.venn()).
#' @param label.col Caracter vector of length 3. Colors of the three area labels.
#' @param alpha Numeric. Transparency of the circle area.
#' 

vennDiagram2 <- function(granges1, granges2, names = c("Peaks1", "Peaks2"), 
                         fill = c("Blue", "Red"), col = c("Darkblue", "Darkred"),
                         cat.col = c("Darkblue", "Darkred"), cat.cex = c(1,1), cat.pos = c(1,11),
                         label.col = c("Black", "Black", "Black"), alpha = 0.2){

  # Load required packages 								 
  library(ChIPpeakAnno)
  library(VennDiagram)

  # Find overlaps between 2 sets of peaks 								 
  overlaps <- findOverlapsOfPeaks(granges1, granges2)
  
  # Draw Pairwise Venn object
  venn <- draw.pairwise.venn(area1 = sum(overlaps$venn_cnt[c(3,4), 3]), 
                             area2 = sum(overlaps$venn_cnt[c(2,4), 3]), 
                             cross.area = overlaps$venn_cnt[4, 3], 
                             category = names, scaled = T, 
                             fill = fill, col = col, alpha = alpha, 
                             cat.col = cat.col, cat.cex = cat.cex, cat.pos = cat.pos,
                             label.col = label.col)

   return(venn)
}