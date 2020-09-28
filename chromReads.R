###############
# CHROM READS #
###############

#' @title chromReads
#' @author amitjavilaventura
#'
#' @usage chromReads(bamfile, text.size = 8, main = NULL, main.size = 13, subtitle = NULL, sub.size = 11, coord.flip = T, axis.x = T, axis.y = T, xlab = "Mapped reads", ylab = "Chromosome", x.size = 9, y.size = 9, legend = F)
#'
#' Function counts the reads mapped to each chromosome and plots a bar graph.
#' It uses the function idxstatsBam() from the Rsamtools R package. Run help("idxstatsBam") for more information.
#'
#' @param bamfile Character. Path and name of the BAM file whose reads are to be mapped.
#' @param main Character. Title of the pie chart. Default: NULL.
#' @param main.size Numerical. Font size of the title. It works only if main is not NULL. Default: 13.
#' @param subtitle Character. Subtitle of the piechart. It works only if main is not NULL. Default: NULL.
#' @param sub.size Numerical. Font size of the subtitle. It works only if main is not NULL. Default: 11.
#' @param coord.flip Logical. If the coordinates must be flipped around. Default: TRUE.
#' @param axis.x Logical. If TRUE, X axis is plotted. Default: TRUE.
#' @param axis.y Logical. If TRUE, Y axis is plotted. Default: TRUE.
#' @param xlab Character. Title of the X axis. Default: "Mapped reads".
#' @param ylab Character. Title of the Y axis. Default: "Chromosome".
#' @param legend Logical. If TRUE, legend is plotted. Default: FALSE.
#'
chromReads <- function(bamfile, text.size = 8, main = NULL, main.size = 13, subtitle = NULL, sub.size = 11,
                       axis.x = T, axis.y = T, xlab = "Mapped reads", ylab = "Chromosome", x.size = 9, y.size = 9,
                       legend = F, percent = T, genome = "mouse" , lab.size = 3){

  # Load required packages
  require(Rsamtools)
  require(ggplot2)
  require(dplyr)

  # Calculate the number of reads mapping to each chromosome with Rsamtools::idxstatsBam()
  chromReads <- idxstatsBam(bamfile)

  # Remove names of strange chromosomes
  chromReads <- chromReads[grep("chr", chromReads$seqnames),]

  # Calculate the total number of mapped reads in "good" chromosomes
  totalReads <- sum(chromReads$mapped)

  # Calculate percentage of mapped reads in each chromosome against all mapped reads
  chromReads$percentage <- chromReads$mapped/totalReads*100

  # Draw a bar graph
  b <- ggplot(data = chromReads, mapping = aes(seqnames, mapped, fill = seqnames)) +

    geom_bar(stat = "identity", show.legend = legend, colour = "Gray15") +

    # General formatting
    theme_bw() +

    # Text formatting
    theme(text = element_text(size = text.size)) +

    # Flip coordinates and set axis limits
    coord_flip(ylim = c(0, max(chromReads$mapped)*1.15))

  # Annotate labels (percentages) in the case we are working on mouse genome
  if(percent == T & genome == "mouse"){
    b <- b +
      annotate("text", label = paste(round(chromReads$percentage[chromReads$seqnames == "chr1"], 2), "%"), y = chromReads$mapped[chromReads$seqnames == "chr1"]+100000, x = 1, hjust = 0, size = lab.size, color = "grey30") +
      annotate("text", label = paste(round(chromReads$percentage[chromReads$seqnames == "chr2"], 2), "%"), y = chromReads$mapped[chromReads$seqnames == "chr2"]+100000, x = 2, hjust = 0, size = lab.size, color = "grey30") +
      annotate("text", label = paste(round(chromReads$percentage[chromReads$seqnames == "chr3"], 2), "%"), y = chromReads$mapped[chromReads$seqnames == "chr3"]+100000, x = 3, hjust = 0, size = lab.size, color = "grey30") +
      annotate("text", label = paste(round(chromReads$percentage[chromReads$seqnames == "chr4"], 2), "%"), y = chromReads$mapped[chromReads$seqnames == "chr4"]+100000, x = 4, hjust = 0, size = lab.size, color = "grey30") +
      annotate("text", label = paste(round(chromReads$percentage[chromReads$seqnames == "chr5"], 2), "%"), y = chromReads$mapped[chromReads$seqnames == "chr5"]+100000, x = 5, hjust = 0, size = lab.size, color = "grey30") +
      annotate("text", label = paste(round(chromReads$percentage[chromReads$seqnames == "chr6"], 2), "%"), y = chromReads$mapped[chromReads$seqnames == "chr6"]+100000, x = 6, hjust = 0, size = lab.size, color = "grey30") +
      annotate("text", label = paste(round(chromReads$percentage[chromReads$seqnames == "chr7"], 2), "%"), y = chromReads$mapped[chromReads$seqnames == "chr7"]+100000, x = 7, hjust = 0, size = lab.size, color = "grey30") +
      annotate("text", label = paste(round(chromReads$percentage[chromReads$seqnames == "chr8"], 2), "%"), y = chromReads$mapped[chromReads$seqnames == "chr8"]+100000, x = 8, hjust = 0, size = lab.size, color = "grey30") +
      annotate("text", label = paste(round(chromReads$percentage[chromReads$seqnames == "chr9"], 2), "%"), y = chromReads$mapped[chromReads$seqnames == "chr9"]+100000, x = 9, hjust = 0, size = lab.size, color = "grey30") +
      annotate("text", label = paste(round(chromReads$percentage[chromReads$seqnames == "chr10"], 2), "%"), y = chromReads$mapped[chromReads$seqnames == "chr10"]+100000, x = 10, hjust = 0, size = lab.size, color = "grey30") +
      annotate("text", label = paste(round(chromReads$percentage[chromReads$seqnames == "chr11"], 2), "%"), y = chromReads$mapped[chromReads$seqnames == "chr11"]+100000, x = 11, hjust = 0, size = lab.size, color = "grey30") +
      annotate("text", label = paste(round(chromReads$percentage[chromReads$seqnames == "chr12"], 2), "%"), y = chromReads$mapped[chromReads$seqnames == "chr12"]+100000, x = 12, hjust = 0, size = lab.size, color = "grey30") +
      annotate("text", label = paste(round(chromReads$percentage[chromReads$seqnames == "chr13"], 2), "%"), y = chromReads$mapped[chromReads$seqnames == "chr13"]+100000, x = 13, hjust = 0, size = lab.size, color = "grey30") +
      annotate("text", label = paste(round(chromReads$percentage[chromReads$seqnames == "chr14"], 2), "%"), y = chromReads$mapped[chromReads$seqnames == "chr14"]+100000, x = 14, hjust = 0, size = lab.size, color = "grey30") +
      annotate("text", label = paste(round(chromReads$percentage[chromReads$seqnames == "chr15"], 2), "%"), y = chromReads$mapped[chromReads$seqnames == "chr15"]+100000, x = 15, hjust = 0, size = lab.size, color = "grey30") +
      annotate("text", label = paste(round(chromReads$percentage[chromReads$seqnames == "chr16"], 2), "%"), y = chromReads$mapped[chromReads$seqnames == "chr16"]+100000, x = 16, hjust = 0, size = lab.size, color = "grey30") +
      annotate("text", label = paste(round(chromReads$percentage[chromReads$seqnames == "chr17"], 2), "%"), y = chromReads$mapped[chromReads$seqnames == "chr17"]+100000, x = 17, hjust = 0, size = lab.size, color = "grey30") +
      annotate("text", label = paste(round(chromReads$percentage[chromReads$seqnames == "chr18"], 2), "%"), y = chromReads$mapped[chromReads$seqnames == "chr18"]+100000, x = 18, hjust = 0, size = lab.size, color = "grey30") +
      annotate("text", label = paste(round(chromReads$percentage[chromReads$seqnames == "chr19"], 2), "%"), y = chromReads$mapped[chromReads$seqnames == "chr19"]+100000, x = 19, hjust = 0, size = lab.size, color = "grey30") +
      annotate("text", label = paste(round(chromReads$percentage[chromReads$seqnames == "chrX"], 2), "%"), y = chromReads$mapped[chromReads$seqnames == "chrX"]+100000, x = 20, hjust = 0, size = lab.size, color = "grey30") +
      annotate("text", label = paste(round(chromReads$percentage[chromReads$seqnames == "chrY"], 2), "%"), y = chromReads$mapped[chromReads$seqnames == "chrY"]+100000, x = 21, hjust = 0, size = lab.size, color = "grey30") +
      annotate("text", label = paste(round(chromReads$percentage[chromReads$seqnames == "chrM"], 2), "%"), y = chromReads$mapped[chromReads$seqnames == "chrM"]+100000, x = 22, hjust = 0, size = lab.size, color = "grey30")
  }

  # Formatting title and subtitle
  if(!is.null(main)){
    b <- b + ggtitle(label = main, subtitle = subtitle) +
      theme(plot.title = element_text(size = main.size, hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(size = sub.size, hjust = 0.5, face = "italic"))
  }

  # Formatting x axis
  if(axis.x){
    b <- b + ylab(label = xlab) + theme(axis.title.x = element_text(size = x.size))
  }
  else{
    b <- b + theme(axis.title.x = element_blank())
  }

  # Formatting y axis
  if(axis.y){
    b <- b + xlab(label = ylab) + theme(axis.title.y = element_text(size = y.size))
  }
  else{
    b <- b + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }

  # Formatting legend
  if(!legend){
    b <- b + theme(legend.position = "none")
  }
  else{
    b <- b + theme(legend.title = element_blank())
  }

  # Return bar graph
  b

}



