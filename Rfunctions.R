######################
# R CUSTOM FUNCTIONS #
######################

# The paths may have to be changed depending on the current working directory
if(!require(here)){ install.packages("here") }; library(here)

# ----- chromReads ----- #
# To count and plot the reads mapped in each chromosome
source(here("R/chromReads.R"))

# ----- pcaplot2 ----- #
# A function to plot the principal components analysis
source(here("R/pcaplot2.R"))

# ----- volcanoplot2 ----- #
# A function to draw volcano plots
source(here("R/volcanoplot2.R"))


# ----- pieAnno3 ----- #
# pieAnno -> A function to draw pie charts with the peaks divided by promoter, gene body (UTRs, introns, exons) and distal (distal intergeninc, downstream).
# filterAnno -> A function that takes the output of Annotate peak and changes the features to promoter or distal.
source(here("R/pieAnno.R"))

# ----- barAnno ----- #
# barAnno -> a function to draw barplots from a list of annotatePeak objects. Can divide peaks in promoter/distal or promoter/gene body/distal
source(here("R/barAnno.R"))

# ----- vennDiagram2 ----- #
# vennDiagram2 -> a function to draw venn diagarams from two GenomicRanges objects.
source(here("R/vennDiagram2.R"))

# ----- read_delim_empty ----- #
# read_delim_empty -> a read.delim()-based function that can read empty files.
source(here("R/read_delim_empty.R"))

# ----- def_enhancers ----- #
# def_enhancers
# def_enhancers_no_k4me3
source(here("R/def_enhancers.R"))

# ----- DEcompare ----- #
# DEcompare
source(here("R/DEcompare.R"))

# ----- chromHMM functions ----- #
# already in the chromHMMviewR package
#source(here("R/chromHMM_functions.R"))






# ----- GOenrichment ----- #
# A function that does enrichGO
# Author: dfernandezperez
goEnrichment <- function(df, ont = "BP", db = org.Mm.eg.db) {
  require(clusterProfiler)
  require(db)
  ego <- enrichGO(gene          = df$ENTREZID,
                  OrgDb         = db,
                  keyType       = 'ENTREZID',
                  ont           = ont,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.1,
                  readable = TRUE)
  return(ego)
}


# ----- peakAnno ----- #
# A function that does annotatePeak
# Author: dfernandezperez
peakAnno <- function(infile, tssRegion = c(-2500, 2500), TxDb, annoDb) {

  # Load packages and set parameters
  require(ChIPseeker)
  require(TxDb, character.only = T)
  require(annoDb, character.only = T)

  # read file and annotate peaks
  peak <- infile %>% as_granges()
  annot <- annotatePeak(peak = peak, TxDb = get(TxDb), annoDb = annoDb, tssRegion = tssRegion)

  # The program classifies promotres in 1kb, 2kb... this line removes that annotation and leaves just "Promoters"
  annot@anno$annotation <- sub(" \\(.*\\)", "", annot@anno$annotation)

  # The program changes the type of object, so go back to GRanges (useful for downstream analysis)
  final <- as.GRanges(annot)
  return(final)
}


# ----- chip_stats ----- #
# funtion that makes a table from a list of df with peak annotation
chip_stats <- function(df, conditions = NULL){

  require(dplyr)

  if(is.null(conditions)){conditions <- names(df)}

  stats <- tibble(Condition = conditions,
                  Peaks     = df %>% purrr::map_dbl(nrow),
                  Targets   = df %>% purrr::map(~dplyr::select(.data = .x, SYMBOL)) %>% purrr::map(~unique(x = .x)) %>%  purrr::map_dbl(nrow))

  return(stats)

}
