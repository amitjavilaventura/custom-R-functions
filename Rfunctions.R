######################
# R CUSTOM FUNCTIONS #
######################

# The paths may have to be changed depending on the current working directory
# I use here, because it sets the working directory in the current project directory where I have the folder R with the Rfunctions.R and other scripts
if(!require(here)){ install.packages("here") }; library(here)

# ----- general helpers ----- #
# AHelper functions
source(here::here("R/general_helpers.R"))

# ----- pcaplot2 ----- #
# A function to plot the principal components analysis
source(here::here("R/pcaplot2.R"))

# ----- deseq2_helpers ----- #
# A helper function for deseq
source(here::here("R/deseq2_helpers.R"))

# ----- count_reads_table.R ----- #
# count reads
source(here::here("R/count_reads_table.R"))

# ----- pieAnno3 ----- #
# pieAnno -> A function to draw pie charts with the peaks divided by promoter, gene body (UTRs, introns, exons) and distal (distal intergeninc, downstream).
# filterAnno -> A function that takes the output of Annotate peak and changes the features to promoter or distal.
source(here::here("R/pieAnno.R"))


# ----- read_delim_empty ----- #
# read_delim_empty -> a read.delim()-based function that can read empty files.
source(here::here("R/read_delim_empty.R"))

# ----- def_enhancers ----- #
# def_enhancers
# def_enhancers_no_k4me3
source(here::here("R/def_enhancers.R"))



# ----- Signals in regions ----- #
# function(s) that take a bigwig(s) and a bed(s) and compute the signal in the desired regions
source(here::here("R/signals_in_regions.R"))


# ----- FindDeNovoTargets ----- #
source(here::here("R/findDeNovoTargets.R"))

# ----- Setdiff ----- #
# Functions from dfernandezperez to find unions and intersections of names
source(here::here("R/setDiff.R"))


# ----- Ggplot Helpers ----- #
# Functions to use with ggplot2 and customize the plots easily
source(here::here("R/ggplot_helpers.R"))

# ----- piC Helpers ----- #
# Functions to use with ggplot2 and customize the plots easily
source(here::here("R/pic_functions.R"))


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
# author: amitjavilaventura
# funtion that makes a table from a list of df with peak annotation
chip_stats <- function(df, conditions = NULL){

  require(dplyr)

  if(is.null(conditions)){conditions <- names(df)}

  stats <- tibble(Condition = conditions,
                  Peaks     = df %>% purrr::map_dbl(nrow),
                  Targets   = df %>% purrr::map(~dplyr::select(.data = .x, SYMBOL)) %>% purrr::map(~unique(x = .x)) %>%  purrr::map_dbl(nrow))

  return(stats)

}


