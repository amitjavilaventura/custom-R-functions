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

# ----- chromHMM functions ----- #
# read_delim_empty -> a read.delim()-based function that can read empty files.
source(here("R/chromHMM_functions.R"))

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
