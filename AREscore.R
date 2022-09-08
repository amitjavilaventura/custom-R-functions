# ============================================= #
# ARE score                                     #
# ============================================= #

#' @title arescore
#'
#' @description
#' An implementation of the algorithm to score ARE's as described here: http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1002433, http://arescore.dkfz.de/info.html
#' I took the function from hpgltools (https://rdrr.io/github/elsayed-lab/hpgltools/man/hpgl_arescore.html), which took the function from SeqTools (https://github.com/lianos/seqtools/blob/master/R/pkg/R/AREScore.R).
#' The function searches ATTTA/AUUUA pentamers in a DNA/RNA sequence and adds a score to the sequence based on the number of pentamers, the distance separating pairs of pentamers and wether they are found inside AU blocks.
#' The higher the score, the more AU-rich elements the sequence will have.
#'
#' @usage arescore(x, basal = 1, overlapping = 1.5, d1.3=0.75, d4.6=0.4, d7.9=0.2, within.AU = 0.3, aub.window.length = 10, aub.p.to.start = 0.8, aub.p.to.end = 0.55)
#'
#' @param x DNAStringSet or RNAStringSet. DNA/RNA StringSet containing the UTR sequences of interest
#' @param basal Numerical of length 1. Basal score for each pentamer found. Default = 1.
#' @param overlapping Numerical of length 1. Score for each pair of overlapping pentamers. Default = 1.5.
#' @param d1.3 Numerical of length 1. Score for each pair of pentamers separated by 1 to 3 nucleotides. Default = 0.75.
#' @param d4.6 Numerical of length 1. Score for each pair of pentamers separated by 4 to 6 nucleotides. Default = 0.4.
#' @param d7.9 Numerical of length 1. Score for each pair of pentamers separated by 7 to 9 nucleotides. Default = 0.2.
#' @param within.AU Numerical of length 1. Score for each pentamer fount in AU-blocks. Default = 0.3.
#' @param aub.window.length Numerical (integer) of length 1. The length of the sliding window to find AU-blocks. An AU-block will be minimum this size. Default = 10.
#' @param aub.p.to.start Numerical of length 1. Minimum percentage of AU in a sliding window to call it the start of an AU-block. Default = 0.8.
#' @param aub.p.to.end Numerical of length 1. Maximum percentage of AU in a sliding window to call it the end of the sliding block. Default = 0.55.
#' @param return_aub Logical of length 1. Whether to return the coordinates of the AU-blocks in each sequence as a IRangesList. Default = FALSE.
#'
#' @return a DataFrame (S4 object) of scores
#' @seealso [IRanges] [Biostrings] [GenomicRanges] [SeqTools]
#'
#' @export
arescore <- function(x, basal = 1, overlapping = 1.5, d1.3=0.75, d4.6=0.4, d7.9=0.2, within.AU = 0.3, aub.window.length = 10, aub.p.to.start = 0.8, aub.p.to.end = 0.55, return_aub = F) {

  # Load required packages
  require(Biostrings)
  require(SeqTools)

  # Define type of the sequence (dna or rna)
  xtype <- match.arg(substr(class(x), 1, 3), c("DNA", "RNA"))

  # Define pentamer for the au-rich elements (attta or auuua) depending on the type of sequence
  # Define also the overmer (two pentamers overlapping, etc)
  if (xtype == "DNA") {
    pentamer <- "ATTTA"
    overmer <- "ATTTATTTA"
  } else {
    pentamer <- "AUUUA"
    overmer <- "AUUUAUUUA"
  }

  # Convert the sequences to DNAstringset
  x <- as(x, "DNAStringSet")

  # Search for matches of the pentamer
  # Search for matches of the overmer
  pmatches <- Biostrings::vmatchPattern(pentamer, x)
  omatches <- Biostrings::vmatchPattern(overmer, x)

  # Add score to the sequence based on the number of pentamers per the basal score
  # Add score to the sequence based on the number of overmers per the overlapping score
  basal.score <- S4Vectors::elementNROWS(pmatches) * basal
  over.score <- S4Vectors::elementNROWS(omatches) * overlapping

  # Create a dataframe of no cluster
  no.cluster <- data.frame(d1.3 = 0, d4.6 = 0, d7.9 = 0)

  # Create a cluster dataframe and travel through the list of pattern matches
  # to stablish whether a sequence is a AU cluster or not:
  # - If there are less than 2 pentamers (i.e., 1 or 0), the sequence is not a au-cluster.
  # - Take the number of pentamers separated by a distance of 1-to-3, 4-to-6 or 7-to-9 nucleotides
  #   (these numbers will be used to compute the score afterwards).
  clust <- lapply(pmatches, function(m) {
    # If the number of p-matches is lower than 2, return no.cluster dataframe
    if (length(m) < 2) { return(no.cluster) }
    # Compute the distance between auuua pentamers in the sequence
    wg <- BiocGenerics::width(IRanges::gaps(m))
    # Create the clust dataframe by taking the number of gaps below distance 3, between 4 and 6, and between 7 and 9.
    data.frame(d1.3=sum(wg <= 3), d4.6=sum(wg >= 4 & wg <= 6), d7.9=sum(wg >= 7 & wg <= 9))
  })

  # Bind the clusters dataframes into one
  clust <- do.call(rbind, clust)

  # Compute the distance scores (distance between pentamers per the stablished score)
  dscores <- clust$d1.3 * d1.3 + clust$d4.6 * d4.6 + clust$d7.9 *  d7.9

  # Identify the AU-blocks using this custom function and create a IRanges object with these AU-blocks
  au.blocks <- identifyAUBlocks(x, window.length = aub.window.length, p.to.start = aub.p.to.start, p.to.end = aub.p.to.end)

  # Count the score of a pentamer overlapping with a AU-block
  aub.score <- sum(IRanges::countOverlaps(pmatches, au.blocks) * within.AU)

  # Sum all the scores: basal (number of pentamers), overlapping (pentamers overlapping), distance (distance between pentamers), au-block (pentamers overlapping with AU-block)
  score <- basal.score + over.score + dscores + aub.score
  # Create a dataframe with the score, number of pentamers, number of overlapping oentamers, the au-blocks and the number of aublocks
  ans <- S4Vectors::DataFrame(score = score,
                              n.pentamer = S4Vectors::elementNROWS(pmatches),
                              n.overmer = S4Vectors::elementNROWS(omatches),
                              au.blocks = au.blocks,
                              n.au.blocks = S4Vectors::elementNROWS(au.blocks))
  # Add this dataframe to the cluster dataframe created previously
  res <- cbind(ans, S4Vectors::DataFrame(clust))

  # Remove the IRangesList with the coordinates of the AU-blocks from the final result
  if (!return_aub) { res$au.blocks <- NULL }

  # Return result
  return(res)
}

# Identify AU blocks -------------------------- #

#' @title identifyAUBlocks
#'
#' @description
#' Function called inside arescore() that finds AU-blocks (AU-rich reagions)
#' arescore() is implementation of the algorithm to score ARE's as described here: http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1002433, http://arescore.dkfz.de/info.html
#' I took the function from hpgltools (https://rdrr.io/github/elsayed-lab/hpgltools/man/hpgl_arescore.html),
#'   which took the function from SeqTools (https://github.com/lianos/seqtools/blob/master/R/pkg/R/AREScore.R).
#'
#' @usage identifyAUBlocks(x, window.length = 20, p.to.start = 0.8, p.to.end = 0.55)
#'
#' @param x DNAStringSet or RNAStringSet. DNA/RNA StringSet containing the UTR sequences of interest
#' @param window.length Numerical (integer) of length 1. The length of the sliding window to find AU-blocks. An AU-block will be minimum this size. Default = 10.
#' @param p.to.start Numerical of length 1. Minimum percentage of AU in a sliding window to call it the start of an AU-block. Default = 0.8.
#' @param p.to.end Numerical of length 1. Maximum percentage of AU in a sliding window to call it the end of the sliding block. Default = 0.55.
#'
#' @return An IRangesList with the coordinates (start, end, width) of the AU-blocks in the sequence
#'
identifyAUBlocks <- function (x, window.length = 20, p.to.start = 0.8, p.to.end = 0.55) {

  # Define type of sequence (RNA or DNA)
  xtype <- match.arg(substr(class(x), 1, 3), c("DNA", "RNA"))

  # Check that inputs are OK.
  stopifnot(S4Vectors::isSingleNumber(window.length) & window.length >= 5 & window.length <= 50)
  stopifnot(S4Vectors::isSingleNumber(p.to.start) & p.to.start >= 0.5 & p.to.start <= 0.95)
  stopifnot(S4Vectors::isSingleNumber(p.to.end) & p.to.end >= 0.2 & p.to.end <= 0.7)
  stopifnot(p.to.start > p.to.end)

  # Define whether to search for AU or AT in AU-blocks (sequence is RNA or DNA)
  if (xtype == "DNA") { AU <- "AT" } else { AU <- "AU" }

  # Take x and convert to the RNA/DNA string set
  y <- as(x, sprintf("%sStringSet", xtype))

  # Compute the widths of the sequences
  widths <- BiocGenerics::width(x)

  # If SeqTools is not installed, install it from github.com/lianos/seqtools.
  # Then load it (if not loaded before)
  if (!"SeqTools" %in% installed.packages()) {
    message("Installing lianos/seqtools from github to get the find_au_start_end() C function.")
    test <- devtools::install_github("lianos/seqtools/R/pkg")
  }
  lib_result <- suppressMessages(requireNamespace("SeqTools"))
  att_result <- suppressMessages(try(attachNamespace("SeqTools"), silent = TRUE))

  # Create a function to iterate through sequences in the DNAstrinSet and
  # compute letter frequency in sliding window of min.length, which slides one nucleotide at a time,
  # and find the AU-blocks of each sequence
  fun <- function(i) {

    # Pick one of the sequence in the DNA/RNA string set
    one_seq <- x[[i]]

    # Check that width of the sequence is longer than lenght of sliding window
    if(length(one_seq) < window.length) { return(IRanges::IRanges()) }

    # Calculate the AU frequency (percentage) using a sliding window
    # If AU is null or has 0 rows, return an empty IRanges object
    # Otherwise convert AU to numeric
    au <- Biostrings::letterFrequencyInSlidingView(x = one_seq, view.width = window.length, letters = AU, as.prob = TRUE)
    if (is.null(au) | nrow(au) == 0) {  return(IRanges::IRanges()) }
    au <- as.numeric(au)

    # Take the windows with higher AU percentage than the specified in p.to.start and marks the first nucleotide
    # Take the windows with lower AU percentage than the specified in p.to.end and marks the last nucleotide
    can.start <- au >= p.to.start
    can.end <- au <= p.to.end

    # Call a C function from the SeqTools package to find the start and end of the AU-blocks
    posts <- .Call("find_au_start_end", au, p.to.start, p.to.end, PACKAGE = "SeqTools")

    # Create a IRanges object with the AU-blocks, then format and reduce it.
    blocks <- IRanges::IRanges(posts$start, posts$end + window.length -  1L)
    IRanges::end(blocks) <- ifelse(IRanges::end(blocks) > widths[i], widths[i], IRanges::end(blocks))
    blocks <- IRanges::reduce(blocks)

    # Return blocks
    return(blocks)
  }

  # Call the created function for all sequences
  au.blocks <- lapply(1:length(x), fun)

  # Create a IRangesList of the coordinates of the AU-blocks and return it
  ret <- IRanges::IRangesList(au.blocks)
  return(ret)
}


# Find start and end of AU-blocks ------------- #

#' @title find_au_start_end
#'
#' @description
#' Function to be called inside  identifyAUBlocks(). Not implemented yet.
#' This function is an implementation of the C-based (RCpp) function find_au_start_end from the package SeqTools (currently used in identifyAUBlocks()).
#'
#' @usage find_au_start_end(au, p.to.start, p.to.end)
#'
find_au_start_end <- function(au, p.to.start, p.to.end){

  in_block = F
  val = 0
  starts = c()
  ends = c()
  for (i in 1:length(au) ){
    val = au.pau
    if( !in_block & val >= p.to.start ){
      in_block = T
      starts <- c(starts, i)
    } else if ( in_block & val < p.to.end & any(i > starts) ){
      in_block = F
      ends <- c(ends, i)
    }
  }

  if(length(ends) == length(starts)-1) { ends = c(ends, length(au)) }

  if(length(ends) != length(starts)) { stop("find_au_start_end : the lengths of the start and end vectors don't match.") }

  res <- matrix(c(starts,ends), ncol = 2)

  return(res)
}


