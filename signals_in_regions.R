###############################
## BW signals in BED regions ##
###############################

#' @title signal_in_regions
#' @author amitjavilaventura
#'
#' @usage signal_in_regions(bigwig, regions, names, names_order = names, operation = "mean")
#'
#' Function that takes a path to a bigwig files (bigwig) and several paths to a BED files (regions) to compute the signal in those regions.
#' It return a dataframe with all the coordinates their scores and the bedfile of origin.
#'
#' @param bigwig Character of length 1. Path to the desired bigwig file
#' @param regions Character of non defined length. Path(s) to the desired BED files.
#' @param names Character of the same length as 'regions'. Names of the regions in the bedfiles.
#' @param names_order Character with the same values as in 'names' but in the desired order for further plotting. Default: names
#' @param operation Caracter of length 1. Operation to calculate the score. One of c("sum", "mean", "max", "min"). Default: "mean"
#'
#'


signal_in_regions <- function(bigwig, regions, names, names_order = names, operation = "mean"){

  require(megadepth)
  require(dplyr)

  if(length(regions) != length(names)){ stop("'regions' and 'names' must have the same length") }

  list <- list()
  for(i in 1:length(regions)){

    coverage <- get_coverage(bigwig, regions[i], op = operation)  %>%
      as_tibble() %>% mutate(condition = names[i])

    list[[names[i]]] <- coverage
  }

  list <- bind_rows(list) %>%
    mutate(condition = factor(condition, levels = names_order))

  return(list)

}
