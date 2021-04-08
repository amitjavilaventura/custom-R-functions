###############################
## BW signals in BED regions ##
###############################

# signal_in_regions() ---------------------------

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
signal_in_regions <- function(bigwig, regions, names, names_order = names, operation = "mean", merge = T){

  require(megadepth)
  require(dplyr)

  if(length(regions) != length(names)){ stop("'regions' and 'names' must have the same length") }

  list <- list()
  for(i in 1:length(regions)){

    coverage <- get_coverage(bigwig, regions[i], op = operation)  %>%
      as_tibble() %>% mutate(condition = names[i])

    list[[names[i]]] <- coverage
  }

  if(merge){ list <- bind_rows(list) %>% mutate(condition = factor(condition, levels = names_order)) }

  return(list)

}


# signals_in_1region() --------------------------

#' @title signals_in_1region -- to test
#' @author amitjavilaventura
#'
#' @usage signals_in_1region(bigwigs, region, names, names_order = names, operation = "mean")
#'
#' Function that takes a path to a bigwig files (bigwig) and several paths to a BED files (regions) to compute the signal in those regions.
#' It return a dataframe with all the coordinates their scores and the bedfile of origin.
#'
#' @param bigwigs Character of length non defined length. Path(s) to the desired bigwig files
#' @param region Character of length 1. Path to the desired BED file.
#' @param names Character of the same length as 'bigwigs'. Names of the regions in the bedfiles.
#' @param names_order Character with the same values as in 'names' but in the desired order for further plotting. Default: names
#' @param operation Caracter of length 1. Operation to calculate the score. One of c("sum", "mean", "max", "min"). Default: "mean"
#' @param bind Character of length 1. Bind the coverages for the different bigwigs, either by rows or by columns. One of c("rows", "cols") or NULL. If NULL, coverages won't be merged and a list with different dataframes will be returned
#'
signals_in_1region <- function(bigwigs, regions, names, names_order = names, operation = "mean", bind = "rows"){

  require(megadepth)
  require(dplyr)
  require(magrittr)

  if(length(bigwigs) != length(names)){ stop("'bigwigs' and 'names' must have the same length") }

  list <- list()
  for(i in 1:length(bigwigs)){

    coverage <- get_coverage(bigwigs[i], regions, op = operation)

    if(bind == "rows"){ coverage <- coverage %>% as_tibble() %>% mutate(signal_from = names[i]) }
    else if(bind == "cols") { coverage <- coverage %>%
      as_tibble() %>%
      dplyr::select(seqnames, start, end, score) %>%
      set_colnames(c("seqnames", "start", "end", names[i]))
    }
    else{ coverage <- coverage %>% as_tibble() }


    list[[names[i]]] <- coverage
  }

  if(bind == "rows"){
    list <- bind_rows(list) %>% mutate(signal_from = factor(signal_from, levels = names_order))

    return(list)

  } else if(bind == "cols"){
    list <- list %>%
      purrr::map(~dplyr::mutate(.data = .x, coords = paste(seqnames,start,end, sep="_"))) %>%
      purrr::map(~dplyr::select(.data = .x, -seqnames,-start, -end))

    join <- list[[1]]
    for(i in 2:(length(list)-1)){
      join <- inner_join(join, list[[i+1]], by = "coords")
    }

    return(join)

  } else{ return(list) }

}


#### some function to boxplot
#### some function to heatmap
