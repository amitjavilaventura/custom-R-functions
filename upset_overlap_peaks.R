######################################
## Upset plot for overlapping peaks ##
######################################

# upset_overlap_peaks() ---------------------------

#' @title upset_overlap_peaks
#' @author amitjavilaventura
#'
#' @usage upset_overlap_peaks(peaks, names = names(peaks), names_order = names)
#'
#' Takes a list of dataframes with regions.
#' Merges these dataframes to one dataframe with non-redundant regions.
#' Add new columns for each condition of the input list
#' It return a dataframe with all the coordinates their scores and the bedfile of origin.
#'
#' @param peaks List of dataframes with the genomic coordinates of the regions to overlap. Dataframes must contain the columns seqnames, start, end.
#' @param names Character with the same length as the 'peaks' list. The names that are given
#' @param names_order Character of the same length as 'regions'. Names of the regions in the bedfiles.
#' @param names_order Character with the same values as in 'names' but in the desired order for further plotting. Default: names
#' @param operation Caracter of length 1. Operation to calculate the score. One of c("sum", "mean", "max", "min"). Default: "mean"
#'


upset_overlap_peaks <- function(peaks, names = names(peaks), names_order = names){

  require(purrr)
  require(dplyr)
  require(tidyr)
  require(magrittr)
  require(plyranges)
  require(UpSetR)

  if(is.null(names)) { stop("'names' must not be NULL")}

  len <- length(peaks)

  peaks <- peaks %>%
    purrr::map(~dplyr::select(.data = .x, seqnames, start, end)) %>%
    purrr::set_names(nm = names) %>%
    purrr::map(~dplyr::mutate(.data = .x, condition = 1)) %>%
    purrr::imap(~set_colnames(x = .x, c("seqnames", "start", "end", .y))) %>%
    purrr::map(~as_granges(.x))


  x <- peaks[[len]]

  for(i in (len-1):1){
    x <- x %>%
      filter_by_non_overlaps(y = peaks[[i]]) %>%
      bind_ranges(peaks[[i]])
  }

  x <- x %>%
    as.data.frame() %>%
    select(seqnames, start, end) %>%
    as_granges()

  for(i in (len):1){
    x <- x %>%
      join_overlap_left(peaks[[i]])
  }

  x <- x %>%
    as.data.frame() %>%
    mutate(peak_id = paste(seqnames,start,end, sep = "_")) %>%
    select(peak_id, everything(), -seqnames, -start, -end, -strand, -width) %>%
    replace(is.na(.), 0) %>% unique()


  plot <- UpSetR::upset(data = x, order.by = "freq")

}
