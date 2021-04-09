######################################
## Upset plot for overlapping peaks ##
######################################

# upset_overlap_peaks() ---------------------------

#' @title upset_overlap_peaks
#' @author amitjavilaventura
#'
#' @seealso plyranges
#' @seealso UpSetR
#'
#' @usage upset_overlap_peaks(peaks, names = names(peaks), names_order = names)
#'
#' Takes a list of dataframes with regions.
#' Merges these dataframes to one dataframe with non-redundant regions.
#' Add new columns for each condition of the input list and writes a 1 if the peak is present in the corresponding condition.
#'
#' @param peaks List of dataframes with the genomic coordinates of the regions to overlap. Dataframes must contain the columns seqnames, start, end.
#' @param names Character with the same length as the 'peaks' list. The names that are given to the diferente objects in the 'peaks' list. Default: 'names(peaks)'
#' @param names_order Character of the same length as 'names'. Same values as in 'names' but with the desired order of priority. Default: 'names'
#' @param draw_plot Logical of length 1. If TRUE, the plot will be returned, if FALSE a dataframe of all the peaks with the conditions where that peak is found will be returned. Default: TRUE
#' @param ... Arguments to be passed through 'UpSetR::upset()'.


upset_overlap_peaks <- function(peaks, names = names(peaks), names_order = names, draw_plot = T, ...){

  require(purrr)
  require(dplyr)
  require(tidyr)
  require(magrittr)
  require(plyranges)
  require(UpSetR)

  if(is.null(names)) { stop("'names' must not be NULL")}

  len <- length(peaks)

  peaks <- peaks[names_order] %>%
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

  if(!draw_plot){
    x <- x %>% as.data.frame() %>% replace(is.na(.), 0) %>% select(-strand, -width) %>% unique()
    return(x)
  }
  else{

    x <- x %>%
      as.data.frame() %>%
      mutate(peak_id = paste(seqnames,start,end, sep = "_")) %>%
      select(peak_id, everything(), -seqnames, -start, -end, -strand, -width) %>%
      replace(is.na(.), 0) %>% unique()


    plot <- UpSetR::upset(data = x, order.by = "freq", ...)

    return(plot)
  }
}
