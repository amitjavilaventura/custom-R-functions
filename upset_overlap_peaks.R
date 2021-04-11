######################################
## Upset plot for overlapping peaks ##
######################################

makeVenn4upSet <- function(peaks, conds = names(peaks), conds_order = conds, plot = F){

  require(pkgcond)
  require(purrr) %>% suppress_messages() %>% suppress_warnings()
  require(dplyr) %>% suppress_messages() %>% suppress_warnings()
  require(plyranges) %>% suppress_messages() %>% suppress_warnings()
  require(magrittr) %>% suppress_messages() %>% suppress_warnings()
  require(ChIPpeakAnno) %>% suppress_messages() %>% suppress_warnings()

  if(is.null(conds)) { stop("'conds' must not be NULL")}

  len <- length(peaks)

  overlaps <- peaks[conds_order] %>%
    purrr::map(~as_granges(.x)) %>%
    makeVennDiagram(plot = plot) %>%
    suppress_messages() %>% suppress_warnings()

  overlaps <- overlaps$vennCounts

  matrix <- matrix(data = rep(0, len), ncol = len, byrow=T) %>% set_colnames(conds_order)
  for(row in 1:nrow(overlaps)){

    counts <- overlaps[row, len+1]

    m <- matrix(rep(overlaps[row, 1:len], counts), ncol = len, byrow = T)

    matrix <- rbind(matrix, m)

  }

    x <- matrix %>%
      na.omit %>%
      as.data.frame() %>%
      dplyr::mutate(rowSum = rowSums(.)) %>%
      dplyr::filter(rowSum != 0) %>%
      dplyr::mutate(peak = paste("peak", 1:nrow(.), sep = "")) %>%
      dplyr::select(peak, everything(), -rowSum)

    return(list("matrix" = x, "vennCounts" = overlaps))

}

# upset_overlap_peaks() -------------------------

#' @title upset_overlap_peaks
#' @author amitjavilaventura
#'
#' @seealso plyranges
#' @seealso UpSetR
#'
#' @usage upset_overlap_peaks(peaks, names = names(peaks), names_order = names)
#'
#' @param peaks List of dataframes with the genomic coordinates of the regions to overlap. Dataframes must contain the columns seqnames, start, end.
#' @param conds Character with the same length as the 'peaks' list. The names that are given to the diferente objects in the 'peaks' list. Default: 'names(peaks)'
#' @param conds_order Character of the same length as 'conds'. Same values as in 'conds' but with the desired order of priority. Default: 'conds'.
#' @param title (Not available yet) Character of length 1 or NULL. Title of the plot shown at the top. Passed through grid.arrange(). If NULL, title is not shown. Default: NULL
#' @param order.by Same as in UpSetR::upset(). One of "freq", "degree" or both.
#' @param mainbar.y.label Same as in UpSetR::upset().
#' @param sets.x.label Same as in UpSetR::upset().
#' @param ... Further arguments to be passed through 'UpSetR::upset()'.

upsetOverlapPeaks <- function(peaks, conds = names(peaks), conds_order = conds,
                              title = NULL,
                              order.by = "freq", mainbar.y.label = "Intersect size", sets.x.label = "Set size",
                              ...){

  require(UpSetR)
  require(grid)
  require(gridExtra)

  # call makeVenn4upSet to retrieve the peak matrix
  x <- makeVenn4upSet(peaks, conds, conds_order, plot = F)

  # draw upset plot with options
  if(is.null(title)) {
    upset <- UpSetR::upset(data = x$matrix, order.by = order.by, sets.x.label = sets.x.label, mainbar.y.label = mainbar.y.label)
    return(upset)

  } else {
    plot.new()
    UpSetR::upset(data = x$matrix, order.by = order.by, sets.x.label = sets.x.label, mainbar.y.label = mainbar.y.label)
    grid.edit('arrange', name="upset")
    vp <- grid.grab()
    grid.arrange( grobs = list( vp ), top=title, cols=1 )
  }
}


### THIS FUNCTION DOES NOT WORK WELL, LOOK AT OTHER WAY TO COMPUTE THE OVERLAPS BECAUSE THIS WAY THERE ARE MANY MORE PEAKS THAN THE SUM.
### IT WORKS WITH PLYRANGES SO IT'S EXTREMELY FAST, BUT GIVES MORE PEAKS THAN EXPECTED.
### THE LIST OF PEAKS IS OK, BUT IN SOME CONDITIONS, EXCEPT THE FIRST 1 IN CONDS_ORDER, THERE ARE MORE OVERLAPS THAN PEAKS.
### ALTERNATIVE IS ChIPpeakAnno::findOverlapsOfPeaks or ChIPpeakAnno::makeVennDiagram, BUT THEY ARE SLOW (makeVennDiagram is less slow)
### KEEP TRYING WITH PLYRANGES.
### JUST SEEN THAT COMPLEX HEATMAP CAN DO THE OVERLAP OF GENOMIC REGIONS AND THE UPSET PLOT ... :-I
### COMPLEXHEATMAP::UPSET HAS MANY MORE OVERLAPS.
overlap4upset <- function(peaks, conds = names(peaks), conds_order = conds){

  require(purrr)
  require(dplyr)
  require(tidyr)
  require(magrittr)
  require(plyranges)

  if(is.null(conds)) { stop("'conds' must not be NULL")}

  len <- length(peaks)

  peaks <- peaks[conds_order] %>%
    purrr::map(~dplyr::select(.data = .x, seqnames, start, end)) %>%
    purrr::set_names(nm = conds) %>%
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

  x <- x %>% as.data.frame() %>% replace(is.na(.), 0) %>% select(-strand, -width) %>% unique()
  return(x)

}
