
### findClosest()
### ===============================================================================================

#' @title findClosest()
#' @author amitjavilaventura
#'
#' @usage findClosest(x, y, k = 2, y_id = "V4")
#'
#' @param x Granges object.
#' @param y Granges object.
#' @param k Numeric of length 1. Integer number showing the number of nearest features in `y` to search for each element in `x`.
#' @param y_id Character of length 1. Name of the column in `y` that will be added to `x`. Usually, the ID of the feature (e.g. gene or peak).
#'

findClosest <- function(x, y, k = 2, y_id = "V4"){

  # Load required packages
  require(GenomicRanges)
  require(dplyr)
  require(tidyr)
  require(magrittr)

  if(!is.integer(k)){ stop("'k' must be an integer number (e.g. 1 or 1L") }

  # Search nearest neighbor elements
  nearest <- nearestKNeighbors(x = x, subject = y, k=k) %>%

    # Remove warnings
    pkgcond::suppress_warnings()

  # Take number of rows in the `x` object
  n_rows <- length(nearest)

  # Append the closests_cols dataframe to the x2 object using bind_cols
  x2 <- x %>% as_tibble()
  y2 <- y %>% as_tibble()

  # Build a dataframe with k columns and n_rows rows
  closests <- matrix(data = rep("", n_rows*k), ncol = k, nrow = n_rows) %>%
    as_tibble() %>%
    set_colnames(value = paste("closest", 1:k, sep=""))

  for(i in 1:k){
    for(j in 1:length(nearest)){
      closests[j, i] <- y2[nearest[[j]][i], y_id]
    }
  }

  x2 <- x2 %>% bind_cols(closests)

  return(x2)
}



### findClosestCoords()
### ===============================================================================================

#' @title findClosestCoords()
#' @author amitjavilaventura
#'
#' @usage findClosestCoords(x, y, k = 2, y_cols = c("seqnames", "start", "end", "strand", "V4")
#'
#' @param x Granges object.
#' @param y Granges object.
#' @param k Numeric of length 1. Integer number showing the number of nearest features in `y` to search for each element in `x`. Default: 2
#' @param y_cols Character of length 6. The names of the columns in `y` to add to `x`. These columns must contain, in this order: seqnames, start, end, width, ID and strand
#'
findClosestCoords <- function(x, y, k = 1,
                              x_cols = c("seqnames", "start", "end", "width", "V4","strand", "annotation"),
                              y_cols = c("seqnames", "start", "end", "width", "V4","strand")){

  # Load required packages
  require(plyranges)
  require(GenomicRanges)
  require(pkgcond)
  require(dplyr)
  require(tidyr)
  require(magrittr)

  # Check if inputs are correct
  if(!(class(x) %in% c("data.frame", "GRanges"))){ stop("'x' must be an object of class 'data.frame' or 'GRanges'")  }
  else if(!(class(y) %in% c("data.frame", "GRanges"))){ stop("'y' must be an object of class 'data.frame' or 'GRanges'")  }
  else if(!equals(k, as.integer(k))){ stop("'k' must be an integer number (e.g. 1 or 1L") }
  else if(is.data.frame(x)) {
    if(all(c("seqnames", "start", "end") %in% colnames(x))) { x <- x %>% as_granges() }
    else { stop("'x' must contain the columns 'seqnames', 'start' and 'end' to convert it to GRanges") }
  }
  else if(is.data.frame(y)) {
    if(all(c("seqnames", "start", "end") %in% colnames(y))) { y <- y %>% as_granges() }
    else { stop("'y' must contain the columns 'seqnames', 'start' and 'end' to convert it to GRanges") }
  }
  else if(!is.character(y_cols)) { stop("'y_cols' must a character vector of length 6.") }
  else if(length(y_cols) > 6) { stop("'y_cols' must a character vector of length 6.") }
  else if(!all(y_cols %in% colnames(as_tibble(y)))) { stop("All the elements in the 'y_cols' input must be colnames of 'y'") }


  # convert x and y to granges
  if(is.data.frame(x)){  x <- as_granges(x) }
  else if(is.data.frame(y)){  y <- as_granges(y) }

  # Search nearest neighbor elements
  nearest  <- nearestKNeighbors(x = x, subject = y, k=k) %>% pkgcond::suppress_warnings()

  # Re convert x and y to tibble
  x <- x %>% as_tibble() %>% dplyr::select(x_cols)
  y <- y %>% as_tibble() %>% dplyr::select(y_cols)

  nearest_list <- nearest %>%
    as.list() %>%
    purrr::map(~matrix(.x, ncol = k)) %>%
    purrr::map(~magrittr::set_colnames(.x, paste("k",1:k, sep = ""))) %>%
    do.call(rbind, .)

  nearest_list <- split(nearest_list, rep(1:ncol(nearest_list), each = nrow(nearest_list))) %>%
    purrr::set_names(paste("k",1:k, sep = "")) %>%
    purrr::map(~tibble(index_y = .x)) %>%
    purrr::map(~dplyr::mutate(.x, index_x = row_number()))

  # Print only the nearest peaks from y in the same order as x
  x_nearest_list <- list()
  for(i in 1:length(nearest_list)){
    y_list <- y[nearest_list[[i]]$index_y %>% na.omit(),]
    x_filt <- bind_cols(x, nearest_list[[i]]) %>% dplyr::filter(!is.na(index_y)) %>% dplyr::select(-contains("index"))
    x_nearest_list[[paste("k",i, sep = "")]] <- bind_cols(x_filt, y_list) %>%
      magrittr::set_colnames(c(paste(x_cols, "x", sep = "_"), paste(y_cols, "y", sep = "_")))
      pkgcond::suppress_warnings() %>% pkgcond::suppress_messages()
  }

  # Return x
  return(x_nearest_list)

}



