
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
findClosestCoords <- function(x, y, k = 2, y_cols = c("seqnames", "start", "end", "width", "V4","strand")){
  
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
  
  # Search nearest neighbor elements
  nearest <- nearestKNeighbors(x = x, subject = y, k=k) %>%
    
    # Remove warnings
    pkgcond::suppress_warnings()
  
  # Take number of rows in the `x` object
  n_rows <- length(nearest)

  # Re convert x and y to tibble
  x <- x %>% as_tibble()
  y <- y %>% as_tibble()

  ### This step could be done without nesting and unnesting, just create the list.
  ### But I leave it like this in case I find a way to do all the code with purrr and not for loops.
  # Build a nested list with n_rows number of lists
  closest <- rep(list(list("chr" = character(1), "start" = numeric(1), "end" = numeric(1), width = numeric(1),
                           "id" = character(1), "strand" = character(1))), times = n_rows) %>%
    
    # Create a tibble column with the nested list
    tibble(k = .) %>% 
    
    # Unnest the nested list to the wide format.
    tidyr::unnest_wider(col = k, names_sep = "_")
  
  
  # For each element nearest feature
  for(i in 1:k){
    
    # Change the colnames of the closest tibble to have the number of the ith closest feature.
    col_names <- names(closest) %>% gsub("k", paste("k", i, sep = ""), .)
    closest_k <- closest %>% set_colnames(value = col_names)
    
    # For each rew in 'x'
    for(j in 1:length(nearest)){
      
      # Take the corresponding columns in y for the ith closest observation in the jth column of x.
      # Append the values to each jth row of the closest_k tibble
      closest_k[j, ] <- y[nearest[[j]][i], y_cols]
      
    }
    
    # Bind the columns of the closest_k tibble to x
    x <- x %>% bind_cols(closest_k)
  }
  
  # Return x
  return(x)
  
}

findClosestCoords(x, y)
