########################
### read_delim_empty ###
########################


# function that will read a file without lines.
# it has been made to be able to list large lists of files using purrr::map
#  and some files may have no lines but are still useful to make plots.

#' @title read_delim_empty
#' @author amitjavilaventura
#'
#' @usage read_delim_empty(file, header = F, sep="\t", dec=".")
#'
#' Function that reads files using read.delim(), with the difference that, if there is a file that have no lines, it reads that file too.
#' This function is useful when reading a list of files with purrr::map (maybe apply()?), when you need to read all the files in order to build plots, etc.
#'
#' @param file character. Path of the file to be read.
#' @param header logical. Whether to read first line as header or not. Default: FALSE
#' @param sep character. Separator that separates the different fields or columns of the file to be read.
#' @param dec character. Symbol to use as decimal separator.
#'

read_delim_empty <- function(file, header = F, sep="\t", dec="."){

  info <- file.info(file)

  if(info$size == 0){
    a <- data.frame(0)
    return(a)
  }
  else{
    a <- read.delim(file = file, header = header, sep = sep, dec = dec)
    return(a)
  }
}
