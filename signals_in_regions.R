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
      as_tibble() %>% mutate(region = names[i])

    list[[names[i]]] <- coverage
  }

  if(merge){ list <- bind_rows(list) %>% mutate(region = factor(region, levels = names_order)) }

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
#' @param bigwigs Character of non defined length. Path(s) to the desired bigwig files
#' @param region Character of length 1. Path to the desired BED file.
#' @param names Character of the same length as 'bigwigs'. Names of the regions in the bedfiles.
#' @param names_order Character with the same values as in 'names' but in the desired order for further plotting. Default: names
#' @param operation Caracter of length 1. Operation to calculate the score. One of c("sum", "mean", "max", "min"). Default: "mean"
#' @param bind Character of length 1. Bind the coverages for the different bigwigs, either by rows or by columns. One of c("rows", "cols") or NULL. If NULL, coverages won't be merged and a list with different dataframes will be returned
#'
signals_in_1region <- function(bigwigs, region, names, names_order = names, operation = "mean", bind = "rows"){

  require(megadepth)
  require(dplyr)
  require(magrittr)

  if(length(bigwigs) != length(names)){ stop("'bigwigs' and 'names' must have the same length") }

  list <- list()
  for(i in 1:length(bigwigs)){

    coverage <- get_coverage(bigwigs[i], region, op = operation)

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

signals_in_regions <- function(bigwigs, regions, 
                               bw_names = names(bigwigs), bw_names_order = bw_names, 
                               bed_names = names(regions), bed_names_order = bed_names,
                               operation = "mean", merge = T){
  
  # Load packages
  require(megadepth)
  require(dplyr)
  
  # Check if inputs are OK
  if(!is.character(bigwigs)){ stop("'bigwigs' and must be a character vector with the path to each BIGWIG.") }
  else if(!is.character(regions) | !is.list(regions)){ stop("'regions' must be a character vector with the path to each BED or a lists of data frames with the seqnames, start and end columns")}
  else if(is.null(bw_names) | is.null(bed_names)){ stop("'bw_names' and 'bed_names' must be a not NULL character vector.") }
  else if(length(bigwigs) != length(bw_names)){ stop("'bigwigs' and 'bw_names' must have the same length.") }
  else if(length(regions) != length(bed_names)){ stop("'regions' and 'bed_names' must have the same length.") }
  
  # Write temporary files in case 'regions' is a list of dataframes.
  if(is.list(regions)){
    
    temp_dir <- tempdir(check = T)
    
    for(i in 1:length(regions)){
      if(!is.data.frame(regions[[i]])){ stop("If 'regions' is a list, each element in 'regions' must be a data frame with the seqnames, start and end columns. ") }
      else{
        name <- paste(bed_names[i], "_temp", sep = "")
        temporary <- tempfile(pattern = name, tmpdir = temp_dir, fileext = ".bed")
        regions[[i]] %>% write.table(file = temporary, quote = F, sep = "\t", row.names = F, col.names = F)
      }
    }
  
  # List the path of the temporary files
  regions <- list.files(path = temp_dir, pattern = ".bed", full.names = T, recursive = T)  
  
  }
  
  # Initialize list to store GRanges objects with the coverage
  signal_list <- list()
  
  # Go through each bigwig
  for(i in 1:length(bigwigs)){ 
    # Go through each bedfile
    for(j in 1:length(regions)){
    
      # Compute the coverage from each bigwig in each region in the bedfiles
      # The operation to compute the coverage in each region is the "mean" by default.
      # 'get_coverage()' returns a GRanges object
      coverage <- megadepth::get_coverage(bigwigs[i], regions[j], op = operation)  %>%
        # Convert GRanges to tibble/dataframe
        dplyr::as_tibble() %>% 
        # Create new columns for the names of condition of each bedfile and each bigwig
        dplyr::mutate(region = bed_names[j], signal_from = bw_names[i])
    
      # Store the coverage dataframe to the list
      signal_list[[paste(bw_names[i], "in", bed_names[j], sep="_")]] <- coverage
    }
  }
  
  # Conditional if merge == T
  # Merge the dataframes by rows
  # Converts the variables 'region' and 'signal_from' to factor
  # Order these variables by levels
  if(merge){ 
    signal_list <- dplyr::bind_rows(signal_list) %>% 
      dplyr::mutate(region = factor(region, levels = bed_names_order),
                    signal_from = factor(signal_from, levels = bw_names_order)) }
  
  # Return list
  return(signal_list)
  
}

#### some function to boxplot
#### some function to heatmap