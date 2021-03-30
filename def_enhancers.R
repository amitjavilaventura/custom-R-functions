
### def_enhancers() DOES NOT WORK
### ===============================================================================================

#' @title def_enhancers
#' @author amitjavilaventura
#'
#' @usage def_enhancers(k4me1, k27ac, k27me3 = NULL, conditions)
#'
#' @returns A list of dataframes with the active, primed and poised enhancers
#'
#' @param k4me1 List of dataframes with the H3K4me1 distal peaks in different conditions. Dataframes must contain the columns seqnames, start, end and condition. Cannot be NULL.
#' @param k27ac List of dataframes with the H3K27ac distal peaks in different conditions. Dataframes must contain the columns seqnames, start, end and condition. Cannot be NULL.
#' @param k27me3 List of dataframes with the H3K27me3 distal peaks in different conditions. Dataframes must contain the columns seqnames, start, end and condition. Can be NULL.
#' @param conditions Character vector with the names of the conditions present in the dataframes present in k4me1, k27ac and k27me3 arguments

def_enhancers  <- function(k4me1, k27ac, k27me3 = NULL, conditions){

  # Load packages
  library(dplyr)
  library(plyranges)

  if(!is.null(k27me3)){

    # define empty lists
    active_enhancers <- list()
    poised_enhancers <- list()
    primed_enhancers <- list()

    for(i in 1:length(conditions)){

      # name of the conditions
      name   <- conditions[i]

      # active enhancers
      active <- plyranges::filter_by_overlaps(x = k4me1 %>% bind_rows() %>%
                                               dplyr::filter(condition == name) %>% as_granges(),
                                              y = k27ac %>% bind_rows() %>%
                                                dplyr::filter(condition == name) %>% as_granges()) %>%

        plyranges::filter_by_non_overlaps(y = k27me3 %>% bind_rows() %>%
                                            dplyr::filter(condition == name) %>% as_granges())


      # primed enhancers
      primed <- plyranges::filter_by_non_overlaps(x = k4me1 %>% bind_rows() %>%
                                                    dplyr::filter(condition == name) %>% as_granges(),
                                                  y = k27ac %>% bind_rows() %>%
                                                    dplyr::filter(condition == name) %>% as_granges()) %>%

        plyranges::filter_by_non_overlaps(y = k27me3 %>% bind_rows() %>%
                                            dplyr::filter(condition == name) %>% as_granges())

      # poised enhancers
      poised <- plyranges::filter_by_overlaps(x = k4me1 %>% bind_rows() %>%
                                                dplyr::filter(condition == name) %>% as_granges(),
                                              y = k27me3 %>% bind_rows() %>%
                                                dplyr::filter(condition == name) %>% as_granges()) %>%

        plyranges::filter_by_non_overlaps(y = k27ac %>% bind_rows() %>%
                                            dplyr::filter(condition == name) %>% as_granges())

      # fill lists
      active_enhancers[[paste(name, "active", sep = "-")]] <- active %>% as.data.frame() %>% dplyr::mutate(enhancer = "Active")
      primed_enhancers[[paste(name, "primed", sep = "-")]] <- primed %>% as.data.frame() %>% dplyr::mutate(enhancer = "Primed")
      poised_enhancers[[paste(name, "poised", sep = "-")]] <- poised %>% as.data.frame() %>% dplyr::mutate(enhancer = "Poised")
    }

    # merge active, primed and poised
    all_enhancers <- c(active_enhancers, primed_enhancers, poised_enhancers)

    # return
    return(all_enhancers)

  } else {

    # define empty lists
    active_enhancers <- list()
    poised_enhancers <- list()
    primed_enhancers <- list()

    for(i in 1:length(conditions)){

      # name of the conditions
      name   <- conditions[i]

      # active enhancers
      active <- plyranges::filter_by_overlaps(x = k4me1 %>% bind_rows() %>%
                                                dplyr::filter(condition == name) %>% as_granges(),
                                              y = k27ac %>% bind_rows() %>%
                                                dplyr::filter(condition == name) %>% as_granges())

      # primed enhancers
      primed <- plyranges::filter_by_non_overlaps(x = k4me1 %>% bind_rows() %>%
                                                    dplyr::filter(condition == name) %>% as_granges(),
                                                  y = k27ac %>% bind_rows() %>%
                                                    dplyr::filter(condition == name) %>% as_granges())


      # fill lists
      active_enhancers[[paste(name, "active", sep = "-")]] <- active %>% as.data.frame() %>% dplyr::mutate(enhancer = "Active")
      primed_enhancers[[paste(name, "primed", sep = "-")]] <- primed %>% as.data.frame() %>% dplyr::mutate(enhancer = "Primed")
    }

    # merge active, primed and poised
    all_enhancers <- c(active_enhancers, primed_enhancers)

    # return
    return(all_enhancers)

  }


}



### def_enhancers_no_k4me3() DOES NOT WORK
### ===============================================================================================

#' @title def_enhancers_no_k4me3
#' @author amitjavilaventura
#'
#' @usage def_enhancers_no_k4me3(k4me1, k27ac, k4me3, k27me3 = NULL, conditions)
#'
#' @returns A list of dataframes with the active, primed and poised enhancers
#'
#' @param k4me1 List of dataframes with the H3K4me1 distal peaks in different conditions. Dataframes must contain the columns seqnames, start, end and condition. Cannot be NULL.
#' @param k27ac List of dataframes with the H3K27ac distal peaks in different conditions. Dataframes must contain the columns seqnames, start, end and condition. Cannot be NULL.
#' @param k4me3 List of dataframes with the H3K27me3 distal peaks in different conditions. Dataframes must contain the columns seqnames, start, end and condition. Cannot be NULL.
#' @param k27me3 List of dataframes with the H3K27me3 distal peaks in different conditions. Dataframes must contain the columns seqnames, start, end and condition. Can be NULL.
#' @param conditions Character vector with the names of the conditions present in the dataframes present in k4me1, k27ac and k27me3 arguments

def_enhancers_no_k4me3  <- function(k4me1, k27ac, k4me3, k27me3 = NULL, conditions){

  # Load packages
  library(dplyr)
  library(plyranges)

  if(!is.null(k27me3)){

    # define empty lists
    active_enhancers <- list()
    poised_enhancers <- list()
    primed_enhancers <- list()

    for(i in 1:length(conditions)){

      # name of the conditions
      name   <- conditions[i]

      # active enhancers
      active <- plyranges::filter_by_overlaps(x = k4me1 %>% bind_rows() %>%
                                                dplyr::filter(condition == name) %>% as_granges(),
                                              y = k27ac %>% bind_rows() %>%
                                                dplyr::filter(condition == name) %>% as_granges()) %>%

        plyranges::filter_by_non_overlaps(y = k27me3 %>% bind_rows() %>%
                                            dplyr::filter(condition == name) %>% as_granges()) %>%

        plyranges::filter_by_non_overlaps(y = k4me3 %>% bind_rows() %>%
                                            dplyr::filter(condition == name) %>% as_granges())


      # primed enhancers
      primed <- plyranges::filter_by_non_overlaps(x = k4me1 %>% bind_rows() %>%
                                                    dplyr::filter(condition == name) %>% as_granges(),
                                                  y = k27ac %>% bind_rows() %>%
                                                    dplyr::filter(condition == name) %>% as_granges()) %>%

        plyranges::filter_by_non_overlaps(y = k27me3 %>% bind_rows() %>%
                                            dplyr::filter(condition == name) %>% as_granges()) %>%

        plyranges::filter_by_non_overlaps(y = k4me3 %>% bind_rows() %>%
                                            dplyr::filter(condition == name) %>% as_granges())

      # poised enhancers
      poised <- plyranges::filter_by_overlaps(x = k4me1 %>% bind_rows() %>%
                                                dplyr::filter(condition == name) %>% as_granges(),
                                              y = k27me3 %>% bind_rows() %>%
                                                dplyr::filter(condition == name) %>% as_granges()) %>%

        plyranges::filter_by_non_overlaps(y = k27ac %>% bind_rows() %>%
                                            dplyr::filter(condition == name) %>% as_granges()) %>%

        plyranges::filter_by_non_overlaps(y = k4me3 %>% bind_rows() %>%
                                            dplyr::filter(condition == name) %>% as_granges())

      # fill lists
      active_enhancers[[paste(name, "active", sep = "-")]] <- active %>% as.data.frame() %>% dplyr::mutate(enhancer = "Active")
      primed_enhancers[[paste(name, "primed", sep = "-")]] <- primed %>% as.data.frame() %>% dplyr::mutate(enhancer = "Primed")
      poised_enhancers[[paste(name, "poised", sep = "-")]] <- poised %>% as.data.frame() %>% dplyr::mutate(enhancer = "Poised")
    }

    # merge active, primed and poised
    all_enhancers <- c(active_enhancers, primed_enhancers, poised_enhancers)

    # return
    return(all_enhancers)

  } else {

    # define empty lists
    active_enhancers <- list()
    poised_enhancers <- list()
    primed_enhancers <- list()

    for(i in 1:length(conditions)){

      # name of the conditions
      name   <- conditions[i]

      # active enhancers
      active <- plyranges::filter_by_overlaps(x = k4me1 %>% bind_rows() %>%
                                                dplyr::filter(condition == name) %>% as_granges(),
                                              y = k27ac %>% bind_rows() %>%
                                                dplyr::filter(condition == name) %>% as_granges()) %>%

        plyranges::filter_by_non_overlaps(y = k4me3 %>% bind_rows() %>%
                                            dplyr::filter(condition == name) %>% as_granges())

      # primed enhancers
      primed <- plyranges::filter_by_non_overlaps(x = k4me1 %>% bind_rows() %>%
                                                      dplyr::filter(condition == name) %>% as_granges(),
                                                    y = k27ac %>% bind_rows() %>%
                                                      dplyr::filter(condition == name) %>% as_granges()) %>%

        plyranges::filter_by_non_overlaps(y = k4me3 %>% bind_rows() %>%
                                            dplyr::filter(condition == name) %>% as_granges())


      # fill lists
      active_enhancers[[paste(name, "active", sep = "-")]] <- active %>% as.data.frame() %>% dplyr::mutate(enhancer = "Active")
      primed_enhancers[[paste(name, "primed", sep = "-")]] <- primed %>% as.data.frame() %>% dplyr::mutate(enhancer = "Primed")
    }

    # merge active, primed and poised
    all_enhancers <- c(active_enhancers, primed_enhancers)

    # return
    return(all_enhancers)

  }


}







