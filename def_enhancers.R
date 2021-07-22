
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



### def_enhancers_no_k4me3()
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


### def_enhancers_k4me3_signal() DOES NOT WORK
### ===============================================================================================
#' @title def_enhancers_no_k4me3
#' @author amitjavilaventura
#'
#' It takes the k4me3 signal from bigwig files, then it computes the coverage score in the k4me3_promo and k4me3_distal peaks.
#' The k4me3_distal peaks that have a score higher than a certain score in k4me3_promo peaks (i.e. mean, median, q1, max...) are saved in a k4me3_distal_out list.
#' The k4me1_distal peaks are filtered and only those that do not overlap with the k4me3_distal_out peaks are retained.
#' Then k4me1_distal_filt peaks are overlapped (or non overlapped) with other mods to define active, primed and poised (poised only if k4me3_distal is not null) enhancers.
#'
#' @usage def_enhancers_k4me3_signal(k4me1_distal, k27ac_distal, k27me3_distal = NULL, conditions, k4me3_beds_promo, k4me3_beds_distal, k4me3_bws, coverage_op = "mean", filter_by_promo_score = "median")
#'
#' @returns A list of dataframes with the active, primed and poised enhancers
#'
#' @param k4me1_distal List of dataframes with the H3K4me1 distal peaks in different conditions. Dataframes must contain the columns seqnames, start, end and condition. Cannot be NULL.
#' @param k27ac_distal List of dataframes with the H3K27ac distal peaks in different conditions. Dataframes must contain the columns seqnames, start, end and condition. The order of the dataframes must be the same as in 'k4me1_distal'. Cannot be NULL.
#' @param k27me3_distal List of dataframes with the H3K27me3 distal peaks in different conditions. Dataframes must contain the columns seqnames, start, end and condition. The order of the dataframes must be the same as in 'k4me1_distal'. Can be NULL.
#' @param conditions Character vector with the names of the conditions present in the dataframes present in k4me1, k27ac and k27me3 arguments. Order of the conditions must be the same as in 'k4me1_distal'.
#' @param k4me3_beds_promo Character vector with the paths to the BED files of the H3K4me3 promoter peaks in each condition. Length must be equal to 'k4me1_distal'. Order of the conditions must be the same as in 'k4me1_distal'. The path must have the name of the conditions.
#' @param k4me3_beds_distal Character vector with the paths to the BED files of the H3K4me3 distal peaks in each condition. Length must be equal to 'k4me1_distal'. Order of the conditions must be the same as in 'k4me1_distal'. The path must have the name of the conditions.
#' @param k4me3_bws Character vector with the paths to the BW files of H3K4me3 in each condition. Order of the conditions must be the same as in 'k4me1_distal'. The path must have the name of the conditions.
#' @param coverage_op Character of length 1 indicating which operation should be used to compute the score from the BW file in each region. Passed through 'megadepth::get_coverage()'. One of c("mean", "sum", "min", "max"). Default: "mean"
#' @param filter_by_promo_score Either a character of lenght 1 with the measure of promoter scores or a numeric of lenght 1 indicating the percentile of the promoter scores to use to filter the distal peaks. If character, one of c("mean", "min", "q0", "q1, "median", "q2", "q3", "max", "q4"). If numeric, a value between 0 and 1. Default: "median"

def_enhancers_k4me3_signal <- function(k4me1_distal, k27ac_distal, k27me3_distal = NULL, conditions,
                                       k4me3_beds_promo, k4me3_beds_distal, k4me3_bws, coverage_op = "mean",
                                       filter_by_promo_score = "median" ) {

  # load packages
  require(dplyr)
  require(megadepth)

  # FILTER THE K4ME3 DISTAL PEAKS USING THE COVERAGE OF THE K4ME3 PROMO PEAKS. Retain those peaks that have to be removed from k4me1 distal peaks
  # =============================================================================================================================================

  # define an empty list for k4me3_distal_out peaks
  k4me3_distal_out_list  <- list()

  # for each condition filter the k4me3 distal peaks
  for (i in 1:length(conditions)) {

    # take the paths to the bed and bw files
    k4me3_distal <- k4me3_beds_distal[str_detect(k4me3_beds_distal, conditions[i])]
    k4me3_promo  <- k4me3_beds_promo[str_detect(k4me3_beds_promo, conditions[i])]
    k4me3_bw     <- k4me3_bws[str_detect(k4me3_bws, conditions[i])]

    # compute coverage in distal and promo peaks
    k4me3_coverage_distal <- get_coverage(bigwig_file = k4me3_bw, annotation = k4me3_distal, op = coverage_op)
    k4me3_coverage_promo  <- get_coverage(bigwig_file = k4me3_bw, annotation = k4me3_promo, op = coverage_op)

    # filter distal peaks and retain those
    # distal peaks with higher score than median score in promo peaks
    if(filter_by_promo_score == "mean"){
      mean_promo  <- k4me3_coverage_promo$score %>% mean
      k4me3_distal_out <- k4me3_coverage_distal %>% as_tibble() %>% dplyr::filter(score >= mean_promo)
    }

    # distal peaks with higher score than quantile x in promo peaks
    else if(filter_by_promo_score %in% c("mean", "min", "q0", "q1", "median", "q2", "q3", "max", "q4")){
      quantiles_promo  <- k4me3_coverage_promo$score %>% quantile()
      quantiles_distal <- k4me3_coverage_distal$score %>% quantile()

      if(filter_by_promo_score %in% c("min", "q0")) { k4me3_distal_out <- k4me3_coverage_distal %>% as_tibble() %>% dplyr::filter(score >= quantiles_promo[1]) }
      else if(filter_by_promo_score == "q1") { k4me3_distal_filt <- k4me3_distal_out %>% as_tibble() %>% dplyr::filter(score >= quantiles_promo[2]) }
      else if(filter_by_promo_score %in% c("median", "q2")) { k4me3_distal_out <- k4me3_coverage_distal %>% as_tibble() %>% dplyr::filter(score >= quantiles_promo[3]) }
      else if(filter_by_promo_score == "q3") { k4me3_distal_filt <- k4me3_distal_out %>% as_tibble() %>% dplyr::filter(score >= quantiles_promo[4]) }
      else if(filter_by_promo_score %in% c("max", "q4")) { k4me3_distal_out <- k4me3_coverage_distal %>% as_tibble() %>% dplyr::filter(score >= quantiles_promo[5]) }
    }

    # distal peaks with higher score than percentile x in promo peaks
    else if(is.numeric(filter_by_promo_score)){
      percentil_promo  <- k4me3_coverage_promo$score %>% quantile(probs = c(filter_by_promo_score))

      k4me3_distal_out <- k4me3_coverage_distal %>% as_tibble() %>% dplyr::filter(score >= percentil_promo) %>%
        dplyr::mutate(condition = conditions[i])

    } # end of if/else (k4me3 filtering out)

    # append each k4me3_distal_filt to the list
    k4me3_distal_out_list[[conditions[i]]] <- k4me3_distal_out %>% mutate(condition = conditions[i])

  } # end of loop for

  # FILTER THE K4me1 DISTAL PEAKS USING filter_by_non_overlaps() WITH THE K4me3_distal_out peaks.
  # =============================================================================================

  # define an empty list for k4me3_distal_out peaks
  k4me1_distal_filt_list  <- list()

  # for each condition filter_by_non_overlaps() with the k4me3_distal_out_peaks
  for(i in 1:length(conditions)){

    name   <- conditions[i]

    k4me1_distal_filt <- plyranges::filter_by_non_overlaps(x = k4me1_distal %>% bind_rows() %>% dplyr::filter(condition == name) %>% as_granges(),
                                                           y = k4me3_distal_out_list %>% bind_rows() %>% dplyr::filter(condition == name) %>% as_granges())

    k4me1_distal_filt_list[[name]] <- k4me1_distal_filt %>% as_tibble()

  }

  # OVERLAP THE MODS IN EACH CONDITION TO DEFINE ACTIVE, PRIMED AND POISED ENHANCERS (k27me3_distal not null)
  # =========================================================================================================

  # overlap the mods of each condition to define active, primed, and poised enhancers
  if(!is.null(k27me3_distal)){

    # define empty lists
    active_enhancers <- list()
    primed_enhancers <- list()
    poised_enhancers <- list()

    # for each condition
    for(i in 1:length(conditions)){

      # name of the conditions
      name   <- conditions[i]

      # active enhancers
      active <- plyranges::filter_by_overlaps(x = k4me1_distal_filt_list %>% bind_rows() %>% dplyr::filter(condition == name) %>% as_granges(),
                                              y = k27ac_distal %>% bind_rows() %>% dplyr::filter(condition == name) %>% as_granges()) %>%

        plyranges::filter_by_non_overlaps(y = k27me3_distal %>% bind_rows() %>% dplyr::filter(condition == name) %>% as_granges())


      # primed enhancers
      primed <- plyranges::filter_by_non_overlaps(x = k4me1_distal_filt_list %>% bind_rows() %>% dplyr::filter(condition == name) %>% as_granges(),
                                                  y = k27ac_distal %>% bind_rows() %>% dplyr::filter(condition == name) %>% as_granges()) %>%

        plyranges::filter_by_non_overlaps(y = k27me3_distal %>% bind_rows() %>% dplyr::filter(condition == name) %>% as_granges())

      # poised enhancers
      poised <- plyranges::filter_by_overlaps(x = k4me1_distal_filt_list %>% bind_rows() %>% dplyr::filter(condition == name) %>% as_granges(),
                                              y = k27me3_distal %>% bind_rows() %>% dplyr::filter(condition == name) %>% as_granges()) %>%

        plyranges::filter_by_non_overlaps(y = k27ac_distal %>% bind_rows() %>% dplyr::filter(condition == name) %>% as_granges())

      # fill lists
      active_enhancers[[paste(name, "active", sep = "-")]] <- active %>% as.data.frame() %>% dplyr::mutate(enhancer = "Active")
      primed_enhancers[[paste(name, "primed", sep = "-")]] <- primed %>% as.data.frame() %>% dplyr::mutate(enhancer = "Primed")
      poised_enhancers[[paste(name, "poised", sep = "-")]] <- poised %>% as.data.frame() %>% dplyr::mutate(enhancer = "Poised")
    }

    # merge active, primed and poised
    all_enhancers <- c(active_enhancers, primed_enhancers, poised_enhancers)

  }

  # OVERLAP THE MODS IN EACH CONDITION TO DEFINE ACTIVE AND PRIMED ENHANCERS (k27me3_distal IS null)
  # ================================================================================================
  else {

    # define empty lists
    active_enhancers <- list()
    primed_enhancers <- list()

    for(i in 1:length(conditions)){

      # name of the conditions
      name   <- conditions[i]

      # active enhancers
      active <- plyranges::filter_by_overlaps(x = k4me1_distal_filt_list %>% bind_rows() %>% dplyr::filter(condition == name) %>% as_granges(),
                                              y = k27ac_distal %>% bind_rows() %>% dplyr::filter(condition == name) %>% as_granges())


      # primed enhancers
      primed <- plyranges::filter_by_non_overlaps(x = k4me1_distal_filt_list %>% bind_rows() %>% dplyr::filter(condition == name) %>% as_granges(),
                                                  y = k27ac_distal %>% bind_rows() %>% dplyr::filter(condition == name) %>% as_granges())


      # fill lists
      active_enhancers[[paste(name, "active", sep = "-")]] <- active %>% as.data.frame() %>% dplyr::mutate(enhancer = "Active")
      primed_enhancers[[paste(name, "primed", sep = "-")]] <- primed %>% as.data.frame() %>% dplyr::mutate(enhancer = "Primed")
    }

    # merge active, primed and poised
    all_enhancers <- c(active_enhancers, primed_enhancers)

  }

  # return
  return(all_enhancers)

}


