### FUNCTION TO COUNT READS FROM FASTQS, BAMFILES AND FCOUNTS OUTPUTS

fqReads <- function(fqcDir, readsCol = "Reads", returnAll = F){

  # Load packages
  require(dplyr)
  require(magrittr)
  require(fastqcr)

  fqc_info <- fastqcr::qc_aggregate(qc.dir = fqcDir, progressbar = F) %>% suppressMessages()
  fqc_stat <- fastqcr::qc_stats(object = fqc_info)
  fqc_num  <- fqc_stat %>% dplyr::select(sample, tot.seq) %>% magrittr::set_colnames(c("Sample", readsCol))

  fqc_num  <- fqc_num %>% dplyr::mutate(Sample = Sample %>% gsub(".trim", "", .) %>% gsub(".filt", "", .))

  if(returnAll){
    output   <- list("BasicStats" = fqc_stat, "NumReads" = fqc_num, "AllInfo" = fqc_info)
    return(output)
  }

  return(fqc_num)

}

bamReads <- function(bamDir, bamPattern = "bam$", bamNames = NULL, readsCol = "align_reads"){

  # Load packages
  require(dplyr)
  require(purrr)
  require(tibble)
  require(magrittr)
  require(Rsamtools)

  bam_paths <- list.files(path = bamDir, pattern = bamPattern, full.names = T, recursive = T)

  if(is.null(bamNames)) { bamNames <- bam_paths %>% purrr::map(~basename(.x) %>% gsub(".bam", "", .)) }
  bam_paths <- bam_paths %>% purrr::set_names(bamNames)

  bam_stats <- bam_paths %>% purrr::map(~Rsamtools::idxstatsBam(.x)) %>%
    purrr::map(~dplyr::select(.x, mapped) %>%
                 sum() %>%
                 tibble::tibble(bamreads = .) %>%
                 magrittr::set_colnames(readsCol)) %>%
    purrr::imap(~dplyr::mutate(.x, Sample = .y)) %>%
    dplyr::bind_rows()

  return(bam_stats)

}


fcountsNum <- function(countsDir, countsPat = "counts.summary$", countsNames = NULL, countsCol = "fcounts"){

  require(dplyr)
  require(purrr)
  require(stringr)
  require(magrittr)
  require(tibble)

  fcounts_paths <- list.files(path = countsDir,  pattern = countsPat, full.names = T, recursive = T)

  if(is.null(countsNames)) { countsNames <- fcounts_paths %>% purrr::map(~basename(.x) %>% gsub("\\..*", "", .)) }
  fcounts_paths <- fcounts_paths %>% purrr::set_names(countsNames)

  fcounts_info  <- purrr::map(fcounts_paths, read.delim, header = T) %>%
    purrr::map(~dplyr::filter(.x, Status == "Assigned")) %>%
    purrr::map(~magrittr::set_colnames(.x, c("Status", countsCol))) %>%
    purrr::imap(~dplyr::mutate(.x, Sample = .y)) %>%
    purrr::map(~dplyr::select(.x, Sample, countsCol)) %>%
    dplyr::bind_rows()

  return(fcounts_info)
}
