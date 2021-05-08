
### findDeNovoTargets()
### ===============================================================================================

#' @title findDeNovoTargets()
#' @author amitjavilaventura
#'
#' @usage findDeNovoTargets(x, y, k = 2, y_id = "V4")
#' 
#' @param x Data frame with columns seqnames, start, end, annotation (Promoter/Distal) and SYMBOL
#' @param y Data frame with columns seqnames, start, end, annotation (Promoter/Distal) and SYMBOL

findDeNovoTargets <- function(x, y){
  
  require(dplyr)
  require(plyranges)

  x <- x %>% dplyr::mutate(anno = if_else(str_detect(annotation, "Promoter"), "Promoter", "Distal"))
  y <- y %>% dplyr::mutate(anno = if_else(str_detect(annotation, "Promoter"), "Promoter", "Distal"))
  
  ## De novo Targeted promoter --> genes with a peak in the promoter only in X
  x.promo <- x %>% dplyr::filter(anno == "Promoter")
  y.promo <- y %>% dplyr::filter(anno == "Promoter")
  genes_target_denovo_prom <- x.promo[!(x.promo$SYMBOL %in%  y.promo$SYMBOL),] %>%
    dplyr::select(seqnames, "start" = geneStart, "end" = geneEnd, SYMBOL, geneLength, "strand" = geneStrand) %>% 
    unique()
  
  ## De novo Targeted by a distal peak --> closest gene to a distal peak in X that do not overlap with a peak in Y.
  x.distal <- x %>% dplyr::filter(anno == "Distal")
  y.distal <- y %>% dplyr::filter(anno == "Distal")
  
  genes_target_denovo_dist <- x.distal %>% as_granges() %>%
    filter_by_non_overlaps(y.distal %>% as_granges) %>%
    as_tibble() %>% 
    dplyr::select(seqnames, "start" = geneStart, "end" = geneEnd, SYMBOL, geneLength, "strand" = geneStrand) %>% 
    unique()
  
  
  ## De novo Targeted by all peaks
  genes_target_denovo_all <- bind_rows(genes_target_denovo_prom %>% dplyr::mutate(targeted = "Promoter"),
                                       genes_target_denovo_dist %>% dplyr::mutate(targeted = "Distal")) %>%
    unique()
  
  
  ## List of de novo targeted in X vs Y.
  target_denovo <- list("genes_target_denovo_all"  = genes_target_denovo_all,
                        "genes_target_denovo_prom" = genes_target_denovo_prom,
                        "genes_target_denovo_dist" = genes_target_denovo_dist) %>% 
    purrr::map(~dplyr::mutate(.x, denovo   = "De novo"))
  
  
  # Genes targeted by in both X and Y
  ## Targeted in the promoter
  genes_target_common_prom <- x.promo %>% 
    dplyr::filter((SYMBOL %in% y.promo$SYMBOL)) %>%
    dplyr::select(seqnames, "start" = geneStart, "end" = geneEnd, SYMBOL, geneLength, "strand" = geneStrand) %>% 
    unique()
  
  ## Targeted by a distal peak
  genes_target_common_dist <- x.distal %>% as_granges() %>%
    filter_by_overlaps(y.distal %>% as_granges) %>%
    as_tibble() %>% 
    dplyr::select(seqnames, "start" = geneStart, "end" = geneEnd, SYMBOL, geneLength, "strand" = geneStrand) %>% 
    unique()
  
  ## Targeted by promoter + distal
  genes_target_common_all <- bind_rows(genes_target_common_prom %>% dplyr::mutate(targeted = "Promoter"),
                                       genes_target_common_dist %>% dplyr::mutate(targeted = "Distal")) %>%
    unique()
  
  ## Lists of common targeted in X and Y
  target_common <- list("genes_target_common_all"  = genes_target_common_all,
                        "genes_target_common_prom" = genes_target_common_prom,
                        "genes_target_common_dist" = genes_target_common_dist) %>% 
    purrr::map(~dplyr::mutate(.x, denovo   = "Common"))
  
  # Genes targeted only in Y and not in X
  ## Targeted in the promoter
  genes_target_notarget_prom <- y.promo %>% 
    dplyr::filter(!(SYMBOL %in% x.promo$SYMBOL)) %>%
    dplyr::select(seqnames, "start" = geneStart, "end" = geneEnd, SYMBOL, geneLength, "strand" = geneStrand) %>% 
    unique()
  
  ## Targeted by a distal peak
  genes_target_notarget_dist <- y.distal %>% as_granges() %>%
    filter_by_non_overlaps(x.distal %>% as_granges) %>%
    as_tibble() %>% 
    dplyr::select(seqnames, "start" = geneStart, "end" = geneEnd, SYMBOL, geneLength, "strand" = geneStrand) %>% 
    unique()
  
  ## Targeted by promoter + distal
  genes_target_notarget_all <- bind_rows(genes_target_notarget_prom %>% dplyr::mutate(targeted = "Promoter"),
                                         genes_target_notarget_dist %>% dplyr::mutate(targeted = "Distal")) %>%
    unique()
  
  ## Lists of common targeted in X and Y
  target_notarget <- list("genes_target_notarget_all"  = genes_target_notarget_all,
                          "genes_target_notarget_prom" = genes_target_notarget_prom,
                          "genes_target_notarget_dist" = genes_target_notarget_dist) %>% 
    purrr::map(~dplyr::mutate(.x, denovo   = "Not target"))
  
  
  all_targeted <- c(target_denovo, target_common, target_notarget) %>% 
    purrr::map(~dplyr::mutate(.x, strand   = if_else(strand == 1, "+", "-")))
  
  # Return
  return(all_targeted)
}


### searchProm()
### ===============================================================================================

#' @title searchProm()
#' @author amitjavilaventura
#'
#' @usage searchProm(x, region = c(2500, 2500))
#' 
#' @param x Data frame with columns seqnames, start, end and strand
#' @param region Numeric of length 2 with the number of bases to take upstream and downstream of the TSS.

searchProm <- function(x, region = c(2500, 2500)){
  
  if(!all(c("+", "-") %in% x$strand)){ stop("The 'strand' variable of 'x' must have '+' and '-' as values.") }
  
  x.tss <- x %>%
    dplyr::mutate(tss = if_else(strand == "+" | strand == 1, start, end)) %>%
    dplyr::mutate(start = tss-region[1], end = tss+region[2]) %>%
    dplyr::select(seqnames, start, end, everything())
    
  return(x.tss)  
  
}
  

apc_wtn_targeted_bcat$genes_target_denovo_all->x
  
