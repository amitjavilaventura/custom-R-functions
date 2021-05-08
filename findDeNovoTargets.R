
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
  
  x <- x %>% dplyr::mutate(anno = if_else(str_detect(annotation, "Promoter"), "Promoter", "Distal"))
  y <- y %>% dplyr::mutate(anno = if_else(str_detect(annotation, "Promoter"), "Promoter", "Distal"))
  
  ## De novo Targeted promoter --> genes with a peak in the promoter only in X
  x.promo <- x %>% dplyr::filter(anno == "Promoter")
  y.promo <- y %>% dplyr::filter(anno == "Promoter")
  genes_target_denovo_prom <- x.promo %>% 
    dplyr::filter(!(SYMBOL %in% y.promo$SYMBOL)) %>%
    dplyr::select(geneChr, geneStart, geneEnd, SYMBOL, geneLength, geneStrand) %>% 
    unique()
  
  ## De novo Targeted by a distal peak --> closest gene to a distal peak in X that do not overlap with a peak in Y.
  x.distal <- x %>% dplyr::filter(anno == "Distal")
  y.distal <- y %>% dplyr::filter(anno == "Distal")
  
  genes_target_denovo_dist <- x.distal %>% as_granges() %>%
    filter_by_non_overlaps(y.distal %>% as_granges) %>%
    as_tibble() %>% 
    dplyr::select(geneChr, geneStart, geneEnd, SYMBOL, geneLength, geneStrand) %>% 
    unique()
  
  
  ## De novo Targeted by all peaks
  genes_target_denovo_all <- bind_rows(genes_target_denovo_prom %>% dplyr::mutate(targeted = "Promoter"),
                                       genes_target_denovo_dist %>% dplyr::mutate(targeted = "Distal")) %>%
    unique()
  
  
  ## List of de novo targeted in X vs Y.
  target_denovo <- list("genes_target_denovo_all"  = genes_target_denovo_all,
                        "genes_target_denovo_prom" = genes_target_denovo_prom,
                        "genes_target_denovo_dist" = genes_target_denovo_dist)
  
  
  # Genes targeted by in both X and Y
  ## Targeted in the promoter
  genes_target_common_prom <- x.promo %>% 
    dplyr::filter((SYMBOL %in% y.promo$SYMBOL)) %>%
    dplyr::select(geneChr, geneStart, geneEnd, SYMBOL, geneLength, geneStrand) %>% 
    unique()
  
  ## Targeted by a distal peak
  genes_target_common_dist <- x.distal %>% as_granges() %>%
    filter_by_overlaps(y.distal %>% as_granges) %>%
    as_tibble() %>% 
    dplyr::select(geneChr, geneStart, geneEnd, SYMBOL, geneLength, geneStrand) %>% 
    unique() 
  
  ## Targeted by promoter + distal
  genes_target_common_all <- bind_rows(genes_target_common_prom %>% dplyr::mutate(targeted = "Promoter"),
                                       genes_target_common_dist %>% dplyr::mutate(targeted = "Distal")) %>%
    unique()
  
  ## Lists of common targeted in X and Y
  target_common <- list("genes_target_common_all"  = genes_target_common_all,
                        "genes_target_common_prom" = genes_target_common_prom,
                        "genes_target_common_dist" = genes_target_common_dist)
  
  
  all_targeted <- c(target_denovo, target_common)
  
  # Return
  return(all_targeted)
}

