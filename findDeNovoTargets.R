
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
  
  # Targeted promoter --> genes with a peak in the promoter only in X
  x.promo <- x %>% dplyr::filter(anno == "Promoter")
  y.promo <- y %>% dplyr::filter(anno == "Promoter")
  genes_target_denovo_prom <- x.promo %>% dplyr::filter(!(SYMBOL %in% y.promo$SYMBOL))
  
  # Targeted by a distal peak --> closest gene to a distal peak in X that do not overlap with a peak in Y.
  x.distal <- x %>% dplyr::filter(anno == "Distal")
  y.distal <- y %>% dplyr::filter(anno == "Distal")
  
  genes_target_denovo_dist <- x.distal %>% as_granges() %>%
    filter_by_non_overlaps(y.distal %>% as_granges) %>%
    as_tibble()
  
  target_denovo <- list("genes_target_denovo_prom" = genes_target_denovo_prom,
                        "genes_target_denovo_dist" = genes_target_denovo_dist)
    
  return(target_denovo)
}

