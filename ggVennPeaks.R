ggVennPeaks <- function(peak_list, peak_names = names(peak_list), percent = T,
                        fill = c("blue", "gold3"), alpha = .4,
                        color = "black", text_color = "black",
                        name_size = 5, label_size = 3){

  require(ChIPpeakAnno)

  x <- makeVenn4upSet(peak_list)
  y <- x$matrix %>%
    as_tibble() %>%
    magrittr::set_colnames(c("Peak", paste("cond", 1:length(peak_list), sep = ""))) %>%
    reshape2::melt() %>% mutate(value = if_else(value == 1, Peak, NULL)) %>%
    tidyr::pivot_wider(names_from = "variable", values_from = "value") %>%
    dplyr::select(-Peak) %>%
    dplyr::as_tibble() %>%
    as.list() %>%
    purrr::set_names(peak_names)

  venn <- ggvenn::ggvenn(data = y, show_percentage = percent,
                         fill_color = fill, fill_alpha = alpha,
                         stroke_color =  color,
                         set_name_color = color, set_name_size = name_size,
                         text_color = text_color, text_size = label_size)


  return(venn)
}
