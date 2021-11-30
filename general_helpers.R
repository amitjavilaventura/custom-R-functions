# =============================================================================================== #
# General helper functions                                                                        #
# =============================================================================================== #

# "Not in" operator, contrary of %in%
"%nin%" <- Negate("%in%")

# Rearrange columns
rearrange_cols <- function(df, col_order = colnames(df), new_colnames = NULL){
  df <- df[col_order]
  if(!is.null(new_colnames)){ colnames(df) <- new_colnames }
  return(df)
}

# Filter list of vectors
filter_vector_list <- function(vec_list, filt_val){
  vec_filt_list <- list()
  for(i in 1:length(vec_list)){
    name <- names(vec_list[i])
    vec_filt_list[[name]] <- vec_list[[i]] %>% magrittr::extract(. %in% filt_val)
  }
  return(vec_filt_list)
}


# bind cols with different length
bind_cols2 <- function(list, name = names(list)){
  N <- list %>% purrr::map(~nrow(.x)) %>% unlist() %>% max()
  for(i in 1:length(list)){
    numrow <- nrow(list[[i]])
    numcol <- ncol(list[[i]])
    n = N-numrow
    if(numrow == N){ list[[i]] > list[[i]] }
    else{ list[[i]][(numrow+1):(numrow+n),] <- NA }
  }
  list_bind <- bind_cols(list) %>% magrittr::set_colnames(name)
  return(list_bind)
}
