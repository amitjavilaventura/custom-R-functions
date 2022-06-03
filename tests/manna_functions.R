### Import MSA ==========
# this has to change because the samples are not always ordered identically
import_msa <- function(msa_file, format = "long", samples, print_info = F){

  # Set column names depending on format
  if(format=="long"){ cnames = factor(c("TEfam", "clu_start", "length",	"div",	"score",	"te_strand:te_start:te_end")) }
  else if(format=="normal"){ cnames = factor(c("TEfam", "length",	"div"))}
  else if(format=="short"){ cnames = factor(c("TEfam"))}
  else{ stop("'format' must be one of 'long', 'normal' or 'short'.")}

  # Set samples
  samples <- factor(samples, levels = samples)

  # Set column names
  cnames = outer(paste(samples, "__", sep = ""), cnames, FUN = "paste0") %>% t()
  dim(cnames) = NULL

  # Checking file size
  finfo <- file.info(msa_file)
  if(finfo$size == 0) {
    warning(paste(msa_file, "the file does not have any available line. Returning empty data.frame", sep = " : "))
    msa <- matrix(rep(NA, length(cnames))) %>% t() %>% as.data.frame() %>% magrittr::set_colnames(cnames)
    return(msa)
  }

  # Import MSA
  suppressMessages(msa <- read.delim(msa_file, header = F))

  # Retrieve info
  msa_info <- msa %>% dplyr::filter(stringr::str_detect(V1, "#"))
  info <- data.frame(score   = paste(c(msa_info[1,]$V1)),
                     samples = paste(c(t(msa_info[2,])), collapse = " "),
                     cluster = paste(c(t(msa_info[3,1:2])), collapse = " "),
                     columns = paste(c(msa_info[4, 1:6]), collapse = " ")) %>% t()

  info       = info %>% as.data.frame(.) %>% tibble::rownames_to_column("info")
  cluster_id = info %>%  dplyr::filter(info == "cluster") %>% dplyr::mutate(cluster = gsub("#ClusterID ", "", V1)) %>% dplyr::pull(cluster)
  score      = info  %>% dplyr::filter(info == "score") %>% dplyr::mutate(score = gsub(".*#Score: ", "", V1)) %>% dplyr::pull(score)
  samples    = info %>% dplyr::filter(info == "samples") %>% dplyr::mutate(samples = gsub("#Samples ", "", V1)) %>% dplyr::pull(samples) %>% strsplit(., " ", F) %>% unlist()
  samples    = samples[which(samples != "" & !is.na(samples) & samples != "NA" & !is.null(samples))]
  vars       = info %>% dplyr::filter(info == "columns") %>% dplyr::mutate(vars = gsub("#", "", V1)) %>% dplyr::pull(vars) %>% strsplit(., " ", F) %>% unlist()
  vars       = vars[which(vars != "")]

  if(print_info){
    print("Information on the multiple alignment:")
    print(paste("Cluster ID:", cluster_id))
    print(paste("Manna Score:", score))
    print(paste("Samples:", paste(samples, collapse = ", ")))
    print(paste("Variables:", paste(vars, collapse = ", ")))
  }

  # if(return_info){ return () }

  # RESET column names. Not all the msa are ordered as specified in the manna command.
  cnames = outer(paste(samples, "__", sep = ""), vars, FUN = "paste0") %>% t()
  dim(cnames) = NULL
  cnames = gsub("-", "", cnames) %>% gsub("'", "", .)

  # Remove information from the data.frame
  msa <- msa %>% dplyr::filter(!stringr::str_detect(V1, "#"))

  if(length(msa) != length(cnames)){
    warning(paste(cluster_id, ": the length of the imported .msa file is not the same as the length of the column names. Returning empty dataframe."))
    msa <- matrix(rep(NA, length(cnames))) %>% t() %>% as.data.frame() %>% magrittr::set_colnames(cnames)
  }

  # Set column names
  colnames(msa) <- cnames

  # Return MSA
  return(msa)
}

# manna_msa_files <- list.files(here::here("output/1general/pirna_clusters_orthologs/yuetal_clusters/manna/"), "mus_strains", full.names = T, recursive = T) %>%
#   purrr::set_names(basename(.) %>% gsub("-mus_strains.msa", "", .))
#
# msa_file=manna_msa_files[1]
#
# msa <- import_msa(msa = manna_msa_files[["pi-Ccrn4l"]], format = "long", samples = c("musculus", "caroli", "pahari"))
# msa <- msa %>% dplyr::mutate(across(everything(), ~gsub("____Mus", "", .)))

### TE dotplots ==========
te_dotplot <- function(msa, color_by = "div", sample_names = c("musculus", "caroli", "pahari")){

  m <- msa %>% tibble() %>%
    dplyr::mutate(across(everything(), ~ifelse(.!="-", ., "")))  %>%
    dplyr::rowwise() %>% dplyr::mutate(id=as.character(paste(across(matches("TEfam")), collapse=",")) %>% gsub("^,*", "", .) %>% gsub(",.*", "", .)) %>% dplyr::ungroup() %>%
    dplyr::group_by(id) %>% dplyr::mutate(id = paste(id, "(", row_number(), ")", sep = "")) %>% dplyr::ungroup() %>%
    dplyr::mutate(id = factor(id, levels = id)) %>%
    dplyr::mutate(across(matches("TEfam"), ~ifelse(.=="", F, T))) %>%
    magrittr::set_colnames(gsub("__TEfam", "", colnames(.))) %>%
    tidyr::pivot_longer(cols = all_of(sample_names), names_to = "samples", values_to = "presence") %>%
    dplyr::rowwise() %>% dplyr::mutate(info=as.character(paste(across(matches(color_by)), collapse=",")) %>% gsub("^,*", "", .) %>% gsub(",.*", "", .)) %>% dplyr::ungroup() %>%
    dplyr::select(id, samples, presence, info) %>%
    dplyr::mutate(info = ifelse(!presence, NA, as.numeric(info))) %>%
    dplyr::filter(presence) %>%
    dplyr::count(id, samples, info) %>%
    dplyr::group_by(id, samples) %>% dplyr::summarise(info = mean(info), n = sum(n)) %>% dplyr::ungroup()

  legend_title <- ifelse(color_by == "div", "Divergence", ifelse(color_by == "len", "Length", "Score"))
  max_scale    <- max(m$info)

  dotplot <- m %>%
    dplyr::mutate(samples = factor(samples, levels = sample_names)) %>%
    ggplot(aes(id, samples)) +
    geom_point(aes(fill = info, color = info), size = 3) +
    ggmitji::theme_custom(legend = "right", hgrid.major = .3, vgrid.major = .3, x.text.angle = 90, x.text.hjust = 1, x.text.vjust = .5, axis.text.size = 9) +
    labs(x = "Transposable element", y = "Samples") +
    scale_color_gradient(aesthetics = c("colour", "fill"), low = "blue", high = "red", limits = c(0, max_scale), breaks = seq(0, max_scale, 10),
                         guide = guide_colorbar(title = legend_title, title.position = "right",
                                                title.theme = element_text(angle = 270, hjust = 0.5, vjust = 0.5, size = 12),
                                                frame.colour = "black", frame.linewidth = 1.5, ticks.colour = NA))

  return(dotplot)

}

summary_te_dotplot <- function(msa, color_by = "div", sample_names = c("musculus", "caroli", "pahari")){

  m <- msa %>% tibble() %>%
    dplyr::mutate(across(everything(), ~ifelse(.!="-", ., "")))  %>%
    dplyr::rowwise() %>% dplyr::mutate(id=as.character(paste(across(matches("TEfam")), collapse=",")) %>% gsub("^,*", "", .) %>% gsub(",.*", "", .)) %>% dplyr::ungroup() %>%
    dplyr::mutate(across(matches("TEfam"), ~ifelse(.=="", F, T))) %>%
    magrittr::set_colnames(gsub("__TEfam", "", colnames(.))) %>%
    tidyr::pivot_longer(cols = all_of(sample_names), names_to = "samples", values_to = "presence") %>%
    dplyr::rowwise() %>% dplyr::mutate(info=as.character(paste(across(matches(color_by)), collapse=",")) %>% gsub("^,*", "", .) %>% gsub(",.*", "", .)) %>% dplyr::ungroup() %>%
    dplyr::select(id, samples, presence, info) %>%
    dplyr::mutate(info = ifelse(!presence, NA, as.numeric(info))) %>%
    dplyr::filter(presence) %>%
    dplyr::count(id, samples, info) %>%
    dplyr::group_by(id, samples) %>% dplyr::summarise(info = mean(info), n = sum(n)) %>% dplyr::ungroup()

  legend_title <- ifelse(color_by == "div", "Mean divergence", ifelse(color_by == "len", "Mean length", "Mean score"))
  max_scale    <- max(m$info)

  dotplot <- m %>%
    dplyr::mutate(samples = factor(samples, levels = sample_names)) %>%
    ggplot(aes(samples, id)) +
    geom_point(aes(fill = info, color = info, size = n)) +
    ggmitji::theme_custom(legend = "right", hgrid.major = .3, vgrid.major = .3, x.text.angle = 0, x.text.hjust = .5, x.text.vjust = .5, axis.text.size = 9) +
    labs(y = "Transposable element", x = "Samples") +
    scale_size(guide = guide_legend(title = "TE insertions", title.position = "right", title.theme = element_text(angle = 270, hjust = 0.5, vjust = 0.5, size = 12))) +
    scale_color_gradient(aesthetics = c("colour", "fill"), low = "blue", high = "red", limits = c(0, max_scale), breaks = seq(0, max_scale, 10),
                         guide = guide_colorbar(title = legend_title, title.position = "right",
                                                title.theme = element_text(angle = 270, hjust = 0.5, vjust = 0.5, size = 12),
                                                frame.colour = "black", frame.linewidth = 1.5, ticks.colour = NA))

  return(dotplot)

}

# te_dotplot(msa)
# summary_te_dotplot(msa)

### Num of annotated TEs in each sample ==========
te_num <- function(msa, cluster_id = NULL){

  if(!is.null(cluster_id)) { print(paste("Cluster:", cluster_id)) }

  te_num <- msa %>% tibble() %>%
    dplyr::select(matches("TEfam")) %>% dplyr::mutate(across(matches("TEfam"), ~ifelse(.!="-", 1, 0))) %>%
    magrittr::set_colnames(gsub("__TEfam", "", colnames(.))) %>%
    colSums()
  return(te_num)
}

# number_of_te <- te_num(msa)

### Common TEs =========
get_common_te <- function(msa, cluster_id = NULL){

  if(!is.null(cluster_id)) { print(paste("Cluster:", cluster_id)) }

  m           <- msa %>% dplyr::mutate(across(everything(), ~ifelse(.=="-", NA, .))) %>% tibble() %>% dplyr::select(matches("TEfam"))
  te_num      <- m %>% dplyr::mutate(across(everything(), ~ifelse(is.na(.), 0, 1))) %>% dplyr::rowwise() %>% dplyr::transmute(TEnum = sum(across(everything()))) %>% dplyr::ungroup() %>% dplyr::mutate(i = row_number())
  samples_num <- m %>% ncol()
  common_idx  <- te_num$i[which(te_num$TEnum == samples_num)]
  common_te   <- m[common_idx,] %>% magrittr::set_colnames(gsub("__.*", "", colnames(.)))

  return(common_te)
}

# common_te <- get_common_te(msa)

### Unique TEs ==========

get_unique_te <- function(msa, cluster_id = NULL){

  if(!is.null(cluster_id)) { print(paste("Cluster:", cluster_id)) }

  m           <- msa %>% dplyr::mutate(across(everything(), ~ifelse(.=="-", NA, .))) %>% tibble() %>% dplyr::select(matches("TEfam"))
  te_num      <- m %>% dplyr::mutate(across(everything(), ~ifelse(is.na(.), 0, 1))) %>% dplyr::rowwise() %>% dplyr::transmute(TEnum = sum(across(everything()))) %>% dplyr::ungroup() %>% dplyr::mutate(i = row_number())
  samples_num <- m %>% ncol()
  unique_idx  <- te_num$i[which(te_num$TEnum == 1)]
  unique_te   <- m[unique_idx,] %>% magrittr::set_colnames(gsub("__.*", "", colnames(.)))

  return(unique_te)
}

# unique_te <- get_unique_te(msa)
### Partial TEs ==========

get_partial_te <- function(msa, cluster_id = NULL){

  if(!is.null(cluster_id)) { print(paste("Cluster:", cluster_id)) }

  m           <- msa %>% dplyr::mutate(across(everything(), ~ifelse(.=="-", NA, .))) %>% tibble() %>% dplyr::select(matches("TEfam"))
  te_num      <- m %>% dplyr::mutate(across(everything(), ~ifelse(is.na(.), 0, 1))) %>% dplyr::rowwise() %>% dplyr::transmute(TEnum = sum(across(everything()))) %>% dplyr::ungroup() %>% dplyr::mutate(i = row_number())
  samples_num <- m %>% ncol()
  some_idx    <- te_num$i[which(te_num$TEnum > 1 & te_num$TEnum < samples_num)]
  some_te     <- m[some_idx,] %>% magrittr::set_colnames(colnames(.) %>% gsub("__TEfam", "", .))

  return(some_te)
}

# unique_te <- get_unique_te(msa)

### INFORMATION ==========

mean_info_te <- function(msa, info = "div", te = "all", cluster_id = NULL){

  if(!is.null(cluster_id)) { print(paste("Cluster:", cluster_id)) }

  if(!info %in% c("div", "len", "score")) { stop("'info' should be one of div, len or score") }
  if(!te %in% c("all", "common", "some", "unique", "partial")) { stop("'te' should be one of all, common, some or unique") }

  if(te == "partial") { te = "some" }

  #print(paste("Getting", info, "of", te, "transposable elements"))

  m <- msa %>% dplyr::mutate(across(everything(), ~ifelse(.=="-", NA, .)))

  te_num <- m %>% tibble() %>% dplyr::select(matches("TEfam")) %>% dplyr::mutate(across(everything(), ~ifelse(is.na(.), 0, 1))) %>%
    dplyr::rowwise() %>% dplyr::transmute(TEnum = sum(across(everything()))) %>% dplyr::ungroup() %>% dplyr::mutate(i = row_number())
  samples_num <-  m %>% tibble() %>% dplyr::select(matches("TEfam")) %>% ncol()

  common_idx <- te_num$i[which(te_num$TEnum == samples_num)]
  some_idx   <- te_num$i[which(te_num$TEnum > 1 & te_num$TEnum < samples_num)]
  unique_idx <- te_num$i[which(te_num$TEnum == 1)]

  all     <- m %>% dplyr::select(matches(info)) %>% magrittr::set_colnames(gsub("__.*", "", colnames(.))) %>% dplyr::mutate(across(everything(), ~as.numeric(.))) %>% colMeans(na.rm = T)
  common  <- m[common_idx,] %>% dplyr::select(matches(info)) %>% magrittr::set_colnames(gsub("__.*", "", colnames(.))) %>% dplyr::mutate(across(everything(), ~as.numeric(.))) %>% colMeans(na.rm = T)
  some    <- m[some_idx,] %>% dplyr::select(matches(info)) %>% magrittr::set_colnames(gsub("__.*", "", colnames(.))) %>% dplyr::mutate(across(everything(), ~as.numeric(.))) %>% colMeans(na.rm = T)
  unique  <- m[unique_idx,] %>% dplyr::select(matches(info)) %>% magrittr::set_colnames(gsub("__.*", "", colnames(.))) %>% dplyr::mutate(across(everything(), ~as.numeric(.))) %>% colMeans(na.rm = T)

  info_list <- list("all" = all, "common" = common, "some" = some, "unique" = unique)

  return(info_list[[te]])
}

# mean_info_te(msa, info = "div", te = "all")
# mean_info_te(msa, info = "div", te = "unique")
# mean_info_te(msa, info = "div", te = "common")
#
# get_info_te(msa, info = "len", te = "all")
# get_info_te(msa, info = "len", te = "unique")
# get_info_te(msa, info = "len", te = "common")
