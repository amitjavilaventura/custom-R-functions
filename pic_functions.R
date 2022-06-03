
# Functions that are useful for piRNA and piRNA cluster analysis #
# ============================================================== #


# Get bidirectional clusters ----------
#  It takes a data frame with the coordinates of piRNA clusters in BED format (min 6 cols).
#  It takes those clusters that are divergently transcribed from bidirectional promoters.
#  It does it by looking at the distance separating the cluster pairs (default distance is 500).
#  Note that this function will consider 2 clusters are bidirectional if they are in different strands, but it does not matter the order (i.e., plus-minus or minus-plus).
#  To set to bidirectional only those clusters that are divergently transcribed (minus-plus), set only_divergent = T (default)
get_bidir_pics <- function(pic_anno,
                           bidir_dist     = 500,
                           only_divergent = T,
                           return_pairs   = T,
                           return_dist    = T){

  # Load required packages
  require(dplyr)
  require(magrittr)

  # Format annotation
  pic_anno <- pic_anno %>%
    # Select the first 6 fields (chr, start, end, id, length, strand)
    dplyr::select(1:6) %>%
    # Set column names
    magrittr::set_colnames(c("seqnames", "start", "end", "id", "length", "strand")) %>%
    # Sort by chromosome and starting position
    dplyr::arrange(seqnames, start)

  # Get the strands of the following and previous clusters
  pic_strands <- pic_anno %>%
    # Group by chromosome
    dplyr::group_by(seqnames) %>%
    # Get the strands of the following and previous clusters
    dplyr::mutate(next_strand = lead(strand), prev_strand = lag(strand)) %>%
    # Ungroup
    dplyr::ungroup()

  # Compute the distance between clusters
  pic_dists <- pic_strands %>%
    # Group by chromosome
    dplyr::group_by(seqnames) %>%
    # Compute the start of the next cluster and the end of the previous
    dplyr::mutate(next_start = lead(start), prev_end = lag(end)) %>%
    # Compute the distance with the next cluster
    dplyr::mutate(dist_next = next_start - end) %>%
    # Compute distance with the previous cluster
    dplyr::mutate(dist_prev = start - prev_end) %>%
    # Ungroup
    dplyr::ungroup() %>%
    # Remove next start and previous end fields
    dplyr::select(-next_start, -prev_end)

  # Get directionality of clusters
  pic_direc <- pic_dists %>%
    # Set to NA, those distances greater than 500 and those prev/next strands equal to the strand of the cluster
    dplyr::mutate(dist_next   = ifelse(dist_next > bidir_dist, NA, dist_next),
                  dist_prev   = ifelse(dist_prev > bidir_dist, NA, dist_prev),
                  next_strand = ifelse(next_strand == strand, NA, next_strand),
                  prev_strand = ifelse(prev_strand == strand, NA, prev_strand)) %>%

    # Annotate the directionality of the clusters
    dplyr::mutate(direc = ifelse(!is.na(dist_next) & !is.na(next_strand), "Bidirectional",
                                 ifelse(!is.na(dist_prev) & !is.na(prev_strand), "Bidirectional", "Monodirectional")))

  # Take only those that are divergently transcribed as bidirectional
  if(only_divergent){
    pic_direc <- pic_direc %>%
      dplyr::mutate(direc = ifelse(!is.na(dist_next) & next_strand == "+", "Bidirectional",
                                   ifelse(!is.na(dist_prev) & prev_strand == "-", "Bidirectional", "Monodirectional")))
  }

  # Fix those that have directionality = NA, change them to Monodirectional
  pic_direc <- pic_direc %>% dplyr::mutate(direc = ifelse(is.na(direc), "Monodirectional", direc))

  # Get pairs
  if(return_pairs){
    pic_direc <- pic_direc %>%
      dplyr::mutate(pair = ifelse(!is.na(dist_next) & dist_next == lead(dist_prev), lead(id), "-"),
                    pair = ifelse(!is.na(dist_prev) & dist_prev == lag(dist_next), lag(id), pair))
  }

  # Compute distance
  if(return_dist){
    pic_direc <- pic_direc %>%
      dplyr::mutate(dist2pair = ifelse(!is.na(dist_next) & dist_next == lead(dist_prev), dist_next, NA),
                    dist2pair = ifelse(!is.na(dist_prev) & dist_prev == lag(dist_next), dist_prev, dist2pair))
  }

  # Remove unnecessary fields
  pic_res <- pic_direc %>% dplyr::select(-dist_next, -dist_prev, -next_strand, -prev_strand)

  # Return
  return(pic_res)

}


# Annotate genic pics ----------
#
pic_annot_genes <- function(pic_anno, three_utr, gene_set, stranded = T, reciprocal_overlap = .33){

  # Load required packages ----------
  require(plyr)
  require(dplyr)
  require(magrittr)
  require(bedtoolsr)

  # Read and write piC annot ----------
  # Read pic annot in order to get only the first 6 columns
  # Then write it to temporary file
  pic_anno.tmp <- read.delim(pic_anno, header = F) %>% dplyr::select(1:6) %>% dplyr::arrange(V1, V2)
  pic_anno.tmp %>% write.table(., file = "pic_annot.tmp", quote = F, sep = "\t", col.names = F, row.names = F)
  pic_anno.tmp <- "pic_annot.tmp"

  # Intersect piCs with 3' UTRs of protein coding genes ----------
  # Bedtools intersect with same strandedness and requiring a 33% of reciprocal overlap
  # Add that these are genic clusters that have 3'UTR piRNAs
  # If more than one gene overlaps with the piRNA cluster, take the one with higher overlap
  print(paste("Intersecting piRNA clusters with 3UTRs of protein coding genes. % of reciprocal overlap:", reciprocal_overlap))
  pic_3utr <- bedtoolsr::bt.intersect(a = pic_anno.tmp, b = three_utr, s = stranded, wo = T, f = reciprocal_overlap, r = T) %>%
    dplyr::mutate(genic = "Genic", type = "3UTR") %>%
    dplyr::select(V1, V2, V3, V4, V5, V6, genic, V10, type, V11) %>%
    magrittr::set_colnames(c("seqnames", "start", "end", "id", "length", "strand", "genic", "gene", "type", "overlap")) %>%
    dplyr::group_by(id) %>% dplyr::filter(overlap == max(overlap)) %>% dplyr::ungroup() %>%
    dplyr::select(-overlap) %>%
    unique()

  # Intersect piCs with protein coding genes ----------
  # Bedtools intersect with same strandedness and requiring a 33% of reciprocal overlap
  # Add that these are genic clusters that have genic piRNAs
  # If they are in the list of 3'UTRs, add also that they have 3'UTR piRNAs
  # If more than one gene overlaps with the piRNA cluster, take the one with higher overlap
  print(paste("Intersecting piRNA clusters with protein coding genes. % of reciprocal overlap:", reciprocal_overlap))
  pic_gene <- bedtoolsr::bt.intersect(a = pic_anno.tmp, b = gene_set, s = stranded, wo = T, f = reciprocal_overlap, r = T) %>%
    dplyr::mutate(genic = "Genic", type = "Gene") %>%
    dplyr::select(V1, V2, V3, V4, V5, V6, genic, V10, type, V13) %>%
    magrittr::set_colnames(c("seqnames", "start", "end", "id", "length", "strand", "genic", "gene", "type", "overlap")) %>%
    dplyr::group_by(id) %>% dplyr::filter(overlap == max(overlap)) %>% dplyr::ungroup() %>%
    dplyr::select(-overlap) %>%
    unique()

  # Remove 3'UTR that are already in gene ----------
  # If they are also in the list of genic piCs, remove them from 3'UTR piCs
  pic_3utr <- pic_3utr %>% dplyr::filter(!id %in% pic_gene$id)

  # Define intergenic piRNA clusters ---------
  # Get the clusters that are not in the genic or 3'UTR clusters
  # Define they are intergenic (type and genic variables = "Intergenic)
  print(paste("Retrieving the closest protein-coding genes up and downstream of intergenic piRNA clusters. % of reciprocal overlap:", reciprocal_overlap))
  pic_inter <- read.delim(pic_anno.tmp, header = F) %>%
    magrittr::set_colnames(c("seqnames", "start", "end", "id", "length", "strand")) %>%
    dplyr::filter(!id %in% c(pic_3utr$id, pic_gene$id)) %>%
    dplyr::mutate(type = "Intergenic", genic = type)

  # Get the closest genes downstream and upstream to each cluster
  #  Bedtools closest with clusters and genes, ignore overlapping and force upstream/downstream, use D = "a"
  up_closest   <- bedtoolsr::bt.closest(a = pic_anno.tmp, b = gene_set, io = T, fu = T, D = "a", f = reciprocal_overlap, r = T, ) %>% dplyr::filter(!duplicated(V4)) %>% dplyr::select(V4, V10)
  down_closest <- bedtoolsr::bt.closest(a = pic_anno.tmp, b = gene_set, io = T, fd = T, D = "a", f = reciprocal_overlap, r = T) %>% dplyr::filter(!duplicated(V4)) %>% dplyr::select(V4, V10)

  # Join up and down closest genes
  #  Paste each pair in alphabetical order
  closest <- dplyr::left_join(up_closest, down_closest, by = "V4") %>%
    dplyr::mutate(gene = ifelse(V10.x < V10.y, paste(V10.x, V10.y, sep = "::"), paste(V10.y, V10.x, sep = "::"))) %>%
    dplyr::select("id" = V4, gene)

  # Join closest genes to the intergenic clusters
  pic_inter <- pic_inter %>% dplyr::inner_join(closest, by = "id")

  # Concatenate all piRNA clusters ----------
  print(paste("Joining all genic and intergenic clusters."))
  pic_all <- dplyr::bind_rows(pic_3utr, pic_gene, pic_inter) %>% dplyr::arrange(id)

  # Number of original clusters vs number of annotated clusters.
  noriginal  = pic_anno.tmp %>% read.delim(header = F) %>% nrow()
  nannotated = pic_all %>% nrow()

  if(noriginal != nannotated){
    warning(paste("The number of original clusters(", noriginal, ") is not the same as the number of annotated clusters (", nannotated, ").", sep = ""))
  } else {
    print(paste("Number of original clusters:", noriginal))
    print(paste("Number of annotated clusters:", nannotated))
  }

  # Remove temporary file for piC annotation ----------
  suppressMessages(file.remove(pic_anno.tmp))

  # Return annotated piCs ----------
  return(pic_all)

}


