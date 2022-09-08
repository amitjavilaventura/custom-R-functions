# bed helper functions
bt.merge.id_as_chr <- function(input, stranded = T, distance = 0, ops = "distinct,sum,distinct", cols = "4,5,6"){

  # load packages
  require(dplyr)
  require(bedtoolsr)

  # read input file
  # set ID (V4) as chr (V1)
  # and arrange by V4 and V2 (start)
  i <- read.delim(input, header = F) %>%
    dplyr::select(V4,V2,V3,V1, everything()) %>%
    dplyr::arrange(V4, V2)

  # write formatted input to temp file
  i %>% write.table(., file = "i.tmp", quote = F, sep = "\t", row.names = F, col.names = F)

  # run bedtoolsr::bt.merge on the created tmp file
  # change V1 with V4 (id and chr)
  # arrange by chr and start
  # fix length (end-start)
  merge <- bedtoolsr::bt.merge(i = "i.tmp", s = stranded, d = distance, c = cols, o = ops) %>%
    dplyr::select(V4,V2,V3,V1, everything()) %>%
    dplyr::arrange(V4, V2) %>%
    dplyr::mutate(V5 = V3-V2)

  # return merged data
  merge
}

# bedtoolsr::bt.intersect, but being able to use 6-column dataframes as input instead of paths to files
bt.intersect2 <- function(a, b, stranded = NULL, by_id = F, u = NULL, v = NULL, c = NULL, f = NULL, r = NULL, wa = NULL, wb = NULL, loj = NULL, wo = NULL, wao = NULL){

  # load packages
  require(dplyr)
  require(magrittr)
  require(bedtoolsr)

  # set temp files
  a.tmp <- "a.tmp"
  b.tmp <- "b.tmp"

  # column names
  bed6_colnames <- c("seqnames","start", "end","id", "length","strand")

  # write temporary files
  if(by_id) {
    a %>% dplyr::select(1:6) %>% magrittr::set_colnames(value = bed6_colnames) %>% dplyr::select(4,2,3,1,5,6) %>% dplyr::arrange(id,start) %>% na.omit() %>%
      write.table(., a.tmp, quote = F, sep = "\t", col.names = F, row.names = F)
    b %>%  dplyr::select(1:6) %>% magrittr::set_colnames(value = bed6_colnames) %>% dplyr::select(4,2,3,1,5,6) %>% dplyr::arrange(id,start) %>% na.omit() %>%
      write.table(., b.tmp, quote = F, sep = "\t", col.names = F, row.names = F)
  } else {
    a %>% dplyr::select(1:6) %>% magrittr::set_colnames(value = bed6_colnames) %>% dplyr::arrange(seqnames,start) %>% na.omit() %>%
      write.table(., a.tmp, quote = F, sep = "\t", col.names = F, row.names = F)
    b %>%  dplyr::select(1:6) %>% magrittr::set_colnames(value = bed6_colnames) %>% dplyr::arrange(seqnames,start) %>% na.omit() %>%
      write.table(., b.tmp, quote = F, sep = "\t", col.names = F, row.names = F)
  }

  # intersect a and b
  intersection <- bedtoolsr::bt.intersect(a = a.tmp, b = b.tmp, s = stranded, u = u, f = f, r = r, v = v, wa = wa, wb = wb, loj = loj,  wo = wo, wao = wao, c = c) %>% dplyr::mutate(V5 = V3-V2)

  # reformat resulting bedfile
  if(by_id) {
    intersection <- intersection %>% dplyr::select(4,2,3,1,5,6) %>% magrittr::set_colnames(bed6_colnames) %>% dplyr::arrange(seqnames,start)
  } else {
    intersection <- intersection %>% magrittr::set_colnames(bed6_colnames) %>% dplyr::arrange(seqnames,start)
  }

  # remove temporary files
  if(file.exists(a.tmp)) { file.remove(a.tmp) }
  if(file.exists(b.tmp)) { file.remove(b.tmp) }

  # return resulting bedfile
  return(intersection)
}

# bedtoolsr::bt.subtract, but being able to use 6-column dataframes as input instead of paths to files
bt.subtract2 <- function(a, b, stranded = NULL, by_id = T){

  # load packages
  require(dplyr)
  require(magrittr)
  require(bedtoolsr)

  # set temp files
  a.tmp <- "a.tmp"
  b.tmp <- "b.tmp"

  # column names
  bed6_colnames <- c("seqnames","start", "end","id", "length","strand")

  # write temporary files
  if(by_id) {
    a %>% dplyr::select(1:6) %>% magrittr::set_colnames(value = bed6_colnames) %>% dplyr::select(4,2,3,1,5,6) %>% dplyr::arrange(id,start) %>% na.omit() %>%
      write.table(., a.tmp, quote = F, sep = "\t", col.names = F, row.names = F)
    b %>%  dplyr::select(1:6) %>% magrittr::set_colnames(value = bed6_colnames) %>% dplyr::select(4,2,3,1,5,6) %>% dplyr::arrange(id,start) %>% na.omit() %>%
      write.table(., b.tmp, quote = F, sep = "\t", col.names = F, row.names = F)
  } else {
    a %>% dplyr::select(1:6) %>% magrittr::set_colnames(value = bed6_colnames) %>% dplyr::arrange(seqnames,start) %>% na.omit() %>%
      write.table(., a.tmp, quote = F, sep = "\t", col.names = F, row.names = F)
    b %>%  dplyr::select(1:6) %>% magrittr::set_colnames(value = bed6_colnames) %>% dplyr::arrange(seqnames,start) %>% na.omit() %>%
      write.table(., b.tmp, quote = F, sep = "\t", col.names = F, row.names = F)
  }

  # subtract b from a
  subtracted <- bedtoolsr::bt.subtract(a = a.tmp, b = b.tmp, s = stranded) %>% dplyr::mutate(V5 = V3-V2)

  # reformat resulting bedfile
  if(by_id) {
    subtracted <- subtracted %>% dplyr::select(4,2,3,1,5,6) %>% magrittr::set_colnames(bed6_colnames) %>% dplyr::arrange(seqnames,start)
  } else {
    subtracted <- subtracted %>% magrittr::set_colnames(bed6_colnames) %>% dplyr::arrange(seqnames,start)
  }

  # remove temporary files
  if(file.exists(a.tmp)) { file.remove(a.tmp) }
  if(file.exists(b.tmp)) { file.remove(b.tmp) }

  # return resulting bedfile
  return(subtracted)
}

# bedtoolsr::bt.merge, but being able to use 6-column dataframes as input instead of paths to files
bt.merge2 <- function(i, stranded = NULL, by_id = T, distance = 0, cols = "4,5,6", ops = "distinct,sum,distinct"){

  # load packages
  require(dplyr)
  require(magrittr)
  require(bedtoolsr)

  # set temp files
  i.tmp <- "i.tmp"

  # column names
  bed6_colnames <- c("seqnames","start", "end","id", "length","strand")

  # write temporary files
  if(by_id) {
    i %>% dplyr::select(1:6) %>% magrittr::set_colnames(value = bed6_colnames) %>% dplyr::select(4,2,3,1,5,6) %>% dplyr::arrange(id,start) %>% na.omit() %>%
      write.table(., i.tmp, quote = F, sep = "\t", col.names = F, row.names = F)
  } else {
    i %>% dplyr::select(1:6) %>% magrittr::set_colnames(value = bed6_colnames) %>% dplyr::arrange(seqnames,start) %>% na.omit() %>%
      write.table(., i.tmp, quote = F, sep = "\t", col.names = F, row.names = F)
  }

  # merge
  merged <- bedtoolsr::bt.merge(i = i.tmp, s = stranded, d = distance, c = cols, o = ops) %>% dplyr::mutate(V5 = V3-V2)

  # reformat resulting bedfile
  if(by_id) {
    merged <- merged %>% dplyr::select(4,2,3,1,5,6) %>% magrittr::set_colnames(bed6_colnames) %>% dplyr::arrange(seqnames,start)
  } else {
    merged <- merged %>% magrittr::set_colnames(bed6_colnames) %>% dplyr::arrange(seqnames,start)
  }

  # remove temporary files
  if(file.exists(i.tmp)) { file.remove(i.tmp) }

  # return resulting bedfile
  return(merged)
}

# Functions with FASTA reference genomes and bedfiles
getfasta_from_bed <- function(fasta_genome, is_gz = T, has_chr = F, bed_regions, name = T, stranded = T){

  require(bedtoolsr)
  require(dplyr)

  tmp_genome = "genome.tmp"
  tmp_region = "regions.tmp"

  if(!has_chr){ bed_regions <- bed_regions %>% dplyr::mutate(across(everything(), ~gsub("chr", "", .))) }

  if(is.data.frame(bed_regions)) { bed_regions[,1:6] %>% write.table(tmp_region, quote = F, sep = "\t", row.names = F, col.names = F) }
  else { tmp_region <- bed_regions  }

  if(is_gz) { system(paste("zcat", fasta_genome, ">", tmp_genome))}
  else { tmp_genome <- fasta_genome }

  tmp_genome
  tmp_region

  bed_fasta <- bedtoolsr::bt.getfasta(fi = tmp_genome, bed = tmp_region, name = name, s = stranded)

  if(tmp_genome == "genome.tmp") { file.remove(tmp_genome) }
  if(tmp_region == "regions.tmp") { file.remove(tmp_region) }

  return(bed_fasta)
}
