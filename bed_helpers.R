# bed helper functions

bt.merge.id_as_chr <- function(input, stranded = T, distance = 0){

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
  merge <- bedtoolsr::bt.merge(i = "i.tmp", s = stranded, d = distance, c = "4,5,6", o = "distinct,sum,distinct") %>%
    dplyr::select(V4,V2,V3,V1, everything()) %>%
    dplyr::arrange(V4, V2) %>%
    dplyr::mutate(V5 = V3-V2)

  # return merged data
  merge

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
