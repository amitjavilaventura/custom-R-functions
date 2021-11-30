## Deseq helpers

deseq_get_results <- function(dds, ref, contrast, pval = 0.05, lfc = 1, lfcShrink = T) {

  require(dplyr)
  require(tibble)
  require(DESeq2)
  require(apeglm)
  require(IHW)

  #DESeqDataSetFromMatrix and DESeq have been already run before this function

  contrast <- contrast %>% gsub("^condition_", "", .)
  contrast_name = paste("condition_", contrast, sep = "")

  # Relevel
  dds$condition <- relevel(dds$condition, ref)

  # Do nbinomialWaldTest to test differences
  dds <- DESeq2::nbinomWaldTest(dds)

  # Make sure that the contrast name is in the results names
  if(!contrast_name %in% resultsNames(dds)){ stop("The contrast name or the reference are wrong") }

  # Get the results
  res <- DESeq2::results(object = dds, filterFun = ihw, name = contrast_name, alpha = 0.05)

  # LFC shrink
  if(lfcShrink) { res <- DESeq2::lfcShrink(dds = dds, coef = contrast_name, res = res, type = "apeglm") }

  # Add degs column
  res <- res %>% as.data.frame() %>%
    tibble::rownames_to_column("Geneid") %>%
    dplyr::mutate(padj = if_else(is.na(padj), 1, padj),
                  DEG  = "NS",
                  DEG  = if_else(padj < pval & log2FoldChange > lfc, "Upregulated", DEG),
                  DEG  = if_else(padj < pval & log2FoldChange < -lfc, "Downregulated", DEG)) %>%
    dplyr::arrange(padj)

  # Return results
  res

}
