# R functions

R custom functions that are useful for different omics analyses.

If you use one of these functions, please, mention me and give a star to this repository.

## TO DOs

Functions and packages to consider:

* `dplyr::`: `add_count()`, `mutate_if()`, `summarise_if()`
* `tidyr::`: `replace_na()`, `drop_na()`, `pivot_longer()`, `pivot_wider()`, `ùnite()`, `join_all()`, `complete()`, `full_seq()`, `unnest_wider()`, `unnest_longer()`, `gather()`
* `forcats::`: `fct_reorder()`, `fct_relevel()`, `fct_rev()`



## `Rfunctions.R`

Script that loads all the functions by doing source of the other R scripts present in this repository.

It use the R package `here` (if you don't have it).

Paths may have to be changed depending on the current directory where `here` works.

## `pcaplot2.R`

### `pcaplot2()`

It is a modification of the function `pcaplot` from R package `pcaExplorer` that allows to change points shapes.

It takes rlog-transformed DEseq2 dataframe as input. 

Packages used:

* `ggplot2`


## `pieAnno.R`

It includes two functions:
	
### `filterAnno()`

Takes the output of `annotatePeak` function from `ChIPseeker` package and changes the annotation features to "Promoter" and "Distal", or "Promoter", "Distal" and "Gene body". 

### `pieAnno()` 

Takes the output of `annotatePeak`, calls `filterAnno()` and plots a ggplot2-based pie chart with only distal and promoter features.

Packages used:

* `ggplot2`

## `signals_in_regions.R`

Functions that compute the coverage (BIGWIG) over several regions (BED).

* Input: 
  + character vector with the path(s) to bed files;
  + character vector with the path(s) to bigwig files.
