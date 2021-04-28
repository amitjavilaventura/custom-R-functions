# Rfunctions

R custom functions that are useful for different omics analyses.

## TO DOs

Functions to consider:

* `tidyr::`,  `replace_na()`, `drop_na()`, `pivot_longer()`, `pivot_wider()`, `ùnite()`, `join_all()`, `complete()`, `full_seq()`, `unnest_wider()`, `unnest_longer()`


## `Rfunctions.R`

Script that loads all the functions by doing source of the other R scripts present in this repository.

It use the R package `here` (if you don't have it).

Paths may have to be changed depending on the current directory where `here` works.

## `chromReads.R`

### `chromReads()`

Takes the `bam` files and counts the aligned reads in each chromosome and then plots a bar graph.
Note that this function won't count nor plot the reads aligned to strange (i.e. *unknown* or *random*) chromosomes.

Packages used:

* `Rsamtools`
* `ggplot2`
* `dplyr`

## `pcaplot2.R`

### `pcaplot2()`

It is a modification of the function `pcaplot` from R package `pcaExplorer` that allows to change points shapes.

It takes rlog-transformed DEseq2 dataframe as input. 

Packages used:

* `ggplot2`


## `pieAnno2.R`

It includes two functions:
	
### `filterAnno2()`

Takes the output of `annotatePeak` function from `ChIPseeker` package and changes the annotation features to "Promoter" and "Distal". 

### `pieAnno2()` 

Takes the output of `annotatePeak`, calls `filterAnno2` and plots a ggplot2-based pie chart with only distal and promoter features.

Packages used:

* `ggplot2`

## `pieAnno3.R`

It includes two functions:
	
### `filterAnno3()`

Takes the output of �`annotatePeak and changes the annotation features to "Promoter", "Gene body" or "Distal". 

### `pieAnno3()`

Takes the output of `annotaePeak`, calls `filterAnno3` and plots a ggplot2-based pie chart with only distal, promoter and gene body regions. 

Packages used:

* `ggplot2`

## `volcanoplot2.R`

### `volcanoPlot2()`

Modification of a function made by dfernandezperez in his RNAseq_snakemake pipeline.

Needs the output of DESeq2 analysis after calling de DEGs in the corresponding contrast. Must have an extra "DEG" column with the genes "Downregulated", "Upregulated" or "NotDE". 

`ggplot2`-based function that draws a volcanoplot with the output of a DESeq analysis. 

Takes the points that are outside of the graph limits and draw them as a triangle at the borders.

Optionally: Draws the labels of the DEGs with lower adjusted p-value and higher using `ggrepl`.

Optionally: Instead of using `geom_point()`, uses `geom_point_rast()` from `ggrastr`.

Packages used:

* `ggplot2`
* `dplyr`
* `ggrastr` (optional)
* `ggrepel` (optional)

## `barAnno.R`

### `barAnno()`

Function that plots a ggplot2-based bargraph with the distribution of the peaks from a list of annotation objects that come from as output of `annotatePeak()` from `ChIPseeker`.

As input needs a list of annotation objects from `annotatePeak()`.

It can plot the annoations as "Promoter", "Gene body" (those that are "Intron", "Exon"...) and "Distal" (those that are "Distal intergenic" and "Downstream), 
or it can plot the annotations as "Promoter" and "Distal" (those that are not "Promoter").

Packages used:

* `dplyr`
* `purrr`
* `ggplot2`
* `ggpubr`
* `magrittr`

## `chromHMM_functions.R`

Functions to visualize `ChromHMM` states data. Already in the [`chromHMMviewR` package](https://github.com/amitjavilaventura/chromHMMviewR).

### `chromHMM_emission2heatmap()`

Function that draws a ggplot2-based heatmap with the feature/mark likelihood in each state.

* **Input**: `emissions_{states}.txt` (where '{states}' is the number of states). Can be read directly from the function or an already read dataframe. 
* **_Alternative_ input**: if stored in a variable (i.e. `a`), the output of the function can be reused by looking at the field `$data` from the variable (i.e. `a$data`).
* Dependencies:

  + `dplyr`
  + `reshape2`
  + `ggplot2`
  + `ggpubr`
  + `magrittr`
 
### `chromHMM_enrich2heatmap()`

Function that draws a ggplot2-based heatmap with the state enrichment in genomic regions.

* **Input**: `{condition}_{states}_overlap.txt` (where '{states}' is the number of states and {condition} is the condition of determined in the input matrix for `chromHMM`). Can be read directly from the function or an already read dataframe. 
* **_Alternative_ input**: if stored in a variable (i.e. `a`), the output of the function can be reused by looking at the field `$data` from the variable (i.e. `a$data`).
* Dependencies:

  + `dplyr`
  + `reshape2`
  + `ggplot2`
  + `ggpubr`
  + `magrittr`

### `chromHMM_neighbor2heatmap()`

Function that draws a ggplot2-based heatmap with the state enrichment in genomic regions.

* **Input**: `{condition}_{states}_RefSeq{site}_neighborhood.txt` (where '{states}' is the number of states, {condition} is the condition of determined in the input matrix for `chromHMM` and {site} is TSS/TES). Can be read directly from the function or an already read dataframe. 
* **_Alternative_ input**: if stored in a variable (i.e. `a`), the output of the function can be reused by looking at the field `$data` from the variable (i.e. `a$data`).
* Dependencies:

  + `tidyverse`
  + `reshape2`
  + `ggpubr`
  + `magrittr`
