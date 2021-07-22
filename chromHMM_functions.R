###############################
### CHROM HMM VISUALIZATION ###
###############################

# Several functions to help in the visualization of chromHMM outputs

# chromHMM_emission2heatmap  ------
# ----------------------------------------------------------------------------------------------- #

chromHMM_emission2heatmap <- function(df = NULL, file = NULL, df_melt = NULL, features = NULL, states = NULL,
                                      title = "", subtitle = "", xlab = "Features", color = "Tomato",
                                      label_size = 2, reverse_y = T, plotly = F){

  # PACKAGES
  require(dplyr)
  require(magrittr)
  require(reshape2)
  require(ggplot2)
  require(ggpubr)

  # Require 'regions' and 'states' to be a character vector
  if(!is.character(features) | !is.character(states)){ return("Both 'features' and 'states' must be a not-null character vector") }

  # If data is raw output from chromHMM, use either file or df
  if(is.null(df_melt)){
    # Read data either directly from a file or from a already read df. ChromHMM file (*_emission.txt).
    if(is.null(df) & !is.null(file)){ df <- read.delim(file) }
    else if(!is.null(df) & is.null(file)) { df <- df }
    else{ return("One of 'df', 'file' or 'df_melt' must not be null and the others must.") }

    # Process data
    df <- df %>% set_colnames(c("State", features)) %>% mutate(State = states) #%>% column_to_rownames("State")
    df.m <- df %>% melt() %>% mutate(State = factor(State, levels = states))
  }
  # If data has already been processed within this function and then filtered for some reason, use df_melt
  else{ df.m <- df_melt }

  # Draw heatmap with ggplot2
  g <- ggplot(df.m, aes(variable, State, fill = value)) + geom_tile(color = "gray", show.legend = F) +
    scale_fill_gradient(low = "White", high = color) +
    coord_equal() +
    geom_text(mapping = aes(variable, State, label = round(value, 2)), size = label_size) +
    ggtitle(title, subtitle) + ylab("State") + xlab(xlab) +
    theme_pubr() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
                         axis.title = element_text(face = "bold", size = 12),
                         axis.text = element_text(size = 10),
                         plot.title = element_text(hjust=.5), plot.subtitle = element_text(hjust=.5))

  # Reverse Y axis
  if(reverse_y){ g <- g + scale_y_reverse() }


  # Plotlify
  if(plotly){
    require(plotly)
    g <- ggplotly(g)
  }

  # Return plot
  return(g)
}



# chromHMM_enrich2heatmap -------
# ----------------------------------------------------------------------------------------------- #

chromHMM_enrich2heatmap <- function(df = NULL, file = NULL, df_melt = NULL, regions = NULL, states = NULL,
                                    title = "", subtitle = "", color = "Cornflowerblue", scale_color = "scale",
                                    legend = F, label_size = 2, reverse_y = T, plotly = F){

  # PACKAGES
  require(dplyr)
  require(magrittr)
  require(reshape2)
  require(ggplot2)
  require(ggpubr)

  # Require 'regions' and 'states' to be a character vector
  if(!is.character(regions) | !is.character(states)){ return("Both 'regions' and 'states' must be a not-null character vector") }

  # If data is raw output from chromHMM, use either file or df
  if(is.null(df_melt)){
    # Read data either directly from a file or from a already read df. ChromHMM file (*_overlap.txt).
    if(is.null(df) & !is.null(file)){ df <- read.delim(file) }
    else if(!is.null(df) & is.null(file)) { df <- df }
    else{ return("One of 'df', 'file' or 'df_melt' must not be null and the others must.") }

    # Process data
    df <- df %>% set_colnames(c("State", regions)) %>% filter(State != "Base") %>% mutate(State = states) #%>% column_to_rownames("State")
    df.m <- df %>% melt() %>% group_by(variable) %>%
      mutate(total = sum(value), min = min(value), max = max(value)) %>% dplyr::ungroup()  %>%
      mutate(scale = (value-min)/max, percent = value/total*100) %>%
      mutate(State = factor(State, levels = states))
  }
  # If data has already been processed within this function and then filtered for some reason, use df_melt
  else{ df.m <- df_melt }

  # Draw heatmap with ggplot2
  if( scale_color == "scale" ) { g <- ggplot(df.m, aes(variable, State, fill = scale)) + geom_tile(show.legend = legend) }
  else if( scale_color == "percent") { g <- ggplot(df.m, aes(variable, State, fill = percent)) + geom_tile(show.legend = legend) }
  else if( scale_color == "value") { g <- ggplot(df.m, aes(variable, State, fill = value)) + geom_tile(show.legend = legend) }
  else { return("The option scale_color must be one of 'scale', 'percent' or 'value'. Most reccommended is 'scale'.")}

  g <- g + scale_fill_gradient(low = "White", high = color) +
    coord_equal() +
    geom_text(mapping = aes(variable, State, label = round(value, 2)), size = label_size) +
    ggtitle(title, subtitle) + ylab("State") + xlab("Genomic regions") +
    theme_pubr() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
                         axis.title = element_text(face = "bold", size = 12),
                         axis.text = element_text(size = 10))

  # Reverse Y axis
  if(reverse_y){ g <- g + scale_y_reverse() }

  # Plotlify
  if(plotly){
    require(plotly)
    g <- ggplotly(g)
  }

  # Return plot
  return(g)
}


# chromHMM_neighbor2heatmap ------
# ----------------------------------------------------------------------------------------------- #
chromHMM_neighbor2heatmap <- function(df = NULL, file = NULL, df_melt = NULL, states = NULL,
                                      title = "", subtitle = "", xlab = "Distance to TSS",
                                      color = "Cornflowerblue", legend = F, label_size = 2,
                                      reverse_y = T, plotly = F){

  # PACKAGES
  require(tidyverse)
  require(magrittr)
  require(reshape2)
  require(ggpubr)

  # Require 'regions' and 'states' to be a character vector
  if(!is.character(states)){ return("'states' must be a not-null character vector") }

  # If data is raw output from chromHMM, use either file or df
  if(is.null(df_melt)){
    # Read data either directly from a file or from a already read df. ChromHMM file (*_neighborhood.txt).
    if(is.null(df) & !is.null(file)){ df <- read.delim(file) }
    else if(!is.null(df) & is.null(file)) { df <- df }
    else{ return("One of 'df', 'file' or 'df_melt' must not be null and the others must.") }

    # Get position names
    positions <- colnames(df) %>%
      purrr::keep(~str_detect(string = .x, pattern = "State", negate = T)) %>%
      purrr::map_chr(~gsub(x = .x, pattern = "X\\.", replacement = "-")) %>%
      purrr::map(~gsub(x = .x, pattern = "X", replacement = ""))
    df <- df %>% set_colnames(c("State", positions)) %>% mutate(State = states)

    df.m <- df %>% melt() %>% group_by(variable) %>%
      mutate(total = sum(value), min = min(value), max = max(value)) %>% dplyr::ungroup()  %>%
      mutate(scale = (value-min)/max) %>%
      mutate(State = factor(State, levels = states))
  }
  # If data has already been processed within this function and then filtered for some reason, use df_melt to redo the plot
  else{ df.m <- df_melt }

  # Draw heatmap with ggplot2
  g <- ggplot(df.m, aes(variable, State, fill = value)) + geom_tile(show.legend = legend) +
    coord_equal() +
    scale_fill_gradient(low = "White", high = color) +
    geom_text(mapping = aes(variable, State, label = round(value, 2)), size = label_size) +
    ggtitle(title, subtitle) + ylab("State") + xlab(xlab) +
    theme_pubr() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
                         axis.title = element_text(face = "bold", size = 12),
                         axis.text = element_text(size = 10))
  # Reverse Y axis
  if(reverse_y){ g <- g + scale_y_reverse() }

  # Plotlify
  if(plotly){
    require(plotly)
    g <- ggplotly(g)
  }

  # Return plot
  return(g)
}


# chromHMM_states2barplot ------
# ----------------------------------------------------------------------------------------------- #

chromHMM_states2barplot <-  function(df = NULL, file = NULL, df_processed = NULL, states = NULL,
                                     title = "", subtitle = "", color = c("Cornflowerblue", "Darkblue"), scale_color = "",
                                     plotly = F){

  # PACKAGES
  require(dplyr)
  require(magrittr)
  require(ggplot2)
  require(ggpubr)

  # Require 'regions' and 'states' to be a character vector
  if(!is.character(states)){ return("The argument 'states' must be a not-null character vector") }

  # If data is raw output from chromHMM, use either file or df
  if(is.null(df_processed)){
    # Read data either directly from a file or from a already read df. ChromHMM file (*_segments.bed).
    if(is.null(df) & !is.null(file)){ df <- read.delim(file, header = F) }
    else if(!is.null(df) & is.null(file)) { df <- df }
    else{ return("One of 'df', 'file' or 'df_melt' must not be null and the others must.") }

    # Process data
    num_states <- df %>% set_colnames(c("seqnames", "start", "end", "state")) %>% count(state)
    bases_states <- df %>% set_colnames(c("seqnames", "start", "end", "state")) %>%
      mutate(width = sqrt((end-start)^2)) %>% select(state, width) %>%
      group_by(state) %>% summarise(bases = sum(width)) %>% ungroup()

    df2 <- inner_join(num_states, bases_states) %>%
      mutate(`Mean length` = bases/n, `Number of elements` = n) %>%
      mutate(state = factor(states, level = states))
  }
  # If data has already been processed within this function and then filtered for some reason, use df_melt
  else{ df2 <- df_processed }

  # Draw heatmap with ggplot2¡
  if(scale_color == "meanLength"){
    g <- ggplot(df2, aes(state, bases, fill = `Mean length`)) + geom_col() +
      scale_fill_gradient(low = color[1], high = color[2], limits = c(0, max(df2$`Mean length`)))
  }
  else if(scale_color == "numElements"){
    g <- ggplot(df2, aes(state, bases, fill = `Number of elements`)) + geom_col() +
      scale_fill_gradient(low = color[1], high = color[2], limits = c(0, max(df2$`Number of elements`)))
  }
  else{
    g <- ggplot(df2, aes(state, bases, fill = state)) + geom_col(show.legend = F)
  }


  g <- g + coord_flip() + scale_y_continuous(expand = expansion(mult = c(0.01, 0.1))) +
    ggtitle(title, subtitle) + xlab("State") + ylab("Total bases in each state") +
    theme_pubr() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                         axis.title = element_text(face = "bold", size = 14),
                         legend.position = "right")

  # Plotlify
  if(plotly){
    require(plotly)
    g <- ggplotly(g)
  }

  # Return plot
  return(g)
}
