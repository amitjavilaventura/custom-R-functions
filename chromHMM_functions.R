###############################
### CHROM HMM VISUALIZATION ###
###############################

# Several functions to help in the visualization of chromHMM outputs

# -------------------------------------------------------------------------------------------------
# chromHMM_emission2heatmap
# -------------------------------------------------------------------------------------------------

chromHMM_emission2heatmap <- function(df = NULL, file = NULL, df_melt = NULL, features = NULL, states = NULL,
                                      title = "", subtitle = "", xlab = "Features", color = "Blue"){

  # PACKAGES
  require(dplyr)
  require(magrittr)
  require(ggplot2)
  require(ggpubr)

  # Require 'regions' and 'states' to be a character vector
  if(!is.character(features) | !is.character(states)){ return("Both 'features' and 'states' must be a not-null character vector") }

  # If data is raw output from chromHMM, use either file or df
  if(is.null(df_melt)){
    # Read data either directly from a file or from a already read df. ChromHMM file (*_overlap.txt).
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
  g <- ggplot(df.m, aes(variable, State, fill = value)) + geom_tile(show.legend = F) +
    scale_fill_gradient(low = "White", high = color) +
    geom_text(mapping = aes(variable, State, label = round(value, 2))) +
    ggtitle(title, subtitle) + ylab("State") + xlab(xlab) +
    theme_pubr() + theme(axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1),
                         axis.title = element_text(face = "bold", size = 14))

  return(g)


}



# -------------------------------------------------------------------------------------------------
# chromHMM_enrich2heatmap
# -------------------------------------------------------------------------------------------------

chromHMM_enrich2heatmap <- function(df = NULL, file = NULL, df_melt = NULL, regions = NULL, states = NULL,
                                    title = "", subtitle = "", color = "Blue"){

  # PACKAGES
  require(dplyr)
  require(magrittr)
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
    df.m <- df %>% melt() %>% group_by(variable) %>% mutate(total = sum(value)) %>% dplyr::ungroup() %>%
      mutate(percent = value/total*100, State = factor(State, levels = states))
  }
  # If data has already been processed within this function and then filtered for some reason, use df_melt
  else{ df.m <- df_melt }

  # Draw heatmap with ggplot2
  g <- ggplot(df.m, aes(variable, State, fill = percent)) + geom_tile(show.legend = F) +
    scale_fill_gradient(low = "White", high = color) +
    geom_text(mapping = aes(variable, State, label = round(value, 2))) +
    ggtitle(title, subtitle) + ylab("State") + xlab("Genomic regions") +
    theme_pubr() + theme(axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1),
                         axis.title = element_text(face = "bold", size = 14))

  return(g)
}
