# GGPLOT HELPERS
# =================================================================================================


## theme_custom() ---------------------------------------------------------------------------------
## Custom theme for ggplot2-based plots, based on ggpubr::theme_pubr()
theme_custom <- function(legend = "none", x.text.angle = 0, margin = T, base_size = 12, border = T){

  theme_pubr(legend = legend, border = border, margin = margin, base_size = base_size) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(face = "italic"),
          axis.title = element_text(face = "bold"),
          legend.title = element_blank(),
          axis.text.x = element_text(angle = x.text.angle, hjust = .5, vjust = .5))
}

## calc_boxplot_stat() ----------------------------------------------------------------------------
## Function to use with ggplot2::stat_summary() to build boxplots
##  stat_summary(fun.data = calc_boxplot_stat, geom="boxplot", size = 0.8, width = .8)
calc_boxplot_stat <- function(x){

  # coef for the outliers
  coef <- 1.5

  # total number
  n <- sum(!is.na(x))

  # calculate quantiles
  stats <- quantile(x, probs = c(0.0, 0.25, 0.5, 0.75, 1.0))
  names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")

  # interquartile range
  iqr <- diff(stats[c(2, 4)])

  # set whiskers
  outliers <- x < (stats[2] - coef * iqr) | x > (stats[4] + coef * iqr)

  if (any(outliers)) {
    stats[c(1, 5)] <- range(c(stats[2:4], x[!outliers]), na.rm = TRUE)
  }

  return(stats)

}



## stat_sum_boxplot() -----------------------------------------------------------------------------
## Function to plot boxplots wihtout oultliers by calling stat_summary(fun.data = calc_boxplot_stat, geom="boxplot", size = 0.8, width = .8)
stat_sum_boxplot <- function( size = .5, width = .5){
  stat_summary(fun.data = calc_boxplot_stat, geom="boxplot", size = size, width = width,
               position = position_dodge(width = width), stat = "identity")
}

## calc_n() ---------------------------------------------------------------------------------------
calc_n_boxplot <- function(x){

  x <- x %>% na.omit()
  n <- c(y = mean(fivenum(x)[3:4]), label = length(x))

  return(n)

}


## stat_sum_n() -----------------------------------------------------------------------------------
## Functiion to write the number of observations in the boxplots
stat_sum_n_boxplot <- function(text_size = 2, text_color = "black", width = .5){
  stat_summary(fun.data = calc_n_boxplot, geom = "text", size = text_size, color = text_color,
               position = position_dodge(width = width))
}



## calc_n() ---------------------------------------------------------------------------------------
calc_mean_boxplot <- function(x){

  x <- x %>% na.omit()
  n <- c(y = mean(fivenum(x)[2:3]), label = round(mean(x), digits = 2))

  return(n)

}


## stat_sum_n() -----------------------------------------------------------------------------------
## Functiion to write the number of observations in the boxplots
stat_sum_mean_boxplot <- function(text_size = 2, text_color = "black", width = .5){
  stat_summary(fun.data = calc_mean_boxplot, geom = "text", size = text_size, color = text_color,
               position = position_dodge(width = width))
}
