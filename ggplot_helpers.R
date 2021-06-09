# GGPLOT HELPERS
# =================================================================================================


## theme_custom() ---------------------------------------------------------------------------------
## Custom theme for ggplot2-based plots, based on ggpubr::theme_pubr()
theme_custom <- function(legend = "none", x.text.angle = 0, margin = T, base_size = 12, border = T){

  theme_pubr(legend = legend, border = border, x.text.angle = x.text.angle, margin = margin, base_size = base_size) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(face = "italic"),
          axis.title = element_text(face = "bold"),
          legend.title = element_blank())
}

## calc_boxplot_stat() ----------------------------------------------------------------------------
## Function to use with ggplot2::stat_summary() to build boxplots
##  stat_summary(fun.data = calc_boxplot_stat, geom="boxplot", size = 0.8, width = .8)
calc_boxplot_stat <- function(x) {

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
  stat_summary(fun.data = calc_boxplot_stat, geom="boxplot", size = size, width = size)
}
