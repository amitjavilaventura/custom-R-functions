# GGPLOT HELPERS
# =================================================================================================

## theme_custom() ---------------------------------------------------------------------------------
## Custom theme for ggplot2-based plots, based on ggpubr::theme_pubr()
theme_ggvenn <- function(title.hjust = .5, title.face = "plain", subtitle.face = "italic",
                         title.family = "sans", subtitle.family = "sans",
                         caption.face = "plain", caption.hjust = 1, caption.family = "sans"){

    theme(plot.title = element_text(face = title.face, hjust = title.hjust, family = title.family),
          plot.subtitle = element_text(face = subtitle.face, hjust = title.hjust, family = subtitle.family),
          plot.caption = element_text(face = caption.face, hjust = caption.hjust, family = caption.family))
}


## theme_custom() ---------------------------------------------------------------------------------
## Custom theme for ggplot2-based plots, based on ggpubr::theme_pubr()
theme_custom <- function(legend = "none", margin = T, base_size = 12, border = T,
                        x.text.angle = 0, x.text.hjust = .5, x.text.vjust = .5,
                        title.hjust = .5, title.face = "bold", subtitle.face = "italic",
                        title.family = "sans", subtitle.family = "sans",
                        title.size = 12, subtitle.size = 11, axis.title.size = 10, caption.size = 7,
                        caption.face = "plain", caption.hjust = 1, caption.family = "sans",
                        axis.title.face = "bold", axis.title.familiy = "sans",
                        axis.text.family = "sans", axis.text.size = 10,
                        hgrid.major = 0, hgrid.minor = 0,
                        vgrid.major = 0, vgrid.minor = 0,
                        grid.color = "gray50", grid.type = "dashed"){

  theme_pubr(legend = legend, border = border, margin = margin, base_size = base_size) +
    theme(plot.title = element_text(face = title.face, hjust = title.hjust, family = title.family, size = title.size),
          plot.subtitle = element_text(face = subtitle.face, hjust = title.hjust, family = subtitle.family, size = subtitle.size),
          axis.title = element_text(face = axis.title.face, family = axis.title.familiy, size = axis.title.size),
          legend.title = element_blank(),
          axis.text.x = element_text(angle = x.text.angle, hjust = x.text.hjust, vjust = x.text.vjust),
          axis.text   = element_text(family = axis.text.family, size = axis.text.size),
          plot.caption = element_text(face = caption.face, hjust = caption.hjust, family = caption.family, size = caption.size),
          panel.grid.major.y = element_line(colour = grid.color, size = hgrid.major, linetype = grid.type),
          panel.grid.minor.y = element_line(colour = grid.color, size = hgrid.minor, linetype = grid.type),
          panel.grid.major.x = element_line(colour = grid.color, size = vgrid.major, linetype = grid.type),
          panel.grid.minor.x = element_line(colour = grid.color, size = vgrid.minor, linetype = grid.type))
}


## stat_sum_boxplot() -----------------------------------------------------------------------------
## Function to plot boxplots wihtout oultliers by calling stat_summary(fun.data = calc_boxplot_stat, geom="boxplot", size = 0.8, width = .8)
stat_summary_boxplot <- function( size = .5, width = .5){

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

  stat_summary(fun.data = calc_boxplot_stat, geom="boxplot", size = size, width = width,
               position = position_dodge(width = width), stat = "identity")
}

stat_sum_boxplot <- function( size = .5, width = .5) { stat_summary_boxplot(size = width, width = width) }


## stat_sum_info_boxplot() -----------------------------------------------------------------------------------
## Functiion to write the number of observations in the boxplots
stat_info_boxplot <- function(text_size = 2, text_color = "black", width = .5, y = "lower", statistic = "n", label = ""){

  if(is.character(y)){ y <- tolower(y) }
  if(is.character(statistic)){ statistic <- tolower(statistic) }
  else { stop("'statistic' must be one of 'n', 'mean', 'median', 'sd' ") }

  ## calc_mean_boxplot() ---------------------------------------------------------------------------------------
  calc_info_boxplot <- function(x){

    if(y == "upper"){ y = mean(fivenum(x)[3:4]) } ## upper part of the box
    else if(y == "lower"){ y = mean(fivenum(x)[2:3]) } ## lower part of the box
    else if(y == "min"){ y = fivenum(x)[1] }  ## min value of the boxplot
    else if(y == "max"){ y = fivenum(x)[5] }  ## max value of the boxplot
    else if(y == "q1"){ y = fivenum(x)[2] }  ## q1 value of the boxplot (value 25%)
    else if(y %in% c("median", "q2")){ y = fivenum(x)[3] }  ## median value of the boxplot (value 50%)
    else if(y == "q3"){ y = fivenum(x)[4] }  ## q3 value of the boxplot (value 75%)
    else if(is.numeric(y)) { y <- y }

    x <- x %>% na.omit()

    if(statistic == "n"){ label = paste0(label, length(x), sep = "") }
    else if(statistic == "mean"){ label = paste0(label, round(mean(x), digits = 2), sep = "")}
    else if(statistic == "median"){ label = paste0(label, round(median(x), digits = 2), sep = "")}
    else if(statistic == "sd"){ label = paste0(label, round(sd(x), digits = 2), sep = "")}

    n <- data.frame(y = y, label = label)

    return(n)

  }

  stat_summary(fun.data = calc_info_boxplot, geom = "text", size = text_size, color = text_color,
               position = position_dodge(width = width))
}


## stat_sum_n() -----------------------------------------------------------------------------------
## Functiion to write the number of observations in the boxplots
stat_sum_n_boxplot <- function(text_size = 2, text_color = "black", width = .5, y = "upper"){

  ## calc_n() ---------------------------------------------------------------------------------------
  calc_n_boxplot <- function(x = x){

    if(is.character(y)){ y <- tolower(y) }

    if(y == "upper"){ y = mean(fivenum(x)[3:4]) } ## upper part of the box
    else if(y == "lower"){ y = mean(fivenum(x)[2:3]) } ## lower part of the box
    else if(y == "min"){ y = fivenum(x)[1] }  ## min value of the boxplot
    else if(y == "max"){ y = fivenum(x)[5] }  ## max value of the boxplot
    else if(y == "q1"){ y = fivenum(x)[2] }  ## q1 value of the boxplot
    else if(y %in% c("median", "q3")){ y = fivenum(x)[3] }  ## median value of the boxplot
    else if(y == "q3"){ y = fivenum(x)[4] }  ## q3 value of the boxplot
    else if(is.numeric(y)) { y <- y }

    x <- x %>% na.omit()
    label = paste0("N =", length(x))

    n <- data.frame(y = y, label = label)

    return(n)

  }

  stat_summary(fun.data = calc_n_boxplot, geom = "text", size = text_size, color = text_color,
               position = position_dodge(width = width))
}



## stat_sum_mean_boxplot() -----------------------------------------------------------------------------------
## Functiion to write the number of observations in the boxplots
stat_sum_mean_boxplot <- function(text_size = 2, text_color = "black", width = .5, y = "lower"){

  ## calc_mean_boxplot() ---------------------------------------------------------------------------------------
  calc_mean_boxplot <- function(x){

    if(is.character(y)){ y <- tolower(y) }

    if(y == "upper"){ y = mean(fivenum(x)[3:4]) } ## upper part of the box
    else if(y == "lower"){ y = mean(fivenum(x)[2:3]) } ## lower part of the box
    else if(y == "min"){ y = fivenum(x)[1] }  ## min value of the boxplot
    else if(y == "max"){ y = fivenum(x)[5] }  ## max value of the boxplot
    else if(y == "q1"){ y = fivenum(x)[2] }  ## q1 value of the boxplot
    else if(y %in% c("median", "q3")){ y = fivenum(x)[3] }  ## median value of the boxplot
    else if(y == "q3"){ y = fivenum(x)[4] }  ## q3 value of the boxplot
    else if(is.numeric(y)) { y <- y }

    x <- x %>% na.omit()
    label = round(mean(x), digits = 2)
    label = paste0("Mean =", label)
    n <- data.frame(y = y, label = label)

    return(n)

  }

  stat_summary(fun.data = calc_mean_boxplot, geom = "text", size = text_size, color = text_color,
               position = position_dodge(width = width))
}



## fill facet strips top
## function to fill facet strips (only top ones),
## the output is not a ggplot object
## Credits: I found the original code in issues section of some github repo, but I don't remember
fill_facet_strips_top <- function(p, colors){
  g <- ggplot_gtable(ggplot_build(p))
  strip <- which(grepl('strip-t', g$layout$name))
  fills <- colors
  k <- 1
  for (i in strip) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  plot <- cowplot::ggdraw(g)
  plot
}

## as grob
ggplot_to_grob <- function(p) {
  grob <- cowplot::as_grob(p)
  cowplot::ggdraw(p)
}
