#########################
###     PCAPLOT 2     ###
#########################



#' Modification of the R function pcaplot() from the R package pcaExplorer in
#' order to plot PCA easier after runing plotPCA(returnData = T) to get the plot data
#' Take just the part where the plot is made, without any data calculations.
#'
#' Added parameters
#'	@param point_shape Integer or vector of integers giving the shapes of the points (ex. 16=circle, 17=triangle, 15=square, 18=rhombus)
#'	@param legend_title Logical. If FALSE, the legend title is not printed. Default: FALSE
#'
#' Sample PCA plot for transformed data
#'
#' Plots the results of PCA on a 2-dimensional space
#'
#' @param x A \code{\link{DESeqTransform}} object, with data in \code{assay(x)},
#' produced for example by either \code{\link{rlog}} or
#' \code{\link{varianceStabilizingTransformation}}
#' @param intgroup Interesting groups: a character vector of
#' names in \code{colData(x)} to use for grouping
#' @param ntop Number of top genes to use for principal components,
#' selected by highest row variance
#' @param returnData logical, if TRUE returns a data.frame for further use, containing the
#' selected principal components and intgroup covariates for custom plotting
#' @param title The plot title
#' @param pcX The principal component to display on the x axis
#' @param pcY The principal component to display on the y axis
#' @param text_labels Logical, whether to display the labels with the sample identifiers
#' @param point_size Integer, the size of the points for the samples
#' @param ellipse Logical, whether to display the confidence ellipse for the selected groups
#' @param ellipse.prob Numeric, a value in the interval [0;1)
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#' @examples
#' dds <- makeExampleDESeqDataSet_multifac(betaSD_condition = 3,betaSD_tissue = 1)
#' rlt <- DESeq2::rlogTransformation(dds)
#' pcaplot(rlt, ntop=200)
#'
#' @export
pcaplot2 <- function(x, title = NULL, pcX = 1, pcY = 2, text_labels = TRUE, legend_title = FALSE,
					           point_size = 3, point_shape = 16, ellipse = TRUE, ellipse.prob = 0.95) # customized principal components
{

  # Load packages
  require(ggplot2)

  g <- ggplot(data = x, aes_string(x = paste0("PC",pcX), y = paste0("PC",pcY), color = "group")) +
    geom_point(size = point_size, shape = point_shape) +
    xlab(paste0("PC",pcX,": ", round(percentVar[pcX] * 100,digits = 2), "% variance")) +
    ylab(paste0("PC",pcY,": ", round(percentVar[pcY] * 100,digits = 2), "% variance"))

  ## plot confidence ellipse
  # credit to vince vu, author of ggbiplot
  if(ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))

    ell <- ddply(d, 'group', function(x) {
      if(nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x[[paste0("PC",pcX)]], x[[paste0("PC",pcY)]]))
      mu <- c(mean(x[[paste0("PC",pcX)]]), mean(x[[paste0("PC",pcY)]]))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'),
                 groups = x$group[1])
    })
    # names(ell)[1:2] <- c('xvar', 'yvar')
    if(nrow(ell)>0) {
      g <- g + geom_path(data = ell, aes_string(x="X1",y="X2",color = "groups", group = "groups"))
    }
  }

  if(text_labels)
    g <- g + geom_label_repel(mapping = aes_string(label="names",fill="group"),
                              color="white", show.legend = TRUE)
  if(!is.null(title)) g <- g + ggtitle(title)
  g <- g + theme_bw()
  # changing the legend title
  if(legend_title == FALSE) {g <- g + theme(legend.title=element_blank())}
  # as in http://www.huber.embl.de/msmb/Chap-Graphics.html
  # "well-made PCA plots usually have a width that’s larger than the height"
  g <- g + coord_fixed()
  g <- g + theme(plot.title = element_text(face = "bold"))
  g
}
