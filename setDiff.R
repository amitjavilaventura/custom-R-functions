## MAKE VENN DIAGRAMS
venn_2genes <- function(x, color = c("Darkred", "Darkblue"), fill = color, ...){

  x <- Combinations(x)

  draw.pairwise.venn(area1 = length(x[[1]])+length(x[[3]]), area2 = length(x[[2]])+length(x[[3]]), cross.area = length(x[[3]]),
                     category = c(names(x[1]), names(x[2])), col = color, fill = fill, ...)

}


### THE FUNCTIONS BELOW ARE RETRIEVED FROM DFERNANDEZPEREZ
Combinations <- function(x){
  combs <-
    unlist(lapply(1:length(x),
                  function(j) combn(names(x), j, simplify = FALSE)),
           recursive = FALSE)
  names(combs) <- sapply(combs, function(i) paste0(i, collapse = "_"))
  str(combs)
  elements <-
    lapply(combs, function(i) Setdiff(x[i], x[setdiff(names(x), i)]))
  return(elements)
}


Intersect <- function (x) {
  # Multiple set version of intersect
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    intersect(x[[1]], x[[2]])
  } else if (length(x) > 2){
    intersect(x[[1]], Intersect(x[-1]))
  }
}

Union <- function (x) {
  # Multiple set version of union
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    union(x[[1]], x[[2]])
  } else if (length(x) > 2) {
    union(x[[1]], Union(x[-1]))
  }
}

Setdiff <- function (x, y) {
  # Remove the union of the y's from the common x's.
  # x and y are lists of characters.
  xx <- Intersect(x)
  yy <- Union(y)
  setdiff(xx, yy)
}



