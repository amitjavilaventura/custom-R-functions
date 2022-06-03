# Plot helpers

# Highlight areas of 2D Venn
venn2_highlight <- function(highlight = c("AB", "AnoB"), color = c("red", "pink"),
                            label = NULL, label.pos = "bottom.right", label.face = "italic", label.size = 13, label.col = "black",
                            line.col = "black"){

  library(polyclip)
  library(VennDiagram)
  library(ggplot2)

  invisible(vp <- VennDiagram::draw.pairwise.venn(area1 = 2, area2 = 2, cross.area = 1, category = rep("", 2), cex = 0, label.col = NA, scaled = F, fill = NA, margin = 0, disable.logging = T))

  A <- list(list(x = as.vector(vp[[1]][[1]]), y = as.vector(vp[[1]][[2]])))
  B <- list(list(x = as.vector(vp[[2]][[1]]), y = as.vector(vp[[2]][[2]])))

  AB <- polyclip::polyclip(A,B)

  ABC     <- polyclip::polyclip(A,B) %>% polyclip::polyclip(., C)

  AnoB <- polyclip::polyclip(A,B, op = "minus")
  noAB <- polyclip::polyclip(B,A, op = "minus")

  sets <- list("AB" = AB,  "AnoB" = AnoB, "noAB" = noAB)

  data <- data.frame(A=A, B=B)

  venn <- ggplot(data) +
    geom_polygon(mapping = aes(A.x, A.y), color ="black", fill = NA, alpha = 1, size = 0.2) +
    geom_polygon(mapping = aes(B.x, B.y), color ="black", fill = NA, alpha = 1, size = 0.2)

  for(i in 1:length(highlight)){
    set <- highlight[i]
    new_data <- data.frame(set=sets[[set]])
    venn <- venn + geom_polygon(data = new_data, aes(set.x,set.y), color = "black", fill = color[i], alpha = 1, size = 0.2)
  }

  venn <- venn + cowplot::theme_nothing()

  # Add label if desired
  if(!is.null(label)){
    if(label.pos %in% c("bottom.right", "right.bottom")) { lab.x = 1; lab.y = .2 }
    else if(label.pos %in% c("bottom.centre", "centre.bottom", "bottom.center", "center.bottom")) { lab.x = .5; lab.y = .15 }
    else if(label.pos %in% c("bottom.left", "left.bottom")) { lab.x = 0; lab.y = .2 }
    else if(label.pos %in% c("top.right", "right.top")) { lab.x = 1; lab.y = .8 }
    else if(label.pos %in% c("top.centre", "centre.top", "top.center", "center.top")) { lab.x = .5; lab.y = .85 }
    else if(label.pos %in% c("top.left", "left.top")) { lab.x = 0; lab.y = .8 }
    venn <- venn + annotate(geom = "text", label = label, x = lab.x, y = lab.y, hjust = .5, size = label.size, fontface = label.face, color = label.col)
  }

  venn
}


# Highlight areas of 3D Venn
venn3_highlight <- function(highlight = c("ABC", "AnoBnoC"), color = c("red", "pink"),
                            label = NULL, label.pos = "bottom.right", label.face = "italic", label.size = 13, label.col = "black",
                            line.col = "black"){

  # load required packages
  library(polyclip)
  library(VennDiagram)
  library(ggplot2)

  # Create base venn
  invisible(vp <- VennDiagram::draw.triple.venn(area1 = 0, area2 = 0, area3 = 0, n12 = 0, n23 = 0, n13 = 0, n123 = 0,
                                      category = rep("", 3), cex = 0, label.col = NA, scaled = F, fill = NA,
                                      margin = 0, disable.logging = T))

  # Take the circles
  A <- list(list(x = as.vector(vp[[1]][[1]]), y = as.vector(vp[[1]][[2]])))
  B <- list(list(x = as.vector(vp[[2]][[1]]), y = as.vector(vp[[2]][[2]])))
  C <- list(list(x = as.vector(vp[[3]][[1]]), y = as.vector(vp[[3]][[2]])))

  # Find the intersections and other areas
  AB <- polyclip::polyclip(A,B)
  AC <- polyclip::polyclip(A,C)
  BC <- polyclip::polyclip(B,C)

  # Intersection of all
  ABC     <- polyclip::polyclip(A,B) %>% polyclip::polyclip(., C)

  # Intersection of two, subtracting one
  ABnoC   <- polyclip::polyclip(A,B) %>% polyclip::polyclip(., C, op = "minus")
  AnoBC   <- polyclip::polyclip(A,C) %>% polyclip::polyclip(., B, op = "minus")
  noABC   <- polyclip::polyclip(C, B) %>% polyclip::polyclip(., A, op = "minus")

  # Specific of one area
  AnoBnoC <- polyclip::polyclip(A,B, op = "minus") %>% polyclip::polyclip(., C, op = "minus")
  noABnoC <- polyclip::polyclip(B,A, op = "minus") %>% polyclip::polyclip(., C, op = "minus")
  noAnoBC <- polyclip::polyclip(C,A, op = "minus") %>% polyclip::polyclip(., B, op = "minus")

  # Put all intersections in a list
  sets <- list("ABC" = ABC, "ABnoC" = ABnoC, "AnoBC" = AnoBC, "noABC" = noABC, "AnoBnoC" = AnoBnoC, "noABnoC" = noABnoC, "noAnoBC" = noAnoBC)

  # If highlight == "all" set to all clusters
  if(highlight == "all") { highlight <- names(sets); color <- rep(color, 7)}

  # Write a dataframe with the data from each circle
  data <- data.frame(A=A, B=B, C=C)

  # Draw the circles. Line color = black. Fill = NA.
  venn <- ggplot(data) +
    geom_polygon(mapping = aes(A.x, A.y), color = line.col, fill = NA, alpha = 1, size = 0.2) +
    geom_polygon(mapping = aes(B.x, B.y), color = line.col, fill = NA, alpha = 1, size = 0.2) +
    geom_polygon(mapping = aes(C.x, C.y), color = line.col, fill = NA, alpha = 1, size = 0.2)

  # Highlight desired areas by drawing them
  for(i in 1:length(highlight)){
    set <- highlight[i]
    new_data <- data.frame(set=sets[[set]])
    venn <- venn + geom_polygon(data = new_data, aes(set.x,set.y), color = line.col, fill = color[i], alpha = 1, size = 0.2)
  }

  # Remove all the elements from the grid, axes...
  venn <- venn + cowplot::theme_nothing()

  # Add label if desired
  if(!is.null(label)){
    if(label.pos %in% c("bottom.right", "right.bottom")) { lab.x = 1; lab.y = .2 }
    else if(label.pos %in% c("bottom.left", "left.bottom")) { lab.x = 0; lab.y = .2 }
    else if(label.pos %in% c("top.right", "right.top")) { lab.x = 1; lab.y = .9 }
    else if(label.pos %in% c("top.left", "left.top")) { lab.x = 0; lab.y = .9 }
    venn <- venn + annotate(geom = "text", label = label, x = lab.x, y = lab.y, hjust = .5, size = label.size, fontface = label.face, color = label.col)
  }

  # Return venn with shaded areas
  venn
}


# venn3_highlight(highlight = "ABC", color = "black")
#
# venn3_highlight(highlight = c("AnoBC","ABnoC"), color = c("black", "black"))
# venn3_highlight(highlight = c("ABnoC","noABC"), color = c("black", "black"))
# venn3_highlight(highlight = c("AnoBC","noABC"), color = c("black", "black"))
#
# venn3_highlight(highlight = c("AnoBnoC"), color = c("black"))
# venn3_highlight(highlight = c("noABnoC"), color = c("black"))
# venn3_highlight(highlight = c("noAnoBC"), color = c("black"))
