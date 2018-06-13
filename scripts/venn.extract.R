# This code was copied from https://stackoverflow.com/questions/23559371/how-to-get-the-list-of-items-in-venn-diagram-in-r

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


get.venn.dat <- function(xx.1) {
  combs <- 
    unlist(lapply(1:length(xx.1), 
                  function(j) combn(names(xx.1), j, simplify = FALSE)),
           recursive = FALSE)
  names(combs) <- sapply(combs, function(i) paste0(i, collapse = ""))
  elements <- 
    lapply(combs, function(i) Setdiff(xx.1[i], xx.1[setdiff(names(xx.1), i)]))
}






