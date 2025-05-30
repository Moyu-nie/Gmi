#' Title
#'
#' @param x
#' @param threshold
#'
#' @returns
#' @export
#'
#' @examples
#'
Indexgroup <- function(x, threshold) {
  if (length(x) == 1) {
    return(1)
  }
  group <- rep(0, length(x))
  x0 <- x
  t <- 1
  while (length(x) != 0) {
    diff <- x[1] - x
    loc <- which(abs(diff) <= threshold)
    set <- which(x0 %in% unique(x[loc]))
    group[set] <- t
    t <- t + 1
    x <- x[-loc]
  }
  return(group)
}
