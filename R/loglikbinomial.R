#' Title
#'
#' @param X
#' @param y
#' @param beta
#'
#' @returns
#' @export
#'
#' @examples
loglikbinomial <- function(X, y, beta) {
  link = as.vector(X %*% beta)
  n = length(y)
  return(2 * sum(log(1 + exp(link)) - y * link))
}
