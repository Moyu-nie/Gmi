
#' Title
#'
#' @param x
#' @param lambda
#'
#' @returns
#' @export
#'
#' @examples
St <- function(x, lambda) {
  # soft threshholding
  val <- sign(x) * pmax(abs(x) - lambda, 0)
  return(val)
}
