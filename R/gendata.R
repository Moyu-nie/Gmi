#' generate data
#'
#' @param seed
#' @param n
#' @param p
#' @param alpha
#' @param gamma
#' @param beta
#' @param rho
#'
#' @returns
#' @export
#'
#' @examples
gendata <- function(seed, n, p , alpha, gamma, beta, rho){
  set.seed(1)
  Sigma = rho^(abs(outer(1:p,1:p,"-")))
  x <- MASS::mvrnorm(n, mu = rep(0,p), Sigma = Sigma)
  x <- scale(x, center = TRUE,scale = TRUE)
  # generate Y
  # generate all possible two-way interaction term (matrix)
  xx <- t(apply(x, 1, combn, 2, prod))
  #xx <- scale(xx, center = TRUE,scale = TRUE)
  # generate theta:intercept

  theta <- beta + x%*%alpha + xx%*%gamma
  set.seed(seed)
  Y <- rbinom(n, 1, Sigmoid(theta))
  return(list(x = x , xx = xx , Y = Y))
}

#' Title
#'
#' @param seed
#' @param n
#' @param p
#' @param alpha
#' @param mu
#' @param rho
#'
#' @returns
#' @export
#'
#' @examples
gendata_simp <- function(seed, n , p, alpha, mu, rho){
  set.seed(1)
  Sigma = rho^(abs(outer(1:p,1:p,"-")))
  x <- MASS::mvrnorm(n, mu = rep(0,p), Sigma = Sigma)
  x <- scale(x, center = TRUE,scale = TRUE)
  theta <- mu + x%*%alpha + 4*x[,1]*x[,3] + 4*x[,1]*x[,4] + 4*x[,2]*x[,3] + 4*x[,2]*x[,4]
  set.seed(seed)
  Y <- rbinom(n, 1, Sigmoid(theta))
  return(list(x = x ,  Y = Y))
}

