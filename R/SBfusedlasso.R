#' Title
#'
#' @param X input matrix, of dimension nobs x nvars; each row is an observation
#' @param y response variable, of dimension nobs x
#' @param a_0
#' @param beta
#' @param lam
#' @param eta
#' @param rho1
#' @param rho2
#' @param penalty.type
#' @param epsilon1
#' @param epsilon2
#' @param maxiter1
#' @param maxiter2
#' @param pf
#'
#' @returns
#' @export
#'
#' @examples
SBfusedlasso <- function(X, y, a_0, beta, lam, eta, rho1 = 1, rho2 = 1,
                                   penalty.type, epsilon1=4e-7, epsilon2=4e-7, maxiter1 = 100, maxiter2 = 100, pf) {
  ### version add intercept
  # a_0 : intercept
  # require(cPCG)
  # require(Rlinsolve)
  # require(pcg)
  # eta must be [0,1]
  lambda1 =lam*eta
  lambda2 = lam*(1-eta)

  n = dim(X)[1]
  p = dim(X)[2]
  # fusion penalty coefficients
  if(p<=2){
    warning("p must greater than 2")
  }

  D <- matrix(0,p-1,p)
  for (i in 1:(p-1)) {
    D[i,i] = -1
    D[i,i+1] = 1
  }
  for (i1 in 1:maxiter1) {
    # IWLS
    beta0 = beta
    eta1 = X%*%beta + a_0
    p_hat = c(Sigmoid(eta1))
    w = p_hat*(1-p_hat)
    z = eta1 + (y-p_hat)/w


    z_weight_bar <- sum(w * z) / sum(w)
    x_weight_bar <- colSums(X * w) / sum(w)


    z_star <- sqrt(w) * (z - z_weight_bar)
    X_star <- sweep(X, 2, x_weight_bar, FUN = "-") * sqrt(w)

    ### inner step (ADMM) ####
    # initialization for ADMM
    u1 = rep(0,p)
    u2 = rep(0,p-1)
    a = beta
    b = D%*%beta
    for (i2 in 1:maxiter2) {
      a0 <- a
      b0 <- b
      ### step I
      # "cPCG package" method
      beta = pcgsolve(A = 1/n*eigen_matmul(t(X_star), X_star) + rho1*diag(p) + rho2*eigen_matmul(t(D), D),
                      b = 1/n*eigen_matmul(t(X_star), z_star) + rho1*(a-u1) + rho2*t(D)%*%(b-u2), preconditioner = "Jacobi")
      beta = as.vector(beta)
      ### step II
      if(penalty.type == "lasso"){
        # pf = 0 means some covirates is not be penalized
        a <- St(beta + u1, lambda1*pf/rho1)
        b <- St(D%*%beta + u2, lambda2/rho2)
      }
      ### MCP
      if(penalty.type == "MCP"){
        c = 3
        a = mapply(beta + u1, pf, FUN = function(x,y){if(abs(x) <= c*lambda1*y){x = St(x,  lambda1*y/rho1)/(1-1/(c*rho1))}else{x = x}})
        b = mapply(D%*%beta + u2, FUN = function(x){if(abs(x) <= c*lambda2){x = St(x, lambda2/rho2)/(1-1/(c*rho2))}else{x = x}})
      }

      ### SCAD
      if(penalty.type == "SCAD"){
        c = 3.7
        a = mapply(beta + u1, pf, FUN = function(x,y){if(abs(x) <= lambda1*y + lambda1*y/rho1){x = St(x, lambda1*y/rho1)}else if(((lambda1*y + lambda1*y/rho1) < abs(x)) && (abs(x) <= c*lambda1 )){ x = St(x,c*lambda1*y/((c-1)*rho1))/(1-1/(c-1)*rho1)}else{x = x}})
        b = mapply(D%*%beta + u2, FUN = function(x){if(abs(x) <= lambda2 + lambda2/rho2){x = St(x, lambda2/rho2)}else if(((lambda2 + lambda2/rho2) < abs(x)) && (abs(x) <= c*lambda2 )){ x = St(x,c*lambda2/((c-1)*rho2))/(1-1/(c-1)*rho2)}else{x = x}})
      }
      ### step II
      u1 <- u1 + beta - a
      u2 <- u2 + D%*%beta - b
      #### stoping criteria
      stop.primal = mean((beta-a)^2)+ mean((D%*%beta-b)^2)
      stop.dual = mean((a-a0)^2) + mean((b-b0)^2)
      if((stop.primal<=epsilon1)&(stop.dual<=epsilon2)){
        break
      }
    }

    if(i2==maxiter2){
      warning(sprintf("ADMM desn't converge in %d iterations.",maxiter2))
    }
    ### estimate a_0:intercept
    a_0 = c(z_weight_bar - x_weight_bar%*%beta)

    if(mean((beta0 - beta)^2) < 1e-6){
      break
    }
  }
  # output
  result = list(beta = a, a_0 = a_0, iters=i1)
  return(result)
}
