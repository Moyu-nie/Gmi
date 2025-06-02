
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Gmi

<!-- badges: start -->

<!-- badges: end -->

Gmi: Group-Level Main Effects and Interactions in High-Dimensional Data.

Gmi is a package for analyzing the high dimensional data with
group-level main and interaction effetcs, developed by the Huazhen Lin’s
lab.

## Installation

You can install the development version of Gmi from
[GitHub](https://github.com/), firstly, install the ‘remotes’ package.

``` r
install.packages("remotes")
remotes::install_github("Moyu-nie/Gmi")
```

## Example

This is a basic example which shows you how to solve a common problem:
\### Generate data

``` r
library(Gmi)
p <- 100
n <- 600
alpha2 <- c(rep(c(4,3,2),each = 2),rep(0,p-3*2))
gamma2 <- rep(0,choose(p,2))
mu2 <- 2 # intercept
rho2 <- 0.3 # correlation 

# interact between group1 and group2
aall = paste("X", combn(1 : p, 2, paste, collapse="X"), sep="")
aa1 = outer(1:2, 3:4, f <- function(x, y) {
  paste("X", pmin(x, y), "X", pmax(x, y), sep = "")
})
aa1 = as.vector(aa1)
gamma2[match(aa1, aall)] <- c(rep(4,4))
datList <- gendata(seed = 123, n ,p , alpha2, gamma2, mu2, rho2)
x2 = datList$x
Y2 = datList$Y
```

### Estimation

``` r
# fit
Gmi.fit = Gmi(x2, Y2, beta = rep(0,p), lambda.min.ratio = 0.1, n.lambda = 50, penalty = "SCAD", eta = 0.5, tune = "EBIC")
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], : p
#> must greater than 2
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X, y = y, a_0 = a_0, beta = beta, lam = lam, : ADMM
#> desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X, y = y, a_0 = a_0, beta = beta, lam = lam, : ADMM
#> desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X, y = y, a_0 = a_0, beta = beta, lam = lam, : ADMM
#> desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X, y = y, a_0 = a_0, beta = beta, lam = lam, : ADMM
#> desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X, y = y, a_0 = a_0, beta = beta, lam = lam, : ADMM
#> desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X, y = y, a_0 = a_0, beta = beta, lam = lam, : ADMM
#> desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X, y = y, a_0 = a_0, beta = beta, lam = lam, : ADMM
#> desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X, y = y, a_0 = a_0, beta = beta2, lam = lam, :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X, y = y, a_0 = a_0, beta = beta2, lam = lam, :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> Warning in SBfusedlasso(X = X[, ind1], y = y, a_0 = a_0, beta = beta[ind1], :
#> ADMM desn't converge in 20 iterations.
#> the method is finished
# estimate main effects
print(Gmi.fit$beta.m)
#>  [1]  1.1480599  0.7656352  0.7653269  0.7650195  0.7650415  0.3182777
#>  [7]  0.3197179  0.2058801  0.2045365 -0.2044563 -0.2037450
# estimate interaction effects
print(Gmi.fit$beta.i)
#> [1] 1.339334 1.339746
```
