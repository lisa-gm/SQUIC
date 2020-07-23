# demo function for using generate sample for SQUIC

# load libraries to generate sample data
library(mvtnorm)

# load squic
library(SQUIC)

#############################################################################
# generate sample data
# sample from p-dimensional y_i ~ N(0,Sigma)
# with inv(Sigma) = tridiag(-0.5, 1.25, -0.5)

# generates n samples of size p (Y \in dim n x p)
generate_sample <- function(p, n){

  ld <- 1.25*rep.int(1, p)
  sd <- -0.5*rep.int(1, p-1)

  # construct sparse matrix
  diags_invS <- list(ld, sd)
  invSigma <- bandSparse(p, k = c(0,1), diag = diags_invS, symm=TRUE)

  # get diagonals of chol factor, i.e. Sigma = L*L^t
  L <- chol(invSigma)
  L_inv <- solve(L)

  Id <- diag(p)
  mu <- integer(p)

  # suggest package to generate sample data
  Z <- rmvnorm(n, mean = mu, sigma = Id)

  # create sample matrix: Y = Z * L^(-1)
  Y = as.matrix(Z %*% t(L_inv))

  return(Y)
}

##### run demo #########
# with preset parameters
run_demo <- function(p = 12, n=10000, lambda = 0.0, max_iter = 1000,
                     drop_tol = 10e-6, term_tol = 10e-6, verbose = 0,
                     del_files = TRUE){

  Y <- generate_sample(p,n)

  squic_res <- SQUIC(p, n, Y, lambda, max_iter, drop_tol, term_tol, verbose, del_files)
  print(" ")
  print("estimated inverse covariance matrix : ")
  print(round(squic_res$X,2))

  return(list(X=squic_res$X, W=squic_res$W))

}

res_demo <- run_demo()
