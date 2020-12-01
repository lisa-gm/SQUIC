# demo function for using generate sample for SQUIC

# load squic
library(SQUIC)

##### run demo #########
# with preset parameters
run_demo <- function(lambda = 0.2, max_iter = 1000,
                     drop_tol = 10e-6, term_tol = 10e-6, verbose = 0,
                     del_files = TRUE){
  
  Y_file <- system.file("data", "sample_data_10x100.mat", package = "SQUIC", mustWork = TRUE)
  Y <- data.table::fread(Y_file, header = FALSE)

  squic_res <- SQUIC(Y = Y, lambda = lambda, X_pattern=NULL, max_iter=max_iter, 
                     drop_tol=drop_tol, term_tol=term_tol, 
                     X_init = NULL, W_init = NULL, 
                     verbose=verbose, del_files = del_files)
  print(squic_res$X)
  
  print(" ")
  print("estimated inverse covariance matrix : ")
  print(round(squic_res$X,2))
  print(squic_res$log)

  return(list(X=squic_res$X, W=squic_res$W, log=squic_res$log))

}

res_demo <- run_demo()

