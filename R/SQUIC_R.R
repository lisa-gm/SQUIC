
# import Matrix package, automatically updates DESCRIPTION file
# is this the right spot?
usethis::use_package("Matrix")
usethis::use_package("data.table")

##########################################################
# generate some sort of description, need to update this later
#' Sparse Quadratic Inverse Covariance Estimation
#'
#' This function computes large Sparse Inverse Covariance Matrices quickly. 
#' It returns an estimate of the sparse covariance matrix W and its inverse X of a data matrix Y. 
#' 
#' @param Y Input data, expects a matrix of dimension p x n, where each sample corresponds to one row of p random variables and a total number of n samples.
#' @param lambda Regularization parameter. Non-negative scalar or even positive? Each entry of X will be penalized in the same way with this value. Can be adapted for certain entries using the X_pattern parameter
#' @param X_pattern Optional matrix parameter. Lambda = Lambda_matrix + lambda_scala x (1-pattern(Lambda_matrix) ) (update later). Needs to be given in coo format. See details for further information.
#' @param X0 Optional parameter. Matrix with a sparse initial guess of the inverse covariance matrix. Needs to be given in coo format. Defaults to the identity matrix.
#' @param W0 Optional parameter. Matrix with a sparse initial guess of the covariance matrix. Needs to be given in coo format.Defaults to the identity matrix.
#' @param max_iter Maximum number of iterations.
#' @param drop_tol Drop out tolerance threshold for inverse of the inverse covariance matrix. Defaults to 1e-6.
#' @param term_tol Convergence threshold for iteration. Defaults to 1e-6.
#' @param verbose Integer value. Defining the verbosity of the print statistics. The higher the value the more detailed the output gets. The value range is from: ... Defaults to zero.
#' @return X Estimated sparse inverse covariance matrix
#' @return W Estimated sparse inverse of the inverse covariance matrix
#' @details R version of SQUIC.
#' The higher the values are for lambda the more sparsity is enforced. Possible source of how to choose lambda?
#' 
#' Penalty terms: (update this later). Each entry of the covariance is penalised by the scalar value lambda. In order to individually penalise a particular entry X_ij of the inverse covariance matrix the 
#' X_pattern matrix can be used. To not penalise an entry X_ij at all, one has to specify -lambda for this entry as the penalty term is internally computed as lambda * 1 * 1^T - X_pattern. 
#' 
#' The optional parameters X0 and W0 are for sparse initial guesses of X and W. 
#' They need to be given in COO-Format. Even though they are symmetric the all non-zero entries have to be specified. If no initial guess if known they will default to identity matrices. 
#' @author 
#' @references Large-scale Sparse Inverse Covariance Matrix Estimation, M.Bollhöfer, A. Eftekhari, S. Scheidegger, O. Schenk \url{https://epubs.siam.org/doi/abs/10.1137/17M1147615?journalCode=sjoce3}
#' @seealso \link{BigQuic}, \link{QUIC}
#' @keywords Sparse Covariance Matrix Estimation
#' @importFrom Matrix sparseMatrix
#' @export
#' @examples SQUIC(Y = NULL, lambda = 0.0, X_pattern=NULL, max_iter=1, drop_tol=1e-6, term_tol=1e-6, X_init = NULL, W_init = NULL, verbose=1, del_files = TRUE). TODO: write an example that compares to 
#' SQUIC()

SQUIC <- function(Y = NULL, lambda = 0.0, X_pattern=NULL, max_iter=1, drop_tol=1e-6, term_tol=1e-6, X_init = NULL, W_init = NULL, verbose=1, del_files = TRUE){

  pkg_path = system.file(package = "SQUIC")
  
  # creating paths to files for all file data, that is input and output : X, W, info.log
  # make new directory w/ timestamp for results in original working directory
  old_wd <- getwd()
  res_folder <- paste(old_wd, "/squic_out_", format(Sys.time(), "%d_%m_%y_%H_%M_%OS2"), sep = "")
  dir.create(res_folder)
  
  #######################################################################
  # SYSTEM REQUIREMENTS
  # check if 64bit version
  if(grepl("64", Sys.info()[["machine"]]) == FALSE){
    warning("It seems like you are not using a 64-bit machine. This might lead to problems.")
  }
  
  # CHOOSE EXECUTABLE DEP. ON PLATFORM
  # OSX
  if(Sys.info()[["sysname"]] == "Darwin" && grepl("19", Sys.info()[["release"]])){
    sys_path <- "libs/osx/darwin19"
    squic_exe <- system.file(sys_path, "SQUIC_CMD", package = "SQUIC", mustWork = TRUE)
  } else if(Sys.info()[["sysname"]] == "Darwin" && grepl("17", Sys.info()[["release"]])){
    sys_path <- "libs/osx/darwin17"
    squic_exe <- system.file(sys_path, "SQUIC_CMD", package = "SQUIC", mustWork = TRUE)
  # LINUX  
  } else if(Sys.info()[["sysname"]] == "Linux"){
    sys_path <- "libs/linux"
    squic_exe <- system.file(sys_path, "SQUIC_CMD", package = "SQUIC", mustWork = TRUE)
  # WINDOWS: ifelse(){}
  } else {
    warning("platform not suitable!")
  }
 

  # CALL SQUIC, first set path
  lib_wd <- paste(pkg_path, "/", sys_path, sep = "")
  setwd(lib_wd)
  
  # WRITING MATRIX PARAMETERS TO FILE FOR CMD EXECUTABLE
  ########################################################################################333
  # FOR PACKAGE BUILDING: automatically load sample Y file if no data is provided
  if(is.null(Y)){
    input_file <- system.file("data", "sample_data_10x100.mat", package = "SQUIC", mustWork = TRUE)
    p = 10
    n = 100
  } else{
    # get p and n from Y
    # p: number of random variables
    p = ncol(Y)
    # n: number of samples
    n = nrow(Y)
    # for console executable: write Y to file in R --> executable will get it from file
    input_file = file.path(res_folder, "Y_file.mat")
    # replace by fwrite to make it faster
    data.table::fwrite(Y, input_file, append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
  }
  
  #######################################################################
  # handling optional parameters
  
  # check if X_pattern/lambda_2 matrix has entries
  if(is.null(X_pattern)){
    LambdaMatrix_file = 'NULL'
  } else {
    # save sparse coo matrix to file
    LambdaMatrix_file = file.path(res_folder, "lambda_matrix.dat")
    # TODO: double check that summary() always gives the right pattern
    data.table::fwrite(summary(X_pattern), LambdaMatrix_file, append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
  }
  
  # initial guess X_init
  if(is.null(X_init)){
    X0_loc = 'NULL'
  } else {
    X0_loc = file.path(res_folder, "X0_loc.dat")
    # expects X_init to be in sparse matrix format
    data.table::fwrite(summary(X_init), X0_loc, append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
  }  
  
  # initial guess W_init
  if(is.null(W_init)){
    W0_loc = 'NULL'
  } else {
    X0_loc = file.path(res_folder, "W0_loc.dat")
    # expects X_init to be in sparse matrix format
    data.table::fwrite(summary(W_init), X0_loc, append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
  }
  
  #######################################################################
  # [integer:index_offset] Offset of indexing
  index_offset=1
  
  verbose_R = ""
  if(verbose == 0){
    verbose_R = NULL
  }
  
  X_loc <- file.path(res_folder, "X.dat")
  file.create(X_loc)
  W_loc <- file.path(res_folder, "W.dat")
  file.create(W_loc)
  log_loc <- file.path(res_folder, "info.log")
  file.create(log_loc)
  
  #######################################################################
  # SYSTEM CALL to squic_cmd
  args=paste(p,n, input_file, lambda, LambdaMatrix_file, max_iter, drop_tol, term_tol, X0_loc, W0_loc, index_offset, verbose, X_loc, W_loc, log_loc)

  system2(squic_exe, args = args, stdout = verbose_R)
  print(squic_exe)
  print(args)
  
  # just for testing, delete later.
  X = 0
  W = 0
  
  k = TRUE
  if(k==TRUE){
    
    # collect data from file
    X <- data.table::fread(X_loc,  sep="\t", header=FALSE)
    W <- data.table::fread(W_loc,  sep="\t", header=FALSE)
    
    # generate sparse matrices X,W
    # offset already included
    X<- Matrix::sparseMatrix(i=X[,1],j=X[,2],x=X[,3],dims=list(p,p))
    W<- Matrix::sparseMatrix(i=W[,1],j=W[,2],x=W[,3],dims=list(p,p))
    
    # read log file
    log_all <- scan(log_loc, sep = ",", what = character())
    # format the log file into a character vector
    log <- log_all[log_all != ""]
    log <- log[3:length(log)-1]
    
    # if del_files == TRUE: delete output folder, otherwise stays.
    if(del_files == TRUE){
      unlink(res_folder, recursive = TRUE)
    } else {
      print("result files are stored in your current working directory: ")
      print(res_folder)
    }
  }
  
  # change wd back to what it was before
  setwd(old_wd)
  
  # return estimated inverse cov matrix X and its estimated inverse W
  return (list(X=X, W=W, log=log[2:length(log)]))
}

