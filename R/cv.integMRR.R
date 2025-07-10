### cv.integMRR
### Note that this function requres the function "integMRR.R".

### INPUT
# x: a list of design matrices
# y: a list of response matrices
# fold: the number of fold - default is 5.
# lambda_candidate: candidate values of regularization parameter in group lasso penalty. optional. 
# nu_candidate: candidate values of regularization parameter in covariance matrix. optional. 
# rho: penalty parameter in ADMM - default is 1.
# tol: tolerance of the algorithm - default is 1e-5.

### OUTPUT
# lambda_CV: lambda that minimizes CV
# nu_CV: nu that minimizes CV
# CV_matrix: matrix of CV 
# lambda_seq: sequence of lambda
# nu_seq: nu of lambda

source("./integMRR.R")

cv.integMRR <- function(x, y, fold = 5, lambda_candidate = NULL, nu_candidate = NULL, rho = 1, tol = 1e-5)
{
  ### the number of datasets
  M <- length(x)

  ### Preparation of CV
  options(warn = -1)
  extract <- as.list(NULL)
  for(i in 1:M) extract[[i]] <- cbind(matrix(1:fold, nrow(y[[i]])), sample(1:nrow(y[[i]])))
  xall <- x
  yall <- y

  ### Candidate values of tuning parameters (lamdba, gamma)
  if( is.null(lambda_candidate) ) lambda_candidate <- c(1e-0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5)
  if( is.null(nu_candidate) ) nu_candidate <- c(1e+1, 1e-0, 1e-1, 1e-2)
  CV_mat <- matrix(0, length(lambda_candidate), length(nu_candidate))
  
  ### "fold"-fold CV
  for(itr_CV in 1:fold)
  {
    x_ori <- y_ori <- x_test_cv <- y_test_cv <- 
      x_centered <- y_centered <- extract_num <- as.list(NULL)
    for(i in 1:M) extract_num[[i]] <- extract[[i]][extract[[i]][ ,1] == itr_CV, -1]

    for(i in 1:M)
    {
      x_ori[[i]] <- xall[[i]][-extract_num[[i]], ]
      y_ori[[i]] <- yall[[i]][-extract_num[[i]], ]
      
      x_test_cv[[i]] <- xall[[i]][extract_num[[i]], ]
      y_test_cv[[i]] <- yall[[i]][extract_num[[i]], ]
    }
    
    for( itr_lambda in 1:length(lambda_candidate) )
    {
      lambda <- lambda_candidate[itr_lambda]
      for( itr_nu in 1:length(nu_candidate) )
      {
        nu <- nu_candidate[itr_nu]
        result <- integMRR(x = x_ori, y = y_ori, lambda = lambda, nu = nu, rho = rho, tol = tol)
        B <- result$B
        alpha <- result$alpha
        CVerror <- 0
        for(i in 1:M) CVerror <- CVerror + 
          mean( ( y_test_cv[[i]] - rep(1, nrow(y_test_cv[[i]]))%*%t(alpha[[i]]) - x_test_cv[[i]]%*%B[[i]])^2 )
        CV_mat[ itr_lambda, itr_nu ] <- CV_mat[ itr_lambda, itr_nu ] + CVerror
      }
    }
  }
  
  ### average of CV-value
  CV_mat <- CV_mat/fold
  minimum_value <- 0
  for(i in 1:nrow(CV_mat)){
    for(j in 1:ncol(CV_mat)){
      if( (i==1) && (j==1) ){
        minimum_value <- CV_mat[i, j]
        stock_row <- stock_col <- 1
      }
      if( minimum_value > CV_mat[i, j] ){
        stock_row <- i
        stock_col <- j
        minimum_value <- CV_mat[i, j]
      }
    }
  }
  return( list(lambda_CV = lambda_candidate[stock_row], nu_CV = nu_candidate[stock_col], CV_mat = CV_mat, lambda_seq = lambda_candidate, nu_seq = nu_candidate) )
}
