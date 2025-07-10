### integMR

### INPUT
# x: a list of design matrices
# y: a list of response matrices
# lambda: regularization parameter in group lasso penalty
# rho: penalty parameter in ADMM - default is 1.
# tol: tolerance of the algorithm - default is 1e-5.

### OUTPUT
# B: estimate of the parameter B
# H: estimate of the parameter eta
# U: estimate of the parameter u
# alpha: estimate of the parameter alpha
# Sigma: estimate of the parameter Sigma

### Soft-threshold operator for group lasso
softth <- function(a, b)
{
  z <- 1 - b/(sqrt(sum(a^2)))
  if(z < 0)  return(rep(0, length(a)))
  else return( z * a )
}
### Soft-threshold operator for lasso
softth_elementwise <- function(z, eta)
{
  return( sign(z)*pmax(abs(z) - eta, 0) )
}

### multivariate regression in integrative analysis
integMR <- function(x, y, lambda, rho = 1, tol = 1e-5)
{
  ### the number of datasets
  M <- length(x)

  ### initialization of parameters
  np <- ncol(x[[1]])
  nq <- ncol(y[[1]])
  Sigma <- B <- H <- U <- alpha <- as.list(NULL)
  b <-  matrix(0, np, nq)
  for(i in 1:M) B[[i]] <- b
  U <- H <- B
  for(i in 1:M) alpha[[i]] <- rep(0, nq)
  for(i in 1:M) Sigma[[i]] <- cor(y[[i]])
  
  ### convergence condition - initialization -
  L2_term <- 0
  for(i in 1:nrow(H[[1]])){
    for(j in 1:ncol(H[[1]])){
      L2_Term <- 0
      for(k in 1:M) L2_Term <- L2_Term + H[[k]][i,j]^2
      L2_term <- L2_term + sqrt( L2_Term )
    }
  }
  loss_new <- 0
  for(i in 1:M) loss_new <- loss_new + 
    0.5*mean(diag( ( ( y[[i]] - rep(1, nrow(y[[i]]))%*%t(alpha[[i]]) ) - x[[i]]%*%B[[i]] ) %*% 
                solve(Sigma[[i]]) %*% t( ( ( y[[i]] - rep(1, nrow(y[[i]]))%*%t(alpha[[i]]) ) - x[[i]]%*%B[[i]] ) ) ) ) + 
    0.5*rho*sum( ( H[[i]] - B[[i]] + U[[i]] )^2 ) + 0.5*log(det(Sigma[[i]]))
  loss_new <- loss_new + lambda*L2_term 
  condition <- 1e+10
  
  ################# ADMM START #################
  while( condition > tol )
  {
    ### compute inverse matrix
    x.inv <- z.inv <- as.list(NULL)
    for(i in 1:M) x.inv[[i]] <- solve( kronecker( solve(Sigma[[i]]), t(x[[i]]) %*% x[[i]] ) + nrow(x[[i]])*rho*diag(np*nq) )

    ### Update of B
    for(i in 1:M){
      B.stock <- x.inv[[i]] %*% c( t(x[[i]]) %*% ( ( y[[i]] - rep(1, nrow(y[[i]]))%*%t(alpha[[i]]) ) ) %*% 
                                     solve(Sigma[[i]]) + nrow(x[[i]])*rho*( H[[i]] + U[[i]] ) )
      B[[i]] <- matrix(B.stock, np, nq)
    }
    
    ### Update of eta
    for(i in 1:nrow(H[[1]])){
      for(j in 1:ncol(H[[1]])){
        a <- rep(0, M)
        for(k in 1:M) a[k] <- B[[k]][i,j] - U[[k]][i,j]
        result_a <- softth( a, lambda/rho )
        for(k in 1:M) H[[k]][i,j] <- result_a[k]
      }
    }
    
    ### Update of alpha (alpha is a vector of intercept)
    for(i in 1:M) alpha[[i]] <- t( y[[i]] - x[[i]]%*%B[[i]] )%*% rep(1, nrow(y[[i]])) / nrow(y[[i]])
    
    ### Update of Sigma
    for(i in 1:M) Sigma[[i]] <- t( ( y[[i]] - rep(1, nrow(y[[i]]))%*%t(alpha[[i]]) ) - x[[i]]%*%B[[i]] ) %*%
        ( ( y[[i]] - rep(1, nrow(y[[i]]))%*%t(alpha[[i]]) ) - x[[i]]%*%B[[i]] )/nrow(y[[i]])

    ### Update of u
    for(i in 1:nrow(U[[1]])){
      for(j in 1:ncol(U[[1]])){
        for(k in 1:M) U[[k]][i,j] <- U[[k]][i,j] + H[[k]][i,j] - B[[k]][i,j]
      }
    }

    ### convergence condition
    L2_term <- 0
    for(i in 1:nrow(H[[1]])){
      for(j in 1:ncol(H[[1]])){
        L2_Term <- 0
        for(k in 1:M) L2_Term <- L2_Term + H[[k]][i,j]^2
        L2_term <- L2_term + sqrt( L2_Term )
      }
    }
    loss_old <- loss_new
    loss_new <- 0
    for(i in 1:M) loss_new <- loss_new + 
      0.5*mean(diag( ( ( y[[i]] - rep(1, nrow(y[[i]]))%*%t(alpha[[i]]) ) - x[[i]]%*%B[[i]] ) %*% 
                       solve(Sigma[[i]]) %*% t( ( ( y[[i]] - rep(1, nrow(y[[i]]))%*%t(alpha[[i]]) ) - x[[i]]%*%B[[i]] ) ) ) ) + 
      0.5*rho*sum( ( H[[i]] - B[[i]] + U[[i]] )^2 ) + 0.5*log(det(Sigma[[i]]))
    loss_new <- loss_new + lambda*L2_term 
    condition <- abs( loss_old - loss_new )/loss_old
  }
  ################# ADMM END #################
  return( list( B = B, H = H, U = U, alpha = alpha, Sigma = Sigma ) )
}
