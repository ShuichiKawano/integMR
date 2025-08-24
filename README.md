# Sparse multivariate regression in integrative analysis (integMR)
Author: Shuichi Kawano, Toshikazu Fukushima, Junichi Nakagawa, Mamoru Oshiki

This is an R source code for performing sparse multivariate regression in integrative analysis (integMR). The directory *R* consists of four files as follows:
- **integMR.R** computes sparse multivariate regression in integrative analysis.
- **cv.integMR.R** computes cross-validation of integMR. 
- **integMRR.R** provides sparse multivariate regression with regularized covariance matrix in integrative analysis.
- **cv.integMRR.R** computes cross-validation of integMRR.
  
integMR is introduced in the paper:
Kawano, S., Fukushima, T., Nakagawa, J. and Oshiki, M. (2025) Multivariate regression modeling in integrative analysis via sparse regularization. _Japanese Journal of Statistics and Data Science (Online Access)_ (doi: [10.1007/s42081-025-00312-2](https://doi.org/10.1007/s42081-025-00312-2)).

## Usage example
Read source files.
```
source("R/integMR.R")
source("R/integMRR.R")
source("R/cv.integMR.R")
source("R/cv.integMRR.R")
```

Setting of simulation.
```
library(MASS)
rho_x <- 0.1; rho_y <- 0.1
num_noise_variable <- 5
nn1 <- nn2 <- 50
set.seed(1)
x <- y  <-  as.list(NULL)
true_param <- as.list(NULL)
### setting coefficients
B <- cbind(c(c(0.5,0.5,0.5,0.25,0.25), rep(0, 5)), c(rep(0, 5), c(0.25,0.25,0.15,0.15,0.15)))
true_param[[2]] <- true_param[[1]] <- rbind(B, matrix(0, num_noise_variable, ncol(B)))
### generating x
Sigma_x <- diag( rep(1, nrow(B)) )
for(i in 1:nrow(Sigma_x)) for(j in 1:ncol(Sigma_x)) Sigma_x[i ,j] = rho_x^abs(i-j)
x[[1]] <- mvrnorm(nn1, rep(0, nrow(Sigma_x)), Sigma_x)
x[[1]] <- cbind(x[[1]], matrix(rnorm(nrow(x[[1]])*num_noise_variable), nrow(x[[1]]), num_noise_variable))
x[[2]] <- mvrnorm(nn2, rep(0, nrow(Sigma_x)), Sigma_x)
x[[2]] <- cbind(x[[2]], matrix(rnorm(nrow(x[[2]])*num_noise_variable), nrow(x[[2]]), num_noise_variable))
### generating y
Sigma_y <- diag( rep(1, ncol(B)) )
for(i in 1:nrow(Sigma_y)) for(j in 1:ncol(Sigma_y)) Sigma_y[i ,j] = rho_y^abs(i-j)
y[[1]] <- x[[1]]%*%true_param[[1]] + mvrnorm(nrow(x[[1]]), rep(0, nrow(Sigma_y)), Sigma_y)
y[[2]] <- x[[2]]%*%true_param[[2]] + mvrnorm(nrow(x[[2]]), rep(0, nrow(Sigma_y)), Sigma_y)
```

Perform integMR
```
# Perform integMR in the file integMR.R
integMR(x = x, y = y, lambda = 0.3)

# Perform cv.integMR with five-fold in the file cv.integMR.R
cv.integMR(x = x, y = y)

# Perform integMRR in the file integMRR.R
integMRR(x = x, y = y, lambda = 0.3, nu = 1e-3)

# Perform cv.integMRR with five-fold in the file cv.integMRR.R
cv.integMRR(x = x, y = y)
```
