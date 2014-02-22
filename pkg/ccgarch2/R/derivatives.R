#####################
## cDCC GARCH model ##
#####################
## functions for computing numerical derivatives
loglik2.cdcc.t <- function(param, data, cDCC, mode){        # data is the standardised residuals
  if(is.zoo(data)) data <- as.matrix(data)
  nobs <- dim(data)[1]
  ndim <- dim(data)[2]
  cdcc <- cdcc.est(data, param)
  DCC <- cdcc$cDCC

  lf1 <- dcc.ll2(DCC, data)    
  if(mode=="gradient"){
      lf1
  } else {
      sum(lf1)
  }
  
#   lf1 <- numeric(ndim)
#   for( i in 1:nobs){
#     R1 <- matrix(DCC[i,], ndim, ndim)
#     invR1 <- solve(R1)
#     # lf1[i] <- -0.5*(log(det(R1)) - sum(data[i,]*crossprod(invR1, data[i,])))  
#     lf1[i] <- -0.5*(log(det(R1)) +sum(data[i,]%*%invR1%*%data[i,]))  
#   }
#   if(mode=="gradient"){
#     lf1 + 0.5*rowSums(data^2)  # the second term is unrelated with the optimization, but is included for computing log-lik value
#   } else {
#     sum(lf1) + 0.5*sum(data^2)
#   }
}

cdcc.hessian <- function(param, data, model){        # data is the original one
   if(is.zoo(data)) data <- as.matrix(data)
   
   nobs <- dim(data)[1]
   ndim <- dim(data)[2]
   In <- diag(ndim)
   mu <- matrix(param[1:ndim], nobs, ndim, byrow = TRUE)    # constant in the mean
   data <- data - mu
   param <- param[-(1:ndim)]

   if(model=="diagonal"){
    a <- param[1:ndim]; param <- param[-(1:ndim)]
    A <- diag(param[1:ndim]); param <- param[-(1:ndim)]
    B <- diag(param[1:ndim]); param <- param[-(1:ndim)]    # now "param" contains cDCC parameters
   } else {
    a <- param[1:ndim]; param <- param[-(1:ndim)]
    A <- matrix(param[1:ndim^2], ndim, ndim); param <- param[-(1:ndim^2)]
    B <- matrix(param[1:ndim^2], ndim, ndim); param <- param[-(1:ndim^2)]    # now "param" contains cDCC parameters
   }
   h <- vgarch(a, A, B, data)    # a call to vgarch function
   z <- data/sqrt(h)             # computing the standardized residuals
   
   cdcc <- cdcc.est(z, param)
   DCC <- cdcc$cDCC
   lf1 <- dcc.ll2(DCC, z)    
   sum(lf1)
#   lf1 <- numeric(ndim)
#   for( i in 1:nobs){
#     R1 <- matrix(DCC[i,], ndim, ndim)
#     invR1 <- solve(R1)
#     # lf1[i] <- -0.5*(log(det(R1)) - sum(z[i,]*crossprod(invR1, z[i,])))  
#     lf1[i] <- -0.5*(log(det(R1)) + sum(z[i,]%*%invR1%*%z[i,]))  
#   }
#     sum(lf1) + 0.5*sum(z^2)
}

# log-likelihood function of the GARCH part for computing Jacobian and Hessian
loglik1.dcc.t <- function(param, data, model, mode){
   if(is.zoo(data)) data <- as.matrix(data)
   nobs <- dim(data)[1]
   ndim <- dim(data)[2]
   In <- diag(ndim)
   mu <- matrix(param[1:ndim], nobs, ndim, byrow = TRUE)    # constant in the mean
   data <- data - mu
   param <- param[-(1:ndim)]

   if(model=="diagonal"){
    a <- param[1:ndim]; param <- param[-(1:ndim)]
    A <- diag(param[1:ndim]); param <- param[-(1:ndim)]
    B <- diag(param[1:ndim])
   } else {
    a <- param[1:ndim]; param <- param[-(1:ndim)]
    A <- matrix(param[1:ndim^2], ndim, ndim); param <- param[-(1:ndim^2)]
    B <- matrix(param[1:ndim^2], ndim, ndim)
   }
   h <- vgarch(a, A, B, data)    # a call to vgarch function
   z <- data/sqrt(h)
   lf <- -0.5*ndim*log(2*pi)-0.5*rowSums(log(h))-0.5*rowSums(z^2)
   
   if(mode=="gradient"){
    lf
   } else {
    sum(lf)
   }
}

#####################
## DCC GARCH model ##
#####################
## functions for computing numerical derivatives
loglik2.dcc.t <- function(param, data, DCC, mode){        # data is the standardised residuals
   if(is.zoo(data)) data <- as.matrix(data)
  nobs <- dim(data)[1]
  ndim <- dim(data)[2]
  dcc <- dcc.est(data, param)
  DCC <- dcc$DCC
  
   lf1 <- dcc.ll2(DCC, data)    
   if(mode=="gradient"){
       lf1
   } else {
       sum(lf1)
   }
   
#    lf1 <- numeric(ndim)
#   for( i in 1:nobs){
#     R1 <- matrix(DCC[i,], ndim, ndim)
#     invR1 <- solve(R1)
#     # lf1[i] <- -0.5*(log(det(R1)) - sum(data[i,]*crossprod(invR1, data[i,])))  
#     lf1[i] <- -0.5*(log(det(R1)) + sum(data[i,]%*%invR1%*%data[i,]))  
#   }
#   if(mode=="gradient"){
#     lf1 + 0.5*rowSums(data^2)  # the second term is unrelated with the optimization, but is included for computing log-lik value
#   } else {
#     sum(lf1) + 0.5*sum(data^2)
#   }
}

dcc.hessian <- function(param, data, model){        # data is the original one
   if(is.zoo(data)) data <- as.matrix(data)
   nobs <- dim(data)[1]
   ndim <- dim(data)[2]
   In <- diag(ndim)
   mu <- matrix(param[1:ndim], nobs, ndim, byrow = TRUE)    # constant in the mean
   data <- data - mu
   param <- param[-(1:ndim)]

   if(model=="diagonal"){
    a <- param[1:ndim]; param <- param[-(1:ndim)]
    A <- diag(param[1:ndim]); param <- param[-(1:ndim)]
    B <- diag(param[1:ndim]); param <- param[-(1:ndim)]    # now "param" contains DCC parameters
   } else {
    a <- param[1:ndim]; param <- param[-(1:ndim)]
    A <- matrix(param[1:ndim^2], ndim, ndim); param <- param[-(1:ndim^2)]
    B <- matrix(param[1:ndim^2], ndim, ndim); param <- param[-(1:ndim^2)]    # now "param" contains cDCC parameters
   }
   h <- vgarch(a, A, B, data)    # a call to vgarch function
   z <- data/sqrt(h)             # computing the standardized residuals
   
   dcc <- dcc.est(z, param)
   DCC <- dcc$DCC
  
   lf1 <- dcc.ll2(DCC, z)    
   sum(lf1)
   
#    lf1 <- numeric(ndim)
#   for( i in 1:nobs){
#     R1 <- matrix(DCC[i,], ndim, ndim)
#     invR1 <- solve(R1)
#     # lf1[i] <- -0.5*(log(det(R1)) - sum(z[i,]*crossprod(invR1, z[i,])))  
#     lf1[i] <- -0.5*(log(det(R1)) - sum(z[i,]%*%invR1%*%z[i,]))  
#   }
#     sum(lf1) + 0.5*sum(z^2)
}
