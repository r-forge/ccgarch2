###################################################################
# simulating data from DCC-GARCH(1, 1) DGP
simulateDCC <- function(a, b, Q, a0, A, B, nobs, ncut=1000){
  nobs <- nobs + ncut             # ncut is the number of observations to be removed
  ndim <- nrow(Q)
  
  const <- (1-a-b)*Q
  
  z <- diag(0, nobs, ndim)
  DCC <- diag(0, nobs, ndim^2)
  vecQ2 <- diag(0, nobs, ndim^2)
  Q2 <- Q         # initial value of Q0
  zz <- Q
  h <- diag(0, nobs, ndim)
  eps <- diag(0, nobs, ndim)
  ht <- rep(0, ndim)
  et2 <- rep(0, ndim)
  for(i in 1:nobs){
    Q2 <- const + a*zz + b*Q2
    invdQ2 <- diag(1/sqrt(diag(Q2)))
    D2 <- invdQ2%*%Q2%*%invdQ2
    DCC[i, ] <- as.vector(D2)
    vecQ2[i, ] <- as.vector(Q2)
    
    # zt <- mvrnorm(1, rep(0, ndim), D2)     # standardized residual sdampling from N(0, R_t) whereR_t = D2 is DCC at t
#############################################    
    cholD2 <- t(chol(D2))
    zt <- drop(cholD2%*%rnorm(ndim))
#############################################    
    zz <- zt%o%zt
    z[i, ] <- zt
    
    ht <- a0 + A%*%et2 + B%*%ht           # conditional variance at t
    h[i, ] <- drop(ht)
    eps[i, ] <- zt*sqrt(ht)               # simulated observation at t
    et2 <- eps[i, ]^2
  }
  
  # A character vector/matrix for naming parameters
  name.id <- as.character(1:ndim)
  namev <- diag(0, ndim, ndim)
  for(i in 1:ndim){
    for(j in 1:ndim){
     namev[i, j] <- paste(name.id[i], name.id[j], sep="")
    }
  }
  # naming rows
  rownames(DCC) <- c(rep(0, ncut), 1:(nobs-ncut))
  rownames(vecQ2) <- c(rep(0, ncut), 1:(nobs-ncut))
  rownames(z) <- c(rep(0, ncut), 1:(nobs-ncut))
  rownames(h) <- c(rep(0, ncut), 1:(nobs-ncut))
  rownames(eps) <- c(rep(0, ncut), 1:(nobs-ncut))
  # naming columns
  colnames(DCC) <- paste("R", namev, sep="")
  colnames(vecQ2) <- paste("R", namev, sep="")
  colnames(z) <- paste("Series", name.id, sep="")
  colnames(h) <- paste("Series", name.id, sep="")
  colnames(eps) <- paste("Series", name.id, sep="")
  
  list(DCC=DCC[-(1:ncut), ], z=z[-(1:ncut), ], Q=vecQ2[-(1:ncut), ], h=h[-(1:ncut), ], eps=eps[-(1:ncut), ])
}

###################################################################
# simulating data from DCC-GARCH(1, 1) DGP with GJR-type asymmetry
simulateDCC.lev <- function(a, b, Q, a0, A, B, Lev, nobs, ncut=1000){
  nobs <- nobs + ncut             # ncut is the number of observations to be removed
  ndim <- nrow(Q)
  
  const <- (1-a-b)*Q
  
  z <- diag(0, nobs, ndim)
  DCC <- diag(0, nobs, ndim^2)
  vecQ2 <- diag(0, nobs, ndim^2)
  Q2 <- Q         # initial value of Q0
  zz <- Q
  h <- diag(0, nobs, ndim)
  eps <- diag(0, nobs, ndim)
  ht <- rep(0, ndim)
  et2 <- rep(0, ndim)
  for(i in 1:nobs){
    Q2 <- const + a*zz + b*Q2
    invdQ2 <- diag(1/sqrt(diag(Q2)))
    D2 <- invdQ2%*%Q2%*%invdQ2
    DCC[i, ] <- as.vector(D2)
    vecQ2[i, ] <- as.vector(Q2)
    
    zt <- mvrnorm(1, rep(0, ndim), D2)     # sdampling from N(0m R_t) where R_t is DCC at t
    zz <- zt%o%zt
    z[i, ] <- zt
    
    ind <- (-sign(zt)+1)/2
    
    ht <- a0 + A%*%et2 + B%*%ht + Lev%*%(ind*et2)           # conditional variange at t
    h[i, ] <- drop(ht)
    eps[i, ] <- zt*sqrt(ht)               # simulated observation at t
    et2 <- eps[i, ]^2
  }
  # A character vector/matrix for naming parameters
  name.id <- as.character(1:ndim)
  namev <- diag(0, ndim, ndim)
  for(i in 1:ndim){
    for(j in 1:ndim){
     namev[i, j] <- paste(name.id[i], name.id[j], sep="")
    }
  }
  # naming rows
  rownames(DCC) <- c(rep(0, ncut), 1:(nobs-ncut))
  rownames(vecQ2) <- c(rep(0, ncut), 1:(nobs-ncut))
  rownames(z) <- c(rep(0, ncut), 1:(nobs-ncut))
  rownames(h) <- c(rep(0, ncut), 1:(nobs-ncut))
  rownames(eps) <- c(rep(0, ncut), 1:(nobs-ncut))
  # naming columns
  colnames(DCC) <- paste("R", namev, sep="")
  colnames(vecQ2) <- paste("R", namev, sep="")
  colnames(z) <- paste("Series", name.id, sep="")
  colnames(h) <- paste("Series", name.id, sep="")
  colnames(eps) <- paste("Series", name.id, sep="")
  
  list(DCC=DCC[-(1:ncut), ], z=z[-(1:ncut), ], Q=vecQ2[-(1:ncut), ], h=h[-(1:ncut), ], eps=eps[-(1:ncut), ])
}
