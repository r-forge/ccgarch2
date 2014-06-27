###################################################################
# simulating data from CCC-GARCH(1, 1) DGP
simulateCCC <- function(R, a0, A, B, nobs, ncut=1000){
  nobs <- nobs + ncut             # ncut is the number of observations to be removed
  ndim <- nrow(R)
  
  if(sum(diag(R)) != ndim)  stop("Correlation matrix must have ones on diagonal.")
  
  z <- matrix(rnorm(nobs*ndim), nobs, ndim)
  cholR <- chol(R)
  z <- z%*%cholR  
  
  h <- diag(0, nobs, ndim)
  eps <- diag(0, nobs, ndim)
  ht <- rep(0, ndim)
  et2 <- rep(0, ndim)
  for(i in 1:nobs){
    ht <- a0 + A%*%et2 + B%*%ht           # conditional variance at t
    h[i, ] <- drop(ht)
    eps[i, ] <- z[i, ]*sqrt(ht)               # simulated observation at t
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
  rownames(z) <- c(rep(0, ncut), 1:(nobs-ncut))
  rownames(h) <- c(rep(0, ncut), 1:(nobs-ncut))
  rownames(eps) <- c(rep(0, ncut), 1:(nobs-ncut))
  # naming columns
  colnames(z) <- paste("Series", name.id, sep="")
  colnames(h) <- paste("Series", name.id, sep="")
  colnames(eps) <- paste("Series", name.id, sep="")
  
  list(z=zoo(z[-(1:ncut), ]), h=zoo(h[-(1:ncut), ]), eps=zoo(eps[-(1:ncut), ]))
}

###################################################################
# simulating data from CCC-GARCH(1, 1) DGP with GJR-type asymmetry
simulateCCC.lev <- function(R, a0, A, B, Lev, nobs, ncut=1000){
  nobs <- nobs + ncut             # ncut is the number of observations to be removed
  ndim <- nrow(R)
  
  z <- matrix(rnorm(nobs*ndim), nobs, ndim)
  cholR <- chol(R)
  z <- z%*%cholR  

  ind <- (-sign(z)+1)/2
  ind <- rbind((-sign(colMeans(z))+1)/2, ind)
  
  h <- diag(0, nobs, ndim)
  eps <- diag(0, nobs, ndim)
  ht <- rep(0, ndim)
  et2 <- rep(0, ndim)
  for(i in 1:nobs){
    ht <- a0 + A%*%et2 + B%*%ht + Lev%*%(ind[i,]*et2)           # conditional variange at t
    h[i, ] <- drop(ht)
    eps[i, ] <- z[i,]*sqrt(ht)               # simulated observation at t
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
  rownames(z) <- c(rep(0, ncut), 1:(nobs-ncut))
  rownames(h) <- c(rep(0, ncut), 1:(nobs-ncut))
  rownames(eps) <- c(rep(0, ncut), 1:(nobs-ncut))
  # naming columns
  colnames(z) <- paste("Series", name.id, sep="")
  colnames(h) <- paste("Series", name.id, sep="")
  colnames(eps) <- paste("Series", name.id, sep="")
  
  list(z=zoo(z[-(1:ncut), ]), h=zoo(h[-(1:ncut), ]), eps=zoo(eps[-(1:ncut), ]))
}
