####################################
loglik.ccc.t <- function(param, data, model, mode){
  if(is.zoo(data)) data <- as.matrix(data)
  nobs <- dim(data)[1]
  ndim <- dim(data)[2]
  mu <- matrix(param[1:ndim], nobs, ndim, byrow = TRUE)    # constant in the mean
  data <- data - mu
  param <- param[-(1:ndim)]
  
  para.mat <- p.mat(param, model, ndim)  
  # check if R is positive definite. If not, making R positive definite
  R <- make.pd(para.mat$R)

  h <- vgarch(para.mat$a, para.mat$A, para.mat$B, data)    # a call to vgarch function
  z <- data/sqrt(h)
  invR <- solve(R)
  lf <- -0.5*ndim*log(2*pi) - 0.5*rowSums(log(h)) - 0.5*log(det(R)) - 0.5*rowSums((z%*%invR)*z)
  
  if(mode=="gradient"){
    lf  
  } else {
    sum(lf)
  }
}

loglik.ccc <- function(param, data, model){
  if(is.zoo(data)) data <- as.matrix(data)
  nobs <- dim(data)[1]
  ndim <- dim(data)[2]
  mu <- matrix(param[1:ndim], nobs, ndim, byrow = TRUE)    # constant in the mean
  data <- data - mu
  param <- param[-(1:ndim)]
  
  para.mat <- p.mat(param, model, ndim)  
  # check if R is positive definite. If not, making R positive definite
  R <- make.pd(para.mat$R)
  
  h <- vgarch(para.mat$a, para.mat$A, para.mat$B, data)    # a call to vgarch function
  z <- data/sqrt(h)
  invR <- solve(R)
  lf <- -0.5*nobs*ndim*log(2*pi) - 0.5*sum(log(h)) - 0.5*nobs*log(det(R)) - 0.5*sum((z%*%invR)*z)
  -lf
}

# making a symmetric matrix positive definite. From r-help 2003.12.27
# except normalization.
make.pd <- function(x, tol=1e-6) {
  eig <- eigen(x, symmetric=TRUE)
  rtol <- tol * eig$values[1]
  if(min(eig$values) < rtol) {
    vals <- eig$values
    vals[vals < rtol] <- rtol
    srev <- eig$vectors %*% (vals * t(eig$vectors))
    dimnames(srev) <- dimnames(x)
    srev <- diag(1/sqrt(diag(srev)))%*%srev%*%diag(1/sqrt(diag(srev)))
    return(srev) 
  } else {return(x)}
}

# Inequality constraints for the stationarity of the vector GARCH equation in the (E)CCC GARCH
inEQccc <- function(param, data, model){
  ndim <- ncol(data)
  param <- param[-(1:ndim)]   # removing constants in the mean 
  pmat <- p.mat(param, model, ndim)
  ret <- max(Mod(eigen(pmat$A + pmat$B)$values))     # common to diagonal/extended
  ret
}

#*****************************************************************************************************************
estimateCCC <- function(inia = NULL, iniA = NULL, iniB = NULL, data, model="diagonal", ...){
  nobs <- dim(data)[1]
  ndim <- dim(data)[2]
  
  if(is.zoo(data)){
    d.ind <- index(data)
  } else {
    d.ind <- 1:nobs
  }
  
  data <- as.matrix(data)
  mu <- colMeans(data)
  R <- cor(data)
  
  # setting upper and lower bounds for the constraints
  if(model == "diagonal"){
    LB <- c(rep(-Inf, ndim), rep(0, 3*ndim), rep(-1, ndim*(ndim-1)/2))
  } else {
    LB <- c(rep(-Inf, ndim), rep(0, ndim + 2*ndim^2), rep(-1, ndim*(ndim-1)/2))
  }
  ineqLB <- 0
  ineqUB <- 1
  
  if(!is.null(inia) & !is.null(iniA) & !is.null(iniB)){ # when the initial values are supplied
    if(model=="diagonal"){
      init <- c(mu, inia, diag(iniA), diag(iniB), R[lower.tri(R)])
    } else {
      init <- c(mu, inia, as.vector(iniA), as.vector(iniB), R[lower.tri(R)])
    }
    tryCatch(
      suppressWarnings(
        results <- solnp(pars = init, fun = loglik.ccc,  
                         ineqfun = inEQccc, ineqLB = ineqLB, ineqUB = ineqUB,
                         LB = LB, 
                         control = list(trace=0), 
                         data = data, model = model                   
        )
      ), 
      error = function(e) conditionMessage(e), 
      finally=cat("Bad initial values. Try again without initial values.")
    )
  } else { # when initial values are not supplied
    cat("Initial values are not supplied. Random values are used.")
    # first stage optimisation
    conv <- 1
    ntry <- 0
    inia <- diag(cov(data))     # initial values for the constants in the GARCH part
    while(conv != 0){ # repeating the first stage optimization until it gives a successful convergence
      if(model != "diagonal"){
        ret <- 2
        while(ret > 1){
          ub <- runif(ndim, min=0.0001, max=0.04)
          iniA <- matrix(runif(ndim^2, min=0, max=ub[sample(1:ndim, 1)]), ndim, ndim)
          iniB <- matrix(runif(ndim^2, min=-0.04, max=ub[sample(1:ndim, 1)]), ndim, ndim)
          diag(iniA) <- round(runif(ndim, min=0.04, max=0.05), 4)
          diag(iniB) <- round(runif(ndim, min=0.8, max=0.9), 4)
          init <- c(mu, inia, as.vector(iniA), as.vector(iniB), R[lower.tri(R)])
          ret <- max(Mod(eigen(iniA + iniB)$values))     # common to diagonal/extended
        }
      } else {
        iniA <- diag(round(runif(ndim, min=0.04, max=0.05), 4))
        iniB <- diag(round(runif(ndim, min=0.8, max=0.9), 4))
        init <- c(mu, inia, diag(iniA), diag(iniB), R[lower.tri(R)])
      }
      suppressWarnings(
        results <- solnp(pars = init, fun = loglik.ccc,  
                         ineqfun = inEQccc, ineqLB = ineqLB, ineqUB = ineqUB,
                         LB = LB, 
                         control = list(trace=0), 
                         data = data, model = model                   
        )
      )
      conv <- results$convergence
      ntry <- ntry + 1
    }
  }
  
  mu <- matrix(results$pars[1:ndim], nobs, ndim, byrow = TRUE)
  eps <- data - mu
  
  param <- results$pars[-(1:ndim)]
  estimates <- p.mat(param, model=model, ndim=ndim)
  estimates$R <- make.pd(estimates$R)
  
  # computing conditional variance and std. residuals
  h <- vgarch(estimates$a, estimates$A, estimates$B, eps)    # estimating conditional variances
  std.resid <- eps/sqrt(h)                    # std. residuals
  
  # A character vector/matrix for naming parameters
  name.id <- as.character(1:ndim)
  
  if(is.null(colnames(data))){
    colnames(std.resid) <- paste("Series", name.id, sep="")
    colnames(h) <- paste("Series", name.id, sep="")
    colnames(data) <- paste("Series", name.id, sep="")
  } else {
    colnames(std.resid) <- colnames(data)
    colnames(h) <- colnames(data)
  }
  
  output <- list(
    results = results,
    model = model,
    method = "SQP by Rsolnp package",
    initial = list(a=inia, A=iniA, B=iniB),
    data = zoo(data, d.ind),
    estimates = estimates, 
    h = zoo(h, d.ind),              # conditional variances (not volatility)
    z = zoo(std.resid, d.ind)       # standardized residuals
  )
  
  class(output) <- "ccc"
  return(output)
}

# functions for summarizing output
summary.ccc <- function(object, ...){
  cat("Summarizing outcomes. This takes a while.")
  ndim <- ncol(object$data)
  nobs <- nrow(object$data)
  object$nobs <- nobs
  
  param <- object$results$pars
  object$mu <- param[1:ndim]             # mean estimates (constant)
  names(object$mu) <- paste("mu", 1:ndim, sep="")
  param <- param[-(1:ndim)]
  
  # re-arranging parameter vector into a list with paramerer matrices
  para.mat <- p.mat(param, object$model, ndim)
  para.mat$R <- make.pd(para.mat$R)
  
  # A character vector/matrix for naming parameters
  name.id <- as.character(1:ndim)
  namev <- diag(0, ndim, ndim)
  for(i in 1:ndim){
    for(j in 1:ndim){
      namev[i, j] <- paste(name.id[i], name.id[j], sep="")
    }
  }
  # naming parameters
  if(object$model=="diagonal"){
    vecA <- diag(para.mat$A)
    vecB <- diag(para.mat$B)
    names(para.mat$a) <- paste("a", 1:ndim, sep="")
    names(vecA) <- paste("A", diag(namev), sep="")
    names(vecB) <- paste("B", diag(namev), sep="")
  } else {
    vecA <- as.vector(para.mat$A)
    vecB <- as.vector(para.mat$B)
    names(para.mat$a) <-paste("a", 1:ndim, sep="")
    names(vecA) <-paste("A", namev, sep="")
    names(vecB) <-paste("B", namev, sep="")
  }
  object$R <- para.mat$R[lower.tri(para.mat$R)]
  names(object$R) <- paste("R", namev[lower.tri(namev)], sep="")
  object$garch.par <- c(para.mat$a, vecA, vecB)  # estimates for conditional variance
  
  # computing Jacobian and Hessian for the mean and GARCH part
  all.par <- c(object$mu, object$garch.par, object$R)
  npar <- length(all.par)       # the number of total parameters
  ja <- jacobian(func=loglik.ccc.t, x=all.par, 
                 data=object$data, model=object$model, mode="gradient")  # using jacobian() in numDeriv
  J <- crossprod(ja)  # information matrix
  H <- hessian(func = loglik.ccc.t, x=all.par, 
               data=object$data, model=object$model, mode="hessian", method.args=list(eps=1e-5, d=1e-5))
  invH <- solve(H)
  S <- invH%*%J%*%invH
  se <- sqrt(diag(S))     
  
  coef <- c(object$mu, object$garch.par, object$R)
  zstat <- all.par/se
  pval <- round(2*pnorm(-abs(zstat)), 6)
  coef <- cbind(all.par, se, zstat, pval)
  colnames(coef) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  
  object$coef <- coef
  
  object$convergence <- object$results$convergence    # convergence status
  names(object$convergence) <- "1st"
  object$counts <- c(object$results$outer.iter, object$results$nfuneval)   # the number of iterations
  
  object$logLik <- -tail(object$results$values, 1)  # the value of the log-like at the estimate
  
  # Information Criteria
  object$AIC <- -2*object$logLik + 2*npar
  object$BIC <- -2*object$logLik + log(nobs)*npar
  object$CAIC <- -2*object$logLik + (1+log(nobs))*npar
  
  class(object) <- "summary.ccc"
  return(object)
}

logLik.ccc <- function(object, ...){
  LL <- -tail(object$results$values, 1)
  cat("Log-likelihood at the estimates: ", formatC(LL, digits = 10), sep = "\n")
}

print.summary.ccc <- function(x, digits = max(3, getOption("digits") - 1), ...){
  cat("Constant Conditional Correlation GARCH Model", "\n")
  cat("Conditional variance equation:", x$model, "\n")
  cat("\nCoefficients:", "\n")
  printCoefmat(x$coef, digits = 4, dig.tst = 4)    # printCoefmat() is defined in stats package
  
  cat("\nNumber of Obs.:", formatC(x$nobs, digits = 0), "\n")
  cat("Log-likelihood:", formatC(x$logLik, format="f", digits = 5), "\n\n")
  
  cat("Information Criteria:", "\n")
  cat("  AIC:", formatC(x$AIC, format="f", digits = 3), "\n")
  cat("  BIC:", formatC(x$BIC, format="f", digits = 3), "\n")
  cat(" CAIC:", formatC(x$CAIC, format="f", digits = 3), "\n")
  
  cat("\nOptimization method:", x$method, "\n")
  cat("Convergence :", x$convergence, "\n")
  cat("Iterations  :", "\n")
  cat("  Major iterations :", x$counts[1], "\n")
  cat(" Func. evaluations :", x$counts[2], "\n")
  
  invisible(x)
}

print.ccc <- function(x, ...){
  cat("Estimated model:", x$model, "\n")
  cat("Use summary() to see the estimates and related statistics", "\n")
  invisible(x)
}

