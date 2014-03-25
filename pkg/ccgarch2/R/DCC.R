####################################
# the objective function in the 2nd stage DCC estimation. This is to be minimised!
loglik2.dcc <- function(param, data){        # data is the standardised residuals
  if(is.zoo(data)) data <- as.matrix(data)
  nobs <- dim(data)[1]
  ndim <- dim(data)[2]
  DCC <- dcc.est(data, param)$DCC
  lf1 <- dcc.ll2(DCC, data)    
#   lf1 <- numeric(ndim)
#   for( i in 1:nobs){
#     R1 <- matrix(DCC[i,], ndim, ndim)
#     invR1 <- solve(R1)
#     lf1[i] <- 0.5*(log(det(R1)) + sum(data[i,]*crossprod(invR1, data[i,])))  # negative of the log-likelihood function
#   }
#   sum(lf1) - 0.5*sum(data^2)  # the second term is unrelated with the optimization, but is included for computing log-lik value
    sum(lf1)
}

####################################
# computing DCC 
dcc.est <- function(data, param){
  if(is.zoo(data)) data <- as.matrix(data)
   uncR <- cov(data)
   out <- .Call("dcc_est", data, uncR, param[1], param[2])
   list(DCC=out[[1]], Q=out[[2]])
#  nobs <- nrow(data)
#  ndim <- ncol(data)
#  uncR1 <- cov(data)*((nobs-1)/nobs)    
#  Q1 <- uncR1         # initial value of Q0
###########################################################  
#  zz1 <- colMeans(data)        # initial value of zz0
#  zz1 <- zz1%o%zz1
###########################################################
#  
#  a <- param[1]
#  b <- param[2]
#  const1 <- (1-a-b)*uncR1
#  
#  DCC <- diag(0, nobs, ndim^2)
#  vecQ1 <- diag(0, nobs, ndim^2)
#  for(i in 1:nobs){
#    Q1 <- const1 + a*zz1 + b*Q1
#    invdQ1 <- diag(1/sqrt(diag(Q1)))
#    D1 <- invdQ1%*%Q1%*%invdQ1
#    DCC[i, ] <- as.vector(D1)
#    vecQ1[i, ] <- as.vector(Q1)
#    zz1 <- data[i, ]%o%data[i, ]
#  }
#  list(DCC=DCC, Q=vecQ1)
}

#*****************************************************************************************************************
# The 2nd stage DCC estimation.
dcc.2stg <- function(data, para){ # data must be standardised residuals
    LB <- rep(0, 2)     # lower bound of the parameter
    ineqLB <- 0         # lower bound of the stationarity
    ineqUB <- 1         # upper bound of the stationarity
    
   suppressWarnings(
    step2 <- solnp(pars=para, fun = loglik2.dcc,
                   ineqfun = inEQdcc2, ineqLB = ineqLB, ineqUB = ineqUB,
                   LB = LB, 
                   control = list(trace=0, tol = 1e-9, delta = 1e-10),
                   data = data
                   )
   )
    step2
}

#*****************************************************************************************************************
estimateDCC <- function(inia = NULL, iniA = NULL, iniB = NULL, ini.dcc = NULL, data, model="diagonal", ...){
   nobs <- dim(data)[1]
   ndim <- dim(data)[2]
   In <- diag(ndim)
   
   if(is.zoo(data)){
       d.ind <- index(data)
   } else {
       d.ind <- 1:nobs
   }
   data <- as.matrix(data)

   if(!is.null(inia) & !is.null(iniA) & !is.null(iniB)){ # when the initial values are supplied
     tryCatch(first.stage <- dcc1.estimation(data=data, a=inia, A=iniA, B=iniB, 
                                             model=model), 
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
           iniB <- matrix(runif(ndim^2, min=-0.004, max=ub[sample(1:ndim, 1)]), ndim, ndim)
           diag(iniA) <- round(runif(ndim, min=0.04, max=0.05), 4)
           diag(iniB) <- round(runif(ndim, min=0.8, max=0.9), 4)
           ret <- max(Mod(eigen(iniA + iniB)$values))     # common to diagonal/extended
         }
       } else {
         iniA <- diag(round(runif(ndim, min=0.04, max=0.05), 4))
         iniB <- diag(round(runif(ndim, min=0.8, max=0.9), 4))
       }
       first.stage <- dcc1.estimation(data=data, a=inia, A=iniA, B=iniB, model=model)
       conv <- first.stage$convergence
       ntry <- ntry + 1
     }
   }
   
   mu <- matrix(first.stage$par[1:ndim], nobs, ndim, byrow = TRUE)
   eps <- data - mu

  # re-aranging parameter estimates
   tmp.para <- c(first.stage$par[-(1:ndim)], In[lower.tri(In)])
   estimates <- p.mat(tmp.para, model=model, ndim=ndim)

  # computing conditional variance and std. residuals
   h <- vgarch(estimates$a, estimates$A, estimates$B, eps)    # estimating conditional variances
   std.resid <- eps/sqrt(h)                    # std. residuals

   # second stage optimisation
  if(is.null(ini.dcc)) ini.dcc <- c(0.1, 0.8)
  second.stage <- dcc.2stg(std.resid, ini.dcc)
   dcc <- dcc.est(std.resid, second.stage$par)  # Q and cDCC estimates

  # A character vector/matrix for naming parameters
  name.id <- as.character(1:ndim)
  namev <- diag(0, ndim, ndim)
  for(i in 1:ndim){
    for(j in 1:ndim){
     namev[i, j] <- paste(name.id[i], name.id[j], sep="")
    }
  }
  colnames(dcc$DCC) <- paste("R", namev, sep="")     # column names of the DCC matrix
  
  if(is.null(colnames(data))){
    colnames(std.resid) <- paste("Series", name.id, sep="")
    colnames(h) <- paste("Series", name.id, sep="")
    colnames(data) <- paste("Series", name.id, sep="")
  } else {
    colnames(std.resid) <- colnames(data)
    colnames(h) <- colnames(data)
  }
  
  output <- list(
    f.stage = first.stage,
    s.stage = second.stage,
    model = model,
    method = "SQP by Rsolnp package",
    initial = list(a=inia, A=iniA, B=iniB, dcc.par=ini.dcc),
    data = zoo(data, d.ind),
    DCC = zoo(dcc$DCC, d.ind),   # conditional correlations
    h = zoo(h, d.ind),              # conditional variances (not volatility)
    z = zoo(std.resid, d.ind)       # standardized residuals
  )
  
  class(output) <- "dcc"
  return(output)
}

# functions for summarizing output
summary.dcc <- function(object, ...){
  cat("Summarizing outcomes. This takes a while.")
  ndim <- ncol(object$data)
  nobs <- nrow(object$data)
  object$nobs <- nobs
  In <- diag(ndim)
  
  object$mu <- object$f.stage$par[1:ndim]             # mean estimates (constant)
  names(object$mu) <- paste("mu", 1:ndim, sep="")

  garch.par <- object$f.stage$par[-(1:ndim)]   # parameter estimates for the Vector GARCH
  # re-arranging parameter vector into a list with paramerer matrices
  para.mat <- p.mat(c(garch.par, In[lower.tri(In)]), object$model, ndim)

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
  
  object$garch.par <- c(para.mat$a, vecA, vecB)  # estimates for conditional variance
  
  # computing Jacobian and Hessian for the mean and GARCH part
  ja <- jacobian(func=loglik1.dcc.t, x=c(object$mu, object$garch.par), 
                    data=object$data, model=object$model, mode="gradient")  # using jacobian() in numDeriv
  J <- crossprod(ja)  # information matrix
  # H <- optimHess(c(object$mu, object$garch.par), fn=loglik1.dcc.t, data=object$data, model=object$model, mode="hessian")
  H <- hessian(func = loglik1.dcc.t, x=c(object$mu, object$garch.par), 
                data=object$data, model=object$model, mode="hessian", method.args=list(eps=1e-5, d=1e-5))



  # computing standard errors for the DCC part
  object$Dcc.par <- object$s.stage$par             # estimates for the DCC
  names(object$Dcc.par) <- c("alpha", "beta")

  ja.dcc <- jacobian(func=loglik2.dcc.t, x=object$Dcc.par, 
                        data=object$z, mode="gradient")
  # H.dcc <- optimHess(object$Dcc.par, fn=loglik2.dcc.t, data=object$z, mode="hessian")
  H.dcc <- hessian(func = loglik2.dcc.t, x = object$Dcc.par, data=object$z, 
                    mode="hessian", method.args=list(eps=1e-5, d=1e-5))
  J.dcc <- crossprod(ja.dcc)
  
  cross <- crossprod(ja, ja.dcc)      # npar.garch x 2
  Omega <- rbind(cbind(J, cross), cbind(t(cross), J.dcc))
  
  all.par <- c(object$mu, object$garch.par, object$Dcc.par)
  # g.dcc <- optimHess(all.par, fn=dcc.hessian, data=object$data, model=object$model)  # this is time consuming!!!
  g.dcc <- hessian(func = dcc.hessian, x = all.par, data=object$data, 
                    model=object$model, method.args=list(eps=1e-5, d=1e-5))


  npar <- length(all.par)       # the number of total parameters
  
  g.dcc <- g.dcc[(npar-1):npar, 1:(npar-2)]   # the last 2 rows and the first (npar-2) columns
  G <- rbind(cbind(H, diag(0, nrow(H), 2)), cbind(g.dcc, H.dcc))
  invG <- solve(G)
  invGt <- t(invG)
  
  S.all <- invG%*%Omega%*%invGt
  se.all <- sqrt(diag(S.all))     

  coef <- c(object$mu, object$garch.par, object$Dcc.par)
  zstat <- coef/se.all
  pval <- round(2*pnorm(-abs(zstat)), 6)
  coef <- cbind(coef, se.all, zstat, pval)
  colnames(coef) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")

  object$coef <- coef
  
  object$convergence <- c(object$f.stage$convergence, object$s.stage$convergence)    # convergence status
  names(object$convergence) <- c("1st", "2nd")
  object$counts <- rbind(c(object$f.stage$outer.iter, object$f.stage$nfuneval),       # the number of iterations in the 1st step
                         c(object$s.stage$outer.iter, object$s.stage$nfuneval))           # the number of iterations in the 2nd step
  colnames(object$counts) <- c("Major iterations", "Func. evaluations")
  rownames(object$counts) <- c("1st", "2nd")
  
  object$logLik <- -(tail(object$f.stage$value, 1) + tail(object$s.stage$value, 1))  # the value of the log-like at the estimate

  # Information Criteria
  object$AIC <- -2*object$logLik + 2*npar
  object$BIC <- -2*object$logLik + log(nobs)*npar
  object$CAIC <- -2*object$logLik + (1+log(nobs))*npar
  
  class(object) <- "summary.dcc"
  return(object)
}

logLik.dcc <- function(object, ...){
  LL <- tail(object$f.stage$value, 1) + tail(object$s.stage$value, 1)
  cat("Log-likelihood at the estimates: ", formatC(-LL, digits = 10), sep = "\n")
}

print.summary.dcc <- function(x, digits = max(3, getOption("digits") - 1), ...){
  cat("Dynamic Conditional Correlation GARCH Model", "\n")
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
  cat("Convergence (1st, 2nd):", x$convergence, "\n")
  cat(  "Iterations (1st)      :", x$counts[1,], "\n")
  cat(  "Iterations (2nd)      :", x$counts[2,], "\n")

  invisible(x)
}

print.dcc <- function(x, ...){
  cat("Estimated model:", x$model, "\n")
  cat("Use summary() to see the estimates and related statistics", "\n")
  #print.default(format(x$f.stage$par, digits = 4), print.gap = 1, quote = FALSE)
  invisible(x)
}

