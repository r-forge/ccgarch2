# computing dcc part of the log-likelihood. 
# outout is a vector of length nobs
dcc.ll2 <- function(R, z){
    if(is.zoo(R)) R <- as.matrix(R)
    if(is.zoo(z)) z <- as.matrix(z)
    .Call("dcc_ll2", R, z)
}

####################################
# vector GARCH equatio: this function is common to all GARCH models
vgarch <- function(a, A, B, data){
    if(is.zoo(data)) data <- as.matrix(data)
    dvar <- data^2           # dvar = eps
     .Call("vector_garch", dvar, a, A, B)
}

# vector GARCH equation with GJR-type asymmetry
vgarch.l <- function(a, A, B, Lev, data){
    if(is.zoo(data)) data <- as.matrix(data)
    z <- data           # a vector of residuals 
  z2 <- z^2           # a vector of squared residuals
  nobs <- nrow(data)   # number of observations
  ndim <- ncol(data)   # number of dimension
  
  # leverage term
  sign.z <- -sign(z)        # equal to 1 if z < 0
  sign.z[sign.z < 0] <- 0   # equal to 0 if z < 1
  lev.z <- sign.z*abs(z)    # lev.z = |z| if z < 0, = 0 otherwise

  if(length(a) != ndim)   stop("a is not")
  ht <- matrix(0, nobs, ndim)   
  l.h <- colMeans(z2)         # initial value
  l.z2 <- colMeans(z2)        # initial value
  l.lev.z <- numeric(ndim)
  for(i in 1:nobs){
    ht[i, ] <- l.h <- a + A%*%l.z2 + B%*%l.h + Lev%*%l.lev.z
    l.z2 <- z2[i, ]
    l.lev.z <- as.vector(lev.z[i, ])
  }
  ht
}

####################################
# rearranging a vector of parameters into parameter matrices in the vector GARCH (1, 1) 
# and constant conditional correlation matrix. This is common to all GARCH models in 
# the package.
p.mat <- function(para, model, ndim){
    npara <- length(para)
   if(model=="diagonal"){                          # for the diagonal vector GARCH equation
      a <- para[1:ndim]                          # constant in variance
      A <- diag(para[(ndim+1):(2*ndim)])         # ARCH parameter
      B <- diag(para[(2*ndim+1):(3*ndim)])       # GARCH parameter
      R <- diag(ndim)                              # Constant Conditional Correlation Matrix
      R[lower.tri(R)] <- para[(3*ndim+1):npara]; R <- (R+t(R)); diag(R) <- 0.5*diag(R)
   } else if(model=="extended"){   # for the extended vector GARCH equation
      a <- para[1:ndim]
      A <- matrix(para[(ndim+1):(ndim^2+ndim)], ndim, ndim)
      B <- matrix(para[(ndim^2+ndim+1):(2*ndim^2+ndim)], ndim, ndim)
      R <- diag(ndim)
      R[lower.tri(R)] <- para[(2*ndim^2+ndim+1):npara]; R <- (R+t(R)); diag(R) <- 0.5*diag(R)
   }
   list(a=a, A=A, B=B, R=R)
}

####################################
# Inequality constraint for the stationarity of the vector GARCH equation for the (E)DCC GARCH
inEQdcc1 <- function(pars, data, model){
    ndim <- ncol(data)
    pars <- pars[-(1:ndim)]   # removing constants in the mean 
    In <- diag(ndim)
    pars <- c(pars, In[lower.tri(In)])
    pmat <- p.mat(pars, model, ndim)
    ret <- max(Mod(eigen(pmat$A + pmat$B)$values))     # common to diagonal/extended
#    para.mat <- p.mat(pars, model, ndim)
#        if(model == "diagonal"){    # for the diagonal model
#            ret <- diag(para.mat$A + para.mat$B)
#        } else {    # for the extended model
#            ret <- max(eigen(para.mat$A + para.mat$B)$values)
#        }
   if(model=="diagonal"){
    ret <- c(ret, diag(pmat$A))
   } else {
    ret <- c(ret, as.vector(pmat$A))
   }
    return(ret)
}

inEQdcc2 <- function(param, data){
    param[1] + param[2]
}



####################################
# computing a likelihood for the 1st stage DCC/cDCC: this function is common to all GARCH models
loglik1.dcc1 <- function(param, data, model){
    if(is.zoo(data)) data <- as.matrix(data)
    nobs <- dim(data)[1]
   ndim <- dim(data)[2]
   In <- diag(ndim)
   mu <- matrix(param[1:ndim], nobs, ndim, byrow = TRUE)    # constant in the mean
   data <- data - mu
   param <- param[-(1:ndim)]
   param <- c(param, In[lower.tri(In)])
   para.mat <- p.mat(param, model, ndim)
   h <- vgarch(para.mat$a, para.mat$A, para.mat$B, data)    # a call to vgarch function
   z <- data/sqrt(h)
   lf <- -0.5*nobs*ndim*log(2*pi) - 0.5*sum(log(h)) - 0.5*sum(z^2)
   
   -lf
}

####################################
# the objective function in the 2nd stage DCC estimation. This is to be minimised!
loglik2.cdcc <- function(param, data){        # data is the standardised residuals
    if(is.zoo(data)) data <- as.matrix(data)
    nobs <- dim(data)[1]
    ndim <- dim(data)[2]
    cDCC <- cdcc.est(data, param)$cDCC
    lf1 <- dcc.ll2(cDCC, data)    
    
#   lf1 <- numeric(ndim)
#   for( i in 1:nobs){
#     R1 <- matrix(cDCC[i,], ndim, ndim)
#     invR1 <- solve(R1)
#     lf1[i] <- 0.5*(log(det(R1)) + sum(data[i,]*crossprod(invR1, data[i,])))    # negative of the log-likelihood function
#   }
#   sum(lf1) - 0.5*sum(data^2)  # the second term is unrelated with the optimization, but is included for computing log-lik value
    sum(lf1)
}

####################################
# computing cDCC with restriction on the diagonal elements to be one
cdcc.est <- function(data, param){
  if(is.zoo(data)) data <- as.matrix(data)
   nobs <- nrow(data)
#   uncR <- cov(data)*((nobs-1)/nobs)
   uncR <- cor(data)
   out <- .Call("cdcc_est", data, uncR, param[1], param[2])
   list(cDCC=out[[1]], Q=out[[2]])
}

#*****************************************************************************************************************
# The 1st stage DCC estimation: this function is common to all GARCH models
dcc1.estimation <- function(data, a, A, B, model){
    if(is.zoo(data)) data <- as.matrix(data)
    nobs <- dim(data)[1]
   ndim <- dim(data)[2]
   mu <- colMeans(data)
   if(model=="diagonal"){
    init <- c(a, diag(A), diag(B))
   } else {
    init <- c(a, as.vector(A), as.vector(B))
   }
    init <- c(mu, init)     # adding constant in the mean to the initial value
    npar <- length(init)    # the number of parameters

   # setting upper and lower bounds for the constraints
   if(model == "diagonal"){
       LB <- c(rep(-500, ndim), rep(0, 3*ndim))
       UB <- c(rep(500, ndim), rep(1, 3*ndim))
       ineqLB <- c(0, rep(0, ndim))
       ineqUB <- c(1, rep(1, ndim))
   } else {
       LB <- c(rep(-500, ndim), rep(0, ndim + 2*ndim^2))
       UB <- c(rep(500, ndim), rep(1, ndim + 2*ndim^2))
       ineqLB <- c(0, rep(0, ndim^2))
       ineqUB <- c(1, rep(1, ndim^2))
   }

   # the first stage optimization  
   suppressWarnings(
       step1 <- solnp(pars = init, fun = loglik1.dcc1,  
                      ineqfun = inEQdcc1, ineqLB = ineqLB, ineqUB = ineqUB,
                      LB = LB, #UB = UB, 
                      control = list(trace=0, tol = 1e-9, delta = 1e-10, rho = 2.5),
                      data = data, model = model
       )
   )
#     if(random){
#         tmp.init <- init[-(1:(2*ndim))]
#         distr <- c(rep(2, 2*ndim), rep(3, length(tmp.init)))
#         # distr <- rep(3, length(init))
#         distr.opt = vector(mode = "list", length = length(init))
#         for(i in 1:length(init)){
#             distr.opt[[i]]$mean = init[i]
#             distr.opt[[i]]$sd = sqrt(init[i]^2)*2 + 1e-2
#         }
#         # the first stage optimization  
#         suppressWarnings(
#             step1 <- gosolnp(pars = NULL, fun = loglik1.dcc1,  
#                              ineqfun = inEQdcc1, ineqLB = ineqLB, ineqUB = ineqUB,
#                              LB = LB, UB = UB, 
#                              distr = distr, 
#                              distr.opt = distr.opt,
#                              n.sim = 300000,
#                              control = list(trace=0, tol = 1e-9, delta = 1e-10),
#                              data = data, model = model
#             )
#         )
#     } else {
#         # the first stage optimization  
#         suppressWarnings(
#             step1 <- solnp(pars = init, fun = loglik1.dcc1,  
#                            ineqfun = inEQdcc1, ineqLB = ineqLB, ineqUB = ineqUB,
#                            LB = LB, #UB = UB, 
#                            #distr = distr, 
#                            #distr.opt = distr.opt,
#                            # n.restarts = 2,
#                            #n.sim = 300000,
#                            control = list(trace=0, tol = 1e-9, delta = 1e-10, rho = 2.5),
#                            data = data, model = model
#             )
#         )
#     }
    #   step1 <- optim(par=init, fn=loglik1.dcc1, method=method, control=list(maxit=10^5, ndeps=rep(1e-7, npar), reltol=1e-15), data=data, model=model)
   step1
}


#*****************************************************************************************************************
# The 2nd stage cDCC estimation.
cdcc.2stg <- function(data, para){ # data must be standardised residuals

    LB <- rep(0, 2)     # lower bound of the parameter
    ineqLB <- 0         # lower bound of the stationarity
    ineqUB <- 1         # upper bound of the stationarity
    
   suppressWarnings(
      step2 <- solnp(pars=para, fun = loglik2.cdcc,
                     ineqfun = inEQdcc2, ineqLB = ineqLB, ineqUB = ineqUB,
                     LB = LB, 
                     control = list(trace=0), 
                     data = data
                     )
   )
#    step2 <- constrOptim(theta=para, f=loglik2.cdcc, grad = NULL, ui=resta, ci=restb, mu=1e-5, control=list(maxit=10^5, ndeps=rep(1e-7, 2), reltol=1e-15), data=data)
   step2
}

#*****************************************************************************************************************
estimateCDCC <- function(inia = NULL, iniA = NULL, iniB = NULL, ini.dcc = NULL, 
                         data, model="diagonal", ...){
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
    
   mu <- matrix(first.stage$pars[1:ndim], nobs, ndim, byrow = TRUE)
   eps <- data - mu

   tmp.para <- c(first.stage$pars[-(1:ndim)], In[lower.tri(In)])
   estimates <- p.mat(tmp.para, model=model, ndim=ndim)

   h <- vgarch(estimates$a, estimates$A, estimates$B, eps)    # estimated conditional variances
   std.resid <- eps/sqrt(h)                    # std. residuals

   # second stage optimisation
  if(is.null(ini.dcc)) ini.dcc <- c(0.1, 0.8)
  second.stage <- cdcc.2stg(std.resid, ini.dcc)
  cdcc <- cdcc.est(std.resid, second.stage$pars)  # Q and cDCC estimates

  # A character vector for naming correlations
  name.id <- as.character(1:ndim)
  namev <- diag(0, ndim, ndim)
  for(i in 1:ndim){
    for(j in 1:ndim){
     namev[i, j] <- paste(name.id[i], name.id[j], sep="")
    }
  }
  colnames(cdcc$cDCC) <- paste("R", namev, sep="")     # column names of the DCC matrix

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
    CDCC = zoo(cdcc$cDCC, d.ind),   # conditional correlations
    h = zoo(h, d.ind),              # conditional variances (not volatility)
    z = zoo(std.resid, d.ind)       # standardized residuals
  )
  
  class(output) <- "cdcc"
  return(output)
}

# functions for summarizing output
summary.cdcc <- function(object, ...){
  cat("Summarizing outcomes. This takes a while.")
  ndim <- ncol(object$data)
  nobs <- nrow(object$data)
  object$nobs <- nobs
  In <- diag(ndim)
  
  object$mu <- object$f.stage$pars[1:ndim]             # mean estimates (constant)
  names(object$mu) <- paste("mu", 1:ndim, sep="")

  object$garch.par <- object$f.stage$pars[-(1:ndim)]   # parameter estimates for the Vector GARCH
  # re-arranging parameter vector into a list with paramerer matrices
  para.mat <- p.mat(c(object$garch.par, In[lower.tri(In)]), object$model, ndim)

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
  
  # computing standard errors for the GARCH part
  ja <- jacobian(func=loglik1.dcc.t, x=c(object$mu, object$garch.par), 
                    data=object$data, model=object$model, mode="gradient")  # using jacobian() in numDeriv
  # H <- optimHess(c(object$mu, object$garch.par), fn=loglik1.dcc.t, data=object$data, model=object$model, mode="hessian")
  H <- hessian(func = loglik1.dcc.t, x=c(object$mu, object$garch.par), 
                data=object$data, model=object$model, mode="hessian", method.args=list(eps=1e-5, d=1e-5))
  J <- crossprod(ja)    # information matrix

  # computing standard errors for the cDCC part
  object$cDcc.par <- object$s.stage$pars             # estimates for the cDCC
  names(object$cDcc.par) <- c("alpha", "beta")

  ja.cdcc <- jacobian(func=loglik2.cdcc.t, x=object$cDcc.par, 
                        data=object$z, mode="gradient")
  # H.cdcc <- optimHess(object$cDcc.par, fn=loglik2.cdcc.t, data=object$z, mode="hessian")
  H.cdcc <- hessian(func = loglik2.cdcc.t, x = object$cDcc.par, data=object$z, 
                    mode="hessian", method.args=list(eps=1e-5, d=1e-5))
  J.cdcc <- crossprod(ja.cdcc)
  
  cross <- crossprod(ja, ja.cdcc)      # npar.garch x 2
  Omega <- rbind(cbind(J, cross), cbind(t(cross), J.cdcc))
  
  all.par <- c(object$mu, object$garch.par, object$cDcc.par)
  # g.cdcc <- optimHess(all.par, fn=cdcc.hessian, data=object$data, model=object$model)  # this is time consuming!!!
  g.cdcc <- hessian(func = cdcc.hessian, x = all.par, data=object$data, 
                    model=object$model, method.args=list(eps=1e-5, d=1e-5))
  npar <- length(all.par)
  
  g.cdcc <- g.cdcc[(npar-1):npar, 1:(npar-2)]   # the last 2 rows and the first (npar-2) columns
  
  G <- rbind(cbind(H, diag(0, nrow(H), 2)), cbind(g.cdcc, H.cdcc))
  invG <- solve(G)
  invGt <- t(invG)
  
  S.all <- invG%*%Omega%*%invGt
  se.all <- sqrt(diag(S.all))     

  coef <- c(object$mu, object$garch.par, object$cDcc.par)
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
  
  object$logLik <- -(tail(object$f.stage$values, 1,) + tail(object$s.stage$values, 1))  # the value of the log-like at the estimate
  
  # Information Criteria
  object$AIC <- -2*object$logLik + 2*npar
  object$BIC <- -2*object$logLik + log(nobs)*npar
  object$CAIC <- -2*object$logLik + (1+log(nobs))*npar
  
  class(object) <- "summary.cdcc"
  return(object)
}

logLik.cdcc <- function(object, ...){
  LL <- tail(object$f.stage$values, 1,) + tail(object$s.stage$values, 1)
  cat("Log-likelihood at the estimates: ", formatC(-LL, digits = 10), sep = "\n")
}

# print.summary.cdcc <- function(x, digits = max(3, getOption("digits") - 1), ...){
print.summary.cdcc <- function(x, digits = max(3, getOption("digits") - 1), ...){
        cat("Corrected Dynamic Conditional Correlation GARCH Model", "\n")
  cat("Conditional variance equation:", x$model, "\n")
  
  cat("\nCoefficients:", "\n")
  printCoefmat(x$coef, digits = 2, dig.tst = 4)    # printCoefmat() is defined in stats package
  
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

print.cdcc <- function(x, ...){
  cat("Estimated model:", x$model, "\n")
  cat("Use summary() to see the estimates and related statistics", "\n")
  #print.default(format(x$f.stage$par, digits = 4), print.gap = 1, quote = FALSE)
  invisible(x)
}

