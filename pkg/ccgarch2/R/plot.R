#####################################
## plot method for the class "ccc" ##
#####################################
plot.summary.ccc <- function(x, item, main = NULL, sub = NULL, xlab = NULL, ylab = NULL, lwd = NULL, ...){
  plot.ccc(x, item, main = NULL, sub = NULL, xlab = NULL, ylab = NULL, lwd = NULL, ...)
}

plot.ccc <- function(x, item, main = NULL, sub = NULL, xlab = NULL, ylab = NULL, lwd = NULL, ...){
  # require(zoo)
  if(attr(x, "class")!="ccc" & attr(x, "class")!="summary.ccc") stop("wrong attribute")
  
  ndim <- ncol(x$data)            # the number of dimensions
  Ind <- diag(0, ndim); Ind[lower.tri(Ind)] <- 1

    # common arguments for plot
    if(is.null(sub)) sub <- ""
    if(is.null(xlab)) xlab <- "Time"
    if(is.null(lwd)) lwd <- 1

  if(item == "correlation"){
    stop("Conditional correlations are assumed constant")
  } else if(item == "volatility"){
    # obj <- zoo(sqrt(x$h), index(x$h))
    obj <- sqrt(x$h)
    
    # arguments for plot
    if(is.null(main)) main <- "Volatilities"
    if(is.null(ylab)) ylab <- colnames(x$h)

    plot.zoo(obj, main=main, sub=sub, xlab=xlab, ylab=ylab, lwd=lwd, ...)
  } else if(item == "std.residuals"){
    # obj <- zoo(x$z, index(x$z))
    obj <- x$z
    
    # arguments for plot
    if(is.null(main)) main <- "Standardized Residuals"
    if(is.null(ylab)) ylab <- colnames(x$h)

    plot(obj, main=main, sub=sub, xlab=xlab, ylab=ylab, lwd=lwd)
  } else if(item == "return"){
    # obj <- zoo(x$data, index(x$data))
    obj <- x$data
    
    # arguments for plot
    if(is.null(main)) main <- "Return Data"
    if(is.null(ylab)) ylab <- colnames(x$data)

    plot(obj, main=main, sub=sub, xlab=xlab, ylab=ylab, lwd=lwd )
  } else {
    stop("The argument `item' is one of 'volatility', 'std.residuals', 'return'")
  }
}
#####################################
## plot method for the class "dcc" ##
#####################################
plot.summary.dcc <- function(x, item, main = NULL, sub = NULL, xlab = NULL, ylab = NULL, lwd = NULL, ...){
  plot.dcc(x, item, main = NULL, sub = NULL, xlab = NULL, ylab = NULL, lwd = NULL, ...)
}

plot.dcc <- function(x, item, main = NULL, sub = NULL, xlab = NULL, ylab = NULL, lwd = NULL, ...){
  # require(zoo)
  if(attr(x, "class")!="dcc" & attr(x, "class")!="summary.dcc") stop("wrong attribute")
  
  ndim <- ncol(x$data)            # the number of dimensions
  Ind <- diag(0, ndim); Ind[lower.tri(Ind)] <- 1

    # common arguments for plot
    if(is.null(sub)) sub <- ""
    if(is.null(xlab)) xlab <- "Time"
    if(is.null(lwd)) lwd <- 1

  if(item == "correlation"){
    obj <- x$DCC
    obj <- obj[, as.vector(Ind == 1)]                       # remove duplicated columns
    obj <- zoo(obj, index(x$h))
    
    # arguments for plot
    if(is.null(main)) main <- "Conditional Correlations"
    if(is.null(ylab)) ylab <- colnames(obj)
    
    plot(obj, ylim=c(-1, 1), main=main, sub=sub, xlab=xlab, ylab=ylab, lwd=lwd, ...)
  } else if(item == "volatility"){
    obj <- zoo(sqrt(x$h), index(x$h))

    # arguments for plot
    if(is.null(main)) main <- "Volatilities"
    if(is.null(ylab)) ylab <- colnames(x$h)

    plot(obj, main=main, sub=sub, xlab=xlab, ylab=ylab, lwd=lwd, ...)
  } else if(item == "std.residuals"){
    obj <- zoo(x$z, index(x$z))

    # arguments for plot
    if(is.null(main)) main <- "Standardized Residuals"
    if(is.null(ylab)) ylab <- colnames(x$h)

    plot(obj, main=main, sub=sub, xlab=xlab, ylab=ylab, lwd=lwd, ...)
  } else if(item == "return"){
    obj <- zoo(x$data, index(x$data))

    # arguments for plot
    if(is.null(main)) main <- "Return Data"
    if(is.null(ylab)) ylab <- colnames(x$data)

    plot(obj, main=main, sub=sub, xlab=xlab, ylab=ylab, lwd=lwd, ...)
  } else {
    stop("The argument `item' is one of 'correlation', 'volatility', 'std.residuals', 'return'")
  }
}

######################################
## plot method for the class "cdcc" ##
######################################
plot.summary.cdcc <- function(x, item, main = NULL, sub = NULL, xlab = NULL, ylab = NULL, lwd = NULL, ...){
  plot.cdcc(x, item, main = NULL, sub = NULL, xlab = NULL, ylab = NULL, lwd = NULL, ...)
}

  
plot.cdcc <- function(x, item, main = NULL, sub = NULL, xlab = NULL, ylab = NULL, lwd = NULL, ...){
  # require(zoo)
  if(attr(x, "class")!="cdcc" & attr(x, "class")!="summary.cdcc") stop("wrong attribute")
  
  ndim <- ncol(x$data)            # the number of dimensions
  Ind <- diag(0, ndim); Ind[lower.tri(Ind)] <- 1

    # common arguments for plot
    if(is.null(sub)) sub <- ""
    if(is.null(xlab)) xlab <- "Time"
    if(is.null(lwd)) lwd <- 1

  if(item == "correlation"){
    obj <- x$CDCC
    obj <- obj[, as.vector(Ind == 1)]                       # remove duplicated columns
    obj <- zoo(obj, index(x$h))
    
    # arguments for plot
    if(is.null(main)) main <- "Conditional Correlations"
    if(is.null(ylab)) ylab <- colnames(obj)
    
    plot(obj, ylim=c(-1, 1), main=main, sub=sub, xlab=xlab, ylab=ylab, lwd=lwd, ...)
  } else if(item == "volatility"){
    obj <- zoo(sqrt(x$h), index(x$h))

    # arguments for plot
    if(is.null(main)) main <- "Volatilities"
    if(is.null(ylab)) ylab <- colnames(x$h)

    plot(obj, main=main, sub=sub, xlab=xlab, ylab=ylab, lwd=lwd, ...)
  } else if(item == "std.residuals"){
    obj <- zoo(x$z, index(x$z))

    # arguments for plot
    if(is.null(main)) main <- "Standardized Residuals"
    if(is.null(ylab)) ylab <- colnames(x$h)

    plot(obj, main=main, sub=sub, xlab=xlab, ylab=ylab, lwd=lwd, ...)
  } else if(item == "return"){
    obj <- zoo(x$data, index(x$data))

    # arguments for plot
    if(is.null(main)) main <- "Return Data"
    if(is.null(ylab)) ylab <- colnames(x$data)

    plot(obj, main=main, sub=sub, xlab=xlab, ylab=ylab, lwd=lwd, ...)
  } else {
    stop("The argument `item' is one of 'correlation', 'volatility', 'std.residuals', 'return'")
  }
}


