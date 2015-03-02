lb.test <- function(obj, ...){
    if(is.null(dim(obj))) stop("obj is a vector. Use `Box.test()'.")
    ndim <- ncol(obj)
    lags <- c(1, seq(5, 50, 5))
    nlag <- length(lags)

#     stat <- matrix(0, nlag, ndim)
#     p.val <- matrix(0, nlag, ndim)
    stat <- matrix(0, nlag, 2*ndim)    
    for(i in 1:ndim){
        for(j in 1:nlag){
            LB.test <- Box.test(obj[,i], lag = lags[j], type = "Ljung")
            stat[j, (2*i-1)] <- LB.test$statistic
            stat[j, 2*i] <- LB.test$p.value
#             stat[j, i] <- LB.test$statistic
#             p.val[j, i] <- LB.test$p.value
        }
    }

    odd.ind <- seq(1, 2*ndim, by = 2)
    evn.ind <- seq(2, 2*ndim, by = 2)
    Snames <- numeric(2*ndim)
    if(is.null(colnames(obj))){
        Snames[odd.ind] <- paste("Series", 1:ndim, sep=" ")
    } else {
        Snames[odd.ind] <- colnames(obj)
    }
    
    Snames[evn.ind] <- "p-value"
    colnames(stat) <- Snames
    # colnames(p.val) <- rep("p-value", ndim)
    rownames(stat) <- paste("Lag", lags)
    # rownames(p.val) <- paste("Lag", lags)
    
    #ans <- list(stat = stat, p.val = p.val)    
    
    class(stat) <- "lbtest"
    return(stat)
}

print.lbtest <- function(x, digits = 3, ...){
    cat("\nLjung-Box Test for Serial Correlations", "\n")
    print.default(round(unclass(x), digits = digits), print.gap = 1, quote = FALSE)
  
}
