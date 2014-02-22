# Lomnicki(1961)-Jarque-Bera (1987) test of normality
jb.test <- function(obj, ...){
  x <- obj
  if(!is.vector(x)){
    nobs <- nrow(x)
    ndim <- ncol(x)
    m <- colMeans(x)
    jb <- numeric(ndim)  
    if(is.null(colnames(x))){
      names(jb) <- paste("Series",1:ndim, sep=" ")
    } else {
      names(jb) <- colnames(x)
    }
       for(i in 1:ndim){
          x. <- x[,i] - m[i]
          v <- mean(x.^2)
          std.x <- x./sqrt(v)
          sk <- mean(std.x^3)
          kr <- mean(std.x^4) - 3
          jb[i] <- nobs/6*(sk^2 + 0.25*kr^2)
       }
  } else {
    nobs <- length(x)
    ndim <- 1
    x. <- x - mean(x)
    v <- mean(x.^2)
    std.x <- x./sqrt(v)
    sk <- mean(std.x^3)
    kr <- mean(std.x^4) - 3
    jb <- nobs/6*(sk^2 + 0.25*kr^2)

    if(is.null(names(x))){
      names(jb) <- "Series 1"
    } else {
      names(jb) <- names(jb)
    }

  }
  p.val <- pchisq(jb, 2, lower.tail=FALSE)
  out <- rbind(jb, p.val)
  rownames(out) <- c("test stat", "p-value")
  class(out) <- "jbtest"
  return(out)
}

print.jbtest <- function(x, digits = 3, ...){
  cat("\nLomnicki-Jarque-Bera Test of Normality", "\n")
# print.default(format(x, digits = digits), print.gap = 1, quote = FALSE)
  print.default(round(unclass(x), digits = digits), print.gap = 1, quote = FALSE)
}




