# residual diagnostics
residDiag <- function(object){

    if(class(object) != "ccc" & class(object) != "cdcc" & class(object) != "dcc" & class(object) != "zoo" & 
       class(object) != "summary.ccc" & class(object) != "summary.cdcc" & class(object) != "summary.dcc"){
        # stop if the object is not an appropreate class
        stop("the object must be 'ccc', 'dcc' or 'cdcc' class")
    }
    
    jb <- jb.test(object$z)
    lb <- lb.test(object$z)
    lb2 <- lb.test(object$z^2)
    
    obj <- list(jb=jb, lb=lb, lb2=lb2)
    
    class(obj) <- "residDiag"
    return(obj)
}

print.residDiag <- function(x, digits = 4, ...){
    
    print.jbtest(x$jb)
    cat("\nFor the standardized residuals:", "")
    print.lbtest(x$lb)
    cat("\nFor the squared standardized residuals:", "")
    print.lbtest(x$lb2)
    
}






