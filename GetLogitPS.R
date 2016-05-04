GetLogitPS <- function(dat, lrm.formula){
    capture.output({
        fit <- tryCatch(
            lrm(lrm.formula, data= dat, tol= my.tol),
            error= function(e) return(NULL)
        )
    })
    if(is.null(fit)) return(NULL)
    
    fit$linear.predictors
}
