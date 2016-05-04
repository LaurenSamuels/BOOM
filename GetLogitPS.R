GetLogitPS <- function(dat, lrm.formula){
    tmp <- capture.output({
        fit <- tryCatch(
            lrm(lrm.formula, data= dat),
            error= function(e) return(NULL)
        )
    })
    if(is.null(fit)) return(NULL)
    
    fit$linear.predictors
}
