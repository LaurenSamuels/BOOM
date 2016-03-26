GetLogitPS <- function(dat, lrm.formula){
    fit <- tryCatch(
        lrm(lrm.formula, data= dat),
        error= function(e) return(NULL)
    )
    if(is.null(fit)) return(NULL)

    fit$linear.predictors
}
