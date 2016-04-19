GetPrognosticScore <- function(dat, lm.formula, is.control){
    # fit a prognostic model on the control subjects,
    #    then return prognostic scores for the whole sample
    # is.control is a boolean vector

    fit <- tryCatch(
        lm(lm.formula, data= dat[is.control, ]),
        error= function(e) return(NULL)
    )
    if(is.null(fit)) return(NULL)

    predict(fit, dat)
}
