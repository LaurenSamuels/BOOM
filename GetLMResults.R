GetLMResults <- function(dat, form, true.TE.cont= NULL, justTxEst= FALSE, wts= NULL, tx.ind= 'treat'){
    # returns the TE estimate, SE estimate, and coverage.
    # wts, if used, is a string giving name of weights column in dat
    
    if(is.null(wts)){
        fit <- lm(form, data= dat)
    } else {
    	svyobj <- svydesign(ids= ~0, weights= dat[[wts]], data= dat)
        fit <- svyglm(form, svyobj)
    }

    if(justTxEst){
        coef(fit)[tx.ind]
    } else {
        c(  coef(fit)[tx.ind], 
            sqrt(vcov(fit)[tx.ind, tx.ind]),
            # from confint.svyglm help: The default is a Wald-type confidence interval, 
            #   adding and subtracting a multiple of the standard error. 
            confint(fit)[tx.ind, 1] <= true.TE.cont &&
                true.TE.cont <= confint(fit)[tx.ind, 2])
    }
}

