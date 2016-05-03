CheckAndFixFormula <- function(dat, form) {
    # modify form as necessary to remove terms with only one unique value
    # See http://stackoverflow.com/questions/18171246/error-in-contrasts-when-defining-a-linear-model-in-r
    
    mod.matrix <- model.matrix(form, data= dat)
    n.terms <- ncol(mod.matrix) # this includes the intercept
    form.terms <- terms(form)
    form.response <- all.vars(form[[2]])


    gt1.unique <- apply(mod.matrix, 2, function(x) length(unique(x)) > 1)
    sum.gt1.unique <- sum(gt1.unique)
    if((sum.gt1.unique == n.terms) | 
        ((sum.gt1.unique == n.terms - 1) & "(Intercept)" %in% 
            colnames(mod.matrix))){
        return(list(form= form, removedTermFlag= 0))
    } else {
        # The reformulated model will automatically have an intercept,
        # so even if the original model had no intercept and there's
        # one constant column, we want to get rid of that column
        terms.to.remove <-  
            setdiff(colnames(mod.matrix)[!gt1.unique], "(Intercept)")
        term.positions.to.remove <- match(terms.to.remove,
            attr(form.terms, "term.labels"))
        tmpterms <- drop.terms(form.terms,
            dropx= term.positions.to.remove)         
        list(form= reformulate(attr(tmpterms, "term.labels"), 
                response= form.response),
            removedTermFlag= 1
        )
    }
}
