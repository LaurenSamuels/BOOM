GetPairs <- function(treatvec, X, exact, caliper){
    # Return a matrix of pairs, with tx rownums in col 1 & ctrl rownums in col 2
    # Return NULL if PS estimation fails or no matches found

    # treatvec is a vector of treatment assignments (0's and 1's)
    # logit.ps is a vector of logit propensity scores

    # uses Austin & Small (2014) method 2: greedy NNM on logit of PS
    #    within calipers of width equal to 0.2 * sd(logit(PS))
    #    and using random ordering of treated subjects

    # data is already in random order bec. that's how it was generated.
    # Match() processes treated obsns in order they are presented, so it's
    # processing treated obsns in random order.

    if(is.null(X[1])) return(NULL)

    mm <- tryCatch( 
        Match(
            Y       = NULL,
            Tr      = treatvec,
            X       = X,
            replace = FALSE,
            M       = 1, # the ratio
            ties    = FALSE, # randomly break ties
            exact   = exact,
            caliper = caliper,
            version = "fast"
        ), 
        error= function(e) return(NULL)
    )
    if(is.null(mm) | length(mm$index.treated) == 0) return(NULL)

    ## return the pairs. Tx indices in col 1, ctrl in col 2
    cbind(mm$index.treated, mm$index.control)
}
