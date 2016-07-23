GetPairs <- function(Tr, X, exact, caliper, replace, Weight) {
    # Return a matrix of pairs, with tx rownums in col 1 & ctrl rownums in col 2
    # Return NULL if X is NULL or no matches found

    # all parameters are as in Matching::Match

    if(is.null(X[1])) return(NULL)

    mm <- tryCatch( 
        Match(
            Y             = NULL,
            Tr            = Tr,
            X             = X,
            M             = 1, # the ratio
            exact         = exact,
            caliper       = caliper,
            replace       = replace,
            ties          = FALSE, # randomly break ties
            Weight        = Weight,
            version       = "fast"
        ), 
        error= function(e) return(NULL)
    )
    if(is.null(mm) | length(mm$index.treated) == 0) return(NULL)

    ## return the pairs. Tx indices in col 1, ctrl in col 2
    cbind(mm$index.treated, mm$index.control)
}
