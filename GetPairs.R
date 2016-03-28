GetPairs <- function(Tr, X, exact, caliper, replace){
    # Return a matrix of pairs, with tx rownums in col 1 & ctrl rownums in col 2
    # Return NULL if PS estimation fails or no matches found

    # treatvec is a vector of treatment assignments (0's and 1's)
    # TR, X, exact, caliper, and replace are as in Matching::Match

    if(is.null(X[1])) return(NULL)

    mm <- tryCatch( 
        Match(
            Y       = NULL,
            Tr      = Tr,
            X       = X,
            replace = replace,
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
