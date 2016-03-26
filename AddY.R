AddY <- function(X, beta.low, beta.med, beta.high, beta.v.high){
    # return X with an additional column, y
    # X must have at least the columns x.4...x.10 and TE.cont (true treatment effect)

    lp.core <- with(X, 
        0.5 * beta.low  * x.4 + 
        0.3 * beta.low  * x.4^2 + 

	0.5 * beta.med  * x.5 + 
        0.5 * beta.high * x.6 + 
	0.3 * beta.high  * x.5 * x.6 + 
        
	0.5 * beta.v.high * x.7 +
	0.3 * beta.v.high * x.7^2 +

        # these three are associated with outcome but not tx assignment
        beta.low * x.8 + 
        beta.med * x.9 + 
        beta.high * x.10
    )

    # Generate continuous outcome for each subject
    y <- with(X, 
        TE.cont * treat + 
        lp.core +
        rnorm(nrow(X), 0, 3))
    
    data.frame(X, y)
}

