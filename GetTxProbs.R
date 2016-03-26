GetTxProbs <- function(X, beta.0, beta.low, beta.med, beta.high, beta.v.high){
    # Return treatment probability for each subject,
    #   according to the underlying prevalence beta.0.
    # X needs to have at least columns x.1...x.7
    logit <- with(X,
        beta.0 + 

        # these 3 are associated with tx assignment but not outcome
	beta.low  * x.1 + 
        beta.med  * x.2 + 
	beta.high * x.3 + 

        0.5 * beta.low  * x.4 + 
        0.3 * beta.low  * x.4^2 + 

	0.5 * beta.med  * x.5 + 
        0.5 * beta.high * x.6 + 
	0.3 * beta.high  * x.5 * x.6 + 

	0.5 * beta.v.high * x.7 +
	0.3 * beta.v.high * x.7^2
    )
    exp(logit) / (1 + exp(logit)) 
}

