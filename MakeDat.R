MakeDat <- function(N, true.avg.TE.cont, Beta.0,
        Beta.low, Beta.med, Beta.high, Beta.v.high){
    # return a data.frame of the covariates, true treatment effects, 
    #   treatment indicators, & outcomes

    x.1  <- rnorm(N, 0, 1) 
    x.2  <- rnorm(N, 0, 1) 
    x.3  <- rnorm(N, 0, 1) 
    x.4  <- rnorm(N, 0, 1) 
    x.5  <- rnorm(N, 0, 1) 
    x.6  <- rnorm(N, 0, 1) 
    x.7  <- rnorm(N, 0, 1) 
    x.8  <- rnorm(N, 0, 1) 
    x.9  <- rnorm(N, 0, 1) 
    x.10 <- rnorm(N, 0, 1)
        
    X <- data.frame(x.1 , x.2 , x.3 , x.4 , x.5 , x.6 , x.7 , x.8 , x.9 , x.10)

    # Generate treatment status for each subject
    p.treat <- GetTxProbs(X, Beta.0, Beta.low, Beta.med, Beta.high, Beta.v.high)
    X$treat <- rbinom(N, 1, p.treat)

    X$TE.cont <- rep(true.avg.TE.cont, N)

    # The AddY func returns a data.frame with Y added
    AddY(X, Beta.low, Beta.med, Beta.high, Beta.v.high)
}
