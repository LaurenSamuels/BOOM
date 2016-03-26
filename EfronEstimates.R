EfronEstimates <- function(count.matrix, est.TEs){
    # returns regular and bias-corrected estimates of SE 
    #   for the estimator that is the mean of the est.TEs, using formulas from 
    #   Efron 2014 JASA: Estimation and Accuracy After Model Selection

    # count.matrix is the matrix of bootstrap inclusion counts, 
    #   with nrow= number of bootstraps and ncol= number of subjects in original sample,
    #   and no missing values
    # est.TEs is the vector of estimated treatment effects, 
    #   with length = nrow(count.matrix)

    # In case there were boot samples where the matching failed:
    indicesToKeep <- rowSums(count.matrix) != 0
    count.matrix <- count.matrix[indicesToKeep, ]
    est.TEs <- est.TEs[indicesToKeep]

    # Number of people
    N <- ncol(count.matrix)
    # Number of bootstrap resamples (that actually worked)
    n.boot <- nrow(count.matrix)

    # Avg count for each person
    Y.star.dot.js <- colMeans(count.matrix)

    est.TE <- mean(est.TEs)
    # the centered means. Vector of length n.boot
    second.term.vec <- est.TEs - est.TE

    est.cov.efron <- est.cov.whe <- rep(NA, N)
    for(j in 1:N){
        # from Efron 2014:
        est.cov.efron[j] <- 
            sum((count.matrix[, j] - Y.star.dot.js[j]) * second.term.vec) / n.boot
        
        # from Wager Hastie Efron (WHE) 2014:
        # (also used in Efron BC estimate)
        est.cov.whe[j] <- 
            sum((count.matrix[, j] - 1) * second.term.vec) / n.boot
    }
    est.var.cont.efron <- sum(est.cov.efron^2)

    # for bias-corrected versions
    Z.matrix <- diff.matrix <- matrix(0, nrow= n.boot, ncol= N)
    #Z.matrix.alt <- diff.matrix.alt <- matrix(0, nrow= n.boot, ncol= N)

    for(i in 1:n.boot){
        for(j in 1:N){
            Z.matrix[i, j] <- 
                (count.matrix[i, j] - 1) * second.term.vec[i]
            diff.matrix[i, j] <- 
                Z.matrix[i, j] - est.cov.whe[j] 
    
            # Does it make a difference if we use mean count
            #   rather than expected count? No, but leaving code here as reference
            #Z.matrix.alt[i, j] <- 
            #    (count.matrix[i, j] - Y.star.dot.js[j]) * second.term.vec[i]
            #diff.matrix.alt[i, j] <- 
            #    Z.matrix.alt[i, j] - est.cov.efron[j] 
        }
    } 
    est.var.cont.efron.bc <- est.var.cont.efron -
        (1 / n.boot^2) * sum(rowSums(diff.matrix^2))
    #est.var.cont.efron.bc.alt <- est.var.cont.efron -
    #    (1 / n.boot^2) * sum(rowSums(diff.matrix.alt^2))
    # From WHE p. 1629 eq (7). This also returns same result; leaving code as reference
    #est.var.cont.efron.bc.whe <- est.var.cont.efron -
    #    (N / n.boot^2) * sum(second.term.vec^2)

    c(est.var.cont.efron, est.var.cont.efron.bc#, 
        #est.var.cont.efron.bc.alt, est.var.cont.efron.bc.whe
    )
}

