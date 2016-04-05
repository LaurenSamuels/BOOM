BOOM <- function(dat, n.boot, ps.formula, lm.formula= NULL, 
    mcCores= 2, tx.indicator= "treat", outcome= "y", seed= 1235,
    return.dat= TRUE, recalcDistance= TRUE,
    exactVarNames= NULL, caliper= 0.2, replace= FALSE,
    conf.level= 0.95
    ){
    # Returns a vector of various summary values from the BOOM procedure

    # dat: the original dataset
    # n.boot: number of bootstrap resamples to use
    # ps.formula: propensity score formula to use w/ lrm
    # lm.formula: optional outcome formula
    # mcCores: number of cores for mclapply()
    # tx.indicator: name of the treatment indicator variable in dat (must be 1/0 for now)
    # y: name of the continuous outcome variable in dat
    # seed: random seed for use by this function. 
    # return.dat (boolean): Return the original dataset?
    # recalcDistance (boolean): re-calculate the PS or other distance measure in each resample?
    # exactVarNames: vector of names of variables on which to match exactly
    # caliper: optional caliper for use with Matching. 
    # replace (boolean): Match with replacement?
    # conf.level: level to use for confidence intervals
    
    N <- nrow(dat) # tot number of subjects
    treat <- dat[[tx.indicator]]

    # TODO: error handling: missing values of tx.indicator
    # TODO: error handling: missing values of outcome
    # TODO: error handling: missing values of vars in PS model
    # TODO: error handling: missing values of vars in outcome model
    #***************
    # TODO: error handling: tx indicator not in models
    #***************
    # TODO: error handling: exactvars not in dset

    # TODO: return CIs (ask for level)
    # TODO: allow chr/factor tx indicator?

    # Argument checking
    # from t.test code: getAnywhere("t.test.default")
    if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) || 
        conf.level < 0 || conf.level > 1)) 
        stop("'conf.level' must be a single number between 0 and 1")


    ###########################################################
    # Take care of the random seed
    # from Cole's code in nbpMatching pkg
    if(exists(".Random.seed", envir = .GlobalEnv)) {
        save.seed <- get(".Random.seed", envir= .GlobalEnv)
        on.exit(assign(".Random.seed", save.seed, envir = .GlobalEnv))
    } else {
        on.exit(rm(".Random.seed", envir = .GlobalEnv))
    }
    set.seed(seed, kind= "L'Ecuyer-CMRG") # for mclapply()
    mc.reset.stream() # Not sure this is necessary
    ###########################################################

    #################
    # first: get some estimates on original data
    logitPS.orig <- GetLogitPS(dat, ps.formula)
    PS.orig <- InvLogit(logitPS.orig)
    att.wts.orig <- treat + (1 - treat) * PS.orig / (1 - PS.orig)
    # from Li and Greene 2013:
    match.wts.orig = pmin(1 - PS.orig, PS.orig) / 
        (treat * PS.orig + (1 - treat) * (1 - PS.orig))
    #################


    isTreated <- treat == 1
    # todo someday: the whole indexing/tracking thing would be easier w/ data.table
    tx.orig.dat <- dat[isTreated, ]
    n.treated.orig <- nrow(tx.orig.dat)
    tx.orig.ids.easy <- 1:n.treated.orig

    isControl <- treat == 0
    ctrl.orig.dat <- dat[isControl, ]
    ctrl.orig.ids.easy <- (n.treated.orig + 1) : N

    all.orig.ids.easy <- 1:N

    # for use w/ non-recalculated PS
    logitPS.tx.orig <- logitPS.orig[isTreated]
    logitPS.ctrl.orig <- logitPS.orig[isControl]

    bootStuff <- mclapply(1:n.boot, function(x) {
        # Modified from Austin & Small (2014) --- they did not condition on
        # observed tx & ctrl group sizes
        tx.sample.indices <- 
            sample(1:n.treated.orig, size= n.treated.orig, replace= TRUE)
        tx.sample <- tx.orig.dat[tx.sample.indices, ]

        ctrl.sample.indices <- 
            sample(1:nrow(ctrl.orig.dat), size= nrow(ctrl.orig.dat), replace= TRUE)
        ctrl.sample <- ctrl.orig.dat[ctrl.sample.indices, ]

        # in boot.sample, the first n.treated.orig rows are treated people,
        #   remaining rows are control people
        boot.sample <- rbind(tx.sample, ctrl.sample, make.row.names= FALSE)

        est.mean.tx.tmp <- est.mean.ctrl.tmp <- est.TE.lm.tmp <- NA
        count.vector.tmp <- rep(0, N)
        num.errs.tmp <- 0
        logitPS.ordered.tmp <- rep(NA, N)
        count.vector.matched.tmp <- rep(0, N)

        if (recalcDistance) {
            logitPS <- GetLogitPS(boot.sample, ps.formula)
        } else {
            logitPS <- c(logitPS.tx.orig[tx.sample.indices], logitPS.ctrl.orig[ctrl.sample.indices])
        }
        if (is.null(exactVarNames)) {
            my.X <- logitPS
            my.exact <- FALSE
        } else {
            my.X <- cbind(logitPS, boot.sample[, exactVarNames])
            my.exact <- c(FALSE, rep(TRUE, length(exactVarNames)))
        }
        pairIndices <- GetPairs(
            Tr      = boot.sample[[tx.indicator]], 
            X       = my.X, 
            exact   = my.exact,
            caliper = caliper,
            replace = replace
        )
        # Tx indices are in 1st col, ctrl in 2nd col.
        #    These indices are indices from boot.sample
        # We are handling PS estimation/matching errors by skipping that resample.  
        #    Maybe not the best way, but it rarely happens.
        if(!is.null(pairIndices)){
            est.mean.tx.tmp <- 
                mean(boot.sample[pairIndices[, 1], outcome])
            est.mean.ctrl.tmp <- 
                mean(boot.sample[pairIndices[, 2], outcome])

            if (!is.null(lm.formula)) {
                fit <- lm(lm.formula, 
                    data= boot.sample[c(pairIndices[, 1], pairIndices[, 2]), ])
                est.TE.lm.tmp <- coef(fit)[tx.indicator]
            }

            # fill in counts for Efron calculations
            #   (unrelated to matching, but we do not want to do if
            #    matching was unsuccessful)
            tx.orig.ids.easy.insample <- tx.sample.indices
            ctrl.orig.ids.easy.insample <-
                ctrl.orig.ids.easy[ctrl.sample.indices]
            all.orig.ids.easy.insample <- c(tx.orig.ids.easy.insample,
                ctrl.orig.ids.easy.insample)
            both.tbl <- table(all.orig.ids.easy.insample) 
            # count.vector.tmp is in `easy' order (all tx, then all ctrl)
            count.vector.tmp[as.numeric(names(both.tbl))] <-
                both.tbl

            # For calculating avg logitPS (may not be a useful quantity)
            # logitPS vector is ordered by easy ID
            logitPS.tmp.matrix <- cbind(logitPS, all.orig.ids.easy.insample)
            logitPS.tmp.matrix <- logitPS.tmp.matrix[!duplicated(all.orig.ids.easy.insample), ]
            logitPS.tmp.vec <- logitPS.tmp.matrix[, 1]
            # this vector is in `easy' order (all tx, then all ctrl)
            logitPS.ordered.tmp[logitPS.tmp.matrix[, 2]] <- logitPS.tmp.vec

            # counts for weights
            all.orig.ids.easy.matched <- all.orig.ids.easy.insample[c(pairIndices[, 1],
                pairIndices[, 2])]
            both.tbl.matched <- table(all.orig.ids.easy.matched) 
            # this vector is in `easy' order (all tx, then all ctrl)
            count.vector.matched.tmp[as.numeric(names(both.tbl.matched))] <-
                both.tbl.matched

        } else { # we did not get a match in this resample
            cat("Hit problem in boot", x, ".\n") 
            num.errs.tmp <- num.errs.tmp + 1
        } # end processing for resamples with errors in PS estimation or matching

        # return from mclapply:
        list(
            # TODO: give these names & use names below.
            # scalars
            est.mean.tx.tmp, #1
            est.mean.ctrl.tmp, #2
            est.TE.lm.tmp, #3

            # vector w/ length N
            count.vector.tmp, #4

            # scalar
            num.errs.tmp, #5

            # vectors w/ length N
            logitPS.ordered.tmp, #6
            count.vector.matched.tmp #7
            )
    },
        mc.cores           = mcCores,
        mc.preschedule     = TRUE,
        mc.set.seed        = TRUE,
        mc.allow.recursive = FALSE
    ) # end bootstrap resampling (mc)lapply

    # vectors of length n.boot
    est.means.tx   <- do.call(c, lapply(bootStuff, function(x) x[[1]]))
    est.means.ctrl <- do.call(c, lapply(bootStuff, function(x) x[[2]]))
    est.TEs <- est.means.tx - est.means.ctrl
    est.TEs.lm <- do.call(c, lapply(bootStuff, function(x) x[[3]]))

    # matrix: row = resample; col = subject
    count.matrix <- 
        do.call(rbind, lapply(bootStuff, function(x) x[[4]]))
    count.matrix.tx <- count.matrix[, tx.orig.ids.easy]
    count.matrix.ctrl <- count.matrix[, ctrl.orig.ids.easy]

    # vectors of length n.boot
    num.errs <- do.call(c, lapply(bootStuff, function(x) x[[5]]))

    # matrix: row = resample; col = subject
    logitPS.matrix <- 
        do.call(rbind, lapply(bootStuff, function(x) x[[6]]))
    logitPS.avg.easy <- colMeans(logitPS.matrix, na.rm= TRUE)
    # now get it in the order of the original dataset
    logitPS.avg <- rep(NA, N)
    logitPS.avg[isTreated] <- logitPS.avg.easy[tx.orig.ids.easy]
    logitPS.avg[isControl] <- logitPS.avg.easy[ctrl.orig.ids.easy]
    PS.avg <- InvLogit(logitPS.avg)
    
    # matrix: row = resample; col = subject
    count.matrix.matched <- 
        do.call(rbind, lapply(bootStuff, function(x) x[[7]]))
    #num.tx.matched <- 
    #   rowSums(count.matrix.matched[, tx.orig.ids.easy])
    #row.multiplier <- n.treated.orig / num.tx.matched
    # 1st el. in row.multiplier X each element in 1st row of c.m.m., etc.

    #BOOM.wts.easy <- colMeans(row.multiplier * count.matrix.matched)

    BOOM.wts.easy <- colSums(count.matrix.matched) / colSums(count.matrix)

    # now get it in the order of the original dataset
    BOOM.wts <- rep(NA, N)
    BOOM.wts[isTreated] <- BOOM.wts.easy[tx.orig.ids.easy]
    BOOM.wts[isControl] <- BOOM.wts.easy[ctrl.orig.ids.easy]


    # Summary stats (bagged statistics)
    est.TE <- mean(est.TEs, na.rm= TRUE)
    # TODO: return NA if no model given
    est.TE.lm <- mean(est.TEs.lm, na.rm= TRUE)

    tot.errs <- sum(num.errs, na.rm= TRUE)

    # Efron calculations: Regular & BC
    NUM.THINGS.RETURNED.BY.EE <- 2
    efron2.tx <- efron2.ctrl <- efron.lm <- 
        rep(NA, NUM.THINGS.RETURNED.BY.EE)
    efron2.tx <- 
        EfronEstimates(count.matrix.tx, est.means.tx)
    efron2.ctrl <- 
        EfronEstimates(count.matrix.ctrl, est.means.ctrl)
    est.SE.efron    <- sqrt(efron2.tx[1] + efron2.ctrl[1])
    est.SE.efron.bc <- sqrt(efron2.tx[2] + efron2.ctrl[2])

    # TODO: return NA if no model given
    efron.lm <- 
        EfronEstimates(count.matrix, est.TEs.lm)
    est.SE.lm.efron= sqrt(efron.lm[1])
    est.SE.lm.efron.bc= sqrt(efron.lm[2])

    # TODO: maybe take these out? on the other hand they make an interesting comparison
    est.SE.naive <- sd(est.TEs, na.rm= TRUE)
    # TODO: return NA if no model given
    est.SE.lm.naive <- sd(est.TEs.lm, na.rm= TRUE)

    # todo: add pvals. This code is from t.test, just here as reminder:
    #pval <- 2 * pt(-abs(tstat), df)


    list(  
        # vectors of length N, in order of original dataset
        treat= treat,
        logitPS.orig = logitPS.orig, 
        PS.orig = PS.orig, 
        att.wts.orig = att.wts.orig,
        # from Li and Greene 2013:
        match.wts.orig = match.wts.orig,

        # scalar results from BOOM w/ NO further covariate adj
        est.TE= est.TE, 
        est.SE.naive= est.SE.naive,  # not recommended for use; just for comparison
        est.SE.efron= est.SE.efron,
        est.SE.efron.bc= est.SE.efron.bc,

        # scalar results from BOOM with further covariate adj
        est.TE.lm= est.TE.lm, 
        est.SE.lm.naive= est.SE.lm.naive,  # not recommended for use; just for comparison
        est.SE.lm.efron= est.SE.lm.efron,
        est.SE.lm.efron.bc= est.SE.lm.efron.bc,

        # confidence intervals
        conf.int.efron = GetConfInt(est.TE, est.SE.efron, conf.level),
        conf.int.efron.bc = GetConfInt(est.TE, est.SE.efron.bc, conf.level),
        conf.int.lm.efron = GetConfInt(est.TE.lm, est.SE.lm.efron, conf.level),
        conf.int.lm.efron.bc = GetConfInt(est.TE.lm, est.SE.lm.efron.bc, conf.level),

        # scalars
        tot.errs= tot.errs,
        conf.level= conf.level,

        # vectors of length n.boot
        est.TEs= est.TEs,
        est.TEs.lm= est.TEs.lm,

        # vectors of length N, in order of original dataset
        logitPS.avg= logitPS.avg, # may or may not be useful
        PS.avg= PS.avg,
        BOOM.wts= BOOM.wts,

        dat= if (return.dat) dat else NA
    )
}

