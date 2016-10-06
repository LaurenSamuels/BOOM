BOOM <- function(dat, n.boot, tx.indicator, outcome, 
    distance.type= "logitPropensity",
    matching.pkg= "Matching",
    recalc.distance= TRUE,
    propensity.formula = NULL, 
    prognostic.formula = NULL, 
    outcome.formula= NULL, 
    MD.vars= NULL,
    mc.cores= 2, seed= 1235,
    return.dat= TRUE, 
    conf.level= 0.95,
    exact.var.names= NULL, 
    caliper= 0.2, replace= FALSE,
    threshold= NA
    ){
    # Returns a vector of various summary values from the BOOM procedure

    # dat: the original dataset
    # n.boot: number of bootstrap resamples to use
    # tx.indicator: name of the treatment indicator variable in dat (must be 1/0 for now)
    # outcome: name of the continuous outcome variable in dat
    # distance.type: one of "propensity" (propensity score), "logitPropensity (logit propensity score), "prognostic" (prognostic score), or "MD" (Mahalanobis) 
    # matching.pkg: one of "Matching" or "nbpMatching"
    # recalc.distance (boolean): re-calculate the PS or other distance measure in each resample?
    # propensity.formula: propensity score formula to use w/ lrm
    # prognostic.formula: prognostic score formula. Must not contain factor variables.
    # outcome.formula: optional outcome formula. Currently must not involve interactions with the treatment indicator.
    # mc.cores: number of cores for mclapply()
    # seed: random seed for use by this function. 
    # return.dat (boolean): Return the original dataset?
    # conf.level: level to use for confidence intervals
    # exact.var.names: vector of names of variables on which to match exactly
    # caliper: as in Matching::Match
    # replace: as in Matching::Match
    # threshold: as in nbpMatching::nonbimatch
    
    N         <- nrow(dat) # tot number of subjects
    treat     <- dat[[tx.indicator]]
    isTreated <- treat == 1
    isControl <- treat == 0

    # TODO: error handling: missing values of tx.indicator
    # TODO: error handling: missing values of outcome
    # TODO: error handling: missing values of vars in PS model
    # TODO: error handling: missing values of vars in outcome model
    #***************
    # TODO: error handling: tx indicator not in models
    #***************
    # TODO: error handling: exactvars not in dset
    # TODO: error handling: illegal distance type
    # TODO: error handling: factor vars in prognostic.formula

    # TODO: allow chr/factor tx indicator?
    # TODO: return the whole call (all args used in call)
    # TODO: error handling: can't have recalc.distance= FALSE & matching.pkg= "Matching"
    #        with distance.type = "MD"
    # TODO: message about not using both threshold and caliper
    # TODO: require n.boot >= 2 OR put in processing for n.boot= 1
    # TODO: am still getting "singular information matrix in lrm.fit (rank= 24 ).  Offending variable(s):" on occasion. CheckAndFix is not taking care of all collinearity, just collinearity w/ the intercept. Keep adding to the function.


    # Argument checking
    # from t.test code: getAnywhere("t.test.default")
    if (!missing(conf.level) && (length(conf.level) != 1 || 
        !is.finite(conf.level) || conf.level < 0 || conf.level > 1)) 
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
    logitPS.orig <- PS.orig <- att.wts.orig <- ate.wts.orig <-
        match.wts.orig <- NA
    progscore.orig <- NA

    # people might provide a propensity score formula even if
    #    they are not using it for BOOM distance measure
    if (!is.null(propensity.formula)) {
        logitPS.orig <- GetLogitPS(dat, propensity.formula)
        # TODO: make this nicer/proper error handling 
        if (is.null(logitPS.orig)) print("Can't fit PS model on original data")
        PS.orig <- InvLogit(logitPS.orig)
        att.wts.orig <- treat + (1 - treat) * PS.orig / (1 - PS.orig)
        ate.wts.orig <- treat / PS.orig + (1 - treat) / (1 - PS.orig)
        # from Li and Greene 2013:
        match.wts.orig = pmin(1 - PS.orig, PS.orig) / 
            (treat * PS.orig + (1 - treat) * (1 - PS.orig))
    }
    #################
    # and they might provide a prognostic score formula even if
    #    they are not using it for BOOM distance measure
    if (!is.null(prognostic.formula)) {
        
        # Preparation for fitting inside boot: in case some factors
        #    have only one level
        prog.form.terms <- terms(prognostic.formula)
        prog.form.response <- all.vars(prognostic.formula[[2]])
        prog.form.othervars <- all.vars(prognostic.formula[[3]])


        progscore.orig <- GetPrognosticScore(dat, 
            prognostic.formula, isControl)
        # TODO: make this nicer/proper error handling 
        if (is.null(progscore.orig)) print("Can't fit prognostic model on original data")
    }
    #################
    
    # Preparation for fitting lm inside boot: in case some factors
    #    have only one level
    out.form.terms <- terms(outcome.formula)
    out.form.response <- all.vars(outcome.formula[[2]])
    out.form.othervars <- all.vars(outcome.formula[[3]])


    # todo someday: the whole indexing/tracking thing might be easier w/ data.table
    tx.orig.dat <- dat[isTreated, ]
    n.treated.orig <- sum(isTreated)
    tx.orig.ids.easy <- 1:n.treated.orig

    ctrl.orig.dat <- dat[isControl, ]
    n.ctrl.orig <- N - n.treated.orig
    ctrl.orig.ids.easy <- (n.treated.orig + 1) : N

    all.orig.ids.easy <- 1:N

    # These will stay the same from boot to boot (because of the way sample is generated)
    isTreated.boot <- c(rep(TRUE, n.treated.orig), rep(FALSE, n.ctrl.orig))
    isControl.boot <- !isTreated.boot

    # for use w/ non-recalculated distances
    if (!recalc.distance) {
        if (distance.type == "logitPropensity") {
            logitPS.tx.orig   <- logitPS.orig[isTreated]
            logitPS.ctrl.orig <- logitPS.orig[isControl]
        } else if (distance.type == "propensity") {
            PS.tx.orig   <- PS.orig[isTreated]
            PS.ctrl.orig <- PS.orig[isControl]
        } else if (distance.type == "prognostic") {
            progscore.tx.orig   <- progscore.orig[isTreated]
            progscore.ctrl.orig <- progscore.orig[isControl]
        } else if (distance.type == "MD" & matching.pkg == "nbpMatching") {
            # reorder the dataset
            datForDist <- rbind(tx.orig.dat, ctrl.orig.dat)
            distmat.tx.ctrl.orig <- as.matrix(gendistance(cbind(isTreated.boot, 
                datForDist[, MD.vars]), prevent = 1)$dist)
            # get rid of phantoms
            distmat.tx.ctrl.orig <- distmat.tx.ctrl.orig[1:N, 1:N]
        }
    }


    bootStuff <- mclapply(1:n.boot, function(x) {
        # Modified from Austin & Small (2014) --- they did not condition on
        # observed tx & ctrl group sizes
        tx.sample.indices <- 
            sample(1:n.treated.orig, size= n.treated.orig, replace= TRUE)
        tx.sample <- tx.orig.dat[tx.sample.indices, ]

        ctrl.sample.indices <- 
            sample(1:n.ctrl.orig, size= n.ctrl.orig, replace= TRUE)
        ctrl.sample <- ctrl.orig.dat[ctrl.sample.indices, ]

        
        # Which of the original "easy id's" are in the boot sample?
        # The "easy id's" come from the rearranged dataset with all tx first,
        #   then all control
        tx.orig.ids.easy.insample <- tx.sample.indices
        ctrl.orig.ids.easy.insample <-
            ctrl.orig.ids.easy[ctrl.sample.indices]
        all.orig.ids.easy.insample <- 
            c(tx.orig.ids.easy.insample,
            ctrl.orig.ids.easy.insample)


        # in boot.sample, the first n.treated.orig rows are treated people,
        #   remaining rows are control people
        boot.sample <- rbind(tx.sample, ctrl.sample, make.row.names= FALSE)

        est.mean.tx.tmp <- est.mean.ctrl.tmp <- est.TE.lm.tmp <- NA
        count.vector.tmp <- rep(0, N)
        num.match.errs.tmp <- 0
        num.times.ps.form.changed.tmp   <- NA
        num.times.ps.form.failed.tmp   <- NA
        num.times.prog.form.changed.tmp <- NA
        num.times.out.form.changed.tmp  <- NA
        
        logitPS.ordered.tmp <- rep(NA, N)
        count.vector.matched.tmp <- rep(0, N)


        if (recalc.distance) {
            if (distance.type %in% c("logitPropensity", "propensity")) {
                # the original propensity formula might not work in this sample
                # TODO: maybe: instead of flag, save list of terms removed??
                # TODO: the CheckAndFix function needs some more steps:
                #      remove term if singular cov. matrix
                ps.check <- 
                    CheckAndFixFormula(boot.sample, propensity.formula)
                num.times.ps.form.changed.tmp <- 
                    ps.check$removedTermFlag
                num.times.ps.form.failed.tmp <- 0
                logitPS <- 
                    GetLogitPS(boot.sample, ps.check$form)
                if (is.null(logitPS)) {
                    num.times.ps.form.failed.tmp <- 1
                    PS <- NULL
                } else {
                    PS <- InvLogit(logitPS)
                }
            } else if (distance.type == "prognostic") {
                # the original prognostic formula might not work in this ctrl sample.
                #   (If it works in ctrl sample, it will work in whole boot.sample)
                prog.check <- CheckAndFixFormula(ctrl.sample, prognostic.formula)
                num.times.prog.form.changed.tmp <- prog.check$removedTermFlag

                progscore <- GetPrognosticScore(boot.sample, 
                    prog.check$form, isControl.boot)
            } else if (distance.type == "MD" & matching.pkg == "nbpMatching") {
                distmat <- distancematrix(gendistance(
                    cbind(isTreated.boot, boot.sample[, MD.vars]), 
                    prevent = 1))
            }
            # Note that if distance.type == "MD" & matching.pkg =="Matching",
            #   recalculation of distance is built-in (required) at this point.
            # I think we could allow the distance to not be recalculated
            #   by supplying a custom covariance matrix to Match(), but
            #   I'm not sure and I don't see a reason to investigate at this point
        } else { # using a fixed distance
            if (distance.type == "logitPropensity") {
                logitPS <- c(logitPS.tx.orig[tx.sample.indices], 
                    logitPS.ctrl.orig[ctrl.sample.indices])
            } else if (distance.type == "propensity") {
                PS <- c(PS.tx.orig[tx.sample.indices], 
                    PS.ctrl.orig[ctrl.sample.indices])
            } else if (distance.type == "prognostic") {
                progscore <- c(progscore.tx.orig[tx.sample.indices], 
                    progscore.ctrl.orig[ctrl.sample.indices])
            } else if (distance.type == "MD" & matching.pkg == "nbpMatching") {
                distmat <- distancematrix(distmat.tx.ctrl.orig[all.orig.ids.easy.insample,
                    all.orig.ids.easy.insample])
            }
        }

        if (distance.type == "logitPropensity") {
            my.X <- logitPS
        } else if (distance.type == "propensity") {
            my.X <- PS
        } else if (distance.type == "prognostic") {
            my.X <- progscore
        } else if (distance.type == "MD" & matching.pkg == "Matching") {
            my.X <- boot.sample[, MD.vars]
        }

        if (is.null(exact.var.names)) {
            my.exact <- FALSE
        } else { # we want to match exactly on some vars
            if (distance.type %in% c("logitPropensity", "propensity", "prognostic")) {
                my.X <- cbind(my.X, boot.sample[, exact.var.names])
                my.exact <- 
                    c(FALSE, rep(TRUE, length(exact.var.names)))
            } else if (distance.type == "MD") {
                my.exact <- vector(FALSE, length(MD.vars))
                my.exact[MD.vars %in% exact.var.names] <- TRUE
            }
        }

        if (matching.pkg == "Matching") {
            my.weight <- 
                ifelse(distance.type %in% c("logitPropensity", "propensity", "prognostic"), 
                1, 2)
            pairIndices <- GetPairs(
                Tr            = isTreated.boot,
                X             = my.X,
                exact         = my.exact,
                caliper       = caliper,
                replace       = replace,
                Weight        = my.weight
            )
        # Tx indices are in 1st col, ctrl in 2nd col.
        # These indices are indices from boot.sample
        } else if (matching.pkg == "nbpMatching") {
            matchinfo <- nonbimatch(distmat, threshold= threshold)$matches[1:N, ]
            # get the boot.sample indices of all units that have a legit match
            toKeep <- matchinfo[with(matchinfo, Group2.Row <= N &
                is.finite(Distance)), "Group1.Row"]
            # TODO: make this NULL if no matches
            pairIndices <- cbind(toKeep[toKeep <= n.treated.orig],
                toKeep[toKeep > n.treated.orig])
        }


        # We are handling PS estimation/matching errors 
        #    by skipping that resample.  
        #    Maybe not the best way, but it rarely happens.
        if(!is.null(pairIndices)){ # TODO: maybe set a min # of matches to proceed?
            est.mean.tx.tmp <- 
                mean(boot.sample[pairIndices[, 1], outcome])
            est.mean.ctrl.tmp <- 
                mean(boot.sample[pairIndices[, 2], outcome])
            matched.indices <- c(pairIndices[, 1], pairIndices[, 2])

            if (!is.null(outcome.formula)) {
                out.check <- 
                    CheckAndFixFormula(boot.sample[matched.indices, ], outcome.formula)
                num.times.out.form.changed.tmp <- out.check$removedTermFlag
                fit <- lm(out.check$form, 
                    data= boot.sample[matched.indices, ])
                est.TE.lm.tmp <- coef(fit)[tx.indicator]
            }

            # fill in counts for Efron calculations
            #   (unrelated to matching, but we do not want to do if
            #    matching was unsuccessful)
            both.tbl <- table(all.orig.ids.easy.insample) 
            # count.vector.tmp is in `easy' order (all tx, then all ctrl)
            count.vector.tmp[as.numeric(names(both.tbl))] <-
                both.tbl

            # For calculating avg logitPS (may not be a useful quantity)
            if (distance.type == "propensity") {
                # logitPS vector is ordered by easy ID
                logitPS.tmp.matrix <- cbind(logitPS, all.orig.ids.easy.insample)
                logitPS.tmp.matrix <- logitPS.tmp.matrix[!duplicated(all.orig.ids.easy.insample), ]
                logitPS.tmp.vec <- logitPS.tmp.matrix[, 1]
                # this vector is in `easy' order (all tx, then all ctrl)
                logitPS.ordered.tmp[logitPS.tmp.matrix[, 2]] <- logitPS.tmp.vec
            }

            # counts for weights
            all.orig.ids.easy.matched <- 
                all.orig.ids.easy.insample[c(pairIndices[, 1],
                pairIndices[, 2])]
            both.tbl.matched <- table(all.orig.ids.easy.matched) 
            # this vector is in `easy' order (all tx, then all ctrl)
            count.vector.matched.tmp[as.numeric(names(both.tbl.matched))] <-
                both.tbl.matched
        } else { # we did not get a match in this resample
            # TODO: I took out the cat statement because it might
            #   interfere w/mclapply. But it would be good to have
            #   a different way to track the errors. Maybe
            #   return the whole 1/0 error vecs in addition to the
            #   total counts??
            #cat("Hit problem in boot", x, ".\n") 
            num.match.errs.tmp <- 1
        } # end processing for resamples with errors in PS estimation or matching

        # return from mclapply:
        list(
            # scalars
            est.mean.tx   = est.mean.tx.tmp,
            est.mean.ctrl = est.mean.ctrl.tmp,
            est.TE.lm     = est.TE.lm.tmp,
            # TODO: these are all called "num," but really they're 0/1 flags
            num.match.errs              = num.match.errs.tmp,
            num.times.ps.form.changed   = num.times.ps.form.changed.tmp, 
            num.times.ps.form.failed   = num.times.ps.form.failed.tmp, 
            num.times.prog.form.changed = num.times.prog.form.changed.tmp, 
            num.times.out.form.changed  = num.times.out.form.changed.tmp, 

            # vectors w/ length N
            logitPS.ordered      = logitPS.ordered.tmp,
            count.vector         = count.vector.tmp,
            count.vector.matched = count.vector.matched.tmp
            )
    },
        mc.cores           = mc.cores,
        mc.preschedule     = TRUE,
        mc.set.seed        = TRUE,
        mc.allow.recursive = FALSE
    ) # end bootstrap resampling (mc)lapply

    # vectors of length n.boot
    est.means.tx   <- do.call(c, lapply(bootStuff, function(x) 
        x[["est.mean.tx"]]))
    est.means.ctrl <- do.call(c, lapply(bootStuff, function(x) 
        x[["est.mean.ctrl"]]))
    est.TEs <- est.means.tx - est.means.ctrl
    est.TEs.lm <- do.call(c, lapply(bootStuff, function(x) 
        x[["est.TE.lm"]]))

    # matrix: row = resample; col = subject
    # columns are in 'easy' order
    count.matrix <- do.call(rbind, lapply(bootStuff, function(x) 
        x[["count.vector"]]))
    count.matrix.tx <- count.matrix[, tx.orig.ids.easy]
    count.matrix.ctrl <- count.matrix[, ctrl.orig.ids.easy]
    # vector of length n, in 'easy' order
    count.vector.easy <- colSums(count.matrix)
    # now get it in the order of the original dataset
    count.vector <- rep(NA, N)
    count.vector[isTreated] <- count.vector.easy[tx.orig.ids.easy]
    count.vector[isControl] <- count.vector.easy[ctrl.orig.ids.easy]

    # vectors of length n.boot
    num.match.errs <- do.call(c, lapply(bootStuff, function(x) 
        x[["num.match.errs"]]))
    num.times.ps.form.changed <- do.call(c, lapply(bootStuff, function(x) 
        x[["num.times.ps.form.changed"]]))
    num.times.ps.form.failed <- do.call(c, lapply(bootStuff, function(x) 
        x[["num.times.ps.form.failed"]]))
    num.times.prog.form.changed <- do.call(c, lapply(bootStuff, function(x) 
        x[["num.times.prog.form.changed"]]))
    num.times.out.form.changed <- do.call(c, lapply(bootStuff, function(x) 
        x[["num.times.out.form.changed"]]))

    # matrix: row = resample; col = subject
    # columns are in 'easy' order
    logitPS.matrix <- do.call(rbind, lapply(bootStuff, function(x) 
        x[["logitPS.ordered"]]))
    logitPS.avg.easy <- colMeans(logitPS.matrix, na.rm= TRUE)
    # now get it in the order of the original dataset
    logitPS.avg <- rep(NA, N)
    logitPS.avg[isTreated] <- logitPS.avg.easy[tx.orig.ids.easy]
    logitPS.avg[isControl] <- logitPS.avg.easy[ctrl.orig.ids.easy]
    
    # matrix: row = resample; col = subject
    # columns are in 'easy' order
    count.matrix.matched <- do.call(rbind, lapply(bootStuff, 
        function(x) x[["count.vector.matched"]]))

    BOOM.wts.easy <- ifelse(count.vector.easy == 0, 0,
        colSums(count.matrix.matched) / count.vector.easy)
    # now get it in the order of the original dataset
    BOOM.wts <- rep(NA, N)
    BOOM.wts[isTreated] <- BOOM.wts.easy[tx.orig.ids.easy]
    BOOM.wts[isControl] <- BOOM.wts.easy[ctrl.orig.ids.easy]

    # number of matched pairs per resample
    num.pairs <- rowSums(count.matrix.matched) / 2
    

    # Summary stats (bagged statistics)
    est.TE <- mean(est.TEs, na.rm= TRUE)
    # TODO: return NA if no model given
    est.TE.lm <- mean(est.TEs.lm, na.rm= TRUE)

    tot.match.errs <- sum(num.match.errs, na.rm= TRUE)
    tot.times.ps.form.changed <- sum(num.times.ps.form.changed)
    tot.times.ps.form.failed <- sum(num.times.ps.form.failed)
    tot.times.prog.form.changed <- sum(num.times.prog.form.changed)
    tot.times.out.form.changed <- 
        if (is.null(outcome.formula)) NA else sum(num.times.out.form.changed, na.rm= TRUE)

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
    p.value.efron    <- GetPValue(est.TE, est.SE.efron)
    p.value.efron.bc <- GetPValue(est.TE, est.SE.efron.bc)

    # TODO: return NA if no model given
    efron.lm <- 
        EfronEstimates(count.matrix, est.TEs.lm)
    est.SE.lm.efron <- sqrt(efron.lm[1])
    est.SE.lm.efron.bc <- sqrt(efron.lm[2])
    p.value.lm.efron    <- GetPValue(est.TE.lm, est.SE.lm.efron)
    p.value.lm.efron.bc <- GetPValue(est.TE.lm, est.SE.lm.efron.bc)

    # TODO: maybe take these out? on the other hand they make an interesting comparison
    est.SE.naive <- sd(est.TEs, na.rm= TRUE)
    # TODO: return NA if no model given
    est.SE.lm.naive <- sd(est.TEs.lm, na.rm= TRUE)

    list(  
        # vectors of length N, in order of original dataset
        # These are all things that can be calculated w/o BOOM
        treat          = treat,
        logitPS.orig   = logitPS.orig,
        PS.orig        = PS.orig,
        att.wts.orig   = att.wts.orig,
        ate.wts.orig   = ate.wts.orig,
        # from Li and Greene 2013:
        match.wts.orig = match.wts.orig,
        progscore.orig = progscore.orig,

        # scalar results from BOOM w/ NO further covariate adj
        est.TE           = est.TE,
        # not recommended for use; just for comparison
        est.SE.naive     = est.SE.naive,
        est.SE.efron     = est.SE.efron,
        est.SE.efron.bc  = est.SE.efron.bc,
        p.value.efron    = p.value.efron,
        p.value.efron.bc = p.value.efron.bc,

        # scalar results from BOOM with further covariate adj
        est.TE.lm           = est.TE.lm,
        # not recommended for use; just for comparison
        est.SE.lm.naive     = est.SE.lm.naive,
        est.SE.lm.efron     = est.SE.lm.efron,
        est.SE.lm.efron.bc  = est.SE.lm.efron.bc,
        p.value.lm.efron    = p.value.lm.efron,
        p.value.lm.efron.bc = p.value.lm.efron.bc,

        # confidence intervals
        conf.int.efron       = GetConfInt(est.TE, est.SE.efron, conf.level),
        conf.int.efron.bc    = GetConfInt(est.TE, est.SE.efron.bc, conf.level),
        conf.int.lm.efron    = GetConfInt(est.TE.lm, est.SE.lm.efron, conf.level),
        conf.int.lm.efron.bc = GetConfInt(est.TE.lm, est.SE.lm.efron.bc, conf.level),

        # scalars
        n.boot         = n.boot,
        tot.match.errs = tot.match.errs,
        tot.times.ps.form.changed   = tot.times.ps.form.changed, 
        tot.times.ps.form.failed    = tot.times.ps.form.failed, 
        tot.times.prog.form.changed = tot.times.prog.form.changed, 
        tot.times.out.form.changed  = tot.times.out.form.changed, 
        conf.level     = conf.level,
        avg.num.pairs  = mean(num.pairs),

        # vectors of length n.boot
        est.TEs    = est.TEs,
        est.TEs.lm = est.TEs.lm,
        num.pairs  = num.pairs,

        # vectors of length N, in order of original dataset
        logitPS.avg  = logitPS.avg, # may or may not be useful
        count.vector = count.vector,
        BOOM.wts     = BOOM.wts,

        dat= if (return.dat) dat else NA
    )
}

