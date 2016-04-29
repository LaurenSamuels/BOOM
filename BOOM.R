BOOM <- function(dat, n.boot, tx.indicator, outcome, 
    distance.type= "propensity",
    recalc.distance= TRUE,
    propensity.formula = NULL, 
    prognostic.formula = NULL, 
    outcome.formula= NULL, 
    MD.vars= NULL,
    mc.cores= 2, seed= 1235,
    return.dat= TRUE, 
    conf.level= 0.95,
    exact.var.names= NULL, 
    caliper= 0.2, replace= FALSE
    ){
    # Returns a vector of various summary values from the BOOM procedure

    # dat: the original dataset
    # n.boot: number of bootstrap resamples to use
    # tx.indicator: name of the treatment indicator variable in dat (must be 1/0 for now)
    # outcome: name of the continuous outcome variable in dat
    # distance.type: one of "propensity" (propensity score), "prognostic" (prognostic score), or "MD" (Mahalanobis) 
    # recalc.distance (boolean): re-calculate the PS or other distance measure in each resample?
    # propensity.formula: propensity score formula to use w/ lrm
    # prognostic.formula: prognostic score formula. Must not contain factor variables.
    # outcome.formula: optional outcome formula. Currently must not involve interactions with the treatment indicator.
    # mc.cores: number of cores for mclapply()
    # seed: random seed for use by this function. 
    # return.dat (boolean): Return the original dataset?
    # conf.level: level to use for confidence intervals
    # exact.var.names: vector of names of variables on which to match exactly
    # caliper through end of list: as in Matching::Match
    
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
    # TODO: for MD matching, switch to optimal nbp?

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
        # todo: make this nicer. 
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
        # todo: make this nicer. 
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

    # for use w/ non-recalculated distances
    if (!recalc.distance) {
        if (distance.type == "propensity") {
            logitPS.tx.orig   <- logitPS.orig[isTreated]
            logitPS.ctrl.orig <- logitPS.orig[isControl]
        } else if (distance.type == "prognostic") {
            progscore.tx.orig   <- progscore.orig[isTreated]
            progscore.ctrl.orig <- progscore.orig[isControl]
        } else if (distance.type == "MD") {
            # todo. Doing this would involve switching to nbpmatching
            #   for the MD matching 
        }
    }

    # This does not vary from boot to boot (because of the way sample is generated)
    isControl.boot <- c(rep(FALSE, n.treated.orig), rep(TRUE, n.ctrl.orig))

    bootStuff <- mclapply(1:n.boot, function(x) {
        # Modified from Austin & Small (2014) --- they did not condition on
        # observed tx & ctrl group sizes
        tx.sample.indices <- 
            sample(1:n.treated.orig, size= n.treated.orig, replace= TRUE)
        tx.sample <- tx.orig.dat[tx.sample.indices, ]

        ctrl.sample.indices <- 
            sample(1:n.ctrl.orig, size= n.ctrl.orig, replace= TRUE)
        ctrl.sample <- ctrl.orig.dat[ctrl.sample.indices, ]

        # in boot.sample, the first n.treated.orig rows are treated people,
        #   remaining rows are control people
        boot.sample <- rbind(tx.sample, ctrl.sample, make.row.names= FALSE)

        est.mean.tx.tmp <- est.mean.ctrl.tmp <- est.TE.lm.tmp <- NA
        count.vector.tmp <- rep(0, N)
        num.errs.tmp <- 0
        logitPS.ordered.tmp <- rep(NA, N)
        count.vector.matched.tmp <- rep(0, N)

        if (recalc.distance) {
            if (distance.type == "propensity") {
                logitPS <- 
                    GetLogitPS(boot.sample, propensity.formula)
            } else if (distance.type == "prognostic") {
                # the original prognostic formula might not work in this sample
                # TODO: shld probably do something like this for PS.
                # TODO: keep record of this/ report back
                prog.matrix <- model.matrix(prognostic.formula, data= ctrl.sample)
                prog.rank <- rankMatrix(prog.matrix)
                prog.terms.to.remove <- NA
                prog.estimate.anyway <- FALSE
                while ((prog.rank < ncol(prog.matrix)) | prog.estimate.anyway) {
                    #http://stats.stackexchange.com/questions/16327/testing-for-linear-dependence-among-the-columns-of-a-matrix
                    
                    rank.if.removed <- 
                        sapply(1:ncol(prog.matrix), function (x) 
                            rankMatrix(prog.matrix[, -x])
                        )
                    max.rank.if.removed <- max(rank.if.removed)
                    if (all(rank.if.removed == max.rank.if.removed)) {
                        prog.estimate.anyway <- TRUE
                        # TODO: handle this better
                        warning("Problem with prognostic score formula")
                    }
                    prog.problem.indices <- 
                        which(rank.if.removed == max.rank.if.removed)
                    prog.terms.to.remove <- c(prog.terms.to.remove, 
                        setdiff(colnames(prog.matrix)[prog.problem.indices], 
                        "(Intercept)")[1])
                    prog.matrix <- 
                        prog.matrix[, !(colnames(prog.matrix) %in% prog.terms.to.remove)]
                    prog.rank <- rankMatrix(prog.matrix)
                }
                prog.terms.to.remove <- na.omit(prog.terms.to.remove)
                if (length(prog.terms.to.remove) >= 1) {
                    # some of this might be redundant. but I'm worried about (Intercept)
                    prog.term.positions.to.remove <- match(prog.terms.to.remove,
                        attr(prog.form.terms, "term.labels"))
                    prog.tmpterms <- drop.terms(prog.form.terms,
                        dropx= prog.term.positions.to.remove)         
                    prog.form.inboot <- 
                        reformulate(attr(prog.tmpterms, "term.labels"), 
                        response= prog.form.response )
                } else {
                    prog.form.inboot <- prognostic.formula
                }
                progscore <- GetPrognosticScore(boot.sample, 
                    prog.form.inboot, isControl.boot)
            } else if (distance.type == "MD") {
                # todo, if using nbpmatching rather than Match
            }
        } else { # using a fixed distance
            if (distance.type == "propensity") {
                logitPS <- c(logitPS.tx.orig[tx.sample.indices], 
                    logitPS.ctrl.orig[ctrl.sample.indices])
            } else if (distance.type == "prognostic") {
                progscore <- c(progscore.tx.orig[tx.sample.indices], 
                    progscore.ctrl.orig[ctrl.sample.indices])
            } else if (distance.type == "MD") {
                # todo
                # todo, if using nbpmatching rather than Match
            }
        }

        if (distance.type == "propensity") {
            my.X <- logitPS
        } else if (distance.type == "prognostic") {
            my.X <- progscore
        } else if (distance.type == "MD") {
            my.X <- boot.sample[, MD.vars]
        }

        if (is.null(exact.var.names)) {
            my.exact <- FALSE
        } else { # we want to match exactly on some vars
            if (distance.type %in% c("propensity", "prognostic")) {
                my.X <- cbind(my.X, boot.sample[, exact.var.names])
                my.exact <- 
                    c(FALSE, rep(TRUE, length(exact.var.names)))
            } else if (distance.type == "MD") {
                my.exact <- vector(FALSE, length(MD.vars))
                my.exact[MD.vars %in% exact.var.names] <- TRUE
            }
        }
        my.weight <- 
            ifelse(distance.type %in% c("propensity", "prognostic"), 
            1, 2)
        pairIndices <- GetPairs(
            Tr            = boot.sample[[tx.indicator]],
            X             = my.X,
            exact         = my.exact,
            caliper       = caliper,
            replace       = replace,
            Weight        = my.weight
        )
        # Tx indices are in 1st col, ctrl in 2nd col.
        # These indices are indices from boot.sample

        # We are handling PS estimation/matching errors 
        #    by skipping that resample.  
        #    Maybe not the best way, but it rarely happens.
        if(!is.null(pairIndices)){ # TODO: maybe set a min # of matches to proceed?
            #print(nrow(pairIndices)) # TODO: this is useful for setting caliper; maybe find a way to officially incorporate?
            est.mean.tx.tmp <- 
                mean(boot.sample[pairIndices[, 1], outcome])
            est.mean.ctrl.tmp <- 
                mean(boot.sample[pairIndices[, 2], outcome])
            matched.indices <- c(pairIndices[, 1], pairIndices[, 2])

            if (!is.null(outcome.formula)) {
                # modify outcome.formula as necessary to remove factors with only one level
                # See http://stackoverflow.com/questions/18171246/error-in-contrasts-when-defining-a-linear-model-in-r
                lm.vars.to.remove <- 
                    out.form.othervars[sapply(out.form.othervars, 
                    function(x) length(unique(boot.sample[matched.indices, x])) == 1)]
                if (length(lm.vars.to.remove) >= 1) {
                    # TODO: this assumes that the term names are the same as the variable
                    #     names. This is probably a safe assumption for variables w/
                    #     only a few categories, but it won't cover interactions, etc.
                    # Could do processing for rank after this (see prognostic scores)
                    term.positions.to.remove <- match(lm.vars.to.remove,
                        attr(out.form.terms, "term.labels"))
                    tmpterms <- drop.terms(out.form.terms,
                        dropx= term.positions.to.remove)         
                    out.form.inboot <- 
                        reformulate(attr(tmpterms, "term.labels"), 
                        response= out.form.response )
                } else {
                    out.form.inboot <- outcome.formula
                }

                fit <- lm(out.form.inboot, 
                    data= boot.sample[matched.indices, ])
                est.TE.lm.tmp <- coef(fit)[tx.indicator]
            }

            # fill in counts for Efron calculations
            #   (unrelated to matching, but we do not want to do if
            #    matching was unsuccessful)
            tx.orig.ids.easy.insample <- tx.sample.indices
            ctrl.orig.ids.easy.insample <-
                ctrl.orig.ids.easy[ctrl.sample.indices]
            all.orig.ids.easy.insample <- 
                c(tx.orig.ids.easy.insample,
                ctrl.orig.ids.easy.insample)
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
            cat("Hit problem in boot", x, ".\n") 
            num.errs.tmp <- num.errs.tmp + 1
        } # end processing for resamples with errors in PS estimation or matching

        # return from mclapply:
        list(
            # scalars
            est.mean.tx   = est.mean.tx.tmp,
            est.mean.ctrl = est.mean.ctrl.tmp,
            est.TE.lm     = est.TE.lm.tmp,
            num.errs      = num.errs.tmp,

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

    # vector of length n.boot
    num.errs <- do.call(c, lapply(bootStuff, function(x) 
        x[["num.errs"]]))

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
        n.boot        = n.boot,
        tot.errs      = tot.errs,
        conf.level    = conf.level,
        avg.num.pairs = mean(num.pairs),

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

