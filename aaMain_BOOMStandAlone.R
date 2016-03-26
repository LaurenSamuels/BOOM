# actually run the sims.  

rm(list= ls())

# libraries
library(Hmisc)
library(rms)
library(Matching)
library(Rcpp)
library(survey)
library(parallel)

main.dir <- file.path("~", "Box Sync", "Estimator_Paper1", "code")
code.dir <- file.path(main.dir,"standAloneCode")

source(file.path(code.dir, "args.R"))
source(file.path(code.dir, "InvLogit.R"))
source(file.path(code.dir, "GetTxProbs.R"))
source(file.path(code.dir, "MakeDat.R"))
source(file.path(code.dir, "GetLogitPS.R"))
source(file.path(code.dir, "AddY.R"))
source(file.path(code.dir, "GetPairs.R"))
source(file.path(code.dir, "EfronEstimates.R"))
source(file.path(code.dir, "BOOM.R"))
source(file.path(code.dir, "GetLMResults.R"))
source(file.path(code.dir, "PlotWeights.R"))

# beta.low, beta.med, etc. are in args.R
# TODO: deal w/ special seed-setting for mclapply
#set.seed(123)
set.seed(304)
n <- 250
dat <- MakeDat(n, 
    true.avg.TE = 1,
    Beta.0           = -3.05,
    Beta.low         = beta.low,
    Beta.med         = beta.med,
    Beta.high        = beta.high,
    Beta.v.high      = beta.v.high)
# w/ seed 123, 30 treated


bb <- BOOM(
    dat, 
    n.boot = 10000,
    ps.formula= rightPSFormula, # this and the next are in args.R
    lm.formula= rightOutcomeFormula)

PlotWeights(bb)
treated <- bb$treat == 1
sum(treated)
sum(bb$att.wts.orig[treated])
sum(bb$att.wts.orig[!treated])
sum(bb$BOOM.wts[treated])
sum(bb$BOOM.wts[!treated])

isCtrlAndHasHighAttWt <- !treated & bb$att.wts.orig > 0.5
dat[isCtrlAndHasHighAttWt, ]
psvars <- paste0("x.", 4:10)
dat[isCtrlAndHasHighAttWt, psvars]
summary(dat[isCtrlAndHasHighAttWt, psvars])
summary(dat[treated, psvars])
summary(dat[!treated, psvars])
dat <- within(dat, {
    highCtrlAttWt <- isCtrlAndHasHighAttWt
    treat.factor <- factor(treat)
})
ggplot(data= dat, 
    mapping= aes(x= x.6, y= x.7, alpha= highCtrlAttWt, shape= highCtrlAttWt, colour= treat.factor)) +
    scale_alpha_manual("highCtrlAttWt", values= c(0.3, 1)) +
    geom_point()
