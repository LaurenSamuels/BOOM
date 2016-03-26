# Values used by the main RunSims function and the functions it calls

# N: number of subjects per simulated dataset
N <- 1000

# beta.0.treat.vec: The intercepts to use in the treatment-selection model. 
# This will determine the prevalence of treatment in the simulated dataset. 
beta.0.treat.vec<- c(-3.95, -3.05, -2.07) # to get (.05, .1, .2) treated

# Values of coefficients for treatment and outcome generation
beta.low        <- 0.25
beta.med        <- 0.50
beta.high       <- 0.75
beta.v.high     <- 0.90

rightOutcomeFormula <- y ~ 
    treat + 
    pol(x.4, 2) + 
    x.5 * x.6 + 
    pol(x.7, 2) + 
    x.8 + x.9 + x.10

# for the "right" PS formula we are doing as Austin & Small did,
#   using the variables that affect outcome
rightPSFormula <- treat ~ 
    pol(x.4, 2) + 
    x.5 * x.6 + 
    pol(x.7, 2) + 
    x.8 + x.9 + x.10

# for the "wrong" models we leave out squared terms and interactions
wrongOutcomeFormula <- y ~ 
    treat + 
    x.4 + x.5 + x.6 + x.7 + x.8 + x.9 + x.10

wrongPSFormula <- treat ~ 
    x.4 + x.5 + x.6 + x.7 + x.8 + x.9 + x.10
