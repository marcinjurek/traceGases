source("../advection.r")
source("../simulate-data.r")
source("../smoother.r")
source("../FFBS.r")
source("../filter.r")
source("params-sampler.r")
library(invgamma)
library(Matrix)

######### set parameters #########
#set.seed(1988)
N <- 34**2
M <- 50
NSAMPLES <- 50
FRAC_OBS <- 0.5
T_MAX <- 30
# 34 x 34 parameters
ADV <- 0.02
DIFF <- 0.00004
evolFun <- function(x) evolDiff(x, adv = ADV, diff = DIFF)



## covariance parameters
SIG_02 <- 1
SIG2 <- 0.02
RANGE <- .15
SMOOTH <- 0.5
COVPARMS <- c(SIG2, RANGE, SMOOTH)
PRIOR_COVPARMS <- c(SIG_02, RANGE, SMOOTH)
covfun <- function(locs) GPvecchia::MaternFun(fields::rdist(locs), COVPARMS)
covfun_d_0 <- function(D) GPvecchia::MaternFun(D, PRIOR_COVPARMS)

## likelihood settings
ME_VAR <- 1e-6
DATA_MODEL <- "gauss"  
LIK_PARMS <- list(data.model = DATA_MODEL, sigma = sqrt(ME_VAR))
ALPHA <- 0.001
BETA <- 0.001
SIG2_A <- 1e-6
MU_A <- ADV
SIG2_D <- 1e-10
MU_D <- DIFF
