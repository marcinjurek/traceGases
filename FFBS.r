FFBS <- function(approx, obs, lik.params, prior_covparams, Qcovparms,
                evolFun, Num_samples = 1, verbose = FALSE) {

    samples <- list()
    n <- nrow(approx$locsord)
    Tmax <- length(obs)

    revord <- order(approx$ord)
    locs <- approx$locsord[revord, ]

    dist_mat <- fields::rdist(locs)
    Q <- GPvecchia::MaternFun(dist_mat, Qcovparms)
    Qc <- t(chol(Q))

    Sig0 <- GPvecchia::MaternFun(dist_mat, prior_covparams)
    Sig0c <- t(chol(Sig0))
    
    #Sig0Model <- RMwhittle(nu = prior_covparams[3], scale = prior_covparams[2], var = prior_covparams[1])
    
    for (sample_no in 1:Num_samples) {

        if (verbose) {
            cat(sprintf("%s FFBS: Working on sample no %d\n", Sys.time(), sample_no))
            cat(sprintf("%s FFBS: Simulating data\n", Sys.time()))
        }

        x0 <- Sig0c %*% matrix(rnorm(n), ncol=1)

        XYplus <- simulate_xy(x0, evolFun, NULL, obs, lik.params, Tmax, sig2 = Qcovparms[1],
                                  smooth = Qcovparms[3], range = Qcovparms[2], locs = locs)

        Y <- mapply("-", obs, XYplus$y, SIMPLIFY = FALSE)

        results <- smoothedMeans(approx, Y, lik.params, prior_covparams, Qcovparms, evolFun)

        smeans <- results[[ "means" ]]
        sample <- mapply("+", XYplus$x, smeans, SIMPLIFY = FALSE)

        if (verbose) {
            cat(sprintf("%s FFBS: done working on sample %d\n", Sys.time(), sample_no))
        }
        samples[[sample_no]] <- sample
        
    }

    return( list(filteringResults = results, samples = samples) )
}
