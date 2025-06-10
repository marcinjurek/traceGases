source("settings.r")

######### initialize #########
## generate grid of pred.locs
grid_oneside <- seq(0, 1, length = round(sqrt(N)))
grid <- expand.grid(grid_oneside, grid_oneside)
locs <- as.matrix(grid) 

## set initial state

dist_mat <- fields::rdist(locs)
Sig0Mat <- GPvecchia::MaternFun(dist_mat, PRIOR_COVPARMS)
Sig0c <- t(chol(Sig0Mat))
x0 <- Sig0c %*% matrix(rnorm(N), ncol = 1)

## define Vecchia aproximation
mra <- GPvecchia::vecchia_specify(locs, M, conditioning = "mra")
#exact <- GPvecchia::vecchia_specify(locs, N - 1, conditioning = "firstm")
lowrank <- GPvecchia::vecchia_specify(locs, M, conditioning = "firstm")
#approximations <- list(lowrank = lowrank, exact = exact, mra = mra)
approximations <- list(mra = mra)
#approximations <- list(exact = exact)


Linvs <- list(truth = solve(t(chol(Sig0Mat))))
for (name in names(approximations)) {
    appr <- approximations[[name]]
    matCov <- GPvecchia::getMatCov(appr, covfun_d_0)
    L <- GPvecchia::createL(appr, matCov)
    Linvs[[name]] <- solve(L, sparse = TRUE)
}


XY <- simulate_xy(x0, evolFun, NULL, FRAC_OBS, LIK_PARMS, T_MAX,
                 sig2 = SIG2, smooth = SMOOTH, range = RANGE, locs = locs)


Corr <- GPvecchia::MaternFun(dist_mat, c(1, COVPARMS[-1]))
Li <- solve(t(chol(Corr)))
sse <- 0
for (t in 2:T_MAX) {
    xs1 <- XY$x[[t - 1]]
    xs2 <- XY$x[[t]]
    w <- xs2 - evolFun(xs1)
    sse <- sse + t(w) %*% t(Li) %*% Li %*% w
}
est <- as.numeric(sse / ((N - 1) * (T_MAX - 1)))
cat(sprintf("%s Estimate from true field %f\n", Sys.time(), est))


par_samples <- list()

for (name in names(approximations)) {

    cat(sprintf("%s Sampling using %s\n", Sys.time(), name))
    app <- approximations[[name]]
    sigmas <- rep(SIG2, NSAMPLES)
    adv_coeffs <- rep(ADV, NSAMPLES)
    diff_coeffs <- rep(DIFF, NSAMPLES)
    covparams <- COVPARMS

    for (sample_no in 2:NSAMPLES) {

        if (sample_no %% 10 == 0) {
            cat(sprintf("%s\t working on sample %d\n", Sys.time(), sample_no))
        }

        a <- adv_coeffs[sample_no - 1]
        d <- diff_coeffs[sample_no - 1]
        evol_fun <- function(x) evolDiff(x, adv = a, diff = d)
        covparams[1] <- sigmas[sample_no - 1]
        
        smoothed <- FFBS(app, XY$y, LIK_PARMS, PRIOR_COVPARMS, covparams, evol_fun)
        states <- smoothed$samples[[1]]
        #states <- XY$x

        sigma2 <- sig2_sample(states, evolFun, ALPHA, BETA, Li)
        sigmas[sample_no] <- sigma2
        sLi <- (1 / sqrt(sigma2)) * Li
        
        a <- adv_sample(states, d, sLi, SIG2_A, MU_A)
        adv_coeffs[sample_no] <- a
        d <- diff_sample(states, a, sLi, SIG2_D, MU_D)
        diff_coeffs[sample_no] <- d

    }

    par_samples[[name]] <- cbind(sigmas, adv_coeffs, diff_coeffs)
}






pdf("param-samples.pdf")
oldpar <- par(mfrow = c(3, 1))
sigma_range <- range(c(SIG2, unlist(lapply(par_samples, function(t) t[, 1]))))
adv_range <- range(c(ADV, unlist(lapply(par_samples, function(t) t[, 2]))))
diff_range <- range(c(DIFF, unlist(lapply(par_samples, function(t) t[, 3]))))

plot(par_samples$mra[, 1], ylim = sigma_range, type = "l", main = "sill", ylab="sill")
lines(par_samples$mra[, 1], col = "red")
lines(par_samples$mra[, 1], col = "blue")
abline(h = SIG2, lty = 2)
legend("topright", c("exact", "mra", "lowrank"), col = c("black", "red", "blue"),
       lty = c(1, 1, 1))

plot(par_samples$mra[, 2], ylim = adv_range, type = "l", main = "advection", ylab="adv")
lines(par_samples$mra[, 2], col = "red")
lines(par_samples$mra[, 2], col = "blue")
abline(h = ADV, lty = 2)
legend("topright", c("exact", "mra", "lowrank"), col = c("black", "red", "blue"),
       lty = c(1, 1, 1))

plot(par_samples$mra[, 3], ylim = diff_range, type = "l", main = "diffusion", ylab="diff")
lines(par_samples$mra[, 3], col = "red")
lines(par_samples$mra[, 3], col = "blue")
abline(h = DIFF, lty = 2)
legend("topright", c("exact", "mra", "lowrank"), col = c("black", "red", "blue"),
       lty = c(1, 1, 1))
par(oldpar)
dev.off()
