sig2_sample <- function(states, evol_fun, alpha, beta, q_mat_inv,
                        num_samples = 1) {

    n <- length(states[[1]])
    t_max <- length(states)
    new_alpha <- alpha + n * (t_max - 1) / 2
    new_beta <- beta
    
    for (t in 2:t_max) {
        w <- states[[t]] - evol_fun(states[[t - 1]])
        q_w <- q_mat_inv %*% w
        new_beta <- new_beta + 0.5 * (t(q_w) %*% q_w)
    }

    cat(sprintf("%s\tExpectation of the posterior of sig2: %f\n", Sys.time(), new_beta / (new_alpha - 1)))
    draw <- invgamma::rinvgamma(num_samples, new_alpha, as.numeric(new_beta) )
    cat(sprintf("%s\tDrawn sig2: %f\n", Sys.time(), draw))
    return(draw)
  
}


## assuming states is a list with sampled states at each time point in its
## own list slot
adv_sample <- function(states, diff, l_mat_inv, sig2_a, mu_a) {

    n <- nrow(states[[1]])
    t_max <- length(states)
    nx <- sqrt(n)
    ny <- nx
    dx <- 1 / nx
    
    omega <- 4 * diff / (dx ** 2)
    c1 <- 1 - omega
    c2 <- 0.25 * omega
    c3 <- 0.25 * omega
    m1 <- advDiffMat(c1, c2, c3, n, nx, ny)
    
    d1 <- - 2
    d2 <- 0
    d3 <- 1
    m2 <- advDiffMat(d1, d2, d3, n, nx, ny)

    s1 <- 0
    s2 <- 0
    q_mat_inv <- l_mat_inv %*% t(l_mat_inv)
    for (t in 2:t_max) {
        m2x <- m2 %*% states[[t - 1]]
        s1 <- s1 + t(m2x) %*% q_mat_inv %*% (states[[t]] - m1 %*% states[[t - 1]])
        s2 <- s2 + 0.5 * t(m2x) %*% q_mat_inv %*% m2x
    }

    sig2 <- 1 / (1 / sig2_a + 2 * as.numeric(s2))
    mu <- sig2 * (mu_a / sig2_a + as.numeric(s1))
    cat(sprintf("%s\tExpectation of the posterior of a = %f\n", Sys.time(), mu))
    
    while (TRUE) {
        adv_prop <- as.numeric(sqrt(sig2) * rnorm(1) + mu)
        stab_cond <- (omega + 2 * adv_prop / dx) > 1
        pos_cond <- adv_prop >= 0
        if (stab_cond && pos_cond) {
            break
        }
    }
    cat(sprintf("%s\tDrawn advection coefficient a = %f\n", Sys.time(), adv_prop))
    adv_prop
}

diff_sample <- function(states, adv, l_mat_inv, sig2_d, mu_d) {

    n <- nrow(states[[1]])
    t_max <- length(states)
    nx <- sqrt(n)
    ny <- nx
    dx <- 1 / nx

    d1 <- 1 - 2 * adv
    d2 <- 0
    d3 <- adv
    m1 <- advDiffMat(d1, d2, d3, n, nx, ny)
    
    c1 <- - 4 / (dx ** 2)
    c2 <- 1 / (dx ** 2)
    c3 <- 1 / (dx ** 2)
    m2 <- advDiffMat(c1, c2, c3, n, nx, ny)
    
    s1 <- 0
    s2 <- 0
    q_mat_inv <- l_mat_inv %*% t(l_mat_inv)
    for (t in 2:t_max) {
        m2x <- m2 %*% states[[t - 1]]
        s1 <- s1 + t(m2x) %*% q_mat_inv %*% (states[[t]] - m1 %*% states[[t - 1]])
        s2 <- s2 + 0.5 * t(m2x) %*% q_mat_inv %*% m2x
    }

    sig2 <- 1 / (1 / sig2_d + 2 * as.numeric(s2))
    mu <- sig2 * (mu_d / sig2_d + as.numeric(s1))

    cat(sprintf("%s\tExpectation of the posterior of d = %f\n", Sys.time(), mu))
    
    while (TRUE) {
        diff_prop <- as.numeric(sqrt(sig2) * rnorm(1) + mu)
        stab_cond <- (4 * diff_prop / (dx ** 2) + 2 * adv / dx) > 1
        pos_cond <- diff_prop >= 0
        if (stab_cond && pos_cond) {
            break
        }
    }
    cat(sprintf("%s\tDrawn diffusion coefficient d = %f\n", Sys.time(), diff_prop))
    diff_prop
}
