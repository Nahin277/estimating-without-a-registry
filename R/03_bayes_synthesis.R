## =============================================================================
## 03_bayes_synthesis.R  --  Bayesian evidence synthesis, written out by hand.
##
## Three sources of information, combined:
##
##   1. Regional prior        what population-based registries in comparable
##                            countries actually observe
##   2. Local likelihood      the Bangladeshi hospital series
##   3. Bias model            the referral tilt, learned from places that have
##                            BOTH a registry and hospital series
##
##   x_j       ~ Binomial(n_j, q_j)
##   logit(q_j) = logit(theta) + delta + u_j
##   u_j       ~ Normal(0, tau^2)
##   logit(theta) ~ Normal(m_reg, s_reg^2)      <- regional prior
##   delta        ~ Normal(mu_delta, sd_delta^2) <- bias model
##   tau          ~ HalfNormal(0.5)
##
## Then the national count:
##   C_<30 = theta * C_all,  with C_all lognormal around GLOBOCAN's 9,480
##
## The sampler is a 40-line random-walk Metropolis. It is deliberately written
## out rather than delegated to Stan or JAGS, so that nothing about the
## inference is hidden behind an API. Base R only.
## =============================================================================

source("02_bounds.R")

set.seed(SEED)

## -----------------------------------------------------------------------------
## Step 1: turn the comparator registries into a prior on logit(theta)
## -----------------------------------------------------------------------------
reg_logit <- logit(registries$share_u30)
m_reg <- mean(reg_logit)
s_reg <- sd(reg_logit)
## Widen it: Bangladesh is not one of these countries, and "exchangeable with"
## is a weaker claim than "identical to". A 1.5x inflation of the spread is a
## deliberate, stated act of humility.
s_reg <- s_reg * 1.5

cat("\n--- Regional prior -----------------------------------------------------\n")
cat(sprintf("Registry shares under 30: %s\n",
            paste(sprintf("%.1f%%", 100 * registries$share_u30), collapse = ", ")))
cat(sprintf("Prior on theta: centred at %.2f%%, 95%% range %.2f%% - %.2f%%\n",
            100 * expit(m_reg),
            100 * expit(m_reg - 1.96 * s_reg),
            100 * expit(m_reg + 1.96 * s_reg)))

## -----------------------------------------------------------------------------
## Step 2: the log posterior
## -----------------------------------------------------------------------------
x <- studies$x; n <- studies$n; K <- nrow(studies)

log_post <- function(par) {
  b  <- par[1]                 # logit(theta)
  d  <- par[2]                 # delta
  lt <- par[3]                 # log(tau)
  u  <- par[4:(3 + K)]
  tau <- exp(lt)

  q <- expit(b + d + u)
  if (any(!is.finite(q)) || any(q <= 0) || any(q >= 1)) return(-Inf)

  ll <- sum(dbinom(x, n, q, log = TRUE))                       # likelihood
  lp <- dnorm(b,  m_reg,    s_reg,    log = TRUE) +            # regional prior
        dnorm(d,  DELTA_MU, DELTA_SD, log = TRUE) +            # bias model
        sum(dnorm(u, 0, tau, log = TRUE)) +                    # study effects
        dnorm(tau, 0, 0.5, log = TRUE) + lt                    # half-normal + Jacobian
  ll + lp
}

## -----------------------------------------------------------------------------
## Step 3: random-walk Metropolis
## -----------------------------------------------------------------------------
metropolis <- function(logf, init, n_iter, step, burn = n_iter / 2) {
  p    <- length(init)
  cur  <- init
  lcur <- logf(cur)
  out  <- matrix(NA_real_, n_iter, p)
  acc  <- 0
  for (i in 1:n_iter) {
    prop  <- cur + rnorm(p, 0, step)
    lprop <- logf(prop)
    if (log(runif(1)) < lprop - lcur) { cur <- prop; lcur <- lprop; acc <- acc + 1 }
    out[i, ] <- cur
  }
  list(draws = out[(burn + 1):n_iter, , drop = FALSE], accept = acc / n_iter)
}

init  <- c(m_reg, DELTA_MU, log(0.25), rep(0, K))
step  <- c(0.12, 0.12, 0.15, rep(0.18, K))

N_ITER <- 400000
chains <- lapply(1:3, function(ch) {
  set.seed(SEED + ch)
  metropolis(log_post, init + rnorm(length(init), 0, 0.1), N_ITER, step)
})

cat(sprintf("\nAcceptance rates: %s\n",
            paste(sprintf("%.2f", sapply(chains, `[[`, "accept")), collapse = " ")))

draws <- do.call(rbind, lapply(chains, `[[`, "draws"))
thin  <- seq(1, nrow(draws), by = 20)
draws <- draws[thin, , drop = FALSE]

theta_post <- expit(draws[, 1])
delta_post <- draws[, 2]
tau_post   <- exp(draws[, 3])

## Crude Gelman-Rubin on logit(theta) across the three chains.
gr <- function(mats, j) {
  m  <- length(mats); nn <- nrow(mats[[1]])
  xs <- sapply(mats, function(M) M[, j])
  W  <- mean(apply(xs, 2, var))
  B  <- nn * var(colMeans(xs))
  sqrt(((nn - 1) / nn * W + B / nn) / W)
}
Rhat_theta <- gr(lapply(chains, `[[`, "draws"), 1)
Rhat_delta <- gr(lapply(chains, `[[`, "draws"), 2)

## -----------------------------------------------------------------------------
## Step 4: propagate into a national count
## -----------------------------------------------------------------------------
M <- length(theta_post)
s_ln    <- sqrt(log(1 + BD_CASES_CV^2))
C_all   <- rlnorm(M, log(BD_CASES_ALL_AGES) - 0.5 * s_ln^2, s_ln)
C_u30   <- theta_post * C_all
rate_post <- C_u30 / BD_WOMEN_20_29 * 1e5

qs <- function(v) quantile(v, c(0.025, 0.5, 0.975))

post <- list(
  theta = qs(theta_post), delta = qs(delta_post), tau = qs(tau_post),
  count = qs(C_u30), rate = qs(rate_post),
  theta_draws = theta_post, count_draws = C_u30, delta_draws = delta_post,
  Rhat_theta = Rhat_theta, Rhat_delta = Rhat_delta
)

cat("\n--- Posterior ----------------------------------------------------------\n")
cat(sprintf("theta (share under 30) : %.2f%%   (95%% CrI %.2f%% - %.2f%%)\n",
            100 * post$theta[2], 100 * post$theta[1], 100 * post$theta[3]))
cat(sprintf("delta (referral tilt)  : OR %.2f  (95%% CrI %.2f - %.2f)\n",
            exp(post$delta[2]), exp(post$delta[1]), exp(post$delta[3])))
cat(sprintf("tau (study spread)     : %.2f    (95%% CrI %.2f - %.2f)\n",
            post$tau[2], post$tau[1], post$tau[3]))
cat(sprintf("National cases < 30    : %.0f     (95%% CrI %.0f - %.0f)\n",
            post$count[2], post$count[1], post$count[3]))
cat(sprintf("Rate per 100,000 (20-29): %.2f   (95%% CrI %.2f - %.2f)\n",
            post$rate[2], post$rate[1], post$rate[3]))
cat(sprintf("Rhat: theta %.3f, delta %.3f\n", Rhat_theta, Rhat_delta))

## -----------------------------------------------------------------------------
## Step 5: the ladder of models -- how much each ingredient moves the answer
## -----------------------------------------------------------------------------
fit_variant <- function(use_regional, use_bias, label) {
  lp <- function(par) {
    b <- par[1]; d <- par[2]; lt <- par[3]; u <- par[4:(3 + K)]; tau <- exp(lt)
    q <- expit(b + d + u)
    if (any(!is.finite(q)) || any(q <= 0) || any(q >= 1)) return(-Inf)
    out <- sum(dbinom(x, n, q, log = TRUE)) +
           sum(dnorm(u, 0, tau, log = TRUE)) +
           dnorm(tau, 0, 0.5, log = TRUE) + lt
    out <- out + if (use_regional) dnorm(b, m_reg, s_reg, log = TRUE)
                 else dnorm(b, 0, 5, log = TRUE)           # vague
    out <- out + if (use_bias) dnorm(d, DELTA_MU, DELTA_SD, log = TRUE)
                 else dnorm(d, 0, 1e-4, log = TRUE)        # delta pinned to 0
    out
  }
  th <- unlist(lapply(1:3, function(ch) {
    set.seed(SEED + 99 + ch)
    r <- metropolis(lp, init + rnorm(length(init), 0, 0.1), N_ITER, step)
    expit(r$draws[seq(1, nrow(r$draws), 20), 1])
  }))
  data.frame(model = label,
             med = median(th), lo = quantile(th, .025), hi = quantile(th, .975),
             cases = median(th) * BD_CASES_ALL_AGES,
             row.names = NULL)
}

ladder <- rbind(
  data.frame(model = "M0  Naive pooled hospital data",
             med = naive$p_re, lo = naive$ci_re[1], hi = naive$ci_re[2],
             cases = naive$p_re * BD_CASES_ALL_AGES),
  fit_variant(FALSE, FALSE, "M1  Hierarchical, no tilt, no regional prior"),
  fit_variant(TRUE,  FALSE, "M2  + regional registry prior"),
  fit_variant(TRUE,  TRUE,  "M3  + referral-tilt correction (full)")
)

cat("\n--- The ladder ---------------------------------------------------------\n")
cat(sprintf("%-46s %8s %18s %8s\n", "model", "share", "95% interval", "cases"))
for (i in 1:nrow(ladder))
  cat(sprintf("%-46s %7.2f%%  %7.2f%% - %5.2f%% %8.0f\n",
              ladder$model[i], 100 * ladder$med[i],
              100 * ladder$lo[i], 100 * ladder$hi[i], ladder$cases[i]))

## -----------------------------------------------------------------------------
## Step 6: leave-one-study-out
## -----------------------------------------------------------------------------
## If deleting one small hospital series moves the national estimate a lot, the
## estimate is fragile and you should say so.
loo <- data.frame(dropped = c("none", studies$study), theta = NA_real_)
for (i in 0:K) {
  keep <- if (i == 0) 1:K else setdiff(1:K, i)
  xk <- x[keep]; nk <- n[keep]; Kk <- length(keep)
  lp <- function(par) {
    b <- par[1]; d <- par[2]; lt <- par[3]; u <- par[4:(3 + Kk)]; tau <- exp(lt)
    q <- expit(b + d + u)
    if (any(!is.finite(q)) || any(q <= 0) || any(q >= 1)) return(-Inf)
    sum(dbinom(xk, nk, q, log = TRUE)) +
      dnorm(b, m_reg, s_reg, log = TRUE) +
      dnorm(d, DELTA_MU, DELTA_SD, log = TRUE) +
      sum(dnorm(u, 0, tau, log = TRUE)) +
      dnorm(tau, 0, 0.5, log = TRUE) + lt
  }
  th <- unlist(lapply(1:2, function(ch) {
    set.seed(SEED + 500 + 10 * i + ch)
    r <- metropolis(lp, c(m_reg, DELTA_MU, log(.25), rep(0, Kk)), N_ITER,
                    c(.12, .12, .15, rep(.18, Kk)))
    expit(r$draws[seq(1, nrow(r$draws), 20), 1])
  }))
  loo$theta[i + 1] <- median(th)
}
loo$shift_pp <- 100 * (loo$theta - loo$theta[1])

cat("\n--- Leave-one-study-out ------------------------------------------------\n")
for (i in 1:nrow(loo))
  cat(sprintf("  drop %-36s theta = %.2f%%  (shift %+.2f pp)\n",
              loo$dropped[i], 100 * loo$theta[i], loo$shift_pp[i]))
cat(sprintf("Largest single-study shift: %.2f percentage points.\n",
            max(abs(loo$shift_pp))))

## -----------------------------------------------------------------------------
## Step 7: prior sensitivity on the tilt
## -----------------------------------------------------------------------------
sens <- do.call(rbind, lapply(c(0, 0.4, DELTA_MU, 1.1, 1.4), function(mu_d) {
  lp <- function(par) {
    b <- par[1]; d <- par[2]; lt <- par[3]; u <- par[4:(3 + K)]; tau <- exp(lt)
    q <- expit(b + d + u)
    if (any(!is.finite(q)) || any(q <= 0) || any(q >= 1)) return(-Inf)
    sum(dbinom(x, n, q, log = TRUE)) +
      dnorm(b, m_reg, s_reg, log = TRUE) +
      dnorm(d, mu_d, DELTA_SD, log = TRUE) +
      sum(dnorm(u, 0, tau, log = TRUE)) +
      dnorm(tau, 0, 0.5, log = TRUE) + lt
  }
  th <- unlist(lapply(1:3, function(ch) {
    set.seed(SEED + 777 + ch)
    r <- metropolis(lp, init + rnorm(length(init), 0, 0.1), N_ITER, step)
    expit(r$draws[seq(1, nrow(r$draws), 20), 1])
  }))
  data.frame(mu_delta = mu_d, or = exp(mu_d), med = median(th),
             lo = quantile(th, .025), hi = quantile(th, .975),
             cases = median(th) * BD_CASES_ALL_AGES, row.names = NULL)
}))

cat("\n--- Prior sensitivity on the referral tilt -----------------------------\n")
cat(sprintf("%10s %10s %18s %8s\n", "prior OR", "share", "95% interval", "cases"))
for (i in 1:nrow(sens))
  cat(sprintf("%10.2f %9.2f%% %7.2f%% - %5.2f%% %8.0f\n",
              sens$or[i], 100 * sens$med[i], 100 * sens$lo[i],
              100 * sens$hi[i], sens$cases[i]))
cat("\nThe answer moves by a factor of ~4 across this range. That is the\n")
cat("honest headline: the tilt assumption, not the data, is in charge.\n")
