## =============================================================================
## 01_naive_pooling.R  --  The tempting shortcut, and a demonstration that it
##                          answers a different question than the one asked.
##
## Two things happen here:
##   (a) we pool the hospital studies the ordinary way (fixed effect and
##       DerSimonian-Laird random effects, on the logit scale), and
##   (b) we show by simulation that adding more hospital studies makes the
##       answer more *precise* without making it more *correct*.
##
## Base R only.
## =============================================================================

source("00_inputs.R")

## -----------------------------------------------------------------------------
## (a) Pooling
## -----------------------------------------------------------------------------

## Simple pooled proportion: everyone in one bucket.
pooled_x <- sum(studies$x)
pooled_n <- sum(studies$n)
p_pooled <- pooled_x / pooled_n
se_pooled_logit <- sqrt(1 / pooled_x + 1 / (pooled_n - pooled_x))
ci_pooled <- expit(logit(p_pooled) + c(-1.96, 1.96) * se_pooled_logit)

## Per-study logit and its standard error, with the usual 0.5 continuity
## correction guard for small counts.
lo <- logit((studies$x + 0.5) / (studies$n + 1))
se <- sqrt(1 / (studies$x + 0.5) + 1 / (studies$n - studies$x + 0.5))
w  <- 1 / se^2

## Fixed-effect (inverse variance) estimate.
mu_fe    <- sum(w * lo) / sum(w)
se_fe    <- sqrt(1 / sum(w))
ci_fe    <- expit(mu_fe + c(-1.96, 1.96) * se_fe)

## DerSimonian-Laird random effects.
Q      <- sum(w * (lo - mu_fe)^2)
df     <- length(lo) - 1
C      <- sum(w) - sum(w^2) / sum(w)
tau2   <- max(0, (Q - df) / C)
w_re   <- 1 / (se^2 + tau2)
mu_re  <- sum(w_re * lo) / sum(w_re)
se_re  <- sqrt(1 / sum(w_re))
ci_re  <- expit(mu_re + c(-1.96, 1.96) * se_re)
I2     <- max(0, (Q - df) / Q)

naive <- list(
  p_pooled = p_pooled, ci_pooled = ci_pooled,
  p_fe     = expit(mu_fe), ci_fe = ci_fe,
  p_re     = expit(mu_re), ci_re = ci_re,
  logit_re = mu_re, se_re = se_re,
  tau2     = tau2, Q = Q, I2 = I2,
  studies_lo = lo, studies_se = se
)

cat("\n--- Naive pooling of the hospital studies -----------------------------\n")
cat(sprintf("Simple pooled      : %.2f%%  (95%% CI %.2f%% - %.2f%%)\n",
            100 * p_pooled, 100 * ci_pooled[1], 100 * ci_pooled[2]))
cat(sprintf("Fixed effect       : %.2f%%  (95%% CI %.2f%% - %.2f%%)\n",
            100 * naive$p_fe, 100 * ci_fe[1], 100 * ci_fe[2]))
cat(sprintf("Random effects (DL): %.2f%%  (95%% CI %.2f%% - %.2f%%)\n",
            100 * naive$p_re, 100 * ci_re[1], 100 * ci_re[2]))
cat(sprintf("Heterogeneity      : Q = %.2f on %d df, tau^2 = %.3f, I^2 = %.0f%%\n",
            Q, df, tau2, 100 * I2))
cat(sprintf("Naive national count under 30: %.0f cases  (= %.2f%% of %d)\n",
            naive$p_re * BD_CASES_ALL_AGES, 100 * naive$p_re, BD_CASES_ALL_AGES))

## This number is not wrong. It is an honest answer to the question
##   "among breast cancer patients in the kind of hospital that publishes
##    case series in Bangladesh, what share are under 30?"
## It is not an answer to
##   "among all Bangladeshi women who develop breast cancer, what share are
##    under 30?"
## Those two questions have the same words in them and different answers.

## -----------------------------------------------------------------------------
## (b) Precision is not accuracy
## -----------------------------------------------------------------------------
## Suppose the TRUE national share is 2%, and every hospital series over-samples
## young patients with the same referral tilt delta. Collect more and more
## studies. Watch the confidence interval shrink -- around the wrong number.

set.seed(SEED)

theta_true  <- 0.02
delta_true  <- DELTA_MU
q_hosp      <- expit(logit(theta_true) + delta_true)   # what hospitals see
tau_sim     <- 0.30                                    # study-to-study noise
n_per_study <- 150

k_grid  <- c(2, 5, 10, 20, 50, 100, 200, 500)
n_rep   <- 400
conv <- data.frame(k = k_grid, est = NA_real_, lo = NA_real_, hi = NA_real_)

for (i in seq_along(k_grid)) {
  k <- k_grid[i]
  ests <- numeric(n_rep); los <- numeric(n_rep); his <- numeric(n_rep)
  for (r in 1:n_rep) {
    qj <- expit(logit(q_hosp) + rnorm(k, 0, tau_sim))
    xj <- rbinom(k, n_per_study, qj)
    X  <- sum(xj); N <- k * n_per_study
    ph <- X / N
    s  <- sqrt(1 / max(X, 0.5) + 1 / max(N - X, 0.5))
    ests[r] <- ph
    los[r]  <- expit(logit(max(ph, 1e-6)) - 1.96 * s)
    his[r]  <- expit(logit(max(ph, 1e-6)) + 1.96 * s)
  }
  conv$est[i] <- mean(ests); conv$lo[i] <- mean(los); conv$hi[i] <- mean(his)
}

cat("\n--- Precision without accuracy ----------------------------------------\n")
cat(sprintf("True national share  : %.2f%%\n", 100 * theta_true))
cat(sprintf("What hospitals see   : %.2f%%\n", 100 * q_hosp))
cat("\n  studies    estimate     95% CI width\n")
for (i in seq_along(k_grid))
  cat(sprintf("  %6d      %.2f%%       %.2f pp\n",
              conv$k[i], 100 * conv$est[i], 100 * (conv$hi[i] - conv$lo[i])))
cat("\nThe interval collapses. It never moves toward the truth.\n")
cat("More data of the same kind buys confidence, not correctness.\n")

## -----------------------------------------------------------------------------
## Coverage: how often does the naive 95% interval contain the truth?
## -----------------------------------------------------------------------------
cover <- sapply(k_grid, function(k) {
  hit <- replicate(n_rep, {
    qj <- expit(logit(q_hosp) + rnorm(k, 0, tau_sim))
    xj <- rbinom(k, n_per_study, qj)
    X  <- sum(xj); N <- k * n_per_study; ph <- X / N
    s  <- sqrt(1 / max(X, 0.5) + 1 / max(N - X, 0.5))
    l  <- expit(logit(max(ph, 1e-6)) - 1.96 * s)
    u  <- expit(logit(max(ph, 1e-6)) + 1.96 * s)
    (theta_true >= l) && (theta_true <= u)
  })
  mean(hit)
})
cat("\nCoverage of the nominal 95% interval for the TRUE national share:\n")
cat(paste0(sprintf("  k=%-4d %.0f%%\n", k_grid, 100 * cover), collapse = ""))

sim_conv <- conv
sim_cover <- data.frame(k = k_grid, coverage = cover)
sim_truth <- list(theta_true = theta_true, q_hosp = q_hosp,
                  delta_true = delta_true)
