## =============================================================================
## 02_bounds.R  --  Partial identification: the honest headline.
##
## We cannot separate "the national share" from "the referral tilt" using the
## hospital data alone. The data pin down their combination and nothing more:
##
##     logit(q_hospital) = logit(theta_national) + delta
##
## One equation, two unknowns. So instead of pretending we can solve it, we
## trace the whole solution curve and report the segment of it that corresponds
## to defensible values of delta.
##
## This is the result the post leads with, because it is the only result here
## that does not depend on a prior nobody can check.
##
## Base R only.
## =============================================================================

source("01_naive_pooling.R")

## The hospital-side quantity, with its sampling uncertainty, from the
## random-effects fit in 01.
q_logit    <- naive$logit_re
q_logit_se <- naive$se_re

## -----------------------------------------------------------------------------
## The identification curve
## -----------------------------------------------------------------------------
delta_grid <- seq(0, log(5), length.out = 200)

theta_of_delta <- function(d) expit(q_logit - d)
theta_lo_of_d  <- function(d) expit(q_logit - d - 1.96 * q_logit_se)
theta_hi_of_d  <- function(d) expit(q_logit - d + 1.96 * q_logit_se)

curve_df <- data.frame(
  delta   = delta_grid,
  or      = exp(delta_grid),
  theta   = theta_of_delta(delta_grid),
  lo      = theta_lo_of_d(delta_grid),
  hi      = theta_hi_of_d(delta_grid)
)
curve_df$cases <- curve_df$theta * BD_CASES_ALL_AGES
curve_df$rate  <- curve_df$cases / BD_WOMEN_20_29 * 1e5

## -----------------------------------------------------------------------------
## The headline band
## -----------------------------------------------------------------------------
## Defensible range for the referral tilt: young patients between 1.5x and 3x
## more likely, in odds, to appear in a readable Bangladeshi case series.
band <- data.frame(
  assumption = c("No referral tilt at all (delta = 0)",
                 "Mild tilt      (OR = 1.5)",
                 "Central tilt   (OR = 2.1, the learned mean)",
                 "Strong tilt    (OR = 3.0)",
                 "Severe tilt    (OR = 4.0)"),
  or         = c(1, 1.5, exp(DELTA_MU), 3, 4)
)
band$delta <- log(band$or)
band$theta <- theta_of_delta(band$delta)
band$lo    <- theta_lo_of_d(band$delta)
band$hi    <- theta_hi_of_d(band$delta)
band$cases <- band$theta * BD_CASES_ALL_AGES
band$rate  <- band$cases / BD_WOMEN_20_29 * 1e5

cat("\n--- Bounds: what the tilt assumption does to the answer ----------------\n")
cat(sprintf("%-46s %8s %10s %10s\n",
            "assumption", "share", "cases/yr", "per 100k"))
for (i in 1:nrow(band))
  cat(sprintf("%-46s %7.2f%% %10.0f %10.2f\n",
              band$assumption[i], 100 * band$theta[i],
              band$cases[i], band$rate[i]))

## The headline interval: sweep delta across [1.5x, 3x] AND carry the sampling
## uncertainty at each end. This is a *bound*, not a credible interval -- it is
## the union of what is compatible with the data across a range of assumptions.
## Two versions, because they answer slightly different questions.
##   (i)  assumption-only band: what the tilt range alone implies, treating the
##        pooled hospital estimate as exact
##   (ii) full band: the same sweep, plus the sampling noise in the studies
band_lo <- theta_of_delta(DELTA_HI)    # strongest defensible tilt
band_hi <- theta_of_delta(DELTA_LO)    # mildest defensible tilt
head_lo <- theta_lo_of_d(DELTA_HI)     # strongest tilt, lower sampling edge
head_hi <- theta_hi_of_d(DELTA_LO)     # mildest tilt,  upper sampling edge

HEADLINE <- list(
  band_lo = band_lo, band_hi = band_hi,
  band_cases_lo = band_lo * BD_CASES_ALL_AGES,
  band_cases_hi = band_hi * BD_CASES_ALL_AGES,
  share_lo = head_lo, share_hi = head_hi,
  cases_lo = head_lo * BD_CASES_ALL_AGES,
  cases_hi = head_hi * BD_CASES_ALL_AGES,
  rate_lo  = head_lo * BD_CASES_ALL_AGES / BD_WOMEN_20_29 * 1e5,
  rate_hi  = head_hi * BD_CASES_ALL_AGES / BD_WOMEN_20_29 * 1e5
)

cat("\n--- HEADLINE -----------------------------------------------------------\n")
cat("Share of Bangladeshi breast cancers diagnosed before 30:\n")
cat(sprintf("   %.1f%% to %.1f%%   (tilt assumption only)\n",
            100 * band_lo, 100 * band_hi))
cat(sprintf("   %.1f%% to %.1f%%   (tilt assumption + study sampling noise)\n",
            100 * head_lo, 100 * head_hi))
cat("National cases per year under 30:\n")
cat(sprintf("   about %.0f to %.0f women   (tilt only)\n",
            HEADLINE$band_cases_lo, HEADLINE$band_cases_hi))
cat(sprintf("   about %.0f to %.0f women   (tilt + sampling noise)\n",
            HEADLINE$cases_lo, HEADLINE$cases_hi))
cat(sprintf("Annual incidence among women aged 20-29:\n"))
cat(sprintf("   about %.1f to %.1f per 100,000  (roughly 1 in %s to 1 in %s)\n",
            HEADLINE$rate_lo, HEADLINE$rate_hi,
            format(round(1e5 / HEADLINE$rate_hi, -3), big.mark = ","),
            format(round(1e5 / HEADLINE$rate_lo, -3), big.mark = ",")))

## -----------------------------------------------------------------------------
## The two-truths calculation
## -----------------------------------------------------------------------------
## Rare and common are not opposites; they are answers to different questions.
cat("\n--- The two truths -----------------------------------------------------\n")
cat(sprintf("An individual woman aged 20-29 : roughly a 1 in %s chance this year.\n",
            format(round(1e5 / mean(c(HEADLINE$rate_lo, HEADLINE$rate_hi)), -3),
                   big.mark = ",")))
cat(sprintf("The country as a whole         : roughly %.0f such women this year.\n",
            mean(c(HEADLINE$cases_lo, HEADLINE$cases_hi))))
cat("Both sentences describe the same number. Only one of them is comforting.\n")
