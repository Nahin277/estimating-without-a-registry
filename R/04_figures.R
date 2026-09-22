## =============================================================================
## 04_figures.R  --  Four figures. Base R graphics, no packages.
## =============================================================================

## Re-running the MCMC takes ~3 minutes. Cache it so that fiddling with the
## figures is instant. Delete ../results/fit.rds to force a refit.
if (file.exists("../results/fit.rds")) {
  cat("Loading cached fit (delete ../results/fit.rds to refit).\n")
  load("../results/fit.rds")
} else {
  source("03_bayes_synthesis.R")
  dir.create("../results", showWarnings = FALSE, recursive = TRUE)
  save(list = ls(), file = "../results/fit.rds")
}

dir.create("../figures", showWarnings = FALSE, recursive = TRUE)
dir.create("../results", showWarnings = FALSE, recursive = TRUE)

INK   <- "#1B1B1F"
GREY  <- "#8A8A93"
LINE  <- "#D8D8DE"
BLUE  <- "#2E6FB7"
RUST  <- "#C0492B"
TEAL  <- "#3F8F7C"
SAND  <- "#E8B25C"

png_open <- function(f, w = 1800, h = 1150) {
  png(file.path("../figures", f), width = w, height = h, res = 200)
  par(family = "sans", col.axis = INK, col.lab = INK, col.main = INK,
      fg = INK, bg = "white")
}

## -----------------------------------------------------------------------------
## Figure 1 -- the hospital studies, and the gap they leave
## -----------------------------------------------------------------------------
png_open("fig1_what_the_studies_say.png", 1800, 1250)
par(mar = c(5, 15, 4.5, 3))

k    <- nrow(studies)
ys   <- (k + 2):3
p    <- studies$p_hat
plo  <- expit(logit((studies$x + .5)/(studies$n + 1)) - 1.96 * naive$studies_se)
phi  <- expit(logit((studies$x + .5)/(studies$n + 1)) + 1.96 * naive$studies_se)

xlim <- c(0, 0.20)
plot(NA, xlim = xlim, ylim = c(0.2, k + 2.9), axes = FALSE,
     xlab = "", ylab = "")

## registry band
rect(min(registries$share_u30), 0.2, max(registries$share_u30), k + 2.9,
     col = "#EAF1F8", border = NA)
text(mean(range(registries$share_u30)), k + 2.75,
     "what registries elsewhere\nin the region see", col = BLUE, cex = 0.66,
     font = 3)

segments(seq(0, 0.16, 0.04), 0.2, seq(0, 0.16, 0.04), k + 2.4,
         col = LINE, lty = 3)
axis(1, at = seq(0, 0.16, 0.04), labels = paste0(seq(0, 16, 4), "%"),
     col = LINE, cex.axis = 0.8, lwd = 0, lwd.ticks = 1)
axis(2, at = ys, labels = studies$study, las = 1, col = NA, cex.axis = 0.78)

segments(plo, ys, pmin(phi, 0.163), ys, col = GREY, lwd = 2)
arrows(pmin(phi, 0.163)[phi > 0.163], ys[phi > 0.163], 0.168, ys[phi > 0.163],
       length = 0.05, col = GREY, lwd = 2)
points(p, ys, pch = 22, bg = "white", col = INK, cex = sqrt(studies$n) / 5.2,
       lwd = 1.4)
text(0.20, ys, sprintf("%d/%d", studies$x, studies$n),
     cex = 0.68, col = GREY, adj = 1)

## pooled
abline(h = 2.1, col = LINE)
segments(naive$ci_re[1], 1.3, naive$ci_re[2], 1.3, col = RUST, lwd = 3)
points(naive$p_re, 1.3, pch = 23, bg = RUST, col = RUST, cex = 1.5)
text(par("usr")[1] - 0.004, 1.3, "Pooled (random effects)", adj = 1,
     xpd = NA, cex = 0.82, font = 2, col = RUST)

title(main = "Six hospital case series, and the number they seem to give",
      cex.main = 1.05, font.main = 2, adj = 0, line = 2.6)
mtext("Share of breast cancer patients diagnosed before age 30",
      side = 3, adj = 0, line = 1.2, cex = 0.8, col = GREY)
mtext("Box size is proportional to study size. Under-30 counts are illustrative (see 00_inputs.R).",
      side = 1, adj = 0, line = 3.2, cex = 0.66, col = GREY)
dev.off()

## -----------------------------------------------------------------------------
## Figure 2 -- precision is not accuracy
## -----------------------------------------------------------------------------
png_open("fig2_precision_is_not_accuracy.png", 1900, 1000)
par(mfrow = c(1, 2), mar = c(5, 5, 4.5, 1.5))

## left: intervals collapsing around the wrong value
plot(NA, xlim = c(1, length(k_grid)), ylim = c(0.005, 0.07), axes = FALSE,
     xlab = "number of hospital studies pooled", ylab = "estimated share under 30")
abline(h = sim_truth$theta_true, col = TEAL, lwd = 2.5)
abline(h = sim_truth$q_hosp, col = RUST, lwd = 1.2, lty = 2)
segments(1:length(k_grid), sim_conv$lo, 1:length(k_grid), sim_conv$hi,
         col = GREY, lwd = 6, lend = 1)
points(1:length(k_grid), sim_conv$est, pch = 19, col = INK, cex = 0.8)
axis(1, at = 1:length(k_grid), labels = k_grid, col = LINE, cex.axis = 0.78)
axis(2, at = seq(0.01, 0.07, 0.02), labels = paste0(seq(1, 7, 2), "%"),
     las = 1, col = LINE, cex.axis = 0.78)
text(5.2, sim_truth$theta_true, "the truth", col = TEAL, pos = 3,
     cex = 0.75, font = 2)
text(5.2, sim_truth$q_hosp, "what hospitals see", col = RUST,
     pos = 1, cex = 0.75, font = 2)
title("Confidence grows", adj = 0, cex.main = 0.98, line = 1.4)

## right: coverage collapsing
plot(NA, xlim = c(1, length(k_grid)), ylim = c(0, 1), axes = FALSE,
     xlab = "number of hospital studies pooled",
     ylab = "chance the 95% interval covers the truth")
abline(h = 0.95, col = LINE, lty = 2)
lines(1:length(k_grid), sim_cover$coverage, col = RUST, lwd = 3)
points(1:length(k_grid), sim_cover$coverage, pch = 19, col = RUST, cex = 0.9)
axis(1, at = 1:length(k_grid), labels = k_grid, col = LINE, cex.axis = 0.78)
axis(2, at = seq(0, 1, 0.25), labels = paste0(seq(0, 100, 25), "%"),
     las = 1, col = LINE, cex.axis = 0.78)
text(2.2, 0.95, "what it claims", col = GREY, pos = 3, cex = 0.72)
title("Correctness does not", adj = 0, cex.main = 0.98, line = 1.4)

mtext("Adding studies of the same biased kind: the interval shrinks, the target never arrives",
      side = 3, outer = TRUE, line = -2, adj = 0.02, cex = 0.92, font = 2,
      col = INK)
dev.off()

## -----------------------------------------------------------------------------
## Figure 3 -- THE figure. What the data can and cannot pin down.
## -----------------------------------------------------------------------------
png_open("fig3_the_honest_answer.png", 1900, 1200)
par(mar = c(5.5, 5.5, 5, 5.5))

plot(NA, xlim = c(1, 5), ylim = c(0, 0.075), axes = FALSE, xlab = "", ylab = "")

## plausible-tilt band
rect(1.5, 0, 3, 0.075, col = "#F4F0E6", border = NA)
text(2.12, 0.0715, "defensible range for the referral tilt", col = "#9A7B33",
     cex = 0.74, font = 3)

polygon(c(curve_df$or, rev(curve_df$or)), c(curve_df$lo, rev(curve_df$hi)),
        col = "#E4EAF2", border = NA)
lines(curve_df$or, curve_df$theta, col = BLUE, lwd = 3)

axis(1, at = 1:5, labels = paste0(1:5, "x"), col = LINE, cex.axis = 0.82)
axis(2, at = seq(0, 0.07, 0.01), labels = paste0(seq(0, 7, 1), "%"), las = 1,
     col = LINE, cex.axis = 0.8)
axis(4, at = seq(0, 0.07, 0.01),
     labels = format(round(seq(0, 0.07, 0.01) * BD_CASES_ALL_AGES, -1),
                     big.mark = ","),
     las = 1, col = LINE, cex.axis = 0.8)
mtext("How much more likely a young patient is to appear in a hospital case series (odds ratio)",
      side = 1, line = 3, cex = 0.82)
mtext("Share of all breast cancers occurring before 30", side = 2, line = 3.8,
      cex = 0.82)
mtext("Women per year, nationally", side = 4, line = 3.9, cex = 0.82)

## markers
pts <- data.frame(or = c(1, 1.5, exp(DELTA_MU), 3),
                  lab = c("assume no tilt\n(what pooling implies)",
                          "mild tilt", "learned tilt", "strong tilt"))
pts$th <- expit(q_logit - log(pts$or))
points(pts$or, pts$th, pch = 21, bg = "white", col = BLUE, cex = 1.3, lwd = 2)
text(pts$or[1] + 0.06, pts$th[1], pts$lab[1], adj = 0, cex = 0.72, col = RUST,
     font = 2)
text(pts$or[4] + 0.08, pts$th[4], sprintf("%.1f%%  (~%.0f women)",
     100 * pts$th[4], pts$th[4] * BD_CASES_ALL_AGES), adj = 0, cex = 0.72,
     col = INK)
text(pts$or[2] + 0.07, pts$th[2] + 0.003, sprintf("%.1f%%  (~%.0f women)",
     100 * pts$th[2], pts$th[2] * BD_CASES_ALL_AGES), adj = 0, cex = 0.72,
     col = INK)

title("One equation, two unknowns", adj = 0, cex.main = 1.15, line = 3)
mtext("The studies fix where you sit on this curve only once you say how tilted they are. Everything on the blue line fits the data equally well.",
      side = 3, adj = 0, line = 1.3, cex = 0.76, col = GREY)
dev.off()

## -----------------------------------------------------------------------------
## Figure 4 -- what the synthesis gives back
## -----------------------------------------------------------------------------
png_open("fig4_posterior.png", 1900, 1000)
par(mfrow = c(1, 2), mar = c(5, 5, 4.5, 1.5))

h <- hist(100 * post$theta_draws, breaks = 70, plot = FALSE)
plot(h, col = "#DCE6F1", border = "white", main = "", xlab = "share under 30",
     ylab = "", axes = FALSE, xlim = c(0, 8))
axis(1, at = 0:8, labels = paste0(0:8, "%"), col = LINE, cex.axis = 0.78)
abline(v = 100 * naive$p_re, col = RUST, lwd = 2.5, lty = 2)
abline(v = 100 * post$theta[2], col = BLUE, lwd = 2.5)
text(100 * naive$p_re + 0.18, max(h$counts) * 0.94,
     sprintf("naive pooling\n%.1f%%", 100 * naive$p_re), col = RUST, adj = 0,
     cex = 0.75, font = 2)
text(100 * post$theta[2] + 0.55, max(h$counts) * 0.62,
     sprintf("after correction\n%.1f%% (%.1f-%.1f%%)", 100 * post$theta[2],
             100 * post$theta[1], 100 * post$theta[3]),
     col = BLUE, adj = 0, cex = 0.75, font = 2)
title("The share", adj = 0, cex.main = 0.98, line = 1.4)

h2 <- hist(post$count_draws[post$count_draws < 900], breaks = 70, plot = FALSE)
plot(h2, col = "#DFEDE8", border = "white", main = "",
     xlab = "women diagnosed before 30, per year", ylab = "", axes = FALSE)
axis(1, col = LINE, cex.axis = 0.78)
abline(v = post$count[2], col = TEAL, lwd = 2.5)
segments(post$count[1], 0, post$count[3], 0, col = TEAL, lwd = 5, lend = 1)
text(post$count[2] + 15, max(h2$counts) * 0.9,
     sprintf("%.0f women\n(%.0f-%.0f)", post$count[2], post$count[1],
             post$count[3]), col = TEAL, adj = 0, cex = 0.78, font = 2)
title("The count", adj = 0, cex.main = 0.98, line = 1.4)

mtext("After correcting for referral tilt and carrying the uncertainty in GLOBOCAN's national total",
      side = 3, outer = TRUE, line = -2, adj = 0.02, cex = 0.9, font = 2,
      col = INK)
dev.off()

## -----------------------------------------------------------------------------
## Save every number the post quotes
## -----------------------------------------------------------------------------
write.csv(studies,  "../results/table_studies.csv",  row.names = FALSE)
write.csv(band,     "../results/table_bounds.csv",   row.names = FALSE)
write.csv(ladder,   "../results/table_ladder.csv",   row.names = FALSE)
write.csv(sens,     "../results/table_sensitivity.csv", row.names = FALSE)
write.csv(loo,      "../results/table_loo.csv",      row.names = FALSE)
write.csv(sim_cover,"../results/table_coverage.csv", row.names = FALSE)
write.csv(registries, "../results/table_registries.csv", row.names = FALSE)

key <- data.frame(
  quantity = c("GLOBOCAN 2024 all-age cases", "naive pooled share",
               "naive implied count", "posterior share (median)",
               "posterior share 95% CrI low", "posterior share 95% CrI high",
               "posterior count (median)", "posterior count 95% CrI low",
               "posterior count 95% CrI high", "posterior rate per 100k",
               "bound share low (tilt only)", "bound share high (tilt only)",
               "bound share low (tilt+noise)", "bound share high (tilt+noise)",
               "bound count low", "bound count high",
               "posterior tilt OR (median)", "largest LOO shift (pp)",
               "Rhat theta"),
  value = c(BD_CASES_ALL_AGES, naive$p_re, naive$p_re * BD_CASES_ALL_AGES,
            post$theta[2], post$theta[1], post$theta[3],
            post$count[2], post$count[1], post$count[3], post$rate[2],
            HEADLINE$band_lo, HEADLINE$band_hi,
            HEADLINE$share_lo, HEADLINE$share_hi,
            HEADLINE$band_cases_lo, HEADLINE$band_cases_hi,
            exp(post$delta[2]), max(abs(loo$shift_pp)), post$Rhat_theta)
)
write.csv(key, "../results/key_numbers.csv", row.names = FALSE)

cat("\nFigures written to ../figures/, tables to ../results/\n")
print(key, digits = 4)
