# Simulations for Model Selection (Lecture 12)
# Primary Reference: Montgomery, Peck, and Vining (5th/6th ed.), Chapter 10
# Secondary: Freedman, D. A. (1983), A note on screening regression equations,
#            The American Statistician 37(2), 152-155.
#
# Generates:
#   1. train-test-error.pdf  - (a) training vs test error along a nested sequence of models
#   2. consistency-sim.pdf   - (b) P(correct / overfit / underfit) vs n for AIC, HQIC, BIC, Cp
#   3. freedman-paradox.pdf  - (c) "significant" regressors after screening pure noise
#
# All three use set.seed(2026). Run time: well under a minute.

library(ggplot2)


pt_col      <- "#1b9e77"    # Teal
pt_flag_col <- "#d95f02"    # Orange-red
pt_lev_col  <- "#7570b3"    # Purple
ref_col     <- "#386cb0"    # Blue

theme_lecture <- function() {
  theme_minimal(base_size = 11) +
    theme(
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background  = element_rect(fill = "transparent", colour = NA),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.box.background = element_rect(fill = "transparent", colour = NA),
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = 10, face = "bold", colour = "grey15"),
      strip.background = element_rect(fill = "#f0f0f0", colour = NA),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9)
    )
}

# Output directory (run from repo root, lecture folder, or figures/)
fig_dir <- if (file.exists("cement.csv")) {
  "."
} else if (dir.exists("12_model-selection/figures")) {
  "12_model-selection/figures"
} else if (dir.exists("figures")) {
  "figures"
} else {
  "."
}
out_file <- function(filename) file.path(fig_dir, filename)

# Residual sum of squares of y on [1, X] (fast, no formula overhead)
ss_res <- function(X, y) sum(.lm.fit(cbind(1, X), y)$residuals^2)



# (a) Training vs test error along a nested sequence of models

# Fixed design: n = 30, K = 10 candidate regressors ~ N(0, 1).
# True model: y = 1 + 1.0 x1 + 0.8 x2 + 0.5 x3 + e, e ~ N(0, 1)  -> true size p = 4.
# Models: intercept only, {x1}, {x1, x2}, ..., {x1, ..., x10}  (p = 1, ..., 11).
# Training error = SSRes / n;  test error = mean squared error on fresh y at the same X.
set.seed(2026)
n_a <- 30; K_a <- 10; sigma_a <- 1; reps_a <- 2000
X_a <- matrix(rnorm(n_a * K_a), n_a, K_a)
beta_a <- c(1, 0.8, 0.5, rep(0, K_a - 3))
mu_a <- 1 + drop(X_a %*% beta_a)

err <- array(NA_real_, c(reps_a, K_a + 1, 2))
for (r in seq_len(reps_a)) {
  y     <- mu_a + rnorm(n_a, sd = sigma_a)
  y_new <- mu_a + rnorm(n_a, sd = sigma_a)
  for (j in 0:K_a) {
    Xj  <- cbind(1, X_a[, seq_len(j), drop = FALSE])
    fit <- .lm.fit(Xj, y)
    yhat <- y - fit$residuals
    err[r, j + 1, 1] <- mean(fit$residuals^2)
    err[r, j + 1, 2] <- mean((y_new - yhat)^2)
  }
}
tt <- data.frame(
  p = rep(1:(K_a + 1), 2),
  err = c(colMeans(err[, , 1]), colMeans(err[, , 2])),
  type = rep(c("Training error (SSRes / n)", "Test error (new y, same X)"), each = K_a + 1)
)
tt$type <- factor(tt$type, levels = c("Test error (new y, same X)", "Training error (SSRes / n)"))
print(reshape(tt, idvar = "p", timevar = "type", direction = "wide"), digits = 3, row.names = FALSE)
# Training error falls steadily (2.81 -> 0.64); test error bottoms out at p = 4 (the true size)
# at 1.14, close to sigma^2 (1 + p/n) = 1.13, and rises again (1.38 at p = 11) as useless
# regressors add variance.
# Underfitting at p = 1..3 costs far more than overfitting at p = 5..11.

tt_min <- tt[tt$type == "Test error (new y, same X)", ]
tt_min <- tt_min[which.min(tt_min$err), ]

p_tt <- ggplot(tt, aes(x = p, y = err, colour = type)) +
  geom_hline(yintercept = sigma_a^2, colour = "grey50", linetype = "dotted", linewidth = 0.5) +
  annotate("text", x = 11, y = sigma_a^2 - 0.07, label = "sigma^2 == 1", parse = TRUE,
           hjust = 1, size = 2.8, colour = "grey35") +
  geom_vline(xintercept = 4, colour = ref_col, linetype = "dashed", linewidth = 0.5) +
  annotate("text", x = 4.2, y = 2.9, label = "true model (p = 4)", hjust = 0, size = 2.8,
           colour = ref_col) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  geom_point(data = tt_min, size = 3.4, shape = 21, fill = NA, stroke = 1, colour = pt_flag_col) +
  scale_colour_manual(values = c("Training error (SSRes / n)" = pt_col,
                                 "Test error (new y, same X)" = pt_flag_col), name = NULL) +
  scale_x_continuous(name = expression("Number of parameters" ~ italic(p) ~
                                         "(regressors entered in order" ~ italic(x)[1] * "," ~
                                         italic(x)[2] * ", ...)"),
                     breaks = 1:11) +
  scale_y_continuous(name = "Average squared error", limits = c(0.5, 3.0)) +
  theme_lecture() +
  theme(legend.position = c(0.68, 0.58), legend.text = element_text(size = 8.5),
        legend.key.height = unit(0.4, "cm"))

ggsave(out_file("train-test-error.pdf"), p_tt, bg = "transparent",
       device = cairo_pdf, width = 4.8, height = 3.4)



# (b) Consistency vs efficiency: how often does each criterion find the true model?

# K = 5 independent N(0, 1) candidates; true model {x1, x2} with beta = (1, 0.5), sigma = 1.
# All 2^5 = 32 subsets are scored by AIC, HQIC, BIC (MPV form) and Cp (sigma^2 from full model).
set.seed(2026)
K_b <- 5; reps_b <- 1000
n_vals <- c(15, 25, 50, 100, 200, 500, 1000)
beta_b <- c(1, 0.5, 0, 0, 0)
true_set <- c(TRUE, TRUE, FALSE, FALSE, FALSE)

inc <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), K_b)))   # 32 x 5 inclusion matrix
p_sub <- rowSums(inc) + 1
is_correct <- apply(inc, 1, function(s) all(s == true_set))
is_overfit <- apply(inc, 1, function(s) all(s[true_set]) && any(s[!true_set]))
is_underfit <- !is_correct & !is_overfit                        # misses at least one true regressor

crits <- c("AIC", "HQIC", "BIC", "Cp")
outcome <- array(0, c(length(n_vals), length(crits), 3),
                 dimnames = list(n_vals, crits, c("correct", "overfit", "underfit")))

for (i in seq_along(n_vals)) {
  n <- n_vals[i]
  for (r in seq_len(reps_b)) {
    X <- matrix(rnorm(n * K_b), n, K_b)
    y <- drop(X %*% beta_b) + rnorm(n)
    SS <- apply(inc, 1, function(s) ss_res(X[, s, drop = FALSE], y))
    s2_full <- SS[nrow(inc)] / (n - K_b - 1)
    score <- cbind(
      AIC  = n * log(SS / n) + 2 * p_sub,
      HQIC = n * log(SS / n) + 2 * p_sub * log(log(n)),
      BIC  = n * log(SS / n) + p_sub * log(n),
      Cp   = SS / s2_full - n + 2 * p_sub
    )
    pick <- apply(score, 2, which.min)
    outcome[i, , "correct"]  <- outcome[i, , "correct"]  + is_correct[pick]
    outcome[i, , "overfit"]  <- outcome[i, , "overfit"]  + is_overfit[pick]
    outcome[i, , "underfit"] <- outcome[i, , "underfit"] + is_underfit[pick]
  }
}
outcome <- outcome / reps_b
print(round(outcome[, , "correct"], 3))
#        n=15  n=50  n=100  n=1000
# AIC    0.30  0.53  0.59   0.60
# HQIC   0.30  0.65  0.75   0.87
# BIC    0.33  0.74  0.89   0.98
# Cp     0.34  0.57  0.60   0.60
# BIC climbs towards 1 (consistent); AIC and Cp level off at about 0.6 = (1 - 0.157)^3: each of
# the three noise regressors still enters with probability about P(chi^2_1 > 2) = 0.157.
# HQIC sits in between and converges slowly. At n = 15 BIC underfits most often (0.43).

cons <- as.data.frame.table(outcome, responseName = "prob")
names(cons) <- c("n", "crit", "outcome", "prob")
cons$n <- as.numeric(as.character(cons$n))
cons$crit <- factor(cons$crit, levels = crits)
cons$outcome <- factor(cons$outcome, levels = c("correct", "overfit", "underfit"),
                       labels = c("Correct model", "Overfit (superset)", "Underfit (misses x1 or x2)"))

crit_cols <- c(AIC = pt_col, HQIC = pt_lev_col, BIC = pt_flag_col, Cp = ref_col)

p_cons <- ggplot(cons, aes(x = n, y = prob, colour = crit, shape = crit)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.9) +
  facet_wrap(~ outcome, nrow = 1) +
  scale_colour_manual(values = crit_cols, name = NULL,
                      labels = c("AIC", "HQIC", "BIC", expression(italic(C)[italic(p)]))) +
  scale_shape_manual(values = c(AIC = 16, HQIC = 17, BIC = 15, Cp = 18), name = NULL,
                     labels = c("AIC", "HQIC", "BIC", expression(italic(C)[italic(p)]))) +
  scale_x_log10(name = expression("Sample size" ~ italic(n) ~ "(log scale)"),
                breaks = n_vals, labels = n_vals) +
  scale_y_continuous(name = "Proportion of 1000 replications", limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.2)) +
  theme_lecture() +
  theme(legend.position = "bottom", legend.margin = margin(t = -5),
        axis.text.x = element_text(size = 7.5, angle = 45, hjust = 1))

ggsave(out_file("consistency-sim.pdf"), p_cons, bg = "transparent",
       device = cairo_pdf, width = 7.2, height = 3.2)



# (c) Freedman's paradox: screening pure noise produces "significant" regressors

# n = 100, 50 regressors and y all independent N(0, 1): every true coefficient is 0.
# Stage 1: fit all 50, keep regressors with p < 0.25.  Stage 2: refit on the kept ones and
# read the usual t-tests and overall F-test as if no selection had happened.
set.seed(2026)
n_c <- 100; K_c <- 50; reps_c <- 1000; alpha_screen <- 0.25

freedman_once <- function() {
  X <- matrix(rnorm(n_c * K_c), n_c, K_c, dimnames = list(NULL, paste0("z", 1:K_c)))
  y <- rnorm(n_c)
  d <- data.frame(y = y, X)
  s1 <- summary(lm(y ~ ., data = d))$coefficients[-1, ]
  keep <- rownames(s1)[s1[, 4] < alpha_screen]
  if (!length(keep)) return(list(kept = 0, sig = 0, F_p = 1, R2 = 0, fit = NULL))
  fit2 <- lm(reformulate(keep, "y"), data = d)
  s2 <- summary(fit2)
  f  <- s2$fstatistic
  list(kept = length(keep),
       sig  = sum(s2$coefficients[-1, 4] < 0.05),
       F_p  = unname(pf(f[1], f[2], f[3], lower.tail = FALSE)),
       R2   = s2$r.squared,
       fit  = fit2)
}

first <- freedman_once()           # one replication to show on a slide
cat("\nOne replication: kept", first$kept, "of 50; significant at 5%:", first$sig,
    "; R^2 =", round(first$R2, 3), "; overall F p-value =", signif(first$F_p, 3), "\n")
print(summary(first$fit))

sims <- replicate(reps_c - 1, freedman_once()[c("kept", "sig", "F_p", "R2")], simplify = FALSE)
sims <- rbind(data.frame(kept = first$kept, sig = first$sig, F_p = first$F_p, R2 = first$R2),
              do.call(rbind, lapply(sims, as.data.frame)))
fr_summary <- c(mean_kept = mean(sims$kept), mean_sig = mean(sims$sig),
                P_F_sig = mean(sims$F_p < 0.05), mean_R2 = mean(sims$R2))
print(round(fr_summary, 3))
# Without selection we would expect 50 x 0.05 = 2.5 false positives among 50 tests.
# After screening, the stage-2 model keeps 12.3 regressors on average, 5.5 of them look
# "significant", the overall F-test rejects in 94% of replications, and mean R^2 = 0.31.

fr_counts <- as.data.frame(table(factor(sims$sig, levels = 0:max(sims$sig))))
names(fr_counts) <- c("sig", "count")
fr_counts$sig <- as.integer(as.character(fr_counts$sig))

p_fr <- ggplot(fr_counts, aes(x = sig, y = count / reps_c)) +
  geom_col(fill = pt_col, width = 0.75) +
  geom_vline(xintercept = mean(sims$sig), colour = pt_flag_col, linetype = "dashed",
             linewidth = 0.6) +
  geom_vline(xintercept = K_c * 0.05, colour = "grey45", linetype = "dotted", linewidth = 0.6) +
  annotate("text", x = K_c * 0.05 - 0.3, y = max(fr_counts$count / reps_c) * 1.02,
           hjust = 1, size = 2.7, colour = "grey35",
           label = "no selection:
50 x 0.05 = 2.5") +
  annotate("text", x = mean(sims$sig) + 0.4, y = max(fr_counts$count / reps_c) * 1.02,
           hjust = 0, size = 2.9, colour = pt_flag_col,
           label = sprintf("mean = %.1f", mean(sims$sig))) +
  annotate("text", x = max(fr_counts$sig), y = max(fr_counts$count / reps_c) * 0.80,
           hjust = 1, vjust = 1, size = 2.8, colour = "grey25", lineheight = 1.05,
           label = sprintf(paste0("n = 100, 50 pure-noise regressors\n",
                                  "screen at p < 0.25, then refit\n",
                                  "regressors kept: %.1f on average\n",
                                  "overall F significant: %.0f%% of runs"),
                           fr_summary["mean_kept"], 100 * fr_summary["P_F_sig"])) +
  scale_x_continuous(name = "Stage-2 regressors with p-value < 0.05",
                     breaks = 0:max(fr_counts$sig)) +
  scale_y_continuous(name = "Proportion of replications",
                     expand = expansion(mult = c(0, 0.08))) +
  theme_lecture() +
  theme(panel.grid.major.x = element_blank())

ggsave(out_file("freedman-paradox.pdf"), p_fr, bg = "transparent",
       device = cairo_pdf, width = 4.8, height = 3.4)

cat("\nFigures written to", normalizePath(fig_dir), "\n")
