# Delivery Time Data: Influence Measures and Outliers (Lecture 11)
# Primary Reference: Montgomery, Peck, and Vining (5th/6th ed.), Chapter 6
# Data: delivery.csv (n = 25, p = 3)
#
# Generates:
#   1. rvl-delivery.pdf      - Residuals vs. Leverage with Cook's distance contours
#   2. cooksd-delivery.pdf   - Index plot of Cook's distances (cutoffs 4/n and 1.0)
#   3. dfbetas-delivery.pdf  - Index plots of DFBETAS for cases (beta1) and distance (beta2)
#   4. covratio-delivery.pdf - Index plot of COVRATIO with precision threshold bounds
#   5. dffits-delivery.pdf   - Index plot of DFFITS (cutoff 2*sqrt(p/n))

library(ggplot2)

# ---- Plot configuration & theme ----------------------------------------------
pt_col      <- "#1b9e77"    # Standard observation color (teal)
pt_flag_col <- "#d95f02"    # Flagged / influential point color (orange-red)
pt_lev_col  <- "#7570b3"    # High leverage point color (purple)
line_col    <- "#d95f02"    # Highlight line color
ref_col     <- "#386cb0"    # Reference threshold color (blue)

fig_w <- 3.4
fig_h <- 2.5

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

# Determine input/output directories
csv_path <- if (file.exists("delivery.csv")) {
  "delivery.csv"
} else if (file.exists("11_influence-outliers/figures/delivery.csv")) {
  "11_influence-outliers/figures/delivery.csv"
} else if (file.exists("figures/delivery.csv")) {
  "figures/delivery.csv"
} else {
  "delivery.csv"
}
fig_dir <- dirname(csv_path)

out_file <- function(filename) file.path(fig_dir, filename)

# ---- 1. Data & Model Fitting -------------------------------------------------
delivery <- read.csv(csv_path)
n <- nrow(delivery)      # 25
fit <- lm(time ~ cases + distance, data = delivery)
p <- length(coef(fit))   # 3 (intercept, cases, distance)


print(summary(fit))

# ---- 2. Diagnostic Calculations from Definitions -----------------------------
e_raw <- unname(resid(fit))
y_hat <- unname(fitted(fit))
h_ii  <- unname(hatvalues(fit))
h_bar <- p / n                          # Average leverage = 3/25 = 0.12
h_cut2 <- 2 * p / n                     # Hoaglin-Welsch cutoff = 0.24
h_cut3 <- 3 * p / n                     # Conservative cutoff = 0.36

msres <- sum(e_raw^2) / (n - p)         # MSRes = 10.624
r_std <- unname(rstandard(fit))         # Internally studentized residuals (r_i)

# Leave-one-out residual variance S^2_(i)
s2_i  <- ((n - p) * msres - e_raw^2 / (1 - h_ii)) / (n - p - 1)
r_stu <- unname(rstudent(fit))          # Studentized deleted residuals / R-student (t_i)

# Bonferroni outlier test cutoff at alpha = 0.05
alpha     <- 0.05
df_t      <- n - p - 1                  # 25 - 3 - 1 = 21
t_bonf    <- qt(1 - alpha / (2 * n), df = df_t)   # qt(1 - 0.001, 21) = 3.5272
outliers  <- which(abs(r_stu) > t_bonf)

# Cook's Distance
cook_d    <- unname(cooks.distance(fit))
d_cut_f   <- qf(0.50, p, n - p)         # F(0.50, 3, 22) = 0.814 ~ 1.0
d_cut_n   <- 4 / n                      # 4/25 = 0.160

# DFBETAS (n x p matrix)
dfb       <- dfbetas(fit)
dfb_cut   <- 2 / sqrt(n)                # 2/sqrt(25) = 0.40

# DFFITS
dff       <- unname(dffits(fit))
dff_cut   <- 2 * sqrt(p / n)            # 2*sqrt(3/25) = 0.6928

# COVRATIO
cov_r     <- unname(covratio(fit))
cov_cut_l <- 1 - 3 * p / n              # 1 - 0.36 = 0.64
cov_cut_u <- 1 + 3 * p / n              # 1 + 0.36 = 1.36

# Verification of Cook's D factorization: D_i = (r_i^2 / p) * (h_ii / (1 - h_ii))
cook_d_fact <- (r_std^2 / p) * (h_ii / (1 - h_ii))
stopifnot(isTRUE(all.equal(cook_d, cook_d_fact)))

# ---- 3. Console Reporting & Summary Table ------------------------------------
diag_df <- data.frame(
  obs      = delivery$obs,
  time     = delivery$time,
  cases    = delivery$cases,
  dist     = delivery$distance,
  h_ii     = round(h_ii, 4),
  e_i      = round(e_raw, 3),
  r_i      = round(r_std, 3),
  t_i      = round(r_stu, 3),
  cooks_d  = round(cook_d, 4),
  dffits   = round(dff, 3),
  dfb_c    = round(dfb[, "cases"], 3),
  dfb_d    = round(dfb[, "distance"], 3),
  covratio = round(cov_r, 3)
)


print(diag_df)

cat("\n--- Screening Guidelines & Cutoffs ---\n")
cat(sprintf("Average leverage h_bar = p/n      : %.3f\n", h_bar))
cat(sprintf("Leverage cutoff 2p/n              : %.3f\n", h_cut2))
cat(sprintf("Leverage cutoff 3p/n              : %.3f\n", h_cut3))
cat(sprintf("Bonferroni outlier t_crit (alpha=0.05): %.4f (df = %d)\n", t_bonf, df_t))
cat(sprintf("Cook's D cutoff 4/n               : %.3f\n", d_cut_n))
cat(sprintf("Cook's D cutoff F_0.50            : %.3f\n", d_cut_f))
cat(sprintf("DFBETAS size-adjusted cutoff      : %.3f\n", dfb_cut))
cat(sprintf("DFFITS cutoff 2*sqrt(p/n)         : %.3f\n", dff_cut))
cat(sprintf("COVRATIO bounds [1 - 3p/n, 1 + 3p/n]: [%.3f, %.3f]\n\n", cov_cut_l, cov_cut_u))

cat("--- Key Observations Profile ---\n")
cat("Observation 9:\n")
cat(sprintf("  Cases = %d, Distance = %d ft, Time = %.2f min\n",
            delivery$cases[9], delivery$distance[9], delivery$time[9]))
cat(sprintf("  Leverage h_99  = %.4f  (exceeds 2p/n = %.2f and 3p/n = %.2f)\n",
            h_ii[9], h_cut2, h_cut3))
cat(sprintf("  Studentized r_9 = %+.3f,  R-student t_9 = %+.3f  (|t_9| > %.3f => Outlier!)\n",
            r_std[9], r_stu[9], t_bonf))
cat(sprintf("  Cook's D_9     = %.4f  (>> 1.0 and >> 4/n => Dominant influence)\n",
            cook_d[9]))
cat(sprintf("  DFFITS_9       = %+.3f  (exceeds cutoff %.3f)\n",
            dff[9], dff_cut))
cat(sprintf("  DFBETAS_9      : Intercept = %+.3f, Cases = %+.3f, Distance = %+.3f\n",
            dfb[9, 1], dfb[9, 2], dfb[9, 3]))
cat(sprintf("  COVRATIO_9     = %.4f  (< %.2f => Degrades overall estimation precision)\n\n",
            cov_r[9], cov_cut_l))

cat("Observation 22:\n")
cat(sprintf("  Cases = %d, Distance = %d ft, Time = %.2f min\n",
            delivery$cases[22], delivery$distance[22], delivery$time[22]))
cat(sprintf("  Leverage h_22  = %.4f  (exceeds 2p/n = %.2f)\n",
            h_ii[22], h_cut2))
cat(sprintf("  Studentized r_22 = %+.3f, R-student t_22 = %+.3f  (not an outlier)\n",
            r_std[22], r_stu[22]))
cat(sprintf("  Cook's D_22    = %.4f  (exceeds 4/n = %.3f)\n",
            cook_d[22], d_cut_n))
cat(sprintf("  DFFITS_22      = %+.3f  (exceeds cutoff %.3f)\n",
            dff[22], dff_cut))
cat(sprintf("  DFBETAS_22     : Intercept = %+.3f, Cases = %+.3f, Distance = %+.3f\n",
            dfb[22, 1], dfb[22, 2], dfb[22, 3]))
cat(sprintf("  COVRATIO_22    = %.4f  (> %.2f => High leverage enhances precision)\n\n",
            cov_r[22], cov_cut_u))

# ---- 4. Side-by-Side Model Comparison ----------------------------------------
fit_no9    <- lm(time ~ cases + distance, data = delivery[-9, ])
fit_no9_22 <- lm(time ~ cases + distance, data = delivery[-c(9, 22), ])

mod_comp <- data.frame(
  Metric = c("Intercept (b0)", "Cases (b1)", "Distance (b2)",
             "R-squared", "Adj R-squared", "MSRes", "Residual SE", "n"),
  Full_Model = c(
    sprintf("%.3f (%.3f)", coef(fit)[1], coef(summary(fit))[1, 2]),
    sprintf("%.3f (%.3f)", coef(fit)[2], coef(summary(fit))[2, 2]),
    sprintf("%.5f (%.5f)", coef(fit)[3], coef(summary(fit))[3, 2]),
    sprintf("%.4f", summary(fit)$r.squared),
    sprintf("%.4f", summary(fit)$adj.r.squared),
    sprintf("%.3f", msres),
    sprintf("%.3f", summary(fit)$sigma),
    "25"
  ),
  Without_Obs9 = c(
    sprintf("%.3f (%.3f)", coef(fit_no9)[1], coef(summary(fit_no9))[1, 2]),
    sprintf("%.3f (%.3f)", coef(fit_no9)[2], coef(summary(fit_no9))[2, 2]),
    sprintf("%.5f (%.5f)", coef(fit_no9)[3], coef(summary(fit_no9))[3, 2]),
    sprintf("%.4f", summary(fit_no9)$r.squared),
    sprintf("%.4f", summary(fit_no9)$adj.r.squared),
    sprintf("%.3f", sum(resid(fit_no9)^2) / fit_no9$df.residual),
    sprintf("%.3f", summary(fit_no9)$sigma),
    "24"
  ),
  Without_Obs9_22 = c(
    sprintf("%.3f (%.3f)", coef(fit_no9_22)[1], coef(summary(fit_no9_22))[1, 2]),
    sprintf("%.3f (%.3f)", coef(fit_no9_22)[2], coef(summary(fit_no9_22))[2, 2]),
    sprintf("%.5f (%.5f)", coef(fit_no9_22)[3], coef(summary(fit_no9_22))[3, 2]),
    sprintf("%.4f", summary(fit_no9_22)$r.squared),
    sprintf("%.4f", summary(fit_no9_22)$adj.r.squared),
    sprintf("%.3f", sum(resid(fit_no9_22)^2) / fit_no9_22$df.residual),
    sprintf("%.3f", summary(fit_no9_22)$sigma),
    "23"
  )
)


print(mod_comp, row.names = FALSE)


# ---- FIGURE 1: Residuals vs Leverage with Cook's Contours (`rvl-delivery.pdf`) ----
# Generate smooth grid for Cook's distance contours: r = +/- sqrt(p * D * (1 - h) / h)
y_max_curve <- 4.15
contour_data <- data.frame()

for (D_val in c(0.5, 1.0)) {
  # Smallest h where r <= y_max_curve: p*D*(1-h)/h <= y_max^2 => h_min = p*D / (y_max^2 + p*D)
  h_min <- (p * D_val) / (y_max_curve^2 + p * D_val)
  h_seq <- seq(h_min, 0.53, length.out = 200)
  r_pos <- sqrt(p * D_val * (1 - h_seq) / h_seq)
  contour_data <- rbind(
    contour_data,
    data.frame(h = h_seq, r = r_pos, D = factor(D_val), side = "pos"),
    data.frame(h = h_seq, r = -r_pos, D = factor(D_val), side = "neg")
  )
}
contour_data <- subset(contour_data, r >= -3.3 & r <= 4.15)

# Point categories for coloring
point_class <- rep("Typical", n)
point_class[cook_d > 1.0] <- "Cook's D > 1.0"
point_class[cook_d <= 1.0 & h_ii > h_cut2] <- "High Leverage"
delivery$class <- factor(point_class, levels = c("Typical", "High Leverage", "Cook's D > 1.0"))
delivery$h_ii  <- h_ii
delivery$r_std <- r_std

p_rvl <- ggplot() +
  # Safe band / reference lines
  geom_hline(yintercept = 0, linetype = "solid", colour = "grey75", linewidth = 0.5) +
  geom_hline(yintercept = c(-2, 2), linetype = "dotted", colour = "grey60", linewidth = 0.5) +
  geom_vline(xintercept = h_cut2, linetype = "dashed", colour = ref_col, linewidth = 0.6) +
  # Cook's distance contours
  geom_line(
    data = contour_data,
    aes(x = h, y = r, group = interaction(D, side)),
    colour = "#e41a1c",
    linetype = "dashed",
    linewidth = 0.65,
    alpha = 0.80
  ) +
  # Contour text annotations, placed just past the right-hand ends of the curves
  annotate("text", x = 0.537, y = c(1, -1) * sqrt(p * 0.5 * (1 - 0.53) / 0.53),
           label = "D == 0.5", parse = TRUE, hjust = 0, size = 3.0, colour = "#e41a1c") +
  annotate("text", x = 0.537, y = c(1, -1) * sqrt(p * 1.0 * (1 - 0.53) / 0.53),
           label = "D == 1", parse = TRUE, hjust = 0, size = 3.0, colour = "#e41a1c") +
  # Leverage cutoff annotation (top of the plot, clear of the contours)
  annotate("text", x = h_cut2 + 0.008, y = 4.05, label = "2*p/n == 0.24",
           parse = TRUE, hjust = 0, size = 3.0, colour = ref_col) +
  # Observations
  geom_point(
    data = delivery,
    aes(x = h_ii, y = r_std, colour = class, shape = class, size = class)
  ) +
  # Observation labels for key points
  annotate("text", x = h_ii[9] - 0.015, y = r_std[9],
           label = "Obs 9\n(D = 3.42)", hjust = 1, vjust = 0.4, size = 3.2,
           fontface = "bold", colour = pt_flag_col) +
  # Obs 22 sits on the D = 0.5 contour, so its label goes above-left, in empty space
  annotate("text", x = h_ii[22] - 0.01, y = r_std[22] + 0.25,
           label = "Obs 22\n(D = 0.45)", hjust = 1, vjust = 0, size = 3.2,
           fontface = "bold", colour = pt_lev_col) +
  annotate("text", x = h_ii[10] + 0.012, y = r_std[10],
           label = "10", hjust = 0, vjust = 0.3, size = 2.8, colour = "grey35") +
  annotate("text", x = h_ii[16] + 0.012, y = r_std[16],
           label = "16", hjust = 0, vjust = 0.3, size = 2.8, colour = "grey35") +
  scale_colour_manual(
    name = "Diagnostic Status:",
    values = c("Typical" = pt_col, "High Leverage" = pt_lev_col, "Cook's D > 1.0" = pt_flag_col)
  ) +
  scale_shape_manual(
    name = "Diagnostic Status:",
    values = c("Typical" = 16, "High Leverage" = 17, "Cook's D > 1.0" = 18)
  ) +
  scale_size_manual(
    name = "Diagnostic Status:",
    values = c("Typical" = 2.0, "High Leverage" = 3.0, "Cook's D > 1.0" = 3.8)
  ) +
  scale_x_continuous(
    name = expression(Leverage ~ (italic(h)[italic(ii)])),
    limits = c(0.02, 0.60),
    breaks = seq(0, 0.5, by = 0.1)
  ) +
  scale_y_continuous(
    name = expression(Studentized ~ Residual ~ (italic(r)[italic(i)])),
    limits = c(-3.4, 4.2),
    breaks = seq(-3, 4, by = 1)
  ) +
  theme_lecture() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8.5),
    legend.margin = margin(t = -5, b = 2)
  )

ggsave(out_file("rvl-delivery.pdf"), p_rvl, bg = "transparent",
       device = cairo_pdf, width = 5.2, height = 3.8)


# ---- FIGURE 2: Cook's Distance Index Plot (`cooksd-delivery.pdf`) -------------
delivery$cook_d <- cook_d

p_cook <- ggplot(delivery, aes(x = obs, y = cook_d)) +
  # Horizontal cutoffs
  geom_hline(yintercept = d_cut_n, linetype = "dashed", colour = ref_col, linewidth = 0.6) +
  geom_hline(yintercept = 1.0, linetype = "dashed", colour = pt_flag_col, linewidth = 0.6) +
  # Lollipop stems
  geom_segment(aes(x = obs, xend = obs, y = 0, yend = cook_d),
               colour = "grey55", linewidth = 0.65) +
  # Points
  geom_point(aes(colour = cook_d > 1.0, size = cook_d > 1.0)) +
  # Cutoff labels
  # Cutoff labels placed where no stem or other label is nearby
  annotate("text", x = 14, y = d_cut_n + 0.10, label = "4/n == 0.16",
           parse = TRUE, hjust = 0.5, size = 2.9, colour = ref_col) +
  annotate("text", x = 25.6, y = 1.0 + 0.12, label = "D == 1",
           parse = TRUE, hjust = 1, size = 2.9, colour = pt_flag_col) +
  # Observation labels
  annotate("text", x = 9, y = cook_d[9] + 0.18,
           label = "Obs 9 (D = 3.42)", hjust = 0.5, size = 3.2,
           fontface = "bold", colour = pt_flag_col) +
  annotate("text", x = 22, y = cook_d[22] + 0.16,
           label = "Obs 22 (D = 0.45)", hjust = 0.5, size = 3.0,
           fontface = "bold", colour = pt_lev_col) +
  scale_colour_manual(values = c("FALSE" = pt_col, "TRUE" = pt_flag_col), guide = "none") +
  scale_size_manual(values = c("FALSE" = 2.0, "TRUE" = 3.5), guide = "none") +
  scale_x_continuous(name = "Observation Index (i)", breaks = seq(1, 25, by = 2)) +
  scale_y_continuous(
    name = expression(Cook * "'" * s ~ Distance ~ (italic(D)[italic(i)])),
    limits = c(0, 3.8),
    breaks = seq(0, 3.5, by = 0.5)
  ) +
  theme_lecture()

ggsave(out_file("cooksd-delivery.pdf"), p_cook, bg = "transparent",
       device = cairo_pdf, width = 4.8, height = 3.4)


# ---- FIGURE 3: DFBETAS Index Plots for beta1 & beta2 (`dfbetas-delivery.pdf`) --
dfb_long <- rbind(
  data.frame(
    obs = delivery$obs,
    dfbeta = dfb[, "cases"],
    term = "beta[1] ~ (Cases)",
    variable = "Cases"
  ),
  data.frame(
    obs = delivery$obs,
    dfbeta = dfb[, "distance"],
    term = "beta[2] ~ (Distance)",
    variable = "Distance"
  )
)

dfb_long$term <- factor(dfb_long$term, levels = c("beta[1] ~ (Cases)", "beta[2] ~ (Distance)"))
dfb_flagged   <- subset(dfb_long, abs(dfbeta) > dfb_cut)

p_dfb <- ggplot(dfb_long, aes(x = obs, y = dfbeta)) +
  # Reference zero line
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.5) +
  # Cutoff lines (+/- 2/sqrt(n))
  geom_hline(yintercept = c(-dfb_cut, dfb_cut), linetype = "dashed",
             colour = ref_col, linewidth = 0.55) +
  # Absolute cutoff lines (+/- 1.0)
  geom_hline(yintercept = c(-1.0, 1.0), linetype = "dotted",
             colour = pt_flag_col, linewidth = 0.55) +
  # Stems
  geom_segment(aes(x = obs, xend = obs, y = 0, yend = dfbeta),
               colour = "grey60", linewidth = 0.6) +
  # Points
  geom_point(aes(colour = abs(dfbeta) > dfb_cut), size = 1.9) +
  # Labels for flagged observations
  geom_text(
    data = subset(dfb_flagged, obs %in% c(1, 9, 22, 24)),
    aes(label = obs, vjust = ifelse(dfbeta > 0, -0.6, 1.4)),
    size = 2.8, fontface = "bold", colour = "grey20"
  ) +
  # Facet by coefficient
  facet_wrap(~ term, nrow = 1, labeller = label_parsed) +
  scale_colour_manual(values = c("FALSE" = pt_col, "TRUE" = pt_flag_col), guide = "none") +
  scale_x_continuous(name = "Observation Index (i)", breaks = seq(2, 24, by = 4)) +
  scale_y_continuous(
    name = expression(DFBETAS[italic(j) * "," * italic(i)]),
    limits = c(-1.3, 1.8),
    breaks = seq(-1.0, 1.5, by = 0.5)
  ) +
  theme_lecture()

ggsave(out_file("dfbetas-delivery.pdf"), p_dfb, bg = "transparent",
       device = cairo_pdf, width = 2 * fig_w, height = 1.25 * fig_h)
# cat("Generated: dfbetas-delivery.pdf\n")

# ---- FIGURE 4: COVRATIO Index Plot (`covratio-delivery.pdf`) ------------------
delivery$cov_r <- cov_r

p_cov <- ggplot(delivery, aes(x = obs, y = cov_r)) +
  # Safe band where observation does not excessively alter precision
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = cov_cut_l, ymax = cov_cut_u,
           fill = "grey90", alpha = 0.45) +
  # Neutral line at 1.0
  geom_hline(yintercept = 1.0, colour = "grey50", linewidth = 0.5) +
  # Cutoff bounds 1 +/- 3p/n
  geom_hline(yintercept = c(cov_cut_l, cov_cut_u), linetype = "dashed",
             colour = ref_col, linewidth = 0.6) +
  # Lollipop stems
  geom_segment(aes(x = obs, xend = obs, y = 1.0, yend = cov_r),
               colour = "grey60", linewidth = 0.6) +
  # Points
  geom_point(aes(colour = cov_r < cov_cut_l | cov_r > cov_cut_u), size = 2.1) +
  # Threshold text annotations
  annotate("text", x = 1.0, y = cov_cut_u + 0.05,
           label = "1 + 3*p/n == 1.36 ~ (Enhances ~ Precision)",
           parse = TRUE, hjust = 0, size = 2.8, colour = ref_col) +
  annotate("text", x = 1.0, y = cov_cut_l - 0.06,
           label = "1 - 3*p/n == 0.64 ~ (Degrades ~ Precision)",
           parse = TRUE, hjust = 0, size = 2.8, colour = pt_flag_col) +
  # Observation labels
  # Two-line labels anchored by their edge (vjust) so they clear the point and the cutoff line
  annotate("text", x = 9, y = cov_r[9] - 0.05,
           label = "Obs 9 (0.34)\nDegrades precision", hjust = 0.5, vjust = 1, size = 2.9,
           fontface = "bold", colour = pt_flag_col) +
  annotate("text", x = 22, y = cov_r[22] + 0.05,
           label = "Obs 22 (1.40)\nEnhances precision", hjust = 0.5, vjust = 0, size = 2.9,
           fontface = "bold", colour = pt_lev_col) +
  annotate("text", x = 16, y = cov_r[16] + 0.06,
           label = "16", hjust = 0.5, size = 2.8, colour = "grey30") +
  annotate("text", x = 20, y = cov_r[20] - 0.06,
           label = "20", hjust = 0.5, size = 2.8, colour = "grey30") +
  scale_colour_manual(values = c("FALSE" = pt_col, "TRUE" = pt_flag_col), guide = "none") +
  scale_x_continuous(name = "Observation Index (i)", breaks = seq(1, 25, by = 2)) +
  scale_y_continuous(
    name = expression(COVRATIO[italic(i)]),
    limits = c(0.10, 1.68),
    breaks = seq(0.2, 1.6, by = 0.2)
  ) +
  theme_lecture()

ggsave(out_file("covratio-delivery.pdf"), p_cov, bg = "transparent",
       device = cairo_pdf, width = 5.0, height = 3.4)


# ---- FIGURE 5 (BONUS): DFFITS Index Plot (`dffits-delivery.pdf`) ---------------
delivery$dff <- dff

p_dff <- ggplot(delivery, aes(x = obs, y = dff)) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.5) +
  geom_hline(yintercept = c(-dff_cut, dff_cut), linetype = "dashed",
             colour = ref_col, linewidth = 0.6) +
  geom_segment(aes(x = obs, xend = obs, y = 0, yend = dff),
               colour = "grey60", linewidth = 0.6) +
  geom_point(aes(colour = abs(dff) > dff_cut), size = 2.1) +
  # Cutoff labels placed where no point or other label is nearby
  annotate("text", x = 25.6, y = dff_cut + 0.25,
           label = "2*sqrt(p/n) == 0.69", parse = TRUE,
           hjust = 1, size = 2.9, colour = ref_col) +
  annotate("text", x = 13, y = -dff_cut - 0.25,
           label = "-2*sqrt(p/n) == -0.69", parse = TRUE,
           hjust = 0.5, size = 2.9, colour = ref_col) +
  annotate("text", x = 9, y = dff[9] + 0.35,
           label = "Obs 9 (+4.30)", hjust = 0.5, size = 3.0,
           fontface = "bold", colour = pt_flag_col) +
  annotate("text", x = 22, y = dff[22] - 0.35,
           label = "Obs 22 (-1.20)", hjust = 0.5, size = 3.0,
           fontface = "bold", colour = pt_lev_col) +
  annotate("text", x = 20, y = dff[20] - 0.35,
           label = "20", hjust = 0.5, size = 2.8, colour = "grey35") +
  annotate("text", x = 1, y = dff[1] - 0.35,
           label = "1", hjust = 0.5, size = 2.8, colour = "grey35") +
  scale_colour_manual(values = c("FALSE" = pt_col, "TRUE" = pt_flag_col), guide = "none") +
  scale_x_continuous(name = "Observation Index (i)", breaks = seq(1, 25, by = 2)) +
  scale_y_continuous(
    name = expression(DFFITS[italic(i)]),
    limits = c(-1.8, 5.0),
    breaks = seq(-1, 5, by = 1)
  ) +
  theme_lecture()

ggsave(out_file("dffits-delivery.pdf"), p_dff, bg = "transparent",
       device = cairo_pdf, width = 4.8, height = 3.4)
