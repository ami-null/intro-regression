# Hald Cement Data: Model Selection (Lecture 12)
# Primary Reference: Montgomery, Peck, and Vining (5th/6th ed.), Chapter 10
# Data: cement.csv (n = 13, K = 4 candidate regressors; also used in Lecture 10)
#   y  = heat evolved (calories per gram of cement)
#   x1 = tricalcium aluminate (%), x2 = tricalcium silicate (%),
#   x3 = tetracalcium alumino ferrite (%), x4 = dicalcium silicate (%)
#
# Notation: p = number of parameters INCLUDING the intercept.
# Information criteria in the MPV "regression form" (same as extractAIC / step):
#   AIC  = n ln(SSRes/n) + 2p
#   BIC  = n ln(SSRes/n) + p ln(n)
#   HQIC = n ln(SSRes/n) + 2p ln(ln(n))
# R's AIC()/BIC() add n ln(2 pi) + n + 2 (and count sigma^2); rankings are identical.
#
# Generates:
#   1. cement-allsubsets.csv           - all 16 subsets with every criterion (for LaTeX tables)
#   2. cp-plot-cement.pdf              - Mallows' Cp vs p with the Cp = p reference line
#   3. criteria-best-subset-cement.pdf - best subset of each size under six criteria
#   4. ic-allsubsets-cement.pdf        - AIC / BIC / HQIC differences for all 15 subsets
# Prints: forward / backward / stepwise partial-F traces, step() results,
#         Akaike weights, BIC posterior probabilities, competing fits and VIFs.

library(ggplot2)
library(patchwork)


pt_col      <- "#1b9e77"    # Standard point colour (teal)
pt_flag_col <- "#d95f02"    # Highlighted point colour (orange-red)
pt_lev_col  <- "#7570b3"    # Secondary highlight (purple)
ref_col     <- "#386cb0"    # Reference line colour (blue)

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

# Determine input/output directories (run from repo root, lecture folder, or figures/)
csv_path <- if (file.exists("cement.csv")) {
  "cement.csv"
} else if (file.exists("12_model-selection/figures/cement.csv")) {
  "12_model-selection/figures/cement.csv"
} else if (file.exists("figures/cement.csv")) {
  "figures/cement.csv"
} else {
  "cement.csv"
}
fig_dir <- dirname(csv_path)
out_file <- function(filename) file.path(fig_dir, filename)

# plotmath label for a subset such as "x1,x2,x4" -> list(x[1], x[2], x[4])
lab_expr <- function(model) {
  vapply(model, function(m) {
    if (m == "(none)") return("'(none)'")
    idx <- sub("x", "", strsplit(m, ",")[[1]])
    paste0("list(", paste0("italic(x)[", idx, "]", collapse = ", "), ")")
  }, character(1), USE.NAMES = FALSE)
}

# ---- 1. Data ----------------------------------------------------------------------
cement <- read.csv(csv_path)
n  <- nrow(cement)                                 # 13
xv <- paste0("x", 1:4)
K  <- length(xv)                                   # 4

full   <- lm(y ~ x1 + x2 + x3 + x4, data = cement)
sigma2 <- summary(full)$sigma^2                    # 5.983 (full-model MSRes, used in Cp)
SST    <- sum((cement$y - mean(cement$y))^2)       # 2715.76

round(cor(cement[, c("y", xv)]), 3)                # r(x2, x4) = -0.973, r(x1, x3) = -0.824

# ---- 2. All possible regressions (2^K = 16 subsets) ----------------------------------
subsets <- c(list(character(0)),
             unlist(lapply(1:K, function(k) combn(xv, k, simplify = FALSE)), recursive = FALSE))

fit_subset <- function(vars) {
  lm(if (length(vars)) reformulate(vars, "y") else y ~ 1, data = cement)
}

allsub <- do.call(rbind, lapply(subsets, function(vars) {
  m  <- fit_subset(vars)
  p  <- length(vars) + 1
  SS <- sum(resid(m)^2)
  h  <- hatvalues(m)
  k  <- p + 1                                      # parameters counted by AIC(): beta and sigma^2
  data.frame(
    model = if (length(vars)) paste(vars, collapse = ",") else "(none)",
    p     = p,
    SSRes = SS,
    R2    = 1 - SS / SST,
    R2adj = 1 - (SS / (n - p)) / (SST / (n - 1)),
    MSRes = SS / (n - p),
    Cp    = SS / sigma2 - n + 2 * p,
    PRESS = sum((resid(m) / (1 - h))^2),
    AIC   = n * log(SS / n) + 2 * p,
    BIC   = n * log(SS / n) + p * log(n),
    HQIC  = n * log(SS / n) + 2 * p * log(log(n)),
    AIC_R = AIC(m),
    BIC_R = BIC(m),
    AICc  = AIC(m) + 2 * k * (k + 1) / (n - k - 1)
  )
}))

print(allsub, digits = 4, row.names = FALSE)
# Key rows:            p  SSRes   R2adj  MSRes   Cp    PRESS   AIC    BIC    HQIC   AICc
#   x1,x2              3  57.90  0.9744  5.790  2.68   93.88  25.42  27.11  25.07  69.31
#   x1,x4              3  74.76  0.9670  7.476  5.50  121.22  28.74  30.44  28.39  72.63
#   x1,x2,x3           4  48.11  0.9764  5.346  3.04   90.00  25.01  27.27  24.55  72.48
#   x1,x2,x4           4  47.97  0.9764  5.330  3.02   85.35  24.97  27.23  24.51  72.44
#   x1,x3,x4           4  50.84  0.9750  5.648  3.50   94.54  25.73  27.99  25.26  73.19
#   x1,x2,x3,x4        5  47.86  0.9736  5.983  5.00  110.35  26.94  29.77  26.36  79.84

stopifnot(isTRUE(all.equal(allsub$Cp[allsub$p == K + 1], K + 1)))   # Cp = K + 1 for the full model
stopifnot(isTRUE(all.equal(allsub$AIC[16], extractAIC(full)[2])))   # MPV form = extractAIC()

best_by <- function(col, fun = which.min) allsub$model[fun(allsub[[col]])]
winners <- c(R2adj = best_by("R2adj", which.max), MSRes = best_by("MSRes"),
             Cp = best_by("Cp"), PRESS = best_by("PRESS"), AIC = best_by("AIC"),
             BIC = best_by("BIC"), HQIC = best_by("HQIC"), AICc = best_by("AICc"))
print(winners)
# R2adj, MSRes, PRESS, AIC, HQIC -> x1,x2,x4      Cp, BIC, AICc -> x1,x2

csv_out <- allsub[, c("model", "p", "SSRes", "R2", "R2adj", "MSRes", "Cp", "PRESS",
                      "AIC", "BIC", "HQIC", "AIC_R", "BIC_R", "AICc")]
num_cols <- sapply(csv_out, is.numeric)
csv_out[num_cols] <- lapply(csv_out[num_cols], round, digits = 4)
write.csv(csv_out, out_file("cement-allsubsets.csv"), row.names = FALSE)

# ---- 3. Sequential algorithms with partial F ------------------------------------------
# Partial F for adding one regressor = t^2 of its coefficient in the larger model (Lecture 09)
partial_F <- function(small_vars, big_vars) {
  a <- anova(fit_subset(small_vars), fit_subset(big_vars))
  a$F[2]
}

forward_selection <- function(F_in = 4, F_in_fun = NULL) {
  cur <- character(0)
  repeat {
    cand <- setdiff(xv, cur)
    if (!length(cand)) break
    Fs <- sapply(cand, function(v) partial_F(cur, c(cur, v)))
    df_res <- n - (length(cur) + 2)               # residual df of the larger model
    thr <- if (is.null(F_in_fun)) F_in else F_in_fun(df_res)
    cat(sprintf("  in: {%s}  candidates: %s  (threshold %.2f)\n",
                paste(cur, collapse = ","),
                paste(sprintf("%s=%.2f", names(Fs), Fs), collapse = "  "), thr))
    if (max(Fs) <= thr) break
    cur <- c(cur, names(which.max(Fs)))
    cat("    -> enter", tail(cur, 1), "\n")
  }
  cur
}

backward_elimination <- function(F_out = 4) {
  cur <- xv
  repeat {
    if (!length(cur)) break
    Fs <- sapply(cur, function(v) partial_F(setdiff(cur, v), cur))
    cat(sprintf("  in: {%s}  partial F: %s\n", paste(cur, collapse = ","),
                paste(sprintf("%s=%.2f", names(Fs), Fs), collapse = "  ")))
    if (min(Fs) >= F_out) break
    cat("    -> remove", names(which.min(Fs)), "\n")
    cur <- setdiff(cur, names(which.min(Fs)))
  }
  cur
}

stepwise_selection <- function(F_in = 4, F_out = 4, max_steps = 20) {
  stopifnot(F_in >= F_out)
  cur <- character(0)
  for (step_i in seq_len(max_steps)) {
    cand <- setdiff(xv, cur)
    entered <- FALSE
    if (length(cand)) {
      Fs <- sapply(cand, function(v) partial_F(cur, c(cur, v)))
      cat(sprintf("  in: {%s}  F-to-enter: %s\n", paste(cur, collapse = ","),
                  paste(sprintf("%s=%.2f", names(Fs), Fs), collapse = "  ")))
      if (max(Fs) > F_in) {
        cur <- c(cur, names(which.max(Fs)))
        entered <- TRUE
        cat("    -> enter", tail(cur, 1), "\n")
      }
    }
    if (length(cur) > 1) {
      Fo <- sapply(cur, function(v) partial_F(setdiff(cur, v), cur))
      cat(sprintf("  in: {%s}  F-to-remove: %s\n", paste(cur, collapse = ","),
                  paste(sprintf("%s=%.2f", names(Fo), Fo), collapse = "  ")))
      if (min(Fo) < F_out) {
        cat("    -> remove", names(which.min(Fo)), "\n")
        cur <- setdiff(cur, names(which.min(Fo)))
        next
      }
    }
    if (!entered) break
  }
  cur
}

cat("\n---- Forward selection (F_in = 4) ----\n")
fwd <- forward_selection(F_in = 4)
# x4 (22.80) -> x1 (108.22) -> x2 (5.03) -> stop (x3: 0.02).  Final {x4, x1, x2}

cat("\n---- Forward selection (exact 5% test, F_in = F_{0.05, 1, df}) ----\n")
fwd_exact <- forward_selection(F_in_fun = function(df) qf(0.95, 1, df))
# Step 3: x2 has F = 5.03 < F_{0.05,1,9} = 5.12 -> stops at {x4, x1}

cat("\n---- Backward elimination (F_out = 4) ----\n")
bwd <- backward_elimination(F_out = 4)
# remove x3 (0.02) -> remove x4 (1.86) -> stop.  Final {x1, x2}

cat("\n---- Stepwise (F_in = F_out = 4) ----\n")
stp <- stepwise_selection(F_in = 4, F_out = 4)
# enter x4, x1, x2; then x4 has F-to-remove 1.86 -> removed.  Final {x1, x2}

sequential <- list(forward = fwd, forward_exact5 = fwd_exact, backward = bwd, stepwise = stp)
print(lapply(sequential, sort))

# ---- 4. Criterion-based automatic search with step() ------------------------------------
null_fit <- lm(y ~ 1, data = cement)
step_aic_both <- step(null_fit, scope = ~ x1 + x2 + x3 + x4, direction = "both", k = 2)
step_aic_back <- step(full, direction = "backward", k = 2, trace = 0)
step_bic_both <- step(null_fit, scope = ~ x1 + x2 + x3 + x4, direction = "both",
                      k = log(n), trace = 0)
print(list(AIC_both = formula(step_aic_both),        # y ~ x4 + x1 + x2
           AIC_back = formula(step_aic_back),        # y ~ x1 + x2 + x4
           BIC_both = formula(step_bic_both)))       # y ~ x1 + x2

# ---- 5. Akaike weights and BIC posterior model probabilities ----------------------------
nonnull <- allsub[allsub$model != "(none)", ]
model_weights <- function(ic) { w <- exp(-(ic - min(ic)) / 2); w / sum(w) }
nonnull$w_AIC  <- model_weights(nonnull$AIC)
nonnull$pr_BIC <- model_weights(nonnull$BIC)
print(head(nonnull[order(-nonnull$w_AIC), c("model", "p", "AIC", "w_AIC", "BIC", "pr_BIC")], 6),
      digits = 3, row.names = FALSE)
# x1,x2,x4 0.247/0.234   x1,x2,x3 0.242/0.230   x1,x2 0.197/0.248   x1,x3,x4 0.169/0.161

# ---- 6. The two competing models ------------------------------------------------------
fit_12  <- lm(y ~ x1 + x2, data = cement)
fit_124 <- lm(y ~ x1 + x2 + x4, data = cement)
print(round(summary(fit_12)$coefficients, 4))    # 52.58 + 1.468 x1 + 0.662 x2, |t| > 12
print(round(summary(fit_124)$coefficients, 4))   # b2 = 0.416 (p = 0.052), b4 = -0.237 (p = 0.205)
print(round(car::vif(fit_12), 2))                # 1.06 1.06
print(round(car::vif(fit_124), 2))               # 1.07 18.78 18.94

# ---- FIGURE 1: Mallows' Cp plot (`cp-plot-cement.pdf`) ----------------------------------
cp_max <- 10
cp_df  <- nonnull[nonnull$Cp <= cp_max, c("model", "p", "Cp")]
n_hidden <- sum(nonnull$Cp > cp_max)               # 8 subsets with Cp > 10

# Hand-placed labels (the three p = 4 subsets sit within 0.5 of each other)
lab_pos <- data.frame(
  model = c("x1,x2", "x1,x4", "x1,x2,x4", "x1,x2,x3", "x1,x3,x4", "x2,x3,x4", "x1,x2,x3,x4"),
  dx    = c(-0.05,   0.10,    -0.10,      0.10,       0.10,       0.10,       0.00),
  ly    = c(1.95,    5.50,    2.75,       3.10,       3.70,       7.34,       5.75),
  hjust = c(1,       0,       1,          0,          0,          0,          0.5)
)
cp_df <- merge(cp_df, lab_pos, by = "model")
cp_df$px  <- cp_df$p + ifelse(cp_df$model == "x1,x2,x4", -0.04,
                       ifelse(cp_df$model == "x1,x2,x3", 0.04, 0))
cp_df$lab <- lab_expr(cp_df$model)
cp_df$hl  <- cp_df$model %in% c("x1,x2", "x1,x2,x4")

p_cp <- ggplot(cp_df, aes(x = px, y = Cp)) +
  geom_abline(intercept = 0, slope = 1, colour = ref_col, linetype = "dashed", linewidth = 0.6) +
  annotate("text", x = 1.75, y = 2.65, label = "italic(C)[italic(p)] == italic(p)",
           parse = TRUE, hjust = 0, size = 3.0, colour = ref_col) +
  geom_point(aes(colour = hl, size = hl)) +
  geom_text(aes(x = px + dx, y = ly, label = lab, hjust = hjust),
            parse = TRUE, size = 2.9, colour = "grey20") +
  annotate("text", x = 1.75, y = 9.6, hjust = 0, vjust = 1, size = 2.7, colour = "grey40",
           label = sprintf("Not shown: %d subsets with Cp > %d\n(every one-regressor model: Cp > 138)",
                           n_hidden, cp_max)) +
  scale_colour_manual(values = c("FALSE" = pt_col, "TRUE" = pt_flag_col), guide = "none") +
  scale_size_manual(values = c("FALSE" = 2.2, "TRUE" = 3.0), guide = "none") +
  scale_x_continuous(name = expression("Number of parameters" ~ italic(p)),
                     limits = c(1.7, 5.5), breaks = 2:5) +
  scale_y_continuous(name = expression("Mallows'" ~ italic(C)[italic(p)]),
                     limits = c(0, cp_max), breaks = seq(0, 10, by = 2)) +
  theme_lecture()

ggsave(out_file("cp-plot-cement.pdf"), p_cp, bg = "transparent",
       device = cairo_pdf, width = 4.8, height = 3.4)

# ---- FIGURE 2: Best subset of each size (`criteria-best-subset-cement.pdf`) -------------
# For fixed p, R2, R2adj, MSRes, Cp, AIC, BIC, HQIC all rank subsets by SSRes, so they share the
# same best subset; PRESS does not, so its best subset is found separately.
sizes <- 2:5
best_ss    <- do.call(rbind, lapply(sizes, function(pp) {
  s <- allsub[allsub$p == pp, ]; s[which.min(s$SSRes), ] }))
best_press <- do.call(rbind, lapply(sizes, function(pp) {
  s <- allsub[allsub$p == pp, ]; s[which.min(s$PRESS), ] }))
stopifnot(identical(best_ss$model, best_press$model))      # here they coincide: x4, x1x2, x1x2x4, full

size_labels <- paste0(sizes, "\n", c("{x4}", "{x1,x2}", "{x1,x2,x4}", "{x1..x4}"))

panel <- function(df, ycol, ylab, opt = c("min", "max"), log_y = FALSE, extra = NULL) {
  opt <- match.arg(opt)
  df$val <- df[[ycol]]
  i_opt <- if (opt == "min") which.min(df$val) else which.max(df$val)
  g <- ggplot(df, aes(x = p, y = val)) + extra +
    geom_line(colour = "grey55", linewidth = 0.6) +
    geom_point(colour = pt_col, size = 2) +
    geom_point(data = df[i_opt, ], colour = pt_flag_col, size = 3.2) +
    scale_x_continuous(name = NULL, breaks = sizes, labels = size_labels,
                       expand = expansion(add = 0.3)) +
    labs(y = NULL, title = ylab) +
    theme_lecture() +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 7.5, lineheight = 0.9))
  if (log_y) g <- g + scale_y_log10()
  g
}

g_r2    <- panel(best_ss, "R2", "R-squared", "max")
g_r2adj <- panel(best_ss, "R2adj", "Adjusted R-squared", "max")
g_ms    <- panel(best_ss, "MSRes", "MSRes (log scale)", "min", log_y = TRUE)
g_press <- panel(best_press, "PRESS", "PRESS (log scale)", "min", log_y = TRUE)
g_cp    <- panel(best_ss, "Cp", "Cp (log scale), with Cp = p", "min", log_y = TRUE,
                 extra = geom_line(aes(y = p), colour = ref_col, linetype = "dashed",
                                   linewidth = 0.6))

ic_long <- rbind(
  data.frame(p = best_ss$p, val = best_ss$AIC,  crit = "AIC"),
  data.frame(p = best_ss$p, val = best_ss$BIC,  crit = "BIC"),
  data.frame(p = best_ss$p, val = best_ss$HQIC, crit = "HQIC")
)
ic_long$crit <- factor(ic_long$crit, levels = c("AIC", "HQIC", "BIC"))
ic_min <- do.call(rbind, lapply(split(ic_long, ic_long$crit), function(d) d[which.min(d$val), ]))
ic_cols <- c(AIC = pt_col, HQIC = pt_lev_col, BIC = pt_flag_col)

g_ic <- ggplot(ic_long, aes(x = p, y = val, colour = crit)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.8) +
  geom_point(data = ic_min, size = 3.4, shape = 21, fill = NA, stroke = 1) +
  scale_colour_manual(values = ic_cols, name = NULL) +
  scale_x_continuous(name = NULL, breaks = sizes, labels = size_labels,
                     expand = expansion(add = 0.3)) +
  labs(y = NULL, title = "AIC, HQIC, BIC") +
  theme_lecture() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 7.5, lineheight = 0.9),
        legend.position = c(0.75, 0.78), legend.text = element_text(size = 8),
        legend.key.height = unit(0.35, "cm"))

p_best <- (g_r2 | g_r2adj | g_ms) / (g_press | g_cp | g_ic) +
  plot_annotation(caption = "Best subset of each size (p, with subset); large point / ring = optimum of the criterion",
                  theme = theme(plot.caption = element_text(size = 8.5, colour = "grey35"),
                                plot.background = element_rect(fill = "transparent", colour = NA)))

ggsave(out_file("criteria-best-subset-cement.pdf"), p_best, bg = "transparent",
       device = cairo_pdf, width = 7.6, height = 4.6)

# ---- FIGURE 3: Information-criterion differences, all subsets (`ic-allsubsets-cement.pdf`) ----
delta <- function(x) x - min(x)
ic_all <- rbind(
  data.frame(model = nonnull$model, crit = "AIC",  d = delta(nonnull$AIC)),
  data.frame(model = nonnull$model, crit = "HQIC", d = delta(nonnull$HQIC)),
  data.frame(model = nonnull$model, crit = "BIC",  d = delta(nonnull$BIC))
)
ic_all$crit  <- factor(ic_all$crit, levels = c("AIC", "HQIC", "BIC"))
ord <- nonnull$model[order(nonnull$AIC, decreasing = TRUE)]
ic_all$model <- factor(ic_all$model, levels = ord)

p_ic <- ggplot(ic_all, aes(x = d, y = model, colour = crit, shape = crit)) +
  annotate("rect", xmin = 0, xmax = 2, ymin = -Inf, ymax = Inf, fill = ref_col, alpha = 0.08) +
  geom_vline(xintercept = 2, colour = ref_col, linetype = "dashed", linewidth = 0.5) +
  geom_point(size = 2.1, position = position_dodge(width = 0.55)) +
  annotate("text", x = 2.15, y = 1.2, label = "Delta < 2", hjust = 0,
           size = 2.8, colour = ref_col, parse = TRUE) +
  scale_colour_manual(values = ic_cols, name = NULL) +
  scale_shape_manual(values = c(AIC = 16, HQIC = 17, BIC = 15), name = NULL) +
  scale_x_sqrt(name = expression(Delta == IC - min ~ IC ~ "(square-root scale)"),
               breaks = c(0, 0.5, 2, 5, 10, 20, 30, 40),
               labels = c("0", "0.5", "2", "5", "10", "20", "30", "40")) +
  scale_y_discrete(name = NULL, labels = function(l) parse(text = lab_expr(l))) +
  theme_lecture() +
  theme(legend.position = "bottom", legend.margin = margin(t = -5),
        panel.grid.major.y = element_line(colour = "grey92"))

ggsave(out_file("ic-allsubsets-cement.pdf"), p_ic, bg = "transparent",
       device = cairo_pdf, width = 4.8, height = 4.2)

cat("\nFigures written to", normalizePath(fig_dir), "\n")
