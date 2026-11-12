library(ggplot2)

gasoline <- read.csv("gasoline.csv", na.strings = "NA")
gas <- na.omit(gasoline[, c("automobile", "y", "x2", "x3", "x5", "x6", "x10")])
nrow(gas)                                # 30

pt_col   <- "#1b9e77"
line_col <- "#d95f02"
fig_w    <- 3.4
fig_h    <- 2.5

# x2 = horsepower, x3 = torque, x5 = rear axle ratio, x6 = carburetor barrels, x10 = weight
fit_full <- lm(y ~ x2 + x3 + x5 + x6 + x10, data = gas)
summary(fit_full)                        # all t-tests insignificant, overall F = 14.82
deviance(fit_full)                       # SSRes(FM) = 278.62, df = 24
sigma(fit_full)^2                        # MSRes = 11.609
cor(gas[, c("x2", "x3", "x5", "x6", "x10")])

# ---- Partial F-tests ----------------------------------------------------------
fit_red3 <- lm(y ~ x5 + x6, data = gas)               # drop hp, tq, wt (r = 3)
coef(fit_red3)                                        # 3.768, 6.718, -1.713
anova(fit_red3, fit_full)                             # extra SS 314.34, F = 9.03, p = 0.0004

fit_red2 <- lm(y ~ x2 + x3 + x10, data = gas)         # drop ax, cb (r = 2)
anova(fit_red2, fit_full)                             # extra SS 15.03, F = 0.65, p = 0.53

# ---- General linear hypothesis ------------------------------------------------
# Matrix form: SSH = (T b - c)' [T (X'X)^-1 T']^-1 (T b - c), F0 = (SSH / r) / MSRes
glh <- function(fit, T, c = rep(0, nrow(T))) {
    b    <- coef(fit)
    XtXi <- summary(fit)$cov.unscaled
    d    <- T %*% b - c
    M    <- T %*% XtXi %*% t(T)
    ssh  <- drop(t(d) %*% solve(M) %*% d)
    r    <- nrow(T)
    F0   <- (ssh / r) / sigma(fit)^2
    c(SSH = ssh, r = r, F0 = F0, p = pf(F0, r, df.residual(fit), lower.tail = FALSE))
}

# parameter order: (Intercept), x2, x3, x5, x6, x10
T_equal <- rbind(c(0, 1, -1, 0, 0, 0))                # beta_hp = beta_tq
glh(fit_full, T_equal)                                # SSH = 0.0105, F0 = 0.0009, p = 0.976
fit_equal <- lm(y ~ I(x2 + x3) + x5 + x6 + x10, data = gas)
deviance(fit_equal)                                   # 278.634
anova(fit_equal, fit_full)

T_two <- rbind(c(0, 1, -1, 0, 0, 0),                  # beta_hp = beta_tq
               c(0, 0, 0, 1, 0, 0))                   # beta_ax = 0
glh(fit_full, T_two)                                  # SSH = 16.62, F0 = 0.72, p = 0.50
fit_two <- lm(y ~ I(x2 + x3) + x6 + x10, data = gas)
deviance(fit_two)                                     # 295.24
anova(fit_two, fit_full)

T_cb <- rbind(c(0, 0, 0, 0, 1, 0))                    # beta_cb = 1
glh(fit_full, T_cb, c = 1)                            # SSH = 4.18, F0 = 0.36, p = 0.55
fit_cb <- lm(I(y - x6) ~ x2 + x3 + x5 + x10, data = gas)
deviance(fit_cb)                                      # 282.81
((coef(fit_full)["x6"] - 1) / coef(summary(fit_full))["x6", "Std. Error"])^2   # t^2 = 0.36

# Cross-check with the car package
car::linearHypothesis(fit_full, "x2 = x3")
car::linearHypothesis(fit_full, c("x2 = x3", "x5 = 0"))
car::linearHypothesis(fit_full, "x6 = 1")

# ---- Sequential vs. partial sums of squares -----------------------------------
anova(fit_full)                                       # order: hp, tq, ax, cb, wt
anova(lm(y ~ x10 + x6 + x5 + x3 + x2, data = gas))    # reversed order
drop1(fit_full, test = "F")                           # partial
car::Anova(fit_full, type = "III")

# ---- Figure: rejection region for the partial F-test --------------------------
F0   <- anova(fit_red3, fit_full)$F[2]
crit <- qf(0.95, 3, 24)
grid <- data.frame(f = seq(0.001, 11, length.out = 600))
grid$d <- df(grid$f, 3, 24)

p_f <- ggplot(grid, aes(x = f, y = d)) +
    geom_area(data = subset(grid, f >= crit), fill = line_col, alpha = 0.35) +
    geom_line(colour = "grey30") +
    geom_vline(xintercept = crit, linetype = "dashed", colour = line_col) +
    geom_vline(xintercept = F0, colour = pt_col, linewidth = 0.8) +
    annotate("text", x = crit + 0.25, y = 0.55, hjust = 0, size = 3, colour = line_col,
           label = sprintf("critical value = %.2f", crit)) +
    annotate("text", x = F0 + 0.1, y = 0.30, hjust = 0, size = 3, colour = pt_col,
           label = sprintf("F0 = %.2f", F0)) +
    labs(x = "F", y = "Density") +
    coord_cartesian(xlim = c(0, 11)) +
    theme_minimal(base_size = 11) +
    theme(
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background  = element_rect(fill = "transparent", colour = NA)
    )

p_f

ggsave("f_dist_partial.pdf", p_f, bg = "transparent", width = fig_w, height = fig_h)
