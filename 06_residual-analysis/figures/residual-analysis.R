if (!requireNamespace("lmtest", quietly = TRUE)) install.packages("lmtest")
library(lmtest)

save_pdf <- function(name, expr, width = 4.5, height = 3.6) {
    pdf(file.path(paste0(name, ".pdf")), width = width, height = height, pointsize = 11)
    on.exit(dev.off())
    par(mar = c(4, 4, 1.5, 1))
    force(expr)
    invisible(NULL)
}


# ---------------------------------------------------------------------------
# 1. Delivery time data and fitted model
# ---------------------------------------------------------------------------
delivery <- data.frame(
    y  = c(16.68, 11.50, 12.03, 14.88, 13.75, 18.11, 8.00, 17.83, 79.24, 21.50,
           40.33, 21.00, 13.50, 19.75, 24.00, 29.00, 15.35, 19.00, 9.50, 35.10,
           17.90, 52.32, 18.75, 19.83, 10.75),
    x1 = c(7, 3, 3, 4, 6, 7, 2, 7, 30, 5, 16, 10, 4, 6, 9, 10, 6, 7, 3, 17, 10, 26, 9, 8, 4),
    x2 = c(560, 220, 340, 80, 150, 330, 110, 210, 1460, 605, 688, 215, 255, 462, 448, 776,
           200, 132, 36, 770, 140, 810, 450, 635, 150)
)

fit <- lm(y ~ x1 + x2, data = delivery)
summary(fit)          # coefficients 2.341, 1.616, 0.0144; sigma-hat = 3.259

n <- nrow(delivery)
p <- length(coef(fit))


# ---------------------------------------------------------------------------
# 2. The five types of residuals, computed from the definitions
# ---------------------------------------------------------------------------
e   <- unname(resid(fit))                       # raw
h   <- unname(hatvalues(fit))                   # leverage h_ii
mse <- sum(e^2) / (n - p)                       # 10.624
d   <- e / sqrt(mse)                            # standardized
r   <- e / sqrt(mse * (1 - h))                  # studentized (R: rstandard)
pe  <- e / (1 - h)                              # PRESS residuals
s2i <- ((n - p) * mse - e^2 / (1 - h)) / (n - p - 1)   # S^2_(i)
ti  <- e / sqrt(s2i * (1 - h))                  # R-student (R: rstudent)

press <- sum(pe^2)                              # 459.04
press

# Agreement with the built-in functions
# (note the naming: R's rstandard() is the studentized residual r_i,
#  R's rstudent() is the R-student residual t_i)
stopifnot(isTRUE(all.equal(r,  unname(rstandard(fit)))))
stopifnot(isTRUE(all.equal(ti, unname(rstudent(fit)))))

resid_table <- data.frame(
    i = seq_len(n), y = delivery$y, yhat = unname(fitted(fit)),
    h = h, e = e, d = d, r = r, PRESS = pe, Rstudent = ti
)
print(round(resid_table, 4), row.names = FALSE)

# Rows shown in the slides: i = 1, 9, 10, 15, 20, 22
round(resid_table[c(1, 9, 10, 15, 20, 22), c("i", "h", "e", "d", "r", "PRESS", "Rstudent")], 3)

# Observation 9: effect of leaving it out of the variance estimate
c(MSE = mse, S2_9 = s2i[9], sqrt_1_minus_h9 = sqrt(1 - h[9]))    # 10.624, 5.905, 0.708


# ---------------------------------------------------------------------------
# 3. Formal tests
# ---------------------------------------------------------------------------
# Normality: Shapiro-Wilk on the studentized residuals (W = 0.923, p = 0.060)
shapiro.test(rstandard(fit))

# Constant variance
# Breusch-Pagan, studentized (Koenker) form, the default of bptest(): LM = 11.99, df = 2, p = 0.0025
bptest(fit)

# White: auxiliary regression on regressors, squares and cross-product: LM = 14.96, df = 5, p = 0.0105
bptest(fit, ~ x1 + x2 + I(x1^2) + I(x2^2) + x1:x2, data = delivery)

# Sensitivity to observation 9 (diagnostic only)
fit_no9 <- lm(y ~ x1 + x2, data = delivery[-9, ])
bptest(fit_no9)                       # LM = 3.22, p = 0.20
shapiro.test(rstandard(fit_no9))      # W = 0.970, p = 0.68


# ---------------------------------------------------------------------------
# 3b. Leverage: geometric intuition (synthetic simple-regression data)
# ---------------------------------------------------------------------------
set.seed(3)
n_lev <- 18
x_lev <- runif(n_lev, 2, 8)
y_lev <- 1 + 1.2 * x_lev + rnorm(n_lev, sd = 0.6)
x_far <- 16
y_far <- 1 + 1.2 * x_far + 0.3    # close to the trend line despite being far in x

x_all <- c(x_lev, x_far)
y_all <- c(y_lev, y_far)
fit_lev <- lm(y_all ~ x_all)
h_far <- hatvalues(fit_lev)[length(y_all)]     # 0.733
cutoff_lev <- 2 * 2 / (n_lev + 1)              # 0.211
c(h_far = h_far, cutoff = cutoff_lev)

save_pdf("leverage-geometry", {
    plot(x_lev, y_lev, pch = 19, xlim = range(x_all) + c(-0.5, 0.5), ylim = range(y_all) + c(-1, 1),
         xlab = "x", ylab = "y")
    points(x_far, y_far, pch = 19, col = "red")
    abline(fit_lev, lty = 2)
    abline(v = mean(x_lev), lty = 3, col = "blue")
    text(mean(x_lev), max(y_all) + 0.7, expression(bar(x)), col = "blue")
    legend("topleft", legend = c("typical points", "high-leverage point"),
           col = c("black", "red"), pch = 19, bty = "n", cex = 0.8)
})


# ---------------------------------------------------------------------------
# 4. Delivery time figures (studentized residuals)
# ---------------------------------------------------------------------------
rs <- rstandard(fit)
label_obs <- function(x, y, obs = c(9, 22)) text(x[obs], y[obs], labels = obs, pos = 2, cex = 0.8)

save_pdf("qq-delivery", {
    q <- qqnorm(rs, main = "", xlab = "Theoretical quantiles", ylab = "Studentized residual", pch = 19)
    qqline(rs)
    label_obs(q$x, q$y)
})

save_pdf("rvf-delivery", {
    plot(fitted(fit), rs, pch = 19, xlab = expression(hat(y)), ylab = "Studentized residual")
    abline(h = 0, lty = 2)
    label_obs(fitted(fit), rs)
})

save_pdf("rvx1-delivery", {
    plot(delivery$x1, rs, pch = 19, xlab = expression(x[1]~"(cases)"), ylab = "Studentized residual")
    abline(h = 0, lty = 2)
    label_obs(delivery$x1, rs)
})

save_pdf("rvx2-delivery", {
    plot(delivery$x2, rs, pch = 19, xlab = expression(x[2]~"(distance, ft)"), ylab = "Studentized residual")
    abline(h = 0, lty = 2)
    label_obs(delivery$x2, rs)
})


# ---------------------------------------------------------------------------
# 5. Simulated patterns
# ---------------------------------------------------------------------------
set.seed(2025)
m <- 80

# Normal probability plot patterns
qq_panel <- function(x, main) {
    qqnorm(x, main = main, xlab = "Theoretical quantiles", ylab = "Sample quantiles", pch = 19, cex = 0.6)
    qqline(x)
}
save_pdf("qq-patterns", {
    par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
    qq_panel(rnorm(100), "Normal")
    qq_panel(rt(100, df = 3), "Heavy-tailed")
    qq_panel(runif(100, -1, 1), "Light-tailed")
    qq_panel(rexp(100) - 1, "Right-skewed")
}, width = 6, height = 5)

# Residuals vs fitted values patterns
res_plot <- function(fv, res) {
    plot(fv, res, pch = 19, cex = 0.7, xlab = "Fitted value", ylab = "Residual")
    abline(h = 0, lty = 2)
}

f <- sort(runif(m, 5, 30))
save_pdf("pat-ok",     res_plot(f, rnorm(m)))
save_pdf("pat-funnel", res_plot(f, rnorm(m, sd = f / 12)))

p_hat <- sort(runif(m, 0.05, 0.95))
save_pdf("pat-bow", res_plot(p_hat, rnorm(m, sd = 0.2 * sqrt(p_hat * (1 - p_hat)))))

x_nl  <- runif(m, 0, 10)
y_nl  <- 2 + 3 * x_nl - 0.3 * x_nl^2 + rnorm(m, sd = 0.7)
fit_nl <- lm(y_nl ~ x_nl)
save_pdf("pat-nonlinear", res_plot(fitted(fit_nl), resid(fit_nl)))

# a second curvature example, bowing the opposite way, for the mean-function exercise
x_nl2  <- runif(m, 0, 10)
y_nl2  <- 2 + 0.5 * x_nl2 + 0.35 * x_nl2^2 + rnorm(m, sd = 1.0)
fit_nl2 <- lm(y_nl2 ~ x_nl2)
save_pdf("pat-nonlinear2", res_plot(fitted(fit_nl2), resid(fit_nl2)))

# Residuals in time order: positive and negative autocorrelation
ar1 <- function(phi, len = 40) {
    a <- rnorm(len)
    z <- numeric(len)
    z[1] <- a[1]
    for (tt in 2:len) z[tt] <- phi * z[tt - 1] + a[tt]
    z
}
time_plot <- function(z) {
    plot(seq_along(z), z, type = "b", pch = 19, cex = 0.6, xlab = "Time order", ylab = "Residual")
    abline(h = 0, lty = 2)
}
save_pdf("time-positive", time_plot(ar1(0.9)))
save_pdf("time-negative", time_plot(ar1(-0.9)))


# ---------------------------------------------------------------------------
# 6. Soft drink sales data: Durbin-Watson test
# ---------------------------------------------------------------------------
softdrink <- data.frame(
    year = 1:20,
    y = c(3083, 3149, 3218, 3239, 3295, 3374, 3475, 3569, 3597, 3725,
          3794, 3959, 4043, 4194, 4318, 4493, 4683, 4850, 5005, 5236),
    x = c(75, 78, 80, 82, 84, 88, 93, 97, 99, 104,
          109, 115, 120, 127, 135, 144, 153, 161, 170, 182)
)

fit_sd <- lm(y ~ x, data = softdrink)
coef(fit_sd)                                   # 1608.5, 20.091

e_sd <- unname(resid(fit_sd))
d_stat <- sum(diff(e_sd)^2) / sum(e_sd^2)      # 1.08
r1 <- sum(e_sd[-1] * e_sd[-length(e_sd)]) / sum(e_sd^2)
c(d = d_stat, two_times_1_minus_r1 = 2 * (1 - r1))

# H1: positive autocorrelation (default alternative = "greater")
dwtest(fit_sd)
# Tabulated bounds (alpha = 0.05, n = 20, k = 1): d_L = 1.20, d_U = 1.41  ->  d < d_L: reject H0

save_pdf("softdrink-time", {
    plot(softdrink$year, e_sd, type = "b", pch = 19, cex = 0.7, xlab = "Year", ylab = "Residual")
    abline(h = 0, lty = 2)
})
