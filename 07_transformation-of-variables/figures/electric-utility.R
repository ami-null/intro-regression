# Electric utility data (MPV Table 5.2, Examples 5.1 and 5.3)
# Variance-stabilizing transformation and Box-Cox procedure
library(MASS)

eu <- read.csv("electric-utility.csv")

plot_pdf <- function(file, expr, w = 5, h = 4) {
    pdf(file, width = w, height = h)
    par(mar = c(4.2, 4.2, 1, 1))
    force(expr)
    dev.off()
}

# ---- Scatter plot of demand vs usage (Figure 5.1) ----
plot_pdf("eu_scatter.pdf", {
    plot(eu$x, eu$y, pch = 19, xlab = "Usage (kWh)", ylab = "Demand (kW)")
})

# ---- Straight-line fit on the original scale ----
fit0 <- lm(y ~ x, data = eu)
print(summary(fit0))

plot_pdf("eu_rstudent_original.pdf", {
    plot(fitted(fit0), rstudent(fit0), pch = 19,
         xlab = expression(hat(y)), ylab = "R-student")
    abline(h = 0, lty = 2)
})

# ---- Square-root transformation ----
fit_sqrt <- lm(sqrt(y) ~ x, data = eu)
print(summary(fit_sqrt))

# R^2 on the original scale: squared correlation of y with the back-transformed fit
cat("R^2 of y vs (fitted sqrt)^2:", cor(eu$y, fitted(fit_sqrt)^2)^2, "\n")

plot_pdf("eu_rstudent_sqrt.pdf", {
    plot(fitted(fit_sqrt), rstudent(fit_sqrt), pch = 19,
         xlab = expression(hat(y)^"*"), ylab = "R-student")
    abline(h = 0, lty = 2)
})

# ---- Log transformation (for comparison) ----
fit_log <- lm(log(y) ~ x, data = eu)
print(summary(fit_log))

plot_pdf("eu_rstudent_log.pdf", {
    plot(fitted(fit_log), rstudent(fit_log), pch = 19,
         xlab = expression(hat(y)^"*"), ylab = "R-student")
    abline(h = 0, lty = 2)
})

# ---- Box-Cox procedure with MASS::boxcox ----
# Fine grid: a coarse grid only keeps grid points inside the interval and
# reports a slightly narrower interval than the true one
lambda_grid <- seq(-0.5, 1.5, by = 0.001)
bc <- boxcox(y ~ x, data = eu, lambda = lambda_grid, plotit = FALSE)

lambda_hat <- bc$x[which.max(bc$y)]
cutoff <- max(bc$y) - qchisq(0.95, 1) / 2
ci <- range(bc$x[bc$y >= cutoff])
cat("lambda_hat =", lambda_hat, "\n")
cat("95% CI for lambda: (", ci[1], ",", ci[2], ")\n")

plot_pdf("eu_boxcox.pdf", {
    boxcox(y ~ x, data = eu, lambda = lambda_grid, plotit = TRUE)
})

# ---- SSRes(lambda) for the scaled response y^(lambda) (MPV, Example 5.3) ----
n <- nrow(eu)
y_dot <- exp(mean(log(eu$y)))

ssres_scaled <- function(lambda) {
    y_l <- if (abs(lambda) < 1e-8) {
        y_dot * log(eu$y)
    } else {
        (eu$y^lambda - 1) / (lambda * y_dot^(lambda - 1))
    }
    sum(resid(lm(y_l ~ eu$x))^2)
}

lambdas <- c(-1, -0.5, 0, 0.25, 0.5, 0.75, 1)
ssres_tab <- data.frame(lambda = lambdas,
                        SSRes = sapply(lambdas, ssres_scaled))
print(ssres_tab, digits = 6)

# Unscaled SSRes: not comparable across lambda (different units)
cat("Unscaled SSRes for y, sqrt(y), log(y):",
    sum(resid(fit0)^2), sum(resid(fit_sqrt)^2), sum(resid(fit_log)^2), "\n")

# MPV cutoff SS* = SSRes(lambda_hat) * (1 + t^2 / nu)
opt <- optimize(ssres_scaled, c(-1, 2))
nu <- n - 2
ss_star <- opt$objective * (1 + qt(0.975, nu)^2 / nu)
cat("lambda_hat =", opt$minimum, " SSRes(lambda_hat) =", opt$objective,
    " SS* =", ss_star, "\n")

# Same estimate and a Wald interval from the car package
if (requireNamespace("car", quietly = TRUE)) {
    print(summary(car::powerTransform(fit0)))
}

# ---- Model with the selected transformation ----
fit_bc <- lm(sqrt(y) ~ x, data = eu)
print(summary(fit_bc))
