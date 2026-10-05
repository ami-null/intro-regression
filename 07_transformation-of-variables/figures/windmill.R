# Windmill data (MPV Table 5.5, Examples 5.2 and 5.5)
# Transformation to linearize the model and the Box-Tidwell procedure

wm <- read.csv("windmill.csv")

plot_pdf <- function(file, expr, w = 5, h = 4) {
    pdf(file, width = w, height = h)
    par(mar = c(4.2, 4.2, 1, 1))
    force(expr)
    dev.off()
}

# ---- Scatter plot of DC output vs wind velocity (Figure 5.5) ----
plot_pdf("wm_scatter.pdf", {
    plot(wm$x, wm$y, pch = 19,
         xlab = "Wind velocity, x (mph)", ylab = "DC output, y")
})

# ---- Straight-line fit ----
fit0 <- lm(y ~ x, data = wm)
print(summary(fit0))

plot_pdf("wm_rstudent_linear.pdf", {
    plot(fitted(fit0), rstudent(fit0), pch = 19,
         xlab = expression(hat(y)), ylab = "R-student")
    abline(h = 0, lty = 2)
})

# ---- Reciprocal transformation of the regressor ----
wm$x_inv <- 1 / wm$x
fit_inv <- lm(y ~ x_inv, data = wm)
print(summary(fit_inv))

plot_pdf("wm_scatter_inverse.pdf", {
    plot(wm$x_inv, wm$y, pch = 19,
         xlab = "1 / wind velocity, x' = 1/x", ylab = "DC output, y")
    abline(fit_inv)
})

plot_pdf("wm_rstudent_inverse.pdf", {
    plot(fitted(fit_inv), rstudent(fit_inv), pch = 19,
         xlab = expression(hat(y)), ylab = "R-student")
    abline(h = 0, lty = 2)
})

plot_pdf("wm_qq_inverse.pdf", {
    qqnorm(rstudent(fit_inv), pch = 19, main = "")
    qqline(rstudent(fit_inv))
})

xg <- seq(min(wm$x), max(wm$x), length.out = 200)

# ---- Box-Tidwell procedure (MPV, Section 5.4.2) ----
# Each iteration: fit y ~ x; fit y ~ x + x*log(x); alpha = gamma / beta1 + 1
box_tidwell <- function(y, x, n_iter = 4) {
    x_cur <- x
    alpha_cum <- 1
    out <- data.frame(
        iter = integer(),
        beta1 = numeric(),
        gamma = numeric(),
        t_gamma = numeric(), alpha = numeric(),
        alpha_cum = numeric()
    )

    for (k in seq_len(n_iter)) {
        fit_a <- lm(y ~ x_cur)
        b1 <- coef(fit_a)[2]
        w <- x_cur * log(x_cur)
        fit_b <- lm(y ~ x_cur + w)
        g <- coef(fit_b)[3]
        t_g <- summary(fit_b)$coefficients[3, 3]
        alpha_k <- g / b1 + 1
        alpha_cum <- alpha_cum * alpha_k
        x_cur <- x_cur^alpha_k
        out[k, ] <- c(k, b1, g, t_g, alpha_k, alpha_cum)
    }
    out
}

bt <- box_tidwell(wm$y, wm$x, n_iter = 4)
print(bt, digits = 5, row.names = FALSE)

alpha_hat <- bt$alpha_cum[nrow(bt)]
cat("Box-Tidwell estimate of alpha:", alpha_hat, "\n")

# ---- Fit with the estimated power and with alpha = -1 ----
wm$x_bt <- wm$x^alpha_hat
fit_bt <- lm(y ~ x_bt, data = wm)
print(summary(fit_bt))
cat("R^2 (estimated alpha):", summary(fit_bt)$r.squared, "\n")
cat("R^2 (alpha = -1):", summary(fit_inv)$r.squared, "\n")

plot_pdf("wm_fit_boxtidwell.pdf", {
    plot(wm$x, wm$y, pch = 19,
         xlab = "Wind velocity, x (mph)", ylab = "DC output, y")
    lines(xg, predict(fit_bt, data.frame(x_bt = xg^alpha_hat)), lwd = 2)
    lines(xg, predict(fit_inv, data.frame(x_inv = 1 / xg)), lwd = 2, lty = 2)
    legend("bottomright",
           legend = c(bquote(alpha == .(round(alpha_hat, 2))), expression(alpha == -1)),
           lty = c(1, 2), lwd = 2, bty = "n")
})
