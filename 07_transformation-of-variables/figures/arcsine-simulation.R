# Simulated binomial proportions: arcsine square-root transformation
# (variance-stabilizing transformation for binomial proportions, MPV Table 5.1)

plot_pdf <- function(file, expr, w = 5, h = 4) {
    pdf(file, width = w, height = h)
    par(mar = c(4.2, 4.2, 1, 1))
    force(expr)
    dev.off()
}

set.seed(2026)
m <- 20                                   # trials per observation
levels_x <- seq(0, 1, length.out = 10)    # 10 design points
reps <- 30                                # replicates per design point
x <- rep(levels_x, each = reps)
p <- 0.05 + 0.90 * x                      # true success probability, linear in x
y <- rbinom(length(x), size = m, prob = p) / m
sim <- data.frame(x = x, p = p, y = y, y_star = asin(sqrt(y)))
write.csv(sim, "arcsine-simulated.csv", row.names = FALSE)

# ---- Sample variance of y within each design point vs. group mean ----
grp <- aggregate(
    cbind(y, y_star) ~ x,
    data = sim,
    FUN = function(v) c(mean = mean(v), var = var(v))
)

grp <- do.call(data.frame, grp)
names(grp) <- c("x", "y_mean", "y_var", "ystar_mean", "ystar_var")
print(grp, digits = 3)

pg <- seq(0, 1, length.out = 200)

plot_pdf("as_variance_original.pdf", {
    plot(grp$y_mean, grp$y_var, pch = 19, ylim = c(0, max(grp$y_var) * 1.1),
         xlab = "Group mean of y", ylab = "Group variance of y")
    lines(pg, pg * (1 - pg) / m, lty = 2)
})


plot_pdf("as_variance_transformed.pdf", {
    plot(grp$ystar_mean, grp$ystar_var, pch = 19,
         ylim = c(0, max(grp$ystar_var) * 1.1),
         xlab = expression("Group mean of " * y^"*"),
         ylab = expression("Group variance of " * y^"*"))
    abline(h = 1 / (4 * m), lty = 2)
})


# ---- Regression fits and residual plots ----
fit0 <- lm(y ~ x, data = sim)
fit_as <- lm(y_star ~ x, data = sim)
print(summary(fit0))
print(summary(fit_as))

plot_pdf("as_rstudent_original.pdf", {
    plot(fitted(fit0), rstudent(fit0), pch = 19,
         xlab = expression(hat(y)), ylab = "R-student")
    abline(h = 0, lty = 2)
})


plot_pdf("as_rstudent_transformed.pdf", {
    plot(fitted(fit_as), rstudent(fit_as), pch = 19,
         xlab = expression(hat(y)^"*"), ylab = "R-student")
    abline(h = 0, lty = 2)
})
