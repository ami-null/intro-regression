# Shapes of the linearizable functions (cf. MPV Figure 5.4)
# Each panel shows one function form for a few parameter values

pdf("lin_shapes.pdf", width = 9, height = 5)
par(mfrow = c(2, 3), mar = c(3.4, 3.4, 2.4, 0.6), mgp = c(2.1, 0.6, 0),
    cex.main = 1.7, cex.axis = 1.3, cex.lab = 1.4)

x <- seq(0.05, 3, length.out = 300)

draw <- function(curves, main, ylim, legend_txt, pos = "topleft") {
    plot(NA, xlim = range(x), ylim = ylim, xlab = "x", ylab = "y", main = main)
    for (k in seq_along(curves)) lines(x, curves[[k]], lwd = 2, lty = k)
    legend(pos, legend = legend_txt, lty = seq_along(curves), lwd = 2,
           bty = "n", cex = 1.35)
}

# Power: y = beta0 * x^beta1
draw(list(x^2, x^0.5, x^-0.5),
     expression(y == beta[0] * x^beta[1]), c(0, 4),
     c(expression(beta[1] > 1), expression(0 < beta[1] * " < 1"),
       expression(beta[1] < 0)), pos = "topright")

# Exponential: y = beta0 * exp(beta1 x)
draw(list(exp(0.6 * x), exp(-1.2 * x) * 4),
     expression(y == beta[0] * e^{beta[1] * x}), c(0, 6.5),
     c(expression(beta[1] > 0), expression(beta[1] < 0)), pos = "top")

# Logarithmic: y = beta0 + beta1 ln x
draw(list(2 + log(x), 2 - log(x)),
     expression(y == beta[0] + beta[1] * ln * x), c(-1.5, 5.5),
     c(expression(beta[1] > 0), expression(beta[1] < 0)), pos = "top")

# Reciprocal in x: y = beta0 + beta1 / x
draw(list(3 - 0.5 / x, 1 + 0.5 / x),
     expression(y == beta[0] + beta[1] / x), c(0, 4.5),
     c(expression(beta[1] < 0), expression(beta[1] > 0)), pos = "topright")

# Hyperbolic: y = x / (beta0 x - beta1)
xh <- x
draw(list(xh / (0.5 * xh + 0.3), xh / (0.5 * xh - 0.02)),
     expression(y == x / (beta[0] * x - beta[1])), c(0, 4.5),
     c(expression(beta[1] < 0), expression(beta[1] > 0)), pos = "topright")

plot.new()
dev.off()
