# Lecture 13: fixed-design sampling distributions with known relative variances.
library(ggplot2)
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
out <- dirname(normalizePath(script))
set.seed(1303)
x <- seq(1, 10, length.out = 50)
design <- model.matrix(~ x)
w <- 1 / x^2
# Each row maps one response vector to its fitted slope.
slope_map <- function(weights) {
    solve(crossprod(design, weights * design), t(weights * design))[2, ]
}
maps <- rbind(OLS = slope_map(rep(1, length(x))), WLS = slope_map(w))
errors <- matrix(rnorm(length(x) * 5000), nrow = length(x)) * x
slopes <- maps %*% (2 + 3 * x + errors)
results <- data.frame(mean_slope = rowMeans(slopes), empirical_sd = apply(slopes, 1, sd),
                      exact_sd = sqrt(rowSums(sweep(maps, 2, x, "*")^2)))
capture.output(print(results), file = file.path(out, "slope-simulation.txt"))
samples <- data.frame(slope = as.vector(t(slopes)), method = rep(rownames(maps), each = ncol(slopes)))
p <- ggplot(samples, aes(slope, colour = method, linetype = method)) +
    geom_density(linewidth = 0.9) +
    geom_vline(xintercept = 3, colour = "grey40", linetype = "dotted") +
    scale_colour_manual(values = c(OLS = "#386cb0", WLS = "#d95f02")) +
    labs(x = "Estimated slope", y = "Density", colour = NULL, linetype = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
ggsave(file.path(out, "slope-simulation.pdf"), p, width = 6, height = 3.5, bg = "transparent")
