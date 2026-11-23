# Lecture 13: run with Rscript; outputs are saved beside this script.
library(ggplot2)

script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
out <- dirname(normalizePath(script))

reading <- c(10, 16)
sd_reading <- c(1, 3)
variance <- sd_reading^2
weight <- 1 / variance
fraction <- weight / sum(weight)
estimate <- weighted.mean(reading, weight)
grid <- seq(8, 18, by = 0.05)
loss1 <- weight[1] * (reading[1] - grid)^2
loss2 <- weight[2] * (reading[2] - grid)^2
curves <- data.frame(mean = rep(grid, 3), loss = c(loss1, loss2, loss1 + loss2),
                     curve = rep(c("Instrument 1", "Instrument 2", "Total"), each = length(grid)))
plot <- ggplot(curves, aes(mean, loss, colour = curve, linetype = curve)) +
    geom_line(linewidth = 0.85) +
    geom_point(data = data.frame(mean = estimate, loss = sum(weight * (reading - estimate)^2)),
               aes(mean, loss), inherit.aes = FALSE, size = 2.5) +
    scale_colour_manual(values = c("#1b9e77", "#d95f02", "#386cb0")) +
    labs(x = "Candidate common mean", y = "Weighted squared error", colour = NULL, linetype = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
ggsave(file.path(out, "precision-loss.pdf"), plot, width = 6.5, height = 3.8, bg = "transparent")
capture.output({
    print(data.frame(reading, sd_reading, variance, weight, fraction), row.names = FALSE)
    print(data.frame(method = c("Ordinary mean", "Precision-weighted mean"),
                     estimate = c(mean(reading), estimate),
                     sampling_variance = c(sum(variance) / length(reading)^2, sum(fraction^2 * variance))))
    cat("Minimum total weighted squared error:", sum(weight * (reading - estimate)^2), "\n")
}, file = file.path(out, "precision-example.txt"))
