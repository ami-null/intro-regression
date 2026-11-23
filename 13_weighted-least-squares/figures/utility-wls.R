# Lecture 13: feasible WLS for electric utility demand (MPV Table 5.2).
# Run with Rscript from any directory. Outputs go beside the script.
library(ggplot2)
library(patchwork)

script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
out <- dirname(normalizePath(script))
customers <- read.csv(file.path(out, "electric-utility.csv"))
ols <- lm(y ~ x, customers)
customers$abs_error <- abs(residuals(ols))
stopifnot(all(customers$x > 0), all(customers$abs_error > 0))
spread <- lm(log(abs_error) ~ log(x), customers)
gamma <- unname(coef(spread)[2])
customers$weight <- (customers$x / max(customers$x))^(-2 * gamma)
wls <- lm(y ~ x, customers, weights = weight)
root_w <- sqrt(customers$weight)
transformed <- lm(I(root_w * y) ~ 0 + root_w + I(root_w * x), customers)
rescaled <- lm(y ~ x, customers, weights = 100 * weight)
normal_equations <- crossprod(model.matrix(wls), customers$weight * residuals(wls))
stopifnot(isTRUE(all.equal(unname(coef(transformed)), unname(coef(wls)))),
          isTRUE(all.equal(coef(wls), coef(rescaled))),
          isTRUE(all.equal(vcov(wls), vcov(rescaled))),
          max(abs(normal_equations)) < 1e-7)

theme_set(theme_minimal(base_size = 11) +
          theme(panel.grid.minor = element_blank(), legend.position = "bottom"))
teal <- "#1b9e77"
blue <- "#386cb0"
orange <- "#d95f02"
save_plot <- function(name, plot, width, height) {
    ggsave(file.path(out, name), plot, width = width, height = height, bg = "transparent")
}
spread_plot <- ggplot(customers, aes(log(x), log(abs_error))) +
    geom_point(colour = teal) +
    geom_abline(intercept = coef(spread)[1], slope = gamma, colour = blue) +
    labs(x = "Log usage (kWh)", y = "Log absolute OLS residual (kW)")
save_plot("utility-spread.pdf", spread_plot, 5.5, 3.7)

grid <- data.frame(x = seq(min(customers$x), max(customers$x), length.out = 150))
lines <- rbind(data.frame(grid, y = predict(ols, grid), model = "OLS"),
               data.frame(grid, y = predict(wls, grid), model = "WLS"))
fits <- ggplot(customers, aes(x, y)) +
    geom_point(colour = teal) +
    geom_line(data = lines, aes(colour = model, linetype = model), linewidth = 0.9) +
    scale_colour_manual(values = c(OLS = blue, WLS = orange)) +
    labs(x = "Usage (kWh)", y = "Demand (kW)", colour = NULL, linetype = NULL)
save_plot("utility-fits.pdf", fits, 6, 3.8)

max_weight <- which.max(customers$weight)
x_limits <- range(c(fitted(ols), fitted(wls)))
residual_panel <- function(fit, error, label) {
    panel <- data.frame(fitted = fitted(fit), error, customer = customers$customer)
    ggplot(panel, aes(fitted, error)) +
        geom_hline(yintercept = 0, colour = blue, linetype = "dashed") +
        geom_point(colour = teal) +
        geom_point(data = panel[max_weight, ], colour = orange, size = 2.5) +
        scale_x_continuous(limits = x_limits) +
        labs(x = "Fitted demand (kW)", y = label)
}
diagnostics <- residual_panel(ols, residuals(ols), "Raw OLS residual (kW)") +
    residual_panel(wls, root_w * residuals(wls), "Weighted WLS residual") +
    plot_annotation(caption = sprintf("Orange point: customer %s (largest weight).", customers$customer[max_weight]))
save_plot("utility-residuals.pdf", diagnostics, 8, 3.5)

at_1000 <- data.frame(x = 1000)
capture.output({
    cat("Demand in kW; usage in kWh. Relative weights normalized to minimum 1.\n")
    cat("Customers:", nrow(customers), " Usage range:", range(customers$x), "\n")
    for (model in list(OLS = ols, spread = spread, WLS = wls)) {
        print(summary(model))
    }
    cat("\nSlope units: kW/kWh; intercept units: kW (OLS and WLS).\n")
    cat("Weights: exponent", -2 * gamma, " range", range(customers$weight), "\n")
    cat("Largest weight: customer", customers$customer[max_weight],
        " WLS leverage", hatvalues(wls)[max_weight], "\n")
    cat("Mean fits at 1000 kWh with model-based 95% confidence intervals (kW)\n")
    print(rbind(OLS = unname(predict(ols, at_1000, interval = "confidence")[1, ]),
                WLS = unname(predict(wls, at_1000, interval = "confidence")[1, ])),
          dimnames = FALSE)
    cat("Columns above: fitted mean, lower limit, upper limit\n")
    cat("Maximum absolute weighted normal equation residual:", max(abs(normal_equations)), "\n")
    cat("Transformed-design and common weight-rescaling checks passed.\n")
    cat("\nCustomer-level diagnostics\n")
    print(data.frame(customers[c("customer", "x", "y", "weight")],
                     ols_fit = fitted(ols), ols_error = residuals(ols), wls_fit = fitted(wls),
                     weighted_error = root_w * residuals(wls), leverage = hatvalues(wls),
                     standardized_error = rstandard(wls)), row.names = FALSE)
}, file = file.path(out, "utility-wls.txt"))
