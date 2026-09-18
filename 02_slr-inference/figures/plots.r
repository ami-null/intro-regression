library(ggplot2)

rocket <- data.frame(
      x = c(15.50,23.75,8.00,17.00,5.50,19.00,24.00,2.50,7.50,11.00,
            13.00,3.75,25.00,9.75,22.00,18.00,6.00,12.50,2.00,21.50),
      y = c(2158.70,1678.15,2316.00,2061.30,2207.50,1708.30,1784.70,
            2575.00,2357.90,2256.70,2165.20,2399.55,1779.80,2336.75,
            1765.30,2053.50,2414.40,2200.50,2654.20,1753.70)
    )

# fit <- lm(y ~ x, data = rocket)
# xbar <- mean(rocket$x); ybar <- mean(rocket$y)
# b1 <- coef(fit)[2]
# se_b1 <- summary(fit)$coefficients["x", "Std. Error"]
# tcrit <- qt(0.975, df = fit$df.residual)
# b1_lo <- b1 - tcrit * se_b1
# b1_hi <- b1 + tcrit * se_b1
#
# xs <- seq(min(rocket$x) - 1, max(rocket$x) + 1, length.out = 100)
# wedge <- data.frame(x = xs, lo = ybar + b1_hi * (xs - xbar), hi = ybar + b1_lo * (xs - xbar))
#
# ggplot(rocket, aes(x, y)) +
#       geom_ribbon(data = wedge, aes(x = x, ymin = lo, ymax = hi),
#                                   inherit.aes = FALSE, fill = "grey40", alpha = 0.15) +
#       geom_line(data = wedge, aes(x = x, y = lo), inherit.aes = FALSE,
#                               linetype = "dashed", colour = "grey40") +
#       geom_line(data = wedge, aes(x = x, y = hi), inherit.aes = FALSE,
#                               linetype = "dashed", colour = "grey40") +
#       geom_smooth(method = "lm", se = FALSE, colour = "#B23A48", linewidth = 1) +
#       geom_point(colour = "#3B6FA0", size = 2) +
#       annotate("point", x = xbar, y = ybar, size = 2) +
#       annotate("text", x = xbar, y = ybar, label = "list(bar(x), bar(y))", parse = TRUE,
#              hjust = -0.15, vjust = -0.8) +
#       labs(x = "Age of propellant, x (weeks)", y = "Shear strength, y (psi)") +
#       theme_minimal(base_size = 12)


xs  <- seq(min(rocket$x) - 0.5, max(rocket$x) + 0.5, length.out = 200)
newdat <- data.frame(x = xs)

ci <- predict(fit, newdata = newdat, interval = "confidence", level = 0.95)
pi <- predict(fit, newdata = newdat, interval = "prediction", level = 0.95)

bands <- data.frame(x = xs, fit = ci[, "fit"],
                                           ci_lo = ci[, "lwr"], ci_hi = ci[, "upr"],
                                           pi_lo = pi[, "lwr"], pi_hi = pi[, "upr"])

ggplot() +
      geom_ribbon(data = bands, aes(x, ymin = pi_lo, ymax = pi_hi, fill = "PI"),
              alpha = 0.20) +
  geom_ribbon(data = bands, aes(x, ymin = ci_lo, ymax = ci_hi, fill = "CI"),
              alpha = 0.35) +
  geom_line(data = bands, aes(x, y = fit), linewidth = 1) +
  geom_point(data = rocket, aes(x, y), size = 2) +
  scale_fill_manual(name = NULL,
                     values = c("CI" = "#3B6FA0", "PI" = "#B23A48")) +
  labs(x = "Age of propellant, x (weeks)", y = "Shear strength, y (psi)") +
      theme_minimal(base_size = 12) +
    theme(
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background  = element_rect(fill = "transparent", colour = NA)
    )

ggsave("fig-ci-pi-band.pdf", bg = "transparent", width = 4.5, height = 3.5)


