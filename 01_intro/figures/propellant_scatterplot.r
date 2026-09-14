library(ggplot2)
propellant <- data.frame(
      age = c(15.50,23.75,8.00,17.00,5.50,19.00,24.00,2.50,7.50,11.00,
                          13.00,3.75,25.00,9.75,22.00,18.00,6.00,12.50,2.00,21.50),
      shear = c(2158.70,1678.15,2316.00,2061.30,2207.50,1708.30,1784.70,
                              2575.00,2357.90,2256.70,2165.20,2399.55,1779.80,2336.75,
                              1765.30,2053.50,2414.40,2200.50,2654.20,1753.70)
    )

fit <- lm(shear ~ age, data = propellant)
propellant$fitted <- fitted(fit)


ggplot(propellant, aes(x = age, y = shear)) +
      geom_point(size = 1, color = "#1b9e77") +
      labs(x = "Age of Propellant (weeks)", y = "Shear Strength (psi)") +
      theme_minimal(base_size = 12) +
      theme(
            panel.background = element_rect(fill = "transparent", colour = NA),
            plot.background  = element_rect(fill = "transparent", colour = NA)
          )
ggsave("propellant_scatter.pdf", bg = "transparent", width = 4.5, height = 3.5)




ggplot(propellant, aes(x = age, y = shear)) +
    geom_point(size = 1, color = "#1b9e77") +
    geom_smooth(method = "lm", se = FALSE, color = "#d95f02", linewidth = 0.8) +
    labs(x = "Age of Propellant (weeks)", y = "Shear Strength (psi)") +
    theme_minimal(base_size = 12) +
    theme(
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background  = element_rect(fill = "transparent", colour = NA)
    )
ggsave("propellant_fitted_line.pdf", bg = "transparent", width = 4.5, height = 3.5)




ggplot(propellant, aes(x = age, y = shear)) +
    geom_segment(aes(xend = age, yend = fitted), linetype = "dotted", color = "#e7298a", linewidth = 0.5) +
    geom_smooth(method = "lm", se = FALSE, color = "#d95f02", linewidth = 0.8) +
    geom_point(size = 1, color = "#1b9e77") +
    labs(x = "Age of Propellant (weeks)", y = "Shear Strength (psi)") +
    theme_minimal(base_size = 12) +
    theme(
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background  = element_rect(fill = "transparent", colour = NA)
    )
ggsave("propellant_residuals.pdf", bg = "transparent", width = 4.5, height = 3.5)
