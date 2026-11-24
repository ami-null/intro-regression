# Run from the figures directory with Rscript distance-contours.R.
library(ggplot2)

delivery <- read.csv("delivery.csv")
z <- scale(delivery[c("cases", "distance")])
points <- data.frame(obs = delivery$obs, z)
selected <- points[match(c(9, 11, 21, 22), points$obs), ]
selected$label_x <- selected$cases + c(-0.35, 0.12, 0.32, 0.33)
selected$label_y <- selected$distance + c(0.12, -0.32, -0.18, 0.12)

angle <- seq(0, 2 * pi, length.out = 361)
circle <- cbind(cos(angle), sin(angle))
# The Cholesky factor maps the unit circle to the unit Mahalanobis ellipse.
shapes <- list(Euclidean = circle, Mahalanobis = circle %*% chol(cov(z)))
contours <- list()
for (metric in names(shapes)) {
    for (radius in 1:3) {
        curve <- radius * shapes[[metric]]
        contours[[length(contours) + 1]] <- data.frame(
            cases = curve[, 1], distance = curve[, 2], metric = metric, radius = radius
        )
    }
}
contours <- do.call(rbind, contours)

p <- ggplot(points, aes(cases, distance)) +
    geom_path(data = contours, aes(group = radius), colour = "#386cb0",
              linewidth = 0.5, linetype = "dashed") +
    geom_point(colour = "#1b9e77", size = 1.9) +
    geom_point(data = selected, colour = "#d95f02", size = 2.4) +
    annotate("point", x = 0, y = 0, shape = 3, size = 3, stroke = 1) +
    geom_text(data = selected, aes(label_x, label_y, label = obs),
              size = 3.2, fontface = "bold", colour = "grey15") +
    facet_wrap(~ metric, nrow = 1) +
    coord_fixed(xlim = c(-3.2, 3.6), ylim = c(-3.2, 3.6), expand = FALSE) +
    labs(x = "Cases (standard deviations from mean)",
         y = "Distance (standard deviations from mean)") +
    theme_minimal(base_size = 11) +
    theme(
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 10, face = "bold", colour = "grey15"),
        strip.background = element_rect(fill = "#f0f0f0", colour = NA),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9)
    )

ggsave("distance-contours.pdf", p, device = cairo_pdf, bg = "transparent",
       width = 6.8, height = 3.4)
