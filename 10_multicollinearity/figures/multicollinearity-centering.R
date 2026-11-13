# Two short examples
# (1) Delivery time data (MPV Example 3.1): mild collinearity, VIF = 3.12
# (2) Windmill data (MPV Table 5.5; source: Joglekar, Schuenemeyer and LaRiccia (1989),
#     The American Statistician 43(3), 135-143): collinearity created by the model itself
#     (x and x^2) and removed by centering

library(ggplot2)

pt_col   <- "#1b9e77"
line_col <- "#d95f02"
fig_w    <- 3.4
fig_h    <- 2.5

theme_lecture <- function(){
    theme_minimal(base_size = 11) +
    theme(
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background  = element_rect(fill = "transparent", colour = NA)
    )
}

# ---- (1) Delivery time data ------------------------------------------------------------
delivery <- read.csv("delivery.csv")
r12 <- cor(delivery$cases, delivery$distance)
r12                                              # 0.824
1 / (1 - r12^2)                                  # VIF = 3.12 for both regressors
fit_d <- lm(time ~ cases + distance, data = delivery)
car::vif(fit_d)                                  # same
eig_d <- eigen(cor(delivery[, c("cases", "distance")]))$values
eig_d                                            # 1 + r, 1 - r = 1.824, 0.176
max(eig_d) / min(eig_d)                          # kappa = 10.4: no problem

# ---- (2) Windmill data: quadratic in wind velocity ------------------------------------
wm <- read.csv("windmill.csv")
wm$xc <- wm$x - mean(wm$x)
mean(wm$x)                                       # 6.13

cor(wm$x, wm$x^2)                                # 0.987
cor(wm$xc, wm$xc^2)                              # 0.178

fit_raw <- lm(y ~ x + I(x^2), data = wm)
fit_cen <- lm(y ~ xc + I(xc^2), data = wm)
car::vif(fit_raw)                                # 38.5
car::vif(fit_cen)                                # 1.03
# Same fitted values, same R^2, same quadratic coefficient; only the parameterization changes
max(abs(fitted(fit_raw) - fitted(fit_cen)))      # ~ 0
rbind(raw = coef(fit_raw), centered = coef(fit_cen))
rbind(raw = coef(summary(fit_raw))[, 2], centered = coef(summary(fit_cen))[, 2])   # se of linear term: 0.0614 -> 0.0101

# ---- Figure: x vs. x^2, before and after centering --------------------------------------
cen_long <- rbind(
    data.frame(panel = "Raw", u = wm$x, u2 = wm$x^2),
    data.frame(panel = "Centered", u = wm$xc, u2 = wm$xc^2)
)
cen_long$panel <- factor(cen_long$panel, levels = c("Raw", "Centered"),
                         labels = c(sprintf("Raw: r = %.3f", cor(wm$x, wm$x^2)),
                                    sprintf("Centered: r = %.3f", cor(wm$xc, wm$xc^2))))

p_cen <- ggplot(cen_long, aes(x = u, y = u2)) +
    geom_point(size = 1.8, colour = pt_col) +
    facet_wrap(~ panel, scales = "free") +
    labs(x = "Wind velocity (raw or centered)", y = "Its square") +
    theme_lecture()

p_cen
ggsave("windmill_centering.pdf", p_cen, bg = "transparent", device = cairo_pdf,
       width = 2 * fig_w, height = 1.1 * fig_h)
