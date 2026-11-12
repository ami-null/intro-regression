library(ggplot2)

delivery <- read.csv("delivery.csv")

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

# ---- Models ------------------------------------------------------------------
fit_full <- lm(time ~ cases + distance, data = delivery)
fit_x1   <- lm(time ~ cases, data = delivery)      # reduced model: test distance
fit_x2   <- lm(time ~ distance, data = delivery)   # reduced model: test cases

summary(fit_full)
deviance(fit_full)                  # SSRes(FM) = 233.73
deviance(fit_x1)                    # SSRes(RM), x1 only = 402.13
deviance(fit_x2)                    # SSRes(RM), x2 only = 1185.39

# ---- Partial F-tests (extra sum of squares) -----------------------------------
anova(fit_x1, fit_full)             # distance given cases: extra SS 168.40, F = 15.85
anova(fit_x2, fit_full)             # cases given distance: extra SS 951.66, F = 89.58
coef(summary(fit_full))[, "t value"]^2   # t^2 equals the partial F for each coefficient

# ---- Sequential vs. partial sums of squares -----------------------------------
anova(fit_full)                                       # sequential, cases first
anova(lm(time ~ distance + cases, data = delivery))   # sequential, distance first
drop1(fit_full, test = "F")                           # partial (each term entered last)
car::Anova(fit_full, type = "III")                  # same partial F-tests

# ---- Orthogonal columns: replace distance by e(distance | cases) --------------
delivery$z <- resid(lm(distance ~ cases, data = delivery))
anova(lm(time ~ cases + z, data = delivery))          # 5382.41, 168.40
anova(lm(time ~ z + cases, data = delivery))          # 168.40, 5382.41 (order irrelevant)

# ---- Added-variable data ------------------------------------------------------
av_distance <- data.frame(
    obs = delivery$obs,
    ex  = resid(lm(distance ~ cases, data = delivery)),
    ey  = resid(lm(time ~ cases, data = delivery))
)
av_cases <- data.frame(
    obs = delivery$obs,
    ex  = resid(lm(cases ~ distance, data = delivery)),
    ey  = resid(lm(time ~ distance, data = delivery))
)

# Properties: slope = full-model coefficient, residuals = full-model residuals
av_fit <- lm(ey ~ ex, data = av_distance)
coef(av_fit)                                  # intercept 0, slope 0.014385
coef(fit_full)["distance"]
max(abs(resid(av_fit) - resid(fit_full)))     # ~ 0
sum(av_distance$ex^2) * coef(fit_full)["distance"]^2    # 168.40 = SSR(beta2 | beta1)
summary(av_fit)$r.squared                     # 0.4188 = 168.40 / 402.13

# cross-check with the car package
car::avPlots(fit_full)

# ---- Figures ------------------------------------------------------------------
# Raw scatter: time vs distance (marginal slope 0.0426)
p_raw <- ggplot(delivery, aes(x = distance, y = time)) +
    geom_point(size = 1.8, colour = pt_col) +
    geom_smooth(method = "lm", se = FALSE, colour = line_col, linewidth = 0.8) +
    labs(x = "Distance (ft)", y = "Delivery time (min)") +
    theme_lecture()
p_raw
ggsave("delivery_raw_distance.pdf", p_raw, bg = "transparent", width = fig_w, height = fig_h)

# Added-variable plot for distance given cases
p_av_distance <- ggplot(av_distance, aes(x = ex, y = ey)) +
    geom_point(size = 1.8, colour = pt_col) +
    geom_smooth(method = "lm", se = FALSE, colour = line_col, linewidth = 0.8) +
    geom_text(data = subset(av_distance, obs %in% c(9, 10, 16)), aes(label = obs),
            hjust = -0.4, vjust = -0.4, size = 3) +
    labs(x = "e(x2 | x1): distance adjusted for cases",
       y = "e(y | x1): time adjusted for cases") +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.12))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
    theme_lecture() +
    theme(axis.title = element_text(size = 9))

p_av_distance
ggsave("delivery_av_distance.pdf", p_av_distance, bg = "transparent", width = fig_w, height = fig_h)

# Added-variable plot for cases given distance
p_av_cases <- ggplot(av_cases, aes(x = ex, y = ey)) +
    geom_point(size = 1.8, colour = pt_col) +
    geom_smooth(method = "lm", se = FALSE, colour = line_col, linewidth = 0.8) +
    geom_text(data = subset(av_cases, obs %in% c(9, 22)), aes(label = obs),
            hjust = -0.4, vjust = -0.4, size = 3) +
    labs(x = "e(x1 | x2): cases adjusted for distance",
       y = "e(y | x2): time adjusted for distance") +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.12))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
    theme_lecture() +
    theme(axis.title = element_text(size = 9))

p_av_cases
ggsave("delivery_av_cases.pdf", p_av_cases, bg = "transparent", width = fig_w, height = fig_h)

# Partial residual plot for distance
delivery$partial_res <- resid(fit_full) + coef(fit_full)["distance"] * delivery$distance
p_pr <- ggplot(delivery, aes(x = distance, y = partial_res)) +
    geom_point(size = 1.8, colour = pt_col) +
    geom_smooth(method = "lm", se = FALSE, colour = line_col, linewidth = 0.8) +
    labs(x = "Distance (ft)", y = "Partial residual (min)") +
    theme_lecture()

p_pr
ggsave("delivery_pr_distance.pdf", p_pr, bg = "transparent", width = fig_w, height = fig_h)
