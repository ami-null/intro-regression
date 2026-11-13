# Simulated data: what multicollinearity does to the sampling distribution of the OLS estimates
# Model: y = 1 + 2 x1 + 2 x2 + e, e ~ N(0, 1), n = 30, X held fixed across replicates

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

set.seed(2026)
n     <- 30
beta  <- c(1, 2, 2)
sigma <- 1
n_rep <- 2000

# ---- Fixed design with sample correlation exactly rho -------------------------
# z1, z2: centered, orthogonal, each with sum of squares n - 1 (sample sd = 1)
Z  <- matrix(rnorm(2 * n), n)
z1 <- Z[, 1] - mean(Z[, 1])
z2 <- resid(lm(Z[, 2] ~ z1))
z1 <- z1 / sqrt(sum(z1^2) / (n - 1))
z2 <- z2 / sqrt(sum(z2^2) / (n - 1))

make_design <- function(rho) data.frame(x1 = z1, x2 = rho * z1 + sqrt(1 - rho^2) * z2)

simulate_betas <- function(rho) {
    X <- make_design(rho)
    mu <- beta[1] + beta[2] * X$x1 + beta[3] * X$x2
    b <- t(replicate(n_rep, {
        y <- mu + rnorm(n, sd = sigma)
        coef(lm(y ~ x1 + x2, data = X))[2:3]
    }))
    data.frame(rho = rho, b1 = b[, 1], b2 = b[, 2])
}

sims <- rbind(simulate_betas(0), simulate_betas(0.95))

# ---- Empirical vs. theoretical spread -----------------------------------------
# Var(b_j) = sigma^2 / [(1 - rho^2) S_jj], S_jj = n - 1;  Corr(b1, b2) = -rho
summary_tab <- do.call(rbind, lapply(split(sims, sims$rho), function(d) {
    rho <- d$rho[1]
    data.frame(
        rho        = rho,
        cor_x      = cor(make_design(rho))[1, 2],
        sd_b1      = sd(d$b1),
        sd_b1_theo = sigma / sqrt((1 - rho^2) * (n - 1)),
        cor_b      = cor(d$b1, d$b2),
        VIF        = 1 / (1 - rho^2)
    )
}))
print(summary_tab, digits = 4)
# rho = 0:    sd(b1) = 0.182 (theory 0.186), corr(b1, b2) =  0.00, VIF = 1
# rho = 0.95: sd(b1) = 0.596 (theory 0.595), corr(b1, b2) = -0.95, VIF = 10.26
# sd ratio 0.595 / 0.186 = 3.2 = sqrt(VIF)

# share of replicates in which an estimate has the wrong sign (true value 2)
aggregate(cbind(neg_b1 = b1 < 0, neg_b2 = b2 < 0) ~ rho, data = sims, FUN = mean)   # 0 and 0.0005

# ---- Figure: joint sampling distribution of (b1, b2) ---------------------------
sims$panel <- factor(sims$rho, labels = c("r[12] == 0", "r[12] == 0.95"))

p_scatter <- ggplot(sims, aes(x = b1, y = b2)) +
    geom_point(size = 0.6, alpha = 0.25, colour = pt_col) +
    annotate("point", x = beta[2], y = beta[3], colour = line_col, size = 2.2, shape = 4, stroke = 1.2) +
    facet_wrap(~ panel, labeller = label_parsed) +
    coord_equal() +
    labs(x = expression(hat(beta)[1]), y = expression(hat(beta)[2])) +
    theme_lecture()

p_scatter
ggsave("sim_beta_scatter.pdf", p_scatter, bg = "transparent", device = cairo_pdf, width = 2 * fig_w, height = 1.2 * fig_h)

# ---- Figure: VIF as a function of the correlation between two regressors -------
vif_curve <- data.frame(r = seq(0, 0.99, by = 0.001))
vif_curve$vif <- 1 / (1 - vif_curve$r^2)

r_delivery <- 0.824                       # cases vs distance, delivery data
p_vif <- ggplot(vif_curve, aes(x = r, y = vif)) +
    geom_hline(yintercept = c(5, 10), linetype = "dashed", colour = "grey55") +
    geom_line(colour = "grey25", linewidth = 0.8) +
    annotate("point", x = r_delivery, y = 1 / (1 - r_delivery^2), colour = line_col, size = 2.2) +
    annotate("text", x = r_delivery + 0.03, y = 1 / (1 - r_delivery^2), hjust = 0, vjust = 1.6,
             size = 2.8, label = "delivery\ndata", lineheight = 0.9) +
    annotate("text", x = 0.02, y = c(5, 10), vjust = -0.4, hjust = 0, size = 2.8, colour = "grey40",
             label = c("VIF = 5", "VIF = 10")) +
    scale_y_log10(breaks = c(1, 2, 5, 10, 20, 50)) +
    labs(x = expression("Correlation " * r[12]), y = "VIF (log scale)") +
    theme_lecture()

p_vif
ggsave("sim_vif_curve.pdf", p_vif, bg = "transparent", device = cairo_pdf, width = fig_w, height = fig_h)
