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

# ---- Figure: the sum is well estimated, the difference is not (r = 0.95) ----------
# Same 2000 replicates. Theory: Var(b1 + b2) = 2 sigma^2 / [(n - 1)(1 + r)],
#                               Var(b1 - b2) = 2 sigma^2 / [(n - 1)(1 - r)].
s95_sd <- subset(sims, rho == 0.95)
comb <- rbind(
    data.frame(what = "Sum", value = s95_sd$b1 + s95_sd$b2, centre = beta[2] + beta[3]),
    data.frame(what = "Difference", value = s95_sd$b1 - s95_sd$b2, centre = beta[2] - beta[3])
)
comb_sd <- tapply(comb$value, comb$what, sd)
comb_sd                                                  # Sum 0.19, Difference 1.17
sqrt(2 * sigma^2 / ((n - 1) * c(1 + 0.95, 1 - 0.95)))    # theory: 0.188, 1.174
comb_sd["Difference"] / comb_sd["Sum"]                    # about 6 times wider

comb$panel <- factor(comb$what, levels = c("Sum", "Difference"),
                     labels = c("hat(beta)[1] + hat(beta)[2]~~'(combined credit)'",
                                "hat(beta)[1] - hat(beta)[2]~~'(who closed the deal)'"))
# equal x-ranges in both panels, so the widths can be compared by eye
span   <- 3.6
frame_ <- data.frame(panel = rep(levels(comb$panel), each = 2),
                     value = c(4 - span, 4 + span, -span, span))
frame_$panel <- factor(frame_$panel, levels = levels(comb$panel))
sd_lab <- data.frame(panel = factor(levels(comb$panel), levels = levels(comb$panel)),
                     x = c(4 + span, span) - 0.1,
                     lab = sprintf("sd = %.2f", comb_sd[c("Sum", "Difference")]))

p_comb <- ggplot(comb, aes(x = value)) +
    geom_blank(data = frame_) +
    geom_histogram(binwidth = 0.1, fill = pt_col, colour = NA) +
    geom_vline(aes(xintercept = centre), linetype = "dashed", colour = line_col) +
    geom_text(data = sd_lab, aes(x = x, y = Inf, label = lab), hjust = 1, vjust = 1.6,
              size = 3, colour = "grey20") +
    facet_wrap(~ panel, scales = "free_x", labeller = label_parsed) +
    labs(x = "Estimate (2000 simulated samples, true value dashed)", y = "Count") +
    theme_lecture()

p_comb
ggsave("sim_sum_difference.pdf", p_comb, bg = "transparent", device = cairo_pdf,
       width = 2 * fig_w, height = 1.1 * fig_h)

# ---- Figure: predictions inside and outside the data corridor ---------------------
# r = 0.95 design (salespeople A and B, hours standardized). The slope part of the
# prediction error at a new point (x01, x02) is x01 (b1_hat - b1) + x02 (b2_hat - b2), with
# variance sigma^2 (x01^2 + x02^2 - 2 r x01 x02) / [(1 - r^2) (n - 1)].
rho_c   <- 0.95
new_pts <- data.frame(x1 = c(1.5, 1.5), x2 = c(1.5, -1.5), what = c("joint visit", "A alone"))
new_pts$sd_theory <- with(new_pts, sigma * sqrt((x1^2 + x2^2 - 2 * rho_c * x1 * x2) / ((1 - rho_c^2) * (n - 1))))
s95 <- subset(sims, rho == rho_c)
new_pts$sd_sim <- sapply(seq_len(nrow(new_pts)), function(i)
    sd(new_pts$x1[i] * s95$b1 + new_pts$x2[i] * s95$b2))
new_pts                      # sd: joint visit 0.28, A alone 1.76 (theory and simulation agree)
new_pts$sd_theory[2] / new_pts$sd_theory[1]          # 6.2 times larger outside the corridor
new_pts$lab <- sprintf("%s\nprediction sd %.2f", new_pts$what, new_pts$sd_theory)

corridor <- make_design(rho_c)
p_corr <- ggplot(corridor, aes(x = x1, y = x2)) +
    geom_point(size = 1.2, alpha = 0.6, colour = pt_col) +
    geom_point(data = new_pts, shape = 17, size = 3, colour = line_col) +
    geom_text(data = new_pts, aes(label = lab), hjust = 0, nudge_x = 0.35, size = 2.6,
              lineheight = 0.9, colour = "grey15") +
    coord_equal(xlim = c(-2.6, 4.2), ylim = c(-2.6, 2.6)) +
    labs(x = "Hours with client, salesperson A (standardized)",
         y = "Hours, salesperson B (standardized)") +
    theme_lecture() +
    theme(axis.title = element_text(size = 8))

p_corr
ggsave("sim_prediction_corridor.pdf", p_corr, bg = "transparent", device = cairo_pdf,
       width = 1.2 * fig_w, height = 1.2 * fig_h)

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

# ---- Figure: eigenvalues as the spread of the regressor cloud in a direction -----
# Standardized regressors (sd 1). The sample correlation matrix is [[1, r], [r, 1]], with
# eigenvectors v1 = (1, 1)/sqrt(2), v2 = (1, -1)/sqrt(2) and eigenvalues 1 + r, 1 - r.
# The variance of the data along a unit vector v is v' R v, so along v_j it is lambda_j:
# each segment below runs from -2 sd to +2 sd of the data in that direction.
cloud <- rbind(cbind(make_design(0), rho = 0), cbind(make_design(0.95), rho = 0.95))
sapply(split(cloud, cloud$rho), function(d) eigen(cor(d[, c("x1", "x2")]))$values)
# rho = 0: 1, 1;  rho = 0.95: 1.95, 0.05

eig_dirs <- do.call(rbind, lapply(c(0, 0.95), function(rho) {
    lam <- c(1 + rho, 1 - rho)
    vx  <- c(1, 1) / sqrt(2)                     # x-components of v1, v2
    vy  <- c(1, -1) / sqrt(2)                    # y-components of v1, v2
    data.frame(rho = rho, lam = lam,
               x = -2 * sqrt(lam) * vx, y = -2 * sqrt(lam) * vy,
               xend = 2 * sqrt(lam) * vx, yend = 2 * sqrt(lam) * vy,
               lab = sprintf("lambda[%d] == %s", 1:2, format(lam, nsmall = 2)))
}))
# label position: at the segment end, pushed below the cloud for a very short segment
eig_dirs$lx <- eig_dirs$xend + 0.1
eig_dirs$ly <- eig_dirs$yend - ifelse(eig_dirs$lam < 0.5, 0.55, 0)
panel_lab <- function(rho) factor(rho, levels = c(0, 0.95), labels = c("r[12] == 0", "r[12] == 0.95"))
cloud$panel    <- panel_lab(cloud$rho)
eig_dirs$panel <- panel_lab(eig_dirs$rho)

p_eig <- ggplot(cloud, aes(x = x1, y = x2)) +
    geom_point(size = 0.9, alpha = 0.6, colour = pt_col) +
    geom_segment(data = eig_dirs, aes(x = x, y = y, xend = xend, yend = yend),
                 colour = line_col, linewidth = 0.9) +
    geom_text(data = eig_dirs, aes(x = lx, y = ly, label = lab), parse = TRUE,
              size = 2.6, hjust = 0, colour = "grey15") +
    facet_wrap(~ panel, labeller = label_parsed) +
    coord_equal(xlim = c(-3.2, 4.2), ylim = c(-3.2, 3.2)) +
    labs(x = expression("Standardized " * x[1]), y = expression("Standardized " * x[2])) +
    theme_lecture()

p_eig
ggsave("sim_eigen_directions.pdf", p_eig, bg = "transparent", device = cairo_pdf,
       width = 2 * fig_w, height = 1.2 * fig_h)

# ---- Figure: bias-variance trade-off of shrinking one coefficient --------------
# Illustration (not data): one eigen-direction with true coefficient alpha = 1 and a
# least-squares estimate with sd 1. Ridge reports c * alpha_hat with c = 0.5:
# mean c * alpha = 0.5 (bias -0.5), sd c * 1 = 0.5.
# MSE = variance + bias^2: least squares 1 + 0 = 1, ridge 0.25 + 0.25 = 0.5.
alpha_true <- 1
sd_ls      <- 1
c_shrink   <- 0.5
bv <- data.frame(a = seq(-2.6, 4.2, length.out = 600))
bv <- rbind(
    data.frame(a = bv$a, dens = dnorm(bv$a, alpha_true, sd_ls), est = "Least squares"),
    data.frame(a = bv$a, dens = dnorm(bv$a, c_shrink * alpha_true, c_shrink * sd_ls), est = "Ridge")
)
bv$est <- factor(bv$est, levels = c("Least squares", "Ridge"))
est_col <- c("Least squares" = pt_col, "Ridge" = line_col)
mse <- c(sd_ls^2, (c_shrink * sd_ls)^2 + ((c_shrink - 1) * alpha_true)^2)    # 1, 0.5

p_bv <- ggplot(bv, aes(x = a, y = dens, colour = est)) +
    geom_vline(xintercept = alpha_true, linetype = "dashed", colour = "grey45") +
    geom_line(linewidth = 0.9) +
    annotate("text", x = alpha_true + 0.08, y = 0.83, hjust = 0, size = 2.8, colour = "grey30",
             label = "true value") +
    annotate("text", x = 2.25, y = 0.30, hjust = 0, size = 2.7, colour = pt_col,
             label = sprintf("least squares\nno bias, sd 1\nMSE = %.2f", mse[1]), lineheight = 0.9) +
    annotate("text", x = -0.55, y = 0.62, hjust = 1, size = 2.7, colour = line_col,
             label = sprintf("ridge\nbias -0.5, sd 0.5\nMSE = %.2f", mse[2]), lineheight = 0.9) +
    scale_colour_manual(values = est_col, guide = "none") +
    labs(x = "Estimate of the coefficient in one direction", y = "Density") +
    theme_lecture()

p_bv
ggsave("ridge_bias_variance.pdf", p_bv, bg = "transparent", device = cairo_pdf,
       width = fig_w, height = fig_h)
