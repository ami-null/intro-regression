# Webster, Gunst and Mason data (MPV Table 9.4)
# Source: Webster, J. T., Gunst, R. F. and Mason, R. L. (1974), Latent root regression analysis,
#         Technometrics 16(4), 513-522.
# n = 12, six regressors. x1 + x2 + x3 + x4 is close to 10 in every row: a dependency
# among four regressors that the pairwise correlations do not reveal.

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

web  <- read.csv("webster.csv")
xvar <- paste0("x", 1:6)
X    <- as.matrix(web[, xvar])
p    <- ncol(X)

# ---- Correlation matrix: nothing alarming ---------------------------------------
R <- cor(X)
round(R, 3)
max(abs(R[upper.tri(R)]))                       # 0.505: the largest pairwise |r|

# ---- VIF: severe for x1-x4 ---------------------------------------------------------
fit <- lm(y ~ x1 + x2 + x3 + x4 + x5 + x6, data = web)
summary(fit)                                    # R2 = 0.946; x1-x4 all insignificant, x5 and x6 significant
hatvalues(fit)[1]                               # 1: observation 1 (the only row with sum 11) is fitted exactly
vif <- diag(solve(R))
round(vif, 2)                                   # 182.05 161.36 266.26 297.71 1.92 1.46
round(car::vif(fit), 2)                         # same
round(1 - 1 / vif, 4)                           # R_j^2: 0.9945 0.9938 0.9962 0.9966 0.4792 0.3128

# ---- Eigensystem: one near-zero eigenvalue ------------------------------------------
eig    <- eigen(R, symmetric = TRUE)
lambda <- eig$values
round(lambda, 4)                                # 2.4288 1.5462 0.9221 0.7940 0.3079 0.0011
max(lambda) / min(lambda)                       # kappa = 2196: severe
round(max(lambda) / lambda, 1)                  # condition indices

v6 <- setNames(eig$vectors[, p], xvar)
round(v6, 3)                                    # 0.448 0.421 0.542 0.573 0.006 0.002: only x1-x4 load
rowSums(X[, 1:4])                               # 11, then 10 in every other row

# The eigenvector loadings are for the unit-length-scaled regressors.
# Divide by the column lengths to express the dependency in original units:
W_len <- sqrt(colSums(sweep(X, 2, colMeans(X))^2))
round(v6 / W_len / (v6[1] / W_len[1]), 3)       # 1.00 1.06 1.05 1.05 0.04 0.03: x1 + x2 + x3 + x4 ~ constant

# ---- Variance-decomposition proportions -----------------------------------------
vdp <- function(W) {
    s   <- svd(W)
    phi <- sweep(s$v^2, 2, s$d^2, "/")
    pi  <- t(phi / rowSums(phi))
    colnames(pi) <- colnames(W)
    data.frame(eigenvalue = s$d^2, eta = max(s$d) / s$d, round(pi, 3), check.names = FALSE)
}
W <- sweep(sweep(X, 2, colMeans(X)), 2, W_len, "/")
vdp(W)                                          # last row: eta = 46.9, proportions 0.995-0.998 for x1-x4

# ---- Figure: what one other regressor explains vs. what all the others explain ----
# For each x_j: largest squared pairwise correlation, max_k r_jk^2, against R_j^2 = 1 - 1/VIF_j.
# Both are proportions of the variation in x_j, so they share one 0-1 axis.
diag(R) <- 0
cmp <- data.frame(var = factor(xvar, levels = rev(xvar)),
                  r2_pair = apply(R^2, 2, max), r2_all = 1 - 1 / vif)
round(cmp[, -1], 3)                             # pairwise r^2 <= 0.255; R_j^2 >= 0.994 for x1-x4

cmp_long <- rbind(
    data.frame(var = cmp$var, what = "Best single other regressor", val = cmp$r2_pair),
    data.frame(var = cmp$var, what = "All other regressors", val = cmp$r2_all)
)
cmp_long$what <- factor(cmp_long$what, levels = c("Best single other regressor", "All other regressors"))

p_cmp <- ggplot(cmp, aes(y = var)) +
    geom_segment(aes(x = r2_pair, xend = r2_all, yend = var), colour = "grey70", linewidth = 0.8) +
    geom_point(data = cmp_long, aes(x = val, colour = what), size = 2.2) +
    scale_colour_manual(values = c(line_col, pt_col), name = NULL) +
    scale_y_discrete(labels = parse(text = paste0("x[", 6:1, "]"))) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(x = expression("Share of the variation in " * x[j] * " explained"), y = NULL) +
    theme_lecture() +
    theme(legend.position = "top", legend.text = element_text(size = 8),
          panel.grid.major.y = element_blank())

p_cmp
ggsave("webster_pairwise_vs_rj2.pdf", p_cmp, bg = "transparent", device = cairo_pdf,
       width = 1.2 * fig_w, height = 1.1 * fig_h)

# ---- Figure: eigenvector of the smallest eigenvalue ----------------------------
# The sign of an eigenvector is arbitrary; flip it so the loadings are positive.
v6_plot <- v6 * sign(sum(v6))
ev <- data.frame(var = factor(xvar, levels = xvar), load = v6_plot, lab = paste0("x[", 1:6, "]"))

p_ev <- ggplot(ev, aes(x = var, y = load)) +
    geom_col(width = 0.6, fill = pt_col) +
    geom_hline(yintercept = 0, colour = "grey40") +
    scale_x_discrete(labels = parse(text = ev$lab)) +
    labs(x = NULL, y = expression("Loading in " * v[6] * "  (" * lambda[6] == 0.0011 * ")")) +
    theme_lecture() +
    theme(panel.grid.major.x = element_blank())

p_ev
ggsave("webster_eigvec.pdf", p_ev, bg = "transparent", device = cairo_pdf,
       width = fig_w, height = fig_h)
