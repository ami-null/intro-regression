# Gasoline mileage data (MPV Table B.3; source: Motor Trend, 1975)
# 32 automobiles, y = miles per gallon. Torque (x3) is missing for two cars: n = 30 complete cases.
# x1 displacement, x2 horsepower, x3 torque, x4 compression ratio, x5 rear axle ratio,
# x6 carburetor barrels, x7 transmission speeds, x8 overall length, x9 width, x10 weight,
# x11 transmission type (A/M, not used here)
# Values follow the 5th edition (Matador x6 = 4; the MPV package differs in rows 2 and 32).

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

gasoline <- read.csv("gasoline.csv", na.strings = "NA")
xvar <- paste0("x", 1:10)
xlab <- c(x1 = "Displacement", x2 = "Horsepower", x3 = "Torque", x4 = "Compression",
          x5 = "Axle ratio", x6 = "Carburetor", x7 = "Speeds", x8 = "Length",
          x9 = "Width", x10 = "Weight")
gas <- na.omit(gasoline[, c("automobile", "y", xvar)])
nrow(gas)                                          # 30

# ---- Lecture 9 puzzle: every t-test insignificant, overall F large ---------------
fit5 <- lm(y ~ x2 + x3 + x5 + x6 + x10, data = gas)
summary(fit5)                                      # overall F = 14.66, all |t| < 1.31
round(car::vif(fit5), 2)                           # hp 29.05, tq 47.34, ax 3.18, cb 3.80, wt 8.98

# ---- Full model with all ten regressors -----------------------------------------
fit10 <- lm(reformulate(xvar, "y"), data = gas)
summary(fit10)                                     # R2 = 0.829, F = 9.18 on (10, 19) df, no |t| > 1.96
round(coef(fit10), 4)
# Signs: displacement -0.100 but torque +0.126, although both correlate about -0.85 with mpg
round(cor(gas$y, gas[, xvar]), 3)                  # x1 -0.872, x3 -0.850

# ---- Correlation matrix -----------------------------------------------------------
R <- cor(gas[, xvar])
round(R, 2)
big <- which(abs(R) > 0.9 & upper.tri(R), arr.ind = TRUE)
data.frame(a = unname(xlab[rownames(R)[big[, 1]]]), b = unname(xlab[colnames(R)[big[, 2]]]), r = round(R[big], 3))
# disp-hp 0.941, disp-tq 0.990, hp-tq 0.964, disp-wt 0.924, tq-wt 0.910, length-wt 0.934

# ---- VIF ------------------------------------------------------------------------
vif <- diag(solve(R))
round(vif, 1)                                      # 99.9 44.0 161.2 2.0 8.1 5.7 9.9 11.4 7.7 21.5
round(car::vif(fit10), 1)                          # same
sum(vif > 10)                                      # 5: displacement, horsepower, torque, length, weight

# ---- Eigensystem ------------------------------------------------------------------
eig    <- eigen(R, symmetric = TRUE)
lambda <- eig$values
round(lambda, 4)                                   # 6.9275 1.3571 0.7562 0.5686 0.1539 0.0999 0.0736 0.0390 0.0204 0.0038
max(lambda) / min(lambda)                          # kappa = 1820: severe
round(max(lambda) / lambda, 1)                     # condition indices: three above 100, one above 1000
V <- eig$vectors
dimnames(V) <- list(xlab[xvar], paste0("v", 1:10))
round(V[, 10], 3)                                  # smallest eigenvalue: displacement 0.55, horsepower 0.31, torque -0.77
round(V[, 9], 3)

vdp <- function(W) {
    s   <- svd(W)
    phi <- sweep(s$v^2, 2, s$d^2, "/")
    pi  <- t(phi / rowSums(phi))
    colnames(pi) <- colnames(W)
    data.frame(eigenvalue = s$d^2, eta = max(s$d) / s$d, round(pi, 3), check.names = FALSE)
}
Xc <- sweep(as.matrix(gas[, xvar]), 2, colMeans(gas[, xvar]))
W  <- sweep(Xc, 2, sqrt(colSums(Xc^2)), "/")
vdp(W)                                             # last row (eta = 42.7): displacement 0.79, horsepower 0.56, torque 0.96

# ---- Figure: correlation heatmap --------------------------------------------------
cor_long <- data.frame(
    row = factor(rep(xlab[xvar], times = 10), levels = rev(xlab[xvar])),
    col = factor(rep(xlab[xvar], each = 10), levels = xlab[xvar]),
    r   = as.vector(R)
)
p_heat <- ggplot(cor_long, aes(x = col, y = row, fill = r)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.2f", r)), size = 2.1, colour = "grey10") +
    scale_fill_gradient2(low = line_col, mid = "grey96", high = pt_col, midpoint = 0,
                         limits = c(-1, 1), name = "r") +
    coord_equal() +
    labs(x = NULL, y = NULL) +
    theme_lecture() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          panel.grid = element_blank(),
          legend.key.height = unit(0.8, "cm"))

p_heat
ggsave("gas_corr_heatmap.pdf", p_heat, bg = "transparent", device = cairo_pdf,
       width = 1.6 * fig_w, height = 1.6 * fig_w)

# ---- Figure: scatterplot matrix of the four size variables ---------------------
sv <- c("x1", "x2", "x3", "x10")
pairs_long <- do.call(rbind, lapply(sv, function(a) do.call(rbind, lapply(sv, function(b) {
    if (a == b) return(NULL)
    data.frame(vx = unname(xlab[a]), vy = unname(xlab[b]), x = gas[[a]], y = gas[[b]])
}))))
pairs_long$vx <- factor(pairs_long$vx, levels = xlab[sv])
pairs_long$vy <- factor(pairs_long$vy, levels = xlab[sv])
r_lab <- unique(pairs_long[, c("vx", "vy")])
r_lab$r <- mapply(function(a, b) cor(gas[[names(xlab)[xlab == a]]], gas[[names(xlab)[xlab == b]]]),
                  as.character(r_lab$vx), as.character(r_lab$vy))

p_pairs <- ggplot(pairs_long, aes(x = x, y = y)) +
    geom_point(size = 0.9, colour = pt_col) +
    geom_text(data = r_lab, aes(x = -Inf, y = Inf, label = sprintf("r = %.2f", r)),
              hjust = -0.1, vjust = 1.4, size = 2.4, colour = "grey20") +
    facet_grid(vy ~ vx, scales = "free", switch = "both") +
    labs(x = NULL, y = NULL) +
    theme_lecture() +
    theme(strip.placement = "outside", strip.text = element_text(size = 8),
          axis.text = element_text(size = 6))

p_pairs
ggsave("gas_pairs.pdf", p_pairs, bg = "transparent", device = cairo_pdf,
       width = 1.6 * fig_w, height = 1.6 * fig_w)

# ---- Figure: VIFs ------------------------------------------------------------------
# Values go in the axis labels so they never collide with the reference lines at 5 and 10
vif_df <- data.frame(var = sprintf("%s  %.1f", xlab[xvar], vif), vif = vif)
vif_df$var <- factor(vif_df$var, levels = vif_df$var[order(vif_df$vif)])

p_vif <- ggplot(vif_df, aes(x = vif, y = var)) +
    geom_vline(xintercept = c(5, 10), linetype = "dashed", colour = line_col) +
    geom_col(width = 0.65, fill = pt_col) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.04))) +
    labs(x = "Variance inflation factor (dashed: 5 and 10)", y = NULL) +
    theme_lecture() +
    theme(panel.grid.major.y = element_blank())

p_vif
ggsave("gas_vif_bar.pdf", p_vif, bg = "transparent", device = cairo_pdf,
       width = 1.2 * fig_w, height = 1.2 * fig_h)

# ---- Figure: condition indices -----------------------------------------------------
ci_df <- data.frame(j = 1:10, kappa_j = max(lambda) / lambda)

p_ci <- ggplot(ci_df, aes(x = j, y = kappa_j)) +
    geom_hline(yintercept = c(100, 1000), linetype = "dashed", colour = line_col) +
    geom_segment(aes(xend = j, y = 1, yend = kappa_j), colour = "grey70") +
    geom_point(size = 2.2, colour = pt_col) +
    annotate("text", x = 1, y = c(100, 1000), vjust = -0.5, hjust = 0, size = 2.6, colour = "grey30",
             label = c("moderate (100)", "severe (1000)")) +
    scale_y_log10(breaks = c(1, 10, 100, 1000)) +
    scale_x_continuous(breaks = 1:10) +
    labs(x = expression("Eigenvalue index " * j), y = expression(kappa[j] == lambda[max] / lambda[j])) +
    theme_lecture() +
    theme(panel.grid.minor = element_blank())

p_ci
ggsave("gas_condition_indices.pdf", p_ci, bg = "transparent", device = cairo_pdf,
       width = fig_w, height = fig_h)
