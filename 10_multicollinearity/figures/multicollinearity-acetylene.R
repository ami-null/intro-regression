# Acetylene data (MPV Table 9.1, Example 9.1)
# Source: Kunugi, Tamura and Naito (1961); analysed by Marquardt and Snee (1975)
# y  = conversion of n-heptane to acetylene (%)
# x1 = reactor temperature (deg C), x2 = ratio of H2 to n-heptane (mole ratio), x3 = contact time (s)
# Diagnostics (correlation matrix, VIF, eigensystem) and remedies (ridge regression, PCR)

library(ggplot2)

pt_col   <- "#1b9e77"
line_col <- "#d95f02"
type_col <- c(Linear = "#1b9e77", Interaction = "#d95f02", Quadratic = "#7570b3")
fig_w    <- 3.4
fig_h    <- 2.5

theme_lecture <- function(){
    theme_minimal(base_size = 11) +
    theme(
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background  = element_rect(fill = "transparent", colour = NA)
    )
}

acet <- read.csv("acetylene.csv")
n <- nrow(acet)                                        # 16

# ---- Correlation matrix of the original variables -----------------------------
round(cor(acet[, c("y", "x1", "x2", "x3")]), 3)        # r(x1, x3) = -0.958

# ---- Full quadratic model in unit-normal-scaled regressors (as in MPV) --------
z <- as.data.frame(scale(acet[, c("x1", "x2", "x3")]))
X <- with(z, data.frame(
    x1 = x1, x2 = x2, x3 = x3,
    x1x2 = x1 * x2, x1x3 = x1 * x3, x2x3 = x2 * x3,
    x1sq = x1^2, x2sq = x2^2, x3sq = x3^2
))
term_lab  <- c("x[1]", "x[2]", "x[3]", "x[1]*x[2]", "x[1]*x[3]", "x[2]*x[3]",
               "x[1]^2", "x[2]^2", "x[3]^2")
term_type <- rep(c("Linear", "Interaction", "Quadratic"), each = 3)
p <- ncol(X)                                           # 9 regressors, 10 parameters

fit_ls <- lm(acet$y ~ ., data = X)
summary(fit_ls)                                        # R2 = 0.9977, MSRes = 0.8126 on 6 df
round(coef(fit_ls), 3)
# 35.896, 4.004, 2.778, -8.042, -6.457, -26.980, -3.768, -12.524, -0.973, -11.593

# ---- Unit-length scaling: W'W is the correlation matrix of the regressors -----
unit_length <- function(M) {
    M <- as.matrix(M)
    Mc <- sweep(M, 2, colMeans(M))
    sweep(Mc, 2, sqrt(colSums(Mc^2)), "/")
}
W   <- unit_length(X)
y_c <- acet$y - mean(acet$y)
SST <- sum(y_c^2)
y_s <- y_c / sqrt(SST)                                 # unit-length response
WtW <- crossprod(W)
Wty <- crossprod(W, y_s)
round(WtW, 3)                                          # correlation matrix of the 9 terms

# ---- VIF: diagonal of (W'W)^{-1} -----------------------------------------------
vif <- diag(solve(WtW))
round(vif, 2)
# 375.25, 1.74, 680.28, 31.04, 6563.35, 35.61, 1762.58, 3.16, 1156.77
round(car::vif(fit_ls), 2)                             # same values
round(1 - 1 / vif, 4)                                  # R_j^2 of each term on the others

# ---- Eigensystem analysis -------------------------------------------------------
eig <- eigen(WtW, symmetric = TRUE)
lambda <- eig$values
round(lambda, 4)                       # 4.2052 2.1620 1.1387 1.0405 0.3852 0.0495 0.0136 0.0051 0.0001
kappa <- max(lambda) / min(lambda)
kappa                                  # 43381: severe (> 1000)
round(max(lambda) / lambda, 1)         # condition indices kappa_j = lambda_max / lambda_j
sum(lambda < 0.05)                     # number of near-dependencies (small eigenvalues)
prod(lambda)                           # det(W'W) = 1.4e-09, close to 0
sum(1 / lambda)                        # E||b - beta||^2 / sigma^2 in the unit-length scale = sum of VIFs
sum(vif)

V <- eig$vectors
dimnames(V) <- list(colnames(W), paste0("v", 1:p))
round(V[, p], 3)                       # eigenvector of the smallest eigenvalue: the dependency W v9 ~ 0
max(abs(W %*% V[, p]))                 # 0.0045: this linear combination of the columns is almost zero
# loadings: x1x3 0.80, x1^2 0.41, x3^2 0.33, x3 0.24, x1 0.17 (temperature and contact time terms)

# ---- Variance-decomposition proportions (Belsley, Kuh and Welsch) ------------
# SVD W = U D V'; Var(b_j) is proportional to sum_k v_jk^2 / mu_k^2.
# pi_kj = share of Var(b_j) that comes from component k. Condition index eta_k = mu_max / mu_k.
vdp <- function(W) {
    s   <- svd(W)
    phi <- sweep(s$v^2, 2, s$d^2, "/")               # rows: terms, columns: components
    pi  <- t(phi / rowSums(phi))                     # rows: components, columns: terms
    colnames(pi) <- colnames(W)
    data.frame(eigenvalue = s$d^2, eta = max(s$d) / s$d, round(pi, 3), check.names = FALSE)
}
vdp(W)
# Belsley: eta > 30 with two or more proportions > 0.5 in the same row flags a dependency
# Row 9 (eta = 208): x1, x3, x1x3, x1^2, x3^2 all > 0.83.  Row 7 (eta = 17.6): x1x2 0.82, x2x3 0.79.

# ---- Figure: how the data were collected (temperature vs. contact time) -------
p_x1x3 <- ggplot(acet, aes(x = x1, y = x3)) +
    geom_point(size = 2, colour = pt_col) +
    labs(x = "Reactor temperature (°C)", y = "Contact time (s)") +
    scale_x_continuous(breaks = c(1100, 1200, 1300)) +
    theme_lecture()

p_x1x3
ggsave("acetylene_x1_x3.pdf", p_x1x3, bg = "transparent", device = cairo_pdf,
       width = fig_w, height = fig_h)

# ---- Ridge regression in the unit-length scale --------------------------------
# b_R(k) = (W'W + kI)^{-1} W'y,   k >= 0
ridge <- function(k) {
    A   <- WtW + k * diag(p)
    Ai  <- solve(A)
    b   <- drop(Ai %*% Wty)
    res <- y_s - drop(W %*% b)
    list(k = k, b = setNames(b, colnames(W)),
         ssres = sum(res^2) * SST,                     # back in original y units
         r2 = 1 - sum(res^2),
         vif = diag(Ai %*% WtW %*% Ai))
}

b_ls     <- ridge(0)$b                                 # standardized LS coefficients
ms_s     <- sum((y_s - W %*% b_ls)^2) / (n - p - 1)    # MSRes in the unit-length scale
k_hkb    <- p * ms_s / sum(b_ls^2)                     # Hoerl, Kennard and Baldwin (1975)
k_hkb                                                  # 0.0005: very little shrinkage

k_grid <- sort(unique(c(round(c(seq(0, 0.002, by = 0.00002), seq(0, 0.05, by = 0.0005)), 6), k_hkb)))
trace  <- do.call(rbind, lapply(k_grid, function(k) {
    r <- ridge(k)
    data.frame(k = k, term = names(r$b), b = r$b, r2 = r$r2, ssres = r$ssres, maxvif = max(r$vif))
}))
trace$type  <- factor(term_type[match(trace$term, colnames(W))], levels = names(type_col))
trace$label <- term_lab[match(trace$term, colnames(W))]

# Value of k read from the ridge trace (the coefficients have settled down)
k_sel <- 0.032
fit_k <- ridge(k_sel)
round(cbind(LS = b_ls, ridge = fit_k$b), 4)            # x1x3: -2.034 -> -0.068
c(r2_ls = ridge(0)$r2, r2_ridge = fit_k$r2)            # 0.9977 -> 0.9944: small loss in fit
c(ssres_ls = ridge(0)$ssres, ssres_ridge = fit_k$ssres)   # 4.88 -> 11.92
round(fit_k$vif, 2)                                    # VIFs of the ridge estimator at k_sel: all < 4.2
k_tab <- unique(trace[, c("k", "r2", "maxvif")])
print(k_tab[round(k_tab$k, 6) %in% c(0, 0.001, 0.005, 0.01, 0.02, 0.032, 0.05), ], row.names = FALSE)

# Cross-check with MASS::lm.ridge: it scales by sd with divisor n, so lambda = n * k
rr <- MASS::lm.ridge(acet$y ~ ., data = X, lambda = n * k_sel)
b_orig <- fit_k$b * sqrt(SST) / sqrt(colSums(sweep(as.matrix(X), 2, colMeans(X))^2))
round(cbind(ours = b_orig, MASS = coef(rr)[-1]), 4)    # identical

# ---- Figure: ridge trace ---------------------------------------------------------
# Square-root k axis: the coefficients change fastest for small k. Lines are labelled at k = 0 (least squares).
lab_at <- subset(trace, k == 0)
p_trace <- ggplot(trace, aes(x = k, y = b, group = term, colour = type)) +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_vline(xintercept = k_sel, linetype = "dashed", colour = "grey40") +
    geom_line(linewidth = 0.7) +
    geom_text(data = lab_at, aes(label = label), parse = TRUE, hjust = 1, nudge_x = -0.004,
              size = 2.6, show.legend = FALSE) +
    scale_colour_manual(values = type_col, name = NULL) +
    scale_x_sqrt(breaks = c(0, 0.001, 0.005, 0.01, 0.02, 0.032, 0.05),
                 labels = c("0", "0.001", "0.005", "0.01", "0.02", "0.032", "0.05"),
                 expand = expansion(mult = c(0.12, 0.02))) +
    labs(x = "Biasing parameter k (square-root scale)", y = "Standardized coefficient") +
    theme_lecture() +
    theme(legend.position = "top")

p_trace
ggsave("acetylene_ridge_trace.pdf", p_trace, bg = "transparent", device = cairo_pdf,
       width = 1.4 * fig_w, height = 1.4 * fig_h)

# ---- Principal-component regression ---------------------------------------------
# Z = W V (orthogonal columns, Z'Z = Lambda), alpha_hat = Lambda^{-1} Z'y
Zpc       <- W %*% V
alpha_hat <- drop(crossprod(Zpc, y_s)) / lambda
round(alpha_hat, 4)
round(alpha_hat / sqrt(ms_s / lambda), 2)              # alpha_hat / se: -36.9, 0.4, 32.8, -12.2, 0.3, 2.5, 2.3, 1.9, -1.2

pcr <- function(m) {                                   # keep the first m components
    b   <- drop(V[, 1:m, drop = FALSE] %*% alpha_hat[1:m])
    res <- y_s - drop(W %*% b)
    list(m = m, b = setNames(b, colnames(W)), r2 = 1 - sum(res^2), ssres = sum(res^2) * SST,
         vif = colSums((t(V[, 1:m, drop = FALSE])^2) / lambda[1:m]))
}
pcr_tab <- sapply(p:4, function(m) { r <- pcr(m); c(r$b, R2 = r$r2, SSRes = r$ssres) })
colnames(pcr_tab) <- paste0("m=", p:4)
round(pcr_tab, 4)                                      # m = 9 reproduces least squares; R2 0.9977 (m=9), 0.9938 (m=6)
round(pcr(6)$vif, 2)                                   # VIFs after dropping lambda_7, lambda_8, lambda_9: all < 8.3

# ---- Figure: PCR coefficients vs. number of components retained ----------------
pcr_long <- do.call(rbind, lapply(p:4, function(m) {
    r <- pcr(m)
    data.frame(m = m, term = names(r$b), b = r$b)
}))
pcr_long$type  <- factor(term_type[match(pcr_long$term, colnames(W))], levels = names(type_col))
pcr_long$label <- term_lab[match(pcr_long$term, colnames(W))]

p_pcr <- ggplot(pcr_long, aes(x = m, y = b, group = term, colour = type)) +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.4) +
    geom_text(data = subset(pcr_long, m == p), aes(label = label), parse = TRUE, hjust = 1,
              nudge_x = -0.15, size = 2.6, show.legend = FALSE) +
    scale_colour_manual(values = type_col, name = NULL) +
    scale_x_reverse(breaks = p:4, expand = expansion(mult = c(0.14, 0.03))) +
    labs(x = "Number of principal components retained", y = "Standardized coefficient") +
    theme_lecture() +
    theme(legend.position = "top")

p_pcr
ggsave("acetylene_pcr.pdf", p_pcr, bg = "transparent", device = cairo_pdf,
       width = 1.4 * fig_w, height = 1.4 * fig_h)
