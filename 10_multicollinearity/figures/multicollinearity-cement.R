# Hald cement data (MPV Table B.21; also MPV::cement)
# Source: Woods, H., Steinour, H. H. and Starke, H. R. (1932), Effect of composition of Portland
#         cement on heat evolved during hardening, Industrial and Engineering Chemistry 24(11),
#         1207-1214; popularised by Hald (1952).
# y  = heat evolved (calories per gram of cement), n = 13
# x1 = tricalcium aluminate (%), x2 = tricalcium silicate (%),
# x3 = tetracalcium alumino ferrite (%), x4 = dicalcium silicate (%)
# The four ingredients make up almost all of each sample: x1 + x2 + x3 + x4 is 95-99.

cement <- read.csv("cement.csv")
xvar <- paste0("x", 1:4)
X <- as.matrix(cement[, xvar])

rowSums(X)                                       # 99 97 95 97 98 97 97 98 96 98 98 98 98

# ---- Correlation matrix --------------------------------------------------------------
R <- cor(X)
round(R, 3)                                      # r(x2, x4) = -0.973, r(x1, x3) = -0.824

# ---- Fit and VIF -----------------------------------------------------------------------
fit <- lm(y ~ x1 + x2 + x3 + x4, data = cement)
summary(fit)                                     # R2 = 0.982, overall F significant, no t-test significant at 5%
vif <- diag(solve(R))
round(vif, 2)                                    # 38.50 254.42 46.87 282.51

# ---- Eigensystem ------------------------------------------------------------------------
eig <- eigen(R, symmetric = TRUE)
round(eig$values, 4)                             # 2.2357 1.5761 0.1866 0.0016
max(eig$values) / min(eig$values)                # kappa = 1377: severe
round(setNames(eig$vectors[, 4], xvar), 3)       # 0.241 0.642 0.268 0.677: all four, same sign

# In original units the dependency is close to x1 + x2 + x3 + x4 = constant
W_len <- sqrt(colSums(sweep(X, 2, colMeans(X))^2))
v4 <- eig$vectors[, 4] / W_len
round(v4 / v4[1], 2)                             # 1.00 1.01 1.02 0.99

# ---- One remedy: drop a regressor involved in the dependency --------------------------
fit_no4 <- lm(y ~ x1 + x2 + x3, data = cement)
round(car::vif(fit_no4), 2)                      # 3.25 1.06 3.14
summary(fit_no4)$coefficients
c(R2_full = summary(fit)$r.squared, R2_no4 = summary(fit_no4)$r.squared)   # 0.9824 vs 0.9823
