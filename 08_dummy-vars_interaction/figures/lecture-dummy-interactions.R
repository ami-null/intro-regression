# Dummy variables and interaction terms: tool life data (MPV, Chapter 8)

# ---- Data ------------------------------------------------------------------
tool <- data.frame(
    y  = c(18.73, 14.52, 17.43, 14.54, 13.44, 24.39, 13.34, 22.71, 12.68, 19.32,
           30.16, 27.09, 25.40, 26.05, 33.49, 35.62, 26.07, 36.78, 34.95, 43.67),
    x1 = c(610, 950, 720, 840, 980, 530, 680, 540, 890, 730,
           670, 770, 880, 1000, 760, 590, 910, 650, 810, 500),
    type = factor(rep(c("A", "B"), each = 10))
)
tool$x2 <- as.numeric(tool$type == "B")   # indicator: 1 if type B, 0 if type A

# ---- Models ----------------------------------------------------------------
m_slr <- lm(y ~ x1, data = tool)           # ignores tool type
m_par <- lm(y ~ x1 + x2, data = tool)      # parallel lines
m_int <- lm(y ~ x1 + x2 + x1:x2, data = tool)  # different slopes (same as y ~ x1 * x2)

summary(m_slr)
summary(m_par)
summary(m_int)

# Using the factor directly gives the same fit (R builds the indicator, reference = first level)
summary(lm(y ~ x1 + type, data = tool))
model.matrix(~ x1 + type, data = tool)[c(1, 11), ]

# ---- Tests (reduced vs full) ----------------------------------------------
anova(m_slr, m_par)   # H0: beta2 = 0 (F = t^2)
anova(m_par, m_int)   # H0: beta3 = 0 (parallel lines)
anova(m_slr, m_int)   # H0: beta2 = beta3 = 0 (one common line)

# ---- Separate regressions: same coefficients as the interaction model ------
coef(lm(y ~ x1, data = tool, subset = type == "A"))
coef(lm(y ~ x1, data = tool, subset = type == "B"))
c(coef(m_int)["(Intercept)"] + coef(m_int)["x2"], coef(m_int)["x1"] + coef(m_int)["x1:x2"])

# ---- Plots (saved as PDFs for the slides) ----------------------------------
dir.create("figures", showWarnings = FALSE)
cols <- c(A = "steelblue4", B = "darkorange2")
grid_x <- seq(500, 1000, length.out = 100)

plot_tool <- function(file, model = NULL, ignore_type = FALSE) {
    pdf(file, width = 4, height = 3.2, pointsize = 11)
    par(mar = c(4, 4, 0.5, 0.5), las = 1)
    plot(y ~ x1, data = tool, col = cols[type], pch = 19,
         xlab = "Lathe speed (rpm)", ylab = "Tool life (hours)",
         xlim = c(500, 1000), ylim = c(10, 46))
    if (ignore_type) abline(m_slr, lty = 2, lwd = 2, col = "gray40")
    if (!is.null(model)) {
        for (g in c("A", "B")) {
            lines(grid_x, predict(model, data.frame(x1 = grid_x, x2 = as.numeric(g == "B"))),
                  col = cols[g], lwd = 2)
        }
    }
    legend("topright", legend = c("Type A", "Type B"), col = cols, pch = 19, bty = "n")
    dev.off()
}

plot_tool("tool-slr.pdf", ignore_type = TRUE)  # simple regression ignoring type
plot_tool("tool-parallel.pdf", model = m_par)  # parallel lines
plot_tool("tool-interaction.pdf", model = m_int)  # different slopes

# ---- Factor with k = 3 levels (simulated, for illustration) ----------------
set.seed(1)
d3 <- data.frame(x1 = runif(60, 500, 1000), g = factor(rep(c("A", "B", "C"), each = 20)))
d3$y <- 35 - 0.025 * d3$x1 + c(A = 0, B = 12, C = 6)[as.character(d3$g)] + rnorm(60, sd = 3)

m3 <- lm(y ~ x1 + g, data = d3)
summary(m3)

# Dummy variable trap: all three indicators + intercept -> one coefficient is NA (aliased)
d3$dA <- as.numeric(d3$g == "A"); d3$dB <- as.numeric(d3$g == "B"); d3$dC <- as.numeric(d3$g == "C")
lm(y ~ x1 + dA + dB + dC, data = d3)
qr(model.matrix(~ x1 + dA + dB + dC, data = d3))$rank   # 4, not 5

# Changing the reference level: coefficients change, fit does not
m3_C <- lm(y ~ x1 + relevel(g, ref = "C"), data = d3)
coef(m3_C)
all.equal(fitted(m3), fitted(m3_C))
c(SSRes_ref_A = deviance(m3), SSRes_ref_C = deviance(m3_C))

# Partial F-test for the whole factor, and for equal slopes across levels
anova(lm(y ~ x1, data = d3), m3)
anova(m3, lm(y ~ x1 * g, data = d3))

# ---- Two quantitative regressors: centering reduces collinearity with the product ----
set.seed(2)
q <- data.frame(x1 = rnorm(100, 50, 10), x2 = rnorm(100, 20, 5))
q$y <- 5 + 0.4 * q$x1 + 0.8 * q$x2 + 0.05 * q$x1 * q$x2 + rnorm(100, sd = 3)
cor(q$x1, q$x1 * q$x2)
q$c1 <- q$x1 - mean(q$x1); q$c2 <- q$x2 - mean(q$x2)
cor(q$c1, q$c1 * q$c2)
summary(lm(y ~ x1 * x2, data = q))
summary(lm(y ~ c1 * c2, data = q))
