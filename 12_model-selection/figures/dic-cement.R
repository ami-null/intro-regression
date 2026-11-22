# Hald Cement Data: Deviance Information Criterion (Lecture 12)
# Primary Reference: Spiegelhalter, Best, Carlin and van der Linde (2002), Bayesian measures of
#                    model complexity and fit, JRSS-B 64(4), 583-639.
# Data: cement.csv (n = 13, K = 4)
#
# Bayesian set-up: reference prior pi(beta, sigma^2) proportional to 1/sigma^2. The posterior is
# known exactly, so no MCMC is needed:
#   sigma^2 | y        ~ SSRes / chi^2_{n-p}
#   beta | sigma^2, y  ~ N(beta_hat, sigma^2 (X'X)^{-1})
# Deviance D(theta) = -2 ln L(beta, sigma^2) (full normal log-likelihood, constants kept).
#   D_bar = posterior mean of D(theta)        (fit)
#   p_D   = D_bar - D(theta_bar)              (effective number of parameters)
#   DIC   = D_bar + p_D
# theta_bar = posterior means of beta and sigma^2 (plug-in convention of Spiegelhalter et al.).
# DIC is compared with R's AIC() (same deviance scale: -2 ln L_hat + 2(p + 1)).
#
# Generates:
#   1. dic-vs-aic-cement.pdf          - DIC against AIC for all 15 non-empty subsets
#   2. deviance-posterior-cement.pdf  - posterior draws of D(theta) for {x1, x2}, showing p_D
# Writes:
#   cement-dic.csv                    - DIC table for LaTeX

library(ggplot2)


pt_col      <- "#1b9e77"    # Teal
pt_flag_col <- "#d95f02"    # Orange-red
pt_lev_col  <- "#7570b3"    # Purple
ref_col     <- "#386cb0"    # Blue

theme_lecture <- function() {
  theme_minimal(base_size = 11) +
    theme(
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background  = element_rect(fill = "transparent", colour = NA),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.box.background = element_rect(fill = "transparent", colour = NA),
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = 10, face = "bold", colour = "grey15"),
      strip.background = element_rect(fill = "#f0f0f0", colour = NA),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9)
    )
}

# Determine input/output directories (run from repo root, lecture folder, or figures/)
csv_path <- if (file.exists("cement.csv")) {
  "cement.csv"
} else if (file.exists("12_model-selection/figures/cement.csv")) {
  "12_model-selection/figures/cement.csv"
} else if (file.exists("figures/cement.csv")) {
  "figures/cement.csv"
} else {
  "cement.csv"
}
fig_dir <- dirname(csv_path)
out_file <- function(filename) file.path(fig_dir, filename)

lab_expr <- function(model) {
  vapply(model, function(m) {
    idx <- sub("x", "", strsplit(m, ",")[[1]])
    paste0("list(", paste0("italic(x)[", idx, "]", collapse = ", "), ")")
  }, character(1), USE.NAMES = FALSE)
}

# ---- 1. Data ----------------------------------------------------------------------
cement <- read.csv(csv_path)
n  <- nrow(cement)
xv <- paste0("x", 1:4)
subsets <- unlist(lapply(1:4, function(k) combn(xv, k, simplify = FALSE)), recursive = FALSE)

# ---- 2. Exact posterior simulation and DIC ---------------------------------------------
deviance_fn <- function(y, mu, s2) -2 * sum(dnorm(y, mu, sqrt(s2), log = TRUE))

posterior_dic <- function(vars, S = 20000) {
  m  <- lm(reformulate(vars, "y"), data = cement)
  X  <- model.matrix(m)
  y  <- cement$y
  p  <- ncol(X)
  b_hat <- coef(m)
  SS <- sum(resid(m)^2)
  L  <- t(chol(solve(crossprod(X))))                 # (X'X)^{-1} = L L'

  s2_draw <- SS / rchisq(S, df = n - p)
  z <- matrix(rnorm(p * S), p, S)
  b_draw <- b_hat + (L %*% z) * rep(sqrt(s2_draw), each = p)   # p x S
  mu_draw <- X %*% b_draw                                       # n x S
  D <- vapply(seq_len(S), function(s) deviance_fn(y, mu_draw[, s], s2_draw[s]), numeric(1))

  D_bar <- mean(D)
  D_hat <- deviance_fn(y, X %*% rowMeans(b_draw), mean(s2_draw))
  list(row = data.frame(model = paste(vars, collapse = ","), p = p, k = p + 1,
                        D_bar = D_bar, D_hat = D_hat, p_D = D_bar - D_hat,
                        DIC = 2 * D_bar - D_hat, AIC = AIC(m)),
       D = D)
}

set.seed(2026)
res <- lapply(subsets, posterior_dic)
dic_tab <- do.call(rbind, lapply(res, `[[`, "row"))
dic_tab$diff <- dic_tab$DIC - dic_tab$AIC
print(dic_tab[order(dic_tab$DIC), ], digits = 4, row.names = FALSE)
# DIC tracks AIC within about 1 unit for every subset. p_D is a little below the nominal
# k = p + 1 because sigma^2 is plugged in at its posterior mean. The top three subsets
# ({x1,x2,x4}, {x1,x2,x3}, {x1,x2}) are within Monte Carlo error of each other.

csv_out <- dic_tab
csv_out[sapply(csv_out, is.numeric)] <- lapply(csv_out[sapply(csv_out, is.numeric)], round, 3)
write.csv(csv_out, out_file("cement-dic.csv"), row.names = FALSE)

# ---- FIGURE 1: DIC vs AIC (`dic-vs-aic-cement.pdf`) -------------------------------------
top <- c("x1,x2", "x1,x2,x4", "x1,x2,x3", "x1,x3,x4", "x1,x2,x3,x4", "x1,x4")
dic_tab$top <- dic_tab$model %in% top

p_main <- ggplot(dic_tab, aes(x = AIC, y = DIC)) +
  geom_abline(intercept = 0, slope = 1, colour = ref_col, linetype = "dashed", linewidth = 0.6) +
  geom_point(aes(colour = top), size = 2.2) +
  scale_colour_manual(values = c("FALSE" = "grey55", "TRUE" = pt_flag_col), guide = "none") +
  annotate("text", x = 94, y = 101, label = "DIC == AIC", parse = TRUE, hjust = 1,
           size = 2.9, colour = ref_col) +
  scale_x_continuous(name = "AIC (R form)", limits = c(60, 112)) +
  scale_y_continuous(name = "DIC (flat prior, 20 000 draws)", limits = c(60, 112)) +
  coord_equal() +
  theme_lecture()

# Right panel: the six best subsets as a dot plot (the top three nearly coincide in a scatter)
zoom <- dic_tab[dic_tab$top, ]
zoom_long <- rbind(data.frame(model = zoom$model, crit = "AIC", val = zoom$AIC),
                   data.frame(model = zoom$model, crit = "DIC", val = zoom$DIC))
zoom_long$model <- factor(zoom_long$model, levels = zoom$model[order(zoom$DIC, decreasing = TRUE)])

p_zoom <- ggplot(zoom_long, aes(x = val, y = model, colour = crit, shape = crit)) +
  geom_line(aes(group = model), colour = "grey70", linewidth = 0.5) +
  geom_point(size = 2.4) +
  scale_colour_manual(values = c(AIC = pt_col, DIC = pt_flag_col), name = NULL) +
  scale_shape_manual(values = c(AIC = 16, DIC = 17), name = NULL) +
  scale_x_continuous(name = "Criterion value") +
  scale_y_discrete(name = NULL, labels = function(l) parse(text = lab_expr(l))) +
  labs(title = "The six best subsets") +
  theme_lecture() +
  theme(plot.title = element_text(size = 9.5, hjust = 0.5),
        legend.position = "bottom", legend.margin = margin(t = -5),
        panel.grid.major.y = element_line(colour = "grey92"))

p_dic <- patchwork::wrap_plots(p_main, p_zoom, nrow = 1, widths = c(1, 1.1))

ggsave(out_file("dic-vs-aic-cement.pdf"), p_dic, bg = "transparent",
       device = cairo_pdf, width = 7.2, height = 3.6)

# ---- FIGURE 2: Posterior deviance for {x1, x2} (`deviance-posterior-cement.pdf`) ---------
i12 <- which(dic_tab$model == "x1,x2")
D12 <- res[[i12]]$D
r12 <- dic_tab[i12, ]
round(unlist(r12[c("D_hat", "D_bar", "p_D", "DIC")]), 2)

p_dev <- ggplot(data.frame(D = D12), aes(x = D)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.5, boundary = 0,
                 fill = pt_col, colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = r12$D_hat, colour = ref_col, linewidth = 0.8) +
  geom_vline(xintercept = r12$D_bar, colour = pt_flag_col, linewidth = 0.8, linetype = "dashed") +
  annotate("segment", x = r12$D_hat, xend = r12$D_bar, y = 0.235, yend = 0.235,
           arrow = arrow(ends = "both", length = unit(0.12, "cm")), colour = "grey20") +
  annotate("text", x = (r12$D_hat + r12$D_bar) / 2, y = 0.25, size = 2.9, colour = "grey20",
           label = sprintf("italic(p)[italic(D)] == %.2f", r12$p_D), parse = TRUE) +
  annotate("text", x = r12$D_hat - 0.3, y = 0.19, hjust = 1, size = 2.8, colour = ref_col,
           label = sprintf("italic(D)(bar(theta)) == %.2f", r12$D_hat), parse = TRUE) +
  annotate("text", x = r12$D_bar + 0.3, y = 0.19, hjust = 0, size = 2.8, colour = pt_flag_col,
           label = sprintf("bar(italic(D)) == %.2f", r12$D_bar), parse = TRUE) +
  scale_x_continuous(name = expression("Posterior deviance" ~ italic(D)(theta) ~
                                         "for" ~ list(italic(x)[1], italic(x)[2])),
                     limits = c(r12$D_hat - 4, r12$D_hat + 22)) +
  scale_y_continuous(name = "Density", expand = expansion(mult = c(0, 0.05))) +
  theme_lecture()

ggsave(out_file("deviance-posterior-cement.pdf"), p_dev, bg = "transparent",
       device = cairo_pdf, width = 4.8, height = 3.2)

cat("\nFigures written to", normalizePath(fig_dir), "\n")
