# Penalty per Parameter of AIC, HQIC and BIC (Lecture 12)
# Primary Reference: Montgomery, Peck, and Vining (5th/6th ed.), Chapter 10
# No data: the three criteria share the lack-of-fit term n ln(SSRes/n) and differ only in
# the penalty charged per parameter:
#   AIC  : 2            (constant)
#   HQIC : 2 ln(ln(n))  (exceeds AIC only when n > e^e  = 15.15)
#   BIC  : ln(n)        (exceeds AIC only when n > e^2  =  7.39)
#
# Generates:
#   1. penalty-comparison.pdf - penalty per parameter vs n (log x-axis), Hald n = 13 marked

library(ggplot2)


pt_col      <- "#1b9e77"    # AIC (teal)
pt_flag_col <- "#d95f02"    # BIC (orange-red)
pt_lev_col  <- "#7570b3"    # HQIC (purple)
ref_col     <- "#386cb0"    # Reference line colour (blue)

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

# Output directory (run from repo root, lecture folder, or figures/)
fig_dir <- if (file.exists("cement.csv")) {
  "."
} else if (dir.exists("12_model-selection/figures")) {
  "12_model-selection/figures"
} else if (dir.exists("figures")) {
  "figures"
} else {
  "."
}
out_file <- function(filename) file.path(fig_dir, filename)

# ---- Penalties ----------------------------------------------------------------------
n_grid <- exp(seq(log(5), log(1e4), length.out = 400))
pen <- rbind(
  data.frame(n = n_grid, penalty = 2,                 crit = "AIC"),
  data.frame(n = n_grid, penalty = 2 * log(log(n_grid)), crit = "HQIC"),
  data.frame(n = n_grid, penalty = log(n_grid),       crit = "BIC")
)
pen$crit <- factor(pen$crit, levels = c("AIC", "HQIC", "BIC"))

n_hald <- 13
c(AIC = 2, HQIC = 2 * log(log(n_hald)), BIC = log(n_hald))   # 2.00 1.88 2.56: HQIC < AIC at n = 13
c(BIC_cross = exp(2), HQIC_cross = exp(exp(1)))              # 7.39 15.15
sapply(c(50, 100, 1000), function(n) c(HQIC = 2 * log(log(n)), BIC = log(n)))
#        n=50  n=100 n=1000
# HQIC   2.73  3.05   3.86
# BIC    3.91  4.61   6.91

crit_cols <- c(AIC = pt_col, HQIC = pt_lev_col, BIC = pt_flag_col)

# ---- FIGURE: penalty-comparison.pdf -------------------------------------------------
p_pen <- ggplot(pen, aes(x = n, y = penalty, colour = crit)) +
  geom_vline(xintercept = n_hald, colour = "grey45", linetype = "dashed", linewidth = 0.5) +
  annotate("text", x = n_hald * 1.12, y = 8.6, hjust = 0, size = 2.9, colour = "grey30",
           label = "Hald cement, n = 13") +
  geom_vline(xintercept = c(exp(2), exp(exp(1))), colour = ref_col,
             linetype = "dotted", linewidth = 0.5) +
  annotate("text", x = exp(2) / 1.08, y = 6.6, hjust = 1, size = 2.7, colour = ref_col,
           label = "italic(n) == italic(e)^2", parse = TRUE) +
  annotate("text", x = exp(exp(1)) * 1.08, y = 6.6, hjust = 0, size = 2.7, colour = ref_col,
           label = "italic(n) == italic(e)^italic(e)", parse = TRUE) +
  geom_line(linewidth = 0.9) +
  geom_point(data = data.frame(n = n_hald, penalty = c(2, 2 * log(log(n_hald)), log(n_hald)),
                               crit = factor(c("AIC", "HQIC", "BIC"), levels = levels(pen$crit))),
             size = 2.4) +
  annotate("text", x = 1e4, y = c(2, 2 * log(log(1e4)), log(1e4)) + 0.35, hjust = 1, size = 3.0,
           colour = crit_cols, fontface = "bold",
           label = c("AIC: 2", "HQIC: 2 ln ln n", "BIC: ln n")) +
  scale_colour_manual(values = crit_cols, guide = "none") +
  scale_x_log10(name = expression("Sample size" ~ italic(n) ~ "(log scale)"),
                breaks = c(5, 10, 20, 50, 100, 1000, 10000),
                labels = c("5", "10", "20", "50", "100", "1000", "10000")) +
  scale_y_continuous(name = "Penalty per parameter", limits = c(0, 9.6),
                     breaks = seq(0, 9, by = 1)) +
  theme_lecture()

ggsave(out_file("penalty-comparison.pdf"), p_pen, bg = "transparent",
       device = cairo_pdf, width = 4.8, height = 3.4)

cat("Figure written to", normalizePath(out_file("penalty-comparison.pdf")), "\n")
