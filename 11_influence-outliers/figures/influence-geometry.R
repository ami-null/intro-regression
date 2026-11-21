# Influence and Outliers: Geometric Taxonomy (Lecture 11)
# Illustrates the 4-quadrant taxonomy of regression points:
#   Quadrant 1: Typical Inlier      (Low leverage, small residual)
#   Quadrant 2: y-Outlier           (Low leverage, large residual)
#   Quadrant 3: Good Leverage Point (High leverage, small residual)
#   Quadrant 4: Bad Leverage Point  (High leverage, large residual - influential)
#
# Generates: influence-geometry.pdf

library(ggplot2)

# ---- Plot configuration & theme ----------------------------------------------
pt_base_col   <- "#2b5c8f"    # Slate blue for baseline points
pt_target_col <- "#d95f02"    # Red-orange for target observation
line_full_col <- "#d95f02"    # Line with target point (solid)
line_base_col <- "#636363"    # Line without target point (dashed)

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

# Determine output directory
fig_dir <- if (file.exists("delivery.csv")) {
  "."
} else if (dir.exists("11_influence-outliers/figures")) {
  "11_influence-outliers/figures"
} else if (dir.exists("figures")) {
  "figures"
} else {
  "."
}

# ---- Generate synthetic baseline data ----------------------------------------
set.seed(42)
n_base <- 18
x_base <- seq(2.5, 7.5, length.out = n_base)
# True relationship: y = 2.0 + 1.5 * x + e, e ~ N(0, 0.4^2)
e_base <- rnorm(n_base, mean = 0, sd = 0.40)
y_base <- 2.0 + 1.5 * x_base + e_base

fit_base <- lm(y_base ~ x_base)
b0_base  <- coef(fit_base)[1]
b1_base  <- coef(fit_base)[2]

# Centroid of baseline
x_bar <- mean(x_base)   # 5.0
y_hat_center <- b0_base + b1_base * x_bar

# ---- Define the 4 target points ----------------------------------------------
# 1. Inlier: low leverage (at x = 5.0), small residual (on regression line)
# 2. y-Outlier: low leverage (at x = 5.0), large residual (+10 units above line, so the upward shift of the line is visible)
# 3. Good Leverage: high leverage (at x = 13.5), small residual (conforms to line)
# 4. Bad Leverage: high leverage (at x = 13.5), discordant y (pulled down to y = 6.8)

x_target_high <- 13.5
y_target_good <- b0_base + b1_base * x_target_high + 0.1   # ~21.5
y_target_bad  <- 6.8                                       # severe drop

scenarios <- list(
  q1 = list(
    panel_title = "(a) Typical Inlier\nLow Leverage, Small Residual",
    x_t = 5.0,
    y_t = y_hat_center + 0.15,
    desc = "Target conforms to local pattern\nLine unaffected"
  ),
  q2 = list(
    panel_title = "(b) Response Outlier (y-Outlier)\nLow Leverage, Large Residual",
    x_t = 5.0,
    y_t = y_hat_center + 10,
    desc = "Shifts intercept upward\nSlope virtually unchanged"
  ),
  q3 = list(
    panel_title = "(c) Good Leverage Point\nHigh Leverage, Small Residual",
    x_t = x_target_high,
    y_t = y_target_good,
    desc = "Conforms to linear trend\nStabilizes slope, precision increases"
  ),
  q4 = list(
    panel_title = "(d) Bad Leverage Point (Influential)\nHigh Leverage, Large Residual",
    x_t = x_target_high,
    y_t = y_target_bad,
    desc = "Extreme leverage + discordance\nPulls regression line toward itself"
  )
)

# ---- Build data frames for plotting ------------------------------------------
plot_points <- data.frame()
plot_lines  <- data.frame()
plot_annot  <- data.frame()

grid_x <- seq(1.5, 14.5, length.out = 150)


for (s_id in names(scenarios)) {
  s <- scenarios[[s_id]]
  
  # Data for this scenario (18 baseline + 1 target = 19 points)
  df_scen <- data.frame(
    x = c(x_base, s$x_t),
    y = c(y_base, s$y_t),
    type = c(rep("Baseline", n_base), "Target"),
    panel = s$panel_title
  )
  
  fit_all <- lm(y ~ x, data = df_scen)
  
  # Diagnostics for target point (row 19)
  h_all    <- hatvalues(fit_all)
  r_all    <- rstandard(fit_all)
  cook_all <- cooks.distance(fit_all)
  
  h_t    <- h_all[19]
  r_t    <- r_all[19]
  cook_t <- cook_all[19]
  b0_all <- coef(fit_all)[1]
  b1_all <- coef(fit_all)[2]
  
cat(sprintf("%s:\n", gsub("\n", " - ", s$panel_title)))
cat(sprintf("  Target Point : (x = %.2f, y = %.2f)\n", s$x_t, s$y_t))
cat(sprintf("  Diagnostics  : h_ii = %.3f (2p/n = %.3f), r_i = %+.2f, Cook's D = %.3f\n",
             h_t, 2 * 2 / 19, r_t, cook_t))
cat(sprintf("  Slope without: b1 = %.3f | Slope with: b1 = %.3f (delta = %+.3f)\n\n",
             b1_base, b1_all, b1_all - b1_base))
  
  # Line predictions across grid_x
  lines_scen <- rbind(
    data.frame(
      panel = s$panel_title,
      x = grid_x,
      y = b0_base + b1_base * grid_x,
      model = "Without Target (baseline)"
    ),
    data.frame(
      panel = s$panel_title,
      x = grid_x,
      y = b0_all + b1_all * grid_x,
      model = "With Target (all points)"
    )
  )
  
  # Annotation label showing key diagnostics inside each panel
  # (plotmath ignores "\n", so all three values go on one line)
  diag_label <- sprintf("h[ii] == %.3f*','~~r[i] == %+.2f*','~~D[i] == %.3f",
                        h_t, r_t, cook_t)
  annot_scen <- data.frame(
    panel = s$panel_title,
    x = 1.8,
    y = 23.5,
    label = diag_label
  )
  
  plot_points <- rbind(plot_points, df_scen)
  plot_lines  <- rbind(plot_lines, lines_scen)
  plot_annot  <- rbind(plot_annot, annot_scen)
}

# Preserve factor order of panels
panel_levels <- sapply(scenarios, function(s) s$panel_title)
plot_points$panel <- factor(plot_points$panel, levels = panel_levels)
plot_lines$panel  <- factor(plot_lines$panel, levels = panel_levels)
plot_annot$panel  <- factor(plot_annot$panel, levels = panel_levels)

# Separate points into baseline and target for clean layer styling
pts_base   <- subset(plot_points, type == "Baseline")
pts_target <- subset(plot_points, type == "Target")

# ---- Construct 4-panel ggplot ------------------------------------------------
p_geom <- ggplot() +
  # Reference line for x_bar centroid in regressor space
  geom_vline(xintercept = x_bar, linetype = "dotted", colour = "grey70", linewidth = 0.5) +
  # Fitted lines
  geom_line(
    data = plot_lines,
    aes(x = x, y = y, colour = model, linetype = model, linewidth = model)
  ) +
  # Baseline points
  geom_point(
    data = pts_base,
    aes(x = x, y = y),
    colour = pt_base_col,
    size = 1.9,
    alpha = 0.85
  ) +
  # Target point with highlight ring
  geom_point(
    data = pts_target,
    aes(x = x, y = y),
    colour = pt_target_col,
    fill = "#fdbf6f",
    shape = 21,
    size = 3.6,
    stroke = 1.4
  ) +
  # Diagnostic metrics annotation (filled box so the x-bar line does not cross the text)
  geom_label(
    data = plot_annot,
    aes(x = x, y = y, label = label),
    parse = TRUE,
    hjust = 0,
    vjust = 1,
    size = 2.9,
    colour = "grey25",
    fill = "white",
    linewidth = 0,
    label.padding = unit(0.12, "lines")
  ) +
  # Faceting: 2x2 grid corresponding to the 4 quadrants
  facet_wrap(~ panel, nrow = 2, ncol = 2) +
  scale_colour_manual(
    name = "Regression Model:",
    values = c(
      "Without Target (baseline)" = line_base_col,
      "With Target (all points)"  = line_full_col
    )
  ) +
  scale_linetype_manual(
    name = "Regression Model:",
    values = c(
      "Without Target (baseline)" = "dashed",
      "With Target (all points)"  = "solid"
    )
  ) +
  scale_linewidth_manual(
    name = "Regression Model:",
    values = c(
      "Without Target (baseline)" = 0.75,
      "With Target (all points)"  = 0.95
    )
  ) +
  scale_x_continuous(
    name = expression(Regressor ~ italic(x)),
    limits = c(1.0, 15.0),
    breaks = seq(2, 14, by = 3)
  ) +
  scale_y_continuous(
    name = expression(Response ~ italic(y)),
    limits = c(2.0, 25.5),
    breaks = seq(5, 25, by = 5)
  ) +
  theme_lecture() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 9.5, face = "bold"),
    legend.text = element_text(size = 9),
    legend.margin = margin(t = -3, b = 2)
  )

# ---- Save PDF ----------------------------------------------------------------
out_path <- file.path(fig_dir, "influence-geometry.pdf")
ggsave(
  filename = out_path,
  plot = p_geom,
  bg = "transparent",
  device = cairo_pdf,
  width = 6.8,
  height = 5.2
)
