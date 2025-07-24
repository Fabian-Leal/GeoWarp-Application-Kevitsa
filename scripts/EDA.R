.libPaths("C:/Users/23478671/R_libs")
library(dplyr)
library(ggplot2)
library(patchwork)

# ─────────────────────────────────────────────────────────────────────────────
# 1. Read & prepar
# ─────────────────────────────────────────────────────────────────────────────
data_file <- "C:/Users/23478671/Github/Bayesian-Mixtures-Domaining/data/raw/cluster_0_data.csv"

raw_full <- read.csv(data_file, stringsAsFactors = FALSE)

eda <- raw_full %>%
  select(x = X, y = Y, z = Z,
         Ni_pct, Cu_pct, Pt_ppb, Au_ppb) %>%
  filter(!is.na(x), !is.na(y), !is.na(z),
         !is.na(Ni_pct), !is.na(Cu_pct),
         !is.na(Pt_ppb), !is.na(Au_ppb))

# ─────────────────────────────────────────────────────────────────────────────
# 2. Summary stats
# ─────────────────────────────────────────────────────────────────────────────
stats <- eda %>%
  summarize(
    Ni_min    = min(Ni_pct),    Ni_1Q = quantile(Ni_pct, 0.25), Ni_med = median(Ni_pct),
    Ni_mean   = mean(Ni_pct),   Ni_3Q = quantile(Ni_pct, 0.75), Ni_max = max(Ni_pct),
    Cu_min    = min(Cu_pct),    Cu_1Q = quantile(Cu_pct, 0.25), Cu_med = median(Cu_pct),
    Cu_mean   = mean(Cu_pct),   Cu_3Q = quantile(Cu_pct, 0.75), Cu_max = max(Cu_pct),
    Pt_min    = min(Pt_ppb),    Pt_1Q = quantile(Pt_ppb, 0.25), Pt_med = median(Pt_ppb),
    Pt_mean   = mean(Pt_ppb),   Pt_3Q = quantile(Pt_ppb, 0.75), Pt_max = max(Pt_ppb),
    Au_min    = min(Au_ppb),    Au_1Q = quantile(Au_ppb, 0.25), Au_med = median(Au_ppb),
    Au_mean   = mean(Au_ppb),   Au_3Q = quantile(Au_ppb, 0.75), Au_max = max(Au_ppb)
  )

print(stats)


# ─────────────────────────────────────────────────────────────────────────────
# 3. Histograms densities
# ─────────────────────────────────────────────────────────────────────────────
p1 <- ggplot(eda, aes(x = Ni_pct)) +
  geom_histogram(bins = 80, fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(title = "Histogram of Ni_pct", x = "Ni_pct", y = "Frequency")

p2 <- ggplot(eda, aes(x = Cu_pct)) +
  geom_histogram(bins = 80, fill = "tomato", color = "black") +
  theme_minimal() +
  labs(title = "Histogram of Cu_pct", x = "Cu_pct", y = "Frequency")

p3 <- ggplot(eda, aes(x = Pt_ppb)) +
  geom_histogram(bins = 80, fill = "forestgreen", color = "black") +
  theme_minimal() +
  labs(title = "Histogram of Pt_ppb", x = "Pt_ppb", y = "Frequency")

p4 <- ggplot(eda, aes(x = Au_ppb)) +
  geom_histogram(bins = 80, fill = "darkorange", color = "black") +
  theme_minimal() +
  labs(title = "Histogram of Au_ppb", x = "Au_ppb", y = "Frequency")

(p1 | p2) / (p3 | p4)

# ─────────────────────────────────────────────────────────────────────────────
# 4. Scatter vs Depth (z) binned
# ─────────────────────────────────────────────────────────────────────────────

eda_binned <- eda %>%
  mutate(
    depth    = -z,
    depth_bin = floor(depth / 50) * 50
  ) %>%
  group_by(depth_bin) %>%
  summarise(
    Ni_avg  = mean(Ni_pct,  na.rm = TRUE),
    Cu_avg  = mean(Cu_pct,  na.rm = TRUE),
    Pt_avg  = mean(Pt_ppb,  na.rm = TRUE),
    Au_avg  = mean(Au_ppb,  na.rm = TRUE)
  ) %>%
  ungroup()


eda_binned <- eda %>%
  mutate(
    depth     = -z,
    depth_bin = floor(depth / 50) * 50
  ) %>%
  group_by(depth_bin) %>%
  summarise(
    count   = n(),
    Ni_avg  = mean(Ni_pct,  na.rm = TRUE),
    Ni_sd   = sd(Ni_pct,    na.rm = TRUE),
    Cu_avg  = mean(Cu_pct,  na.rm = TRUE),
    Cu_sd   = sd(Cu_pct,    na.rm = TRUE),
    Pt_avg  = mean(Pt_ppb,  na.rm = TRUE),
    Pt_sd   = sd(Pt_ppb,    na.rm = TRUE),
    Au_avg  = mean(Au_ppb,  na.rm = TRUE),
    Au_sd   = sd(Au_ppb,    na.rm = TRUE)
  ) %>%
  ungroup()

# Function to build a profile plot with:
#    • Horizontal “count bars” at each depth_bin (proportional to count)
#    • ±1 sd error bars and mean points
#    • Labels showing the actual sample count next to each bar
#    • Inverted y-axis (0 at top, more negative depths downward)
make_profile_plot <- function(avg_col, sd_col, colour, title_text, x_label) {
  df_local   <- eda_binned
  max_count  <- max(df_local$count, na.rm = TRUE)
  x_vals     <- df_local[[avg_col]]
  sd_vals    <- df_local[[sd_col]]
  
  # Determine data-range for mean ± sd
  x_min      <- min(x_vals - sd_vals, na.rm = TRUE)
  x_max      <- max(x_vals + sd_vals, na.rm = TRUE)
  data_range <- x_max - x_min
  
  # Use 20% of that data-range as the “bar‐region” width
  bar_region_width <- data_range * 0.20
  bar_start        <- x_min - bar_region_width
  
  # Compute bar_end ∶ proportional to count, and label_x (a bit to the right of bar_end)
  df_local <- df_local %>%
    mutate(
      bar_start = bar_start,
      bar_end   = bar_start + (count / max_count) * bar_region_width,
      label_x   = bar_end + (0.02 * data_range)  # offset 2% of data range to the right
    )
  
  # Determine top (shallowest, typically 0) and bottom (most negative) bins
  shallowest_bin <- max(df_local$depth_bin, na.rm = TRUE)
  deepest_bin    <- min(df_local$depth_bin, na.rm = TRUE)
  
  ggplot(df_local, aes(x = .data[[avg_col]], y = depth_bin)) +
    # 2a. “Count bars” (grey thick segments) behind the points
    geom_segment(
      aes(
        x    = bar_start,
        xend = bar_end,
        y    = depth_bin,
        yend = depth_bin
      ),
      color = "grey80",
      size  = 3
    ) +
    # 2b. ±1 sd error bars (horizontal)
    geom_errorbarh(
      aes(
        xmin = .data[[avg_col]] - .data[[sd_col]],
        xmax = .data[[avg_col]] + .data[[sd_col]]
      ),
      height = 20,    # vertical thickness of the error bar (m)
      color  = colour,
      alpha  = 0.5
    ) +
    # 2c. Mean point at each bin
    geom_point(size = 2, color = colour) +
    # 2d. Sample count labels just to the right of each bar
    geom_text(
      aes(x = label_x, label = count),
      color    = "grey20",
      size     = 3,
      hjust    = 0  # left-justified so text starts at label_x
    ) +
    # 2e. Invert y-axis: 0 at top → most negative at bottom
    scale_y_continuous(
      trans  = "reverse",
      limits = c(shallowest_bin, deepest_bin),
      breaks = seq(shallowest_bin, deepest_bin, by = -50),
      name   = "Depth (m)"
    ) +
    labs(
      title = title_text,
      x     = x_label
    ) +
    theme_minimal()
}

# 3. Build each of the four profiles
q1 <- make_profile_plot(
  avg_col    = "Ni_avg",
  sd_col     = "Ni_sd",
  colour     = "steelblue",
  title_text = "Ni_pct (bin mean ± 1 sd) vs. Depth\nwith Sample‐Count Bars",
  x_label    = "Ni_pct (bin average)"
)

q2 <- make_profile_plot(
  avg_col    = "Cu_avg",
  sd_col     = "Cu_sd",
  colour     = "tomato",
  title_text = "Cu_pct (bin mean ± 1 sd) vs. Depth\nwith Sample‐Count Bars",
  x_label    = "Cu_pct (bin average)"
)

q3 <- make_profile_plot(
  avg_col    = "Pt_avg",
  sd_col     = "Pt_sd",
  colour     = "forestgreen",
  title_text = "Pt_ppb (bin mean ± 1 sd) vs. Depth\nwith Sample‐Count Bars",
  x_label    = "Pt_ppb (bin average)"
)

q4 <- make_profile_plot(
  avg_col    = "Au_avg",
  sd_col     = "Au_sd",
  colour     = "darkorange",
  title_text = "Au_ppb (bin mean ± 1 sd) vs. Depth\nwith Sample‐Count Bars",
  x_label    = "Au_ppb (bin average)"
)

# 4. Arrange the four panels in a 2×2 grid
(q1 | q2) / (q3 | q4)






# ─────────────────────────────────────────────────────────────────────────────
# 5. Spatial scatter (x vs y colored by value)
# ─────────────────────────────────────────────────────────────────────────────
#  Compute 95th‐percentile cutoffs for each variable
Ni_thresh <- quantile(eda$Ni_pct,  0.95, na.rm = TRUE)
Cu_thresh <- quantile(eda$Cu_pct,  0.95, na.rm = TRUE)
Pt_thresh <- quantile(eda$Pt_ppb,  0.95, na.rm = TRUE)
Au_thresh <- quantile(eda$Au_ppb,  0.95, na.rm = TRUE)

r1 <- eda %>%
  filter(Ni_pct <= Ni_thresh) %>%
  ggplot(aes(x = x, y = y, color = Ni_pct)) +
  geom_point(alpha = 0.7, size = 1) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  labs(
    title = paste0("Spatial Ni_pct (≤ 95th%≈", round(Ni_thresh, 2), ")"),
    x     = "X",
    y     = "Y",
    color = "Ni_pct"
  )

r2 <- eda %>%
  filter(Cu_pct <= Cu_thresh) %>%
  ggplot(aes(x = x, y = y, color = Cu_pct)) +
  geom_point(alpha = 0.7, size = 1) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  labs(
    title = paste0("Spatial Cu_pct (≤ 95th%≈", round(Cu_thresh, 2), ")"),
    x     = "X",
    y     = "Y",
    color = "Cu_pct"
  )

r3 <- eda %>%
  filter(Pt_ppb <= Pt_thresh) %>%
  ggplot(aes(x = x, y = y, color = Pt_ppb)) +
  geom_point(alpha = 0.7, size = 1) +
  scale_color_viridis_c(option = "cividis") +
  theme_minimal() +
  labs(
    title = paste0("Spatial Pt_ppb (≤ 95th%≈", round(Pt_thresh, 0), ")"),
    x     = "X",
    y     = "Y",
    color = "Pt_ppb"
  )

r4 <- eda %>%
  filter(Au_ppb <= Au_thresh) %>%
  ggplot(aes(x = x, y = y, color = Au_ppb)) +
  geom_point(alpha = 0.7, size = 1) +
  scale_color_viridis_c(option = "cividis") +
  theme_minimal() +
  labs(
    title = paste0("Spatial Au_ppb (≤ 95th%≈", round(Au_thresh, 0), ")"),
    x     = "X",
    y     = "Y",
    color = "Au_ppb"
  )

# 3. 2×2 grid 
(r1 | r2) / (r3 | r4)










################################################################################################################################################################################################

q1 <- ggplot(eda, aes(x = z, y = Ni_pct)) +
  geom_point(alpha = 0.5, size = 1, color = "steelblue") +
  theme_minimal() +
  labs(title = "Ni_pct vs Depth (z)", x = "Depth (z)", y = "Ni_pct")

q2 <- ggplot(eda, aes(x = z, y = Cu_pct)) +
  geom_point(alpha = 0.5, size = 1, color = "tomato") +
  theme_minimal() +
  labs(title = "Cu_pct vs Depth (z)", x = "Depth (z)", y = "Cu_pct")

q3 <- ggplot(eda, aes(x = z, y = Pt_ppb)) +
  geom_point(alpha = 0.5, size = 1, color = "forestgreen") +
  theme_minimal() +
  labs(title = "Pt_ppb vs Depth (z)", x = "Depth (z)", y = "Pt_ppb")

q4 <- ggplot(eda, aes(x = z, y = Au_ppb)) +
  geom_point(alpha = 0.5, size = 1, color = "darkorange") +
  theme_minimal() +
  labs(title = "Au_ppb vs Depth (z)", x = "Depth (z)", y = "Au_ppb")

(q1 | q2) / (q3 | q4)








# Create a binned summary (50 m bins) with mean ± 1 sd for each variable
eda_binned <- eda %>%
  mutate(
    depth     = -z,
    depth_bin = floor(depth / 50) * 50
  ) %>%
  group_by(depth_bin) %>%
  summarise(
    Ni_avg  = mean(Ni_pct, na.rm = TRUE),
    Ni_sd   = sd(Ni_pct,   na.rm = TRUE),
    Cu_avg  = mean(Cu_pct, na.rm = TRUE),
    Cu_sd   = sd(Cu_pct,   na.rm = TRUE),
    Pt_avg  = mean(Pt_ppb, na.rm = TRUE),
    Pt_sd   = sd(Pt_ppb,   na.rm = TRUE),
    Au_avg  = mean(Au_ppb, na.rm = TRUE),
    Au_sd   = sd(Au_ppb,   na.rm = TRUE)
  ) %>%
  ungroup()

# Determine the shallowest bin (should be 0 if there are samples at z = 0)
shallowest_bin <- max(eda_binned$depth_bin, na.rm = TRUE)    # typically 0
# Determine the deepest (most negative) bin
deepest_bin   <- min(eda_binned$depth_bin, na.rm = TRUE)    # e.g. -250, -300, etc.

# Function to build a plot with ±1 sd horizontal error bars
make_profile_plot <- function(avg_col, sd_col, colour, title_text, x_label) {
  ggplot(eda_binned, aes(x = .data[[avg_col]], y = depth_bin)) +
    # ±1 sd bands (horizontal error bars)
    geom_errorbarh(
      aes(
        xmin = .data[[avg_col]] - .data[[sd_col]],
        xmax = .data[[avg_col]] + .data[[sd_col]]
      ),
      height = 20,         # vertical thickness (20 m) so it doesn't overlap adjacent bins
      color = colour,
      alpha = 0.5
    ) +
    # Mean point at each bin
    geom_point(size = 2, color = colour) +
    scale_y_continuous(
      trans  = "reverse",
      limits = c(shallowest_bin, deepest_bin),
      breaks = seq(shallowest_bin, deepest_bin, by = -50),
      name   = "Depth (m)"
    ) +
    labs(
      title = title_text,
      x     = x_label
    ) +
    theme_minimal()
}

# 2. Build each of the four profiles
q1 <- make_profile_plot(
  avg_col   = "Ni_avg",
  sd_col    = "Ni_sd",
  colour    = "steelblue",
  title_text = "Ni_pct (bin mean ± 1 sd) vs. Depth",
  x_label   = "Ni_pct (bin average)"
)

q2 <- make_profile_plot(
  avg_col   = "Cu_avg",
  sd_col    = "Cu_sd",
  colour    = "tomato",
  title_text = "Cu_pct (bin mean ± 1 sd) vs. Depth",
  x_label   = "Cu_pct (bin average)"
)

q3 <- make_profile_plot(
  avg_col   = "Pt_avg",
  sd_col    = "Pt_sd",
  colour    = "forestgreen",
  title_text = "Pt_ppb (bin mean ± 1 sd) vs. Depth",
  x_label   = "Pt_ppb (bin average)"
)

q4 <- make_profile_plot(
  avg_col   = "Au_avg",
  sd_col    = "Au_sd",
  colour    = "darkorange",
  title_text = "Au_ppb (bin mean ± 1 sd) vs. Depth",
  x_label   = "Au_ppb (bin average)"
)

# 3. Arrange the four plots in a 2×2 grid
(q1 | q2) / (q3 | q4)




