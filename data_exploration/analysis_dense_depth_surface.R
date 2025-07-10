# Analysis of Densely Sampled Depths and Surface Points in XY Grid
# Author: [Your Name]
# Date: 2025-04-30
.libPaths("C:/Users/23478671/R_libs")

# Load libraries
library(dplyr)
library(ggplot2)
library(viridis)

# Path to data file
# (edit as needed)
data_file <- "C:/Users/23478671/Github/Bayesian-Mixtures-Domaining/data/raw/cluster_0_data.csv"

# Read and prepare the data
data <- read.csv(data_file, stringsAsFactors = FALSE) %>%
  rename(
    name   = Name,
    x      = X,
    y      = Y,
    z      = Z,
    target = Cu_pct
  )

# --- 1. Analyze the most densely populated z's (depths) ---

# Define depth bins (adjust bin width as needed)
bin_width <- 10
data$z_bin <- cut(data$z, breaks = seq(floor(min(data$z)), ceiling(max(data$z)), by = bin_width), include.lowest = TRUE)

# Count samples per bin
z_density <- data %>%
  group_by(z_bin) %>%
  summarize(n = n(), mean_z = mean(z, na.rm = TRUE)) %>%
  arrange(desc(n))

# Plot sample density by depth bin
p_z_density <- ggplot(z_density, aes(x = mean_z, y = n)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Sample Density by Depth Bin", x = "Mean Depth (z)", y = "Sample Count") +
  theme_minimal()
print(p_z_density)

# Optionally, get the top N most densely sampled bins
top_z_bins <- head(z_density, 5)
print(top_z_bins)

# --- 2. Analyze the most densely sampled surface points in the XY grid ---

# For each hole, get the surface point (highest z)
surface_points <- data %>%
  group_by(name) %>%
  filter(z == max(z, na.rm = TRUE)) %>%
  ungroup()

# Bin XY surface locations by a specified radius (bin size)
xy_bin_size <- 40  # Adjust as needed (same units as x/y, e.g. meters)

surface_points <- surface_points %>%
  mutate(
    x_bin = floor(x / xy_bin_size) * xy_bin_size,
    y_bin = floor(y / xy_bin_size) * xy_bin_size
  )

xy_density <- surface_points %>%
  group_by(x_bin, y_bin) %>%
  summarize(n = n(), .groups = "drop") %>%
  arrange(desc(n))

# Plot density of binned surface points on XY grid
p_xy_density <- ggplot(xy_density, aes(x = x_bin, y = y_bin, fill = n)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Surface Sample Density (Binned XY Grid)", x = "X bin", y = "Y bin", fill = "Count") +
  theme_minimal()
print(p_xy_density)

# Optionally, print the densest bins
top_xy <- head(xy_density, 5)
print(top_xy)
