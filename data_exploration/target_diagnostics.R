.libPaths("C:/Users/23478671/R_libs")
library(dplyr)
library(ggplot2)
library(patchwork)
library(geowarp)

data_file <- "C:/Users/23478671/Github/Bayesian-Mixtures-Domaining/data/raw/cluster_0_data.csv"

# Read and prepare the data
data <- read.csv(data_file, stringsAsFactors = FALSE) %>%
  rename(
    name   = Name,
    x      = X,
    y      = Y,
    z      = Z,
    target = Pd_ppb
  ) %>%
  mutate(depth = z) %>%
  select(-z)

# Ensure the "plots" directory exists
if (!dir.exists("plots")) dir.create("plots")

# ---- Automated diagnostics for candidate variables ----
library(viridis)
library(gridExtra)

# List of candidate variables (edit as needed)
candidate_vars <- c(
  "Au_ppb", "Pd_ppb", "Pt_ppb", "Ag_ppm", "Cu_pct", "Ni_pct", "Zn_ppm", "Pb_ppm", "Co_ppm", "Cr_ppm", "Mn_ppm",
  "S_pct", "Fe_pct", "Mg_pct", "Al_pct", "Si_pct", "Ca_pct", "Density_gcm3"
)

# Only keep candidate variables that exist in the data
data_vars <- intersect(candidate_vars, names(data))

# Function to create all diagnostic plots for a given variable
diagnostic_plots <- function(var, data, outdir = "plots") {
  
  # --- Function body ---
  
  message("Processing ", var)
  # Remove NAs
  df <- data %>% filter(!is.na(.data[[var]]), !is.na(depth), !is.na(name), .data[[var]] >= 0, depth >= -850)
  if (nrow(df) < 10) return(NULL)
  # For outlier trimming 
  lower <- quantile(df[[var]], 0.25, na.rm = TRUE) - 7 * IQR(df[[var]], na.rm = TRUE)
  upper <- quantile(df[[var]], 0.75, na.rm = TRUE) + 7 * IQR(df[[var]], na.rm = TRUE)
  df_trim <- df %>% filter(.data[[var]] >= lower, .data[[var]] <= upper)
  
  # Filter out holes with too few observations
  min_obs_per_hole <- 50
  hole_counts <- df_trim %>% group_by(name) %>% summarize(n = n())
  good_holes <- hole_counts %>% filter(n >= min_obs_per_hole) %>% pull(name)
  df_trim <- df_trim %>% filter(name %in% good_holes)
  
  # Get min/max for x-axis (target variable) and y-axis (depth) for use in all plots
  x_min <- min(df_trim[[var]], na.rm = TRUE)
  x_max <- max(df_trim[[var]], na.rm = TRUE)
  y_min <- min(df_trim$depth, na.rm = TRUE)
  y_max <- max(df_trim$depth, na.rm = TRUE)
  
  # 1. Mean ± 2*stdev across all holes, binned by depth
  num_bins <- 30
  df_trim$depth_bin <- cut(df_trim$depth, breaks = num_bins)
  summary_by_bin <- df_trim %>%
    group_by(depth_bin) %>%
    summarize(
      mid_depth = mean(depth, na.rm = TRUE),
      mean_val = mean(.data[[var]], na.rm = TRUE),
      sd_val = sd(.data[[var]], na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    arrange(mid_depth)
  p1 <- ggplot(summary_by_bin, aes(x = mean_val, y = mid_depth)) +
    geom_point(color = "blue", size = 1) +
    geom_ribbon(aes(xmin = mean_val - 2 * sd_val, xmax = mean_val + 2 * sd_val), fill = "blue", alpha = 0.2) +
    xlim(x_min, x_max) +
    labs(x = var, y = 'Depth', title = paste(var, 'Mean ± 2*Stdev by Depth Bin')) +
    theme_minimal()
  
  
  # 3. Example holes: 18 plots, each with 1 random hole
  set.seed(123)
  unique_holes <- unique(df_trim$name)
  n_holes <- length(unique_holes)
  n_plots <- min(18, n_holes)
  sample_holes <- sample(unique_holes, n_plots)
  example_plots <- lapply(seq_along(sample_holes), function(i) {
    hole <- sample_holes[i]
    df_sub <- df_trim %>% filter(name == hole)
    ggplot(df_sub, aes(x = .data[[var]], y = depth)) +
      geom_point(size = 0.7, alpha = 0.8, color = "#1b9e77") +
      labs(title = paste("Hole:", hole), x = var, y = "-Depth (z)") +
      xlim(x_min, x_max) +
      ylim(y_min, y_max) +
      theme_minimal(base_size = 8) +
      theme(plot.title = element_text(size = 9),
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 7),
            plot.margin = margin(2,2,2,2,"pt"))
  })
  
  # For arranging: 9 per page (3x3 grid)
  n_per_page <- 9
  n_pages <- ceiling(length(example_plots) / n_per_page)
  example_pages <- split(example_plots, ceiling(seq_along(example_plots) / n_per_page))
  
  
  # 3. 1D Empirical Variogram (spatial autocorrelation) plot -- memory safe
  df_vario <- df_trim %>% select(depth, value = !!sym(var))
  n <- nrow(df_vario)
  max_pairs <- 2000
  
  if (n < 2) {
    # Not enough data for variogram
    p3 <- ggplot() + labs(title = "Not enough data for variogram") + theme_minimal()
  } else {
    # Sample pairs for large datasets
    set.seed(42)
    pair_idx <- if (choose(n, 2) > max_pairs) {
      as.data.frame(t(replicate(max_pairs, sort(sample(1:n, 2)))))
    } else {
      which_mat <- which(upper.tri(matrix(1, n, n)), arr.ind = TRUE)
      data.frame(V1 = which_mat[, 1], V2 = which_mat[, 2])
    }
    lag <- abs(df_vario$depth[pair_idx$V1] - df_vario$depth[pair_idx$V2])
    gamma <- 0.5 * (df_vario$value[pair_idx$V1] - df_vario$value[pair_idx$V2])^2
    lag_bins <- 300
    max_lag <- diff(range(df_vario$depth))
    lag_width <- max_lag / lag_bins
    lag_bin <- cut(lag, breaks = seq(0, max_lag, by = lag_width), include.lowest = TRUE)
    vg <- data.frame(lag = lag, gamma = gamma, lag_bin = lag_bin) %>%
      group_by(lag_bin) %>%
      summarize(
        lag = mean(lag, na.rm = TRUE),
        gamma = mean(gamma, na.rm = TRUE),
        n = n(),
        .groups = "drop"
      ) %>%
      filter(n > 5)
    vg <- vg %>% filter(lag <= 50)
    p3 <- ggplot(vg, aes(x = lag, y = gamma)) +
      geom_point() +
      geom_line() +
      labs(
        x = "Lag (depth difference)",
        y = "Semivariance",
        title = paste("Empirical variogram for", var, "(<= 300m)")
      ) +
      theme_minimal()
  }
  

  # Output to PDF (only the three requested plot types)
  pdf(file.path(outdir, paste0("diagnostics_", var, ".pdf")), width = 8, height = 10)
  print(p1)
  print(p3)
  for (pg in example_pages) {
    print(patchwork::wrap_plots(pg, nrow = 3, ncol = 3))
  }
  dev.off()
}

# --- QUICK CHECK: Generate plots for just Density_gcm3 (remove or comment out for batch runs)
diagnostic_plots("Density_gcm3", data)

# Run diagnostics for all candidate variables
for (var in data_vars) {
  diagnostic_plots(var, data)
}







# Plot sample locations
png(filename = "plots/sample_locations.png", width = 800, height = 600)
ggplot(data, aes(x, y, color = target)) +
  geom_point() +
  scale_color_viridis_c() +
  labs(title = "Sample Locations", x = "X", y = "Y", color = "Density")
dev.off()




