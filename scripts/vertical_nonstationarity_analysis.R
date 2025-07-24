# main.R -------------------------------------------------------------------
setwd("C:/Users/23478671/Github/GeoWarp-Application-Kevitsa")
source("helper_functions/00_utils_io.R")
source("helper_functions/01_preprocess.R")
source("helper_functions/02_split.R")
source("helper_functions/03_models.R")
source("helper_functions/04_visuals.R")
set.seed(1)
library(purrr)      # map()
library(patchwork)  # wrap_plots()
library(ggplot2)
library(dplyr)
library(stats)

csv  <- "cluster_0_data.csv"
grade_vars <- c("Ag_ppm", "Au_ppb", "Cu_pct", "Ni_pct")


nice_label <- function(v) {
  switch(v,
         "Ag_ppm" = "log Ag (ppm)",
         "Au_ppb" = "log Au (ppb)",
         "Cu_pct" = "log Cu (%)",
         "Ni_pct" = "log Ni (%)",
         paste0("log ", v))          # fallback
}





#-------------------------------------------------------------------------------
# 1) Function: fit variogram by depth bin (|Δz| only)
#-------------------------------------------------------------------------------
fit_vario_by_bin <- function(df, bin_width = 50, cutoff = NULL, min_pts = 10) {
  # Create depth bins
  z_breaks <- seq(floor(min(df$z)), ceiling(max(df$z)), by = bin_width)
  
  df %>%
    mutate(z_bin = cut(z, breaks = z_breaks, include.lowest = TRUE)) %>%
    filter(!is.na(resid)) %>%
    group_by(z_bin) %>%
    filter(n() >= min_pts) %>%   # require enough samples
    do({
      bin_df <- as.data.frame(.)
      bin_df$dummy <- 0  # dummy coordinate for 1D variogram
      
      # compute vertical-only empirical variogram
      emp <- variogram(
        resid ~ 1,
        data      = bin_df,
        locations = ~ z + dummy,
        cutoff    = cutoff
      )
      # fit spherical model
      fit <- tryCatch(fit.variogram(emp, vgm("Sph")), error = function(e) NULL)
      
      # extract range and sill
      tibble(
        range = if (!is.null(fit)) fit$range[2] else NA_real_,
        sill  = if (!is.null(fit)) fit$psill[2] else NA_real_
      )
    }) %>%
    ungroup() %>%
    # compute mid-depth of each bin for plotting
    mutate(
      z_mid = (
        readr::parse_number(sub("\\((.+),.*", "\\1", z_bin)) +
          readr::parse_number(sub(".*,(.+)\\]", "\\1", z_bin))
      ) / 2
    )
}

#-------------------------------------------------------------------------------
# 2) Function: compute diagnostics for one metal
#-------------------------------------------------------------------------------
compute_vario_diagnostics <- function(csv_path, var, bin_width = 50, cutoff = NULL) {
  # Load and filter raw data
  df <- load_raw(csv_path, var = var, sample_frac = 1) %>%
    depth_filter(min_cnt = 800) %>%
    mutate(
      log_val = target,  # already log-scale
      x = x, y = y, z = z
    )
  
  # detrend & standardise within depth bins
  z_breaks <- seq(floor(min(df$z)), ceiling(max(df$z)), by = bin_width)
  df_std <- df %>%
    mutate(z_bin = cut(z, breaks = z_breaks, include.lowest = TRUE)) %>%
    group_by(z_bin) %>%
    mutate(
      mean_bin = mean(log_val, na.rm = TRUE),
      sd_bin   = sd(log_val,  na.rm = TRUE),
      resid    = (log_val - mean_bin) / sd_bin
    ) %>%
    ungroup()
  
  # fit variogram per bin
  fit_vario_by_bin(df_std, bin_width = bin_width, cutoff = cutoff)
}

#-------------------------------------------------------------------------------
# 3) Run diagnostics for each metal and plot
#-------------------------------------------------------------------------------
# Settings
set.seed(123)

bin_width <- 30

max_lag   <- 100
grade_vars <- c("Au_ppb", "Ag_ppm", "Cu_pct", "Ni_pct")

metal_cols <- c(
  Au_ppb = "#DAA520",   # gold
  Ag_ppm = "#C0C0C0",   # silver
  Cu_pct = "#B87333",   # copper
  Ni_pct = "#4F9D9D"    # nickel-green/teal
)

metal_labels <- c(
  Au_ppb = "log Au (ppb)",
  Ag_ppm = "log Ag (ppm)",
  Cu_pct = "log Cu (%)",
  Ni_pct = "log Ni (%)"
)

# Compute vario stats
vario_stats <- map_dfr(grade_vars, function(v) {
  compute_vario_diagnostics(
    csv_path = csv,
    var       = v,
    bin_width = bin_width,
    cutoff    = max_lag
  ) %>%
    mutate(variable = v)
})

# Reshape for plotting
vario_long <- vario_stats %>%
  pivot_longer(
    cols      = c(range, sill),
    names_to  = "metric",
    values_to = "value"
  )

# Facet labels
metric_labels <- c(
  range = "Variogram Range (m)",
  sill  = "Variogram Sill"
)

# Plot both metrics vs depth
ggplot(vario_long, aes(x = value, y = z_mid, colour = variable)) +
  geom_path(size = 0.6) +
  facet_wrap(~ metric, scales = "free_x", labeller = as_labeller(metric_labels)) +
  scale_colour_manual(values = metal_cols, labels = metal_labels) +
  scale_y_reverse() +
  labs(x = NULL, y = "Depth (m)", colour = "Metal") +
  theme_minimal()

#-------------------------------------------------------------------------------
# Explanation
# This script:
# 1. Defines a helper to compute vertical-only variograms in contiguous depth bins.
# 2. Detrends each bin's log-grade series to isolate local variability.
# 3. Fits a spherical variogram model per bin, extracting range and sill.
# 4. Repeats for each metal, reshapes results, and plots both metrics against depth.
# The resulting curves reveal how correlation length and local variance change vertically,
# indicating non-stationarity that GeoWarp DepthDeform can address.





# ================================================
# Non‐Stationarity Cycle Analysis for Metals
# ================================================


period_results <- map_df(grade_vars, function(v) {
  # 1) grab & sort your sill series for this metal
  vs <- vario_stats %>%
    filter(variable == v) %>%
    arrange(z_mid)
  
  # 2) turn it into a ts sampled every bin_w metres
  ts_sill <- ts(vs$sill, start = vs$z_mid[1], deltat = bin_width)
  
  
  # 4) find the dominant period from the spectrum
  spec_res <- spectrum(ts_sill, log = "no", plot = FALSE)
  peak_f   <- spec_res$freq[which.max(spec_res$spec)]
  period_m <- (1/peak_f) * bin_width
  
  # 5) return a one‐row tibble
  tibble(variable = v, dominant_period_m = period_m)
})

print(period_results)




# ================================================
# Non‐Stationarity Cycle Analysis  for Your Metals
# ================================================




# ── Prepare a results container ──────────────────
period_results <- tibble(
  variable         = character(),
  acf_lag1         = numeric(),
  median_spacing_m = numeric()
)

# ── Loop over each metal ─────────────────────────
for (v in unique(vario_stats$variable)) {
  
  # 1) Subset & sort by depth
  vs       <- vario_stats %>% 
    filter(variable == v) %>% 
    arrange(z_mid)
  sill_vec <- vs$sill
  z_mid    <- vs$z_mid
  
  # 2) Detrend with loess
  trend_vec <- loess(sill_vec ~ z_mid)$fitted
  resid_vec <- sill_vec - trend_vec
  
  # 3) ACF of residuals (lag‐1 tells you the 15 m correlation)
  acf_res <- acf(resid_vec,
                 lag.max = 5,
                 plot    = FALSE)
  acf_lag1 <- acf_res$acf[2]  # the lag‐1 autocorrelation
  
  
  # 6) Peak‐to‐peak spacing in metres
  #    find local maxima of resid_vec
  peaks       <- which(diff(sign(diff(resid_vec))) == -2) + 1
  peak_depths <- z_mid[peaks]
  spacings    <- diff(peak_depths)          # metres between peaks
  median_sp   <- median(spacings, na.rm = TRUE)
  
  # 7) Store results
  period_results <- period_results %>%
    add_row(
      variable         = v,
      acf_lag1         = acf_lag1,
      median_spacing_m = median_sp
    )
}

# ── Print your summary table ─────────────────────
print(period_results)
