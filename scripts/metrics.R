.libPaths("C:/Users/23478671/R_libs")
library(dplyr)
library(purrr)
library(geowarp)

setwd("C:/Users/23478671/Github/GeoWarp-Application-Kevitsa")
source("helper_functions/00_utils_io.R")   # load_raw()
source("helper_functions/01_preprocess.R") # depth_filter(), normalise_xyz_target()
source("helper_functions/02_split.R")      # km_split()
source("helper_functions/03_models.R")     # fit_kriging(), fit_idw()



# ─────────────────────────────────────────────────────────────
# 0.  Setup & load everything
# ─────────────────────────────────────────────────────────────

fits_dir   <- "model_fits"           # sub-directory holding the *.rds files

model_files <- list.files(
  path       = fits_dir,
  pattern    = "\\.rds$",         # anything ending in .rds
  full.names = TRUE
)

get_target <- function(fn) sub(".*_([^_]+_[^_]+)\\.rds$", "\\1", basename(fn))
targets    <- unique(map_chr(model_files, get_target))

data_dir <- "train_test_data"   # sub-folder inside the repo

train_files <- list.files(
  path       = data_dir,
  pattern    = "^train_data_.*\\.rds$",
  full.names = TRUE
)

test_files  <- list.files(
  path       = data_dir,
  pattern    = "^test_data_.*\\.rds$",
  full.names = TRUE
)

train_list <- set_names(
  train_files,
  sub("^train_data_(.*)\\.rds$", "\\1", basename(train_files))
) %>% map(readRDS)

test_list  <- set_names(
  test_files,
  sub("^test_data_(.*)\\.rds$", "\\1", basename(test_files))
) %>% map(readRDS)
# helper to un-normalise a vector given its stored min/max


csv  <- "cluster_0_data.csv"

# ------------------------------------------------------------------
# simple helper -----------------------------------------------------------
get_ranges_from_csv <- function(targ) {
  raw   <- load_raw(csv, var = targ, sample_frac = 1)
  clean <- depth_filter(raw, min_cnt = 800)
  
  list(
    x      = list(min = min(clean$x,      na.rm = TRUE),
                  max = max(clean$x,      na.rm = TRUE)),
    y      = list(min = min(clean$y,      na.rm = TRUE),
                  max = max(clean$y,      na.rm = TRUE)),
    z      = list(min = min(clean$z,      na.rm = TRUE),
                  max = max(clean$z,      na.rm = TRUE)),
    target = list(min = min(clean$target, na.rm = TRUE),
                  max = max(clean$target, na.rm = TRUE))
  )
}


ranges_list <- setNames(lapply(targets, get_ranges_from_csv), targets)

# sanity check
print(ranges_list[["Cu_pct"]])


unnorm_vec <- function(x_norm, rng) x_norm * (rng$max - rng$min) + rng$min

# 1.  GeoWarp metrics (log space, plain names)
# ─────────────────────────────────────────────────────────────
metrics_on_geowarp <- function(fit, df_norm, rng_target) {
  
  pm_norm  <- mean_profile(fit, df = df_norm)                 # pred (normalised)
  psd_norm <- sqrt(marginal_variance_profile(fit, df = df_norm))
  
  # back-transform target, mean, and SD
  actual <- unnorm_vec(df_norm$target, rng_target)
  pred   <- unnorm_vec(pm_norm,       rng_target)
  psd    <- psd_norm * (rng_target$max - rng_target$min)    # scale only, no shift
  
  resid <- actual - pred
  
  tibble(
    ME             = mean(-resid),
    MAE            = mean(abs(resid), na.rm = TRUE),
    RMSE           = sqrt(mean(resid^2,              na.rm = TRUE)),
    MAPE           = mean(abs(resid / actual),       na.rm = TRUE) * 100,
    Pct_within_1SD = mean(abs(resid) <=  psd,         na.rm = TRUE) * 100,
    Pct_within_2SD = mean(abs(resid) <= 2 * psd,      na.rm = TRUE) * 100
  )
}


metrics_on_gstat <- function(pred_df, df_norm, rng_tg) {
  
  # un-normalise
  pred <- unnorm_vec(pred_df$var1.pred, rng_tg)
  sd   <- sqrt(pred_df$var1.var) * (rng_tg$max - rng_tg$min)
  obs  <- unnorm_vec(df_norm$target,  rng_tg)
  
  tibble(
    MAE            = mean(abs(obs - pred),              na.rm = TRUE),
    RMSE           = sqrt(mean((obs - pred)^2,          na.rm = TRUE)),
    Pct_within_1SD = mean(abs(obs - pred) <= sd,        na.rm = TRUE),
    Pct_within_2SD = mean(abs(obs - pred) <= 2*sd,    na.rm = TRUE)
  )
}

# ── loop over every .rds file ────────────────────────────────────────────────
gw_results <- map_dfr(model_files, function(mf) {
  
  target <- get_target(mf)
  rngs   <- ranges_list[[target]]
  model  <- sub(paste0("_", target, "\\.rds$"), "", basename(mf))
  
  df_norm <- test_list[[target]]
  if (is.null(df_norm)) {
    warning("No test data for ", target, "; skipping ", basename(mf))
    return(NULL)
  }
  
  if (grepl("^gw_", model, ignore.case = TRUE)) {
    # GeoWarp fit
    fit <- readRDS(mf)
    out <- metrics_on_geowarp(fit, df_norm, rngs$target)
    
  } else {
    # Classical prediction frame
    fit_df <- as.data.frame(readRDS(mf))
    out    <- metrics_on_gstat(fit_df, df_norm, rngs$target)
  }
  
  out %>% mutate(model = model, target = target, .before = 1)
})






all_results_rounded <- gw_results %>%
  # round each metric to the desired precision
  mutate(
    ME             = round(ME,             2),
    MAE            = round(MAE,            2),
    RMSE           = round(RMSE,           2),
    MAPE           = round(MAPE,           0),
    Pct_within_1SD = round(Pct_within_1SD, 3),
    Pct_within_2SD = round(Pct_within_2SD, 3)
  ) %>%
  group_by(target)



filtered_results <- all_results_rounded %>%
  select(model, target, MAE, RMSE, Pct_within_1SD, Pct_within_2SD) %>% 
  # keep only the four targets
  filter(target %in% c("Ag_ppm", "Au_ppb", "Cu_pct","Ni_pct")) %>%
  filter(!model %in% c("gw_fit_cv", "gw_fit_noWarp_cv"))

filtered_results


