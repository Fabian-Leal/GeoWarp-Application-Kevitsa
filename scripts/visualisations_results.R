.libPaths("C:/Users/23478671/R_libs")
setwd("C:/Users/23478671/Github/Geowarp-Resource-Estimation/GeoWarp_resource_estimation/application_kevitsa")
source("helper_functions/00_utils_io.R")   # load_raw()
source("helper_functions/01_preprocess.R") # depth_filter(), normalise_xyz_target()
source("helper_functions/02_split.R")      # km_split()
source("helper_functions/03_models.R")     # fit_kriging(), fit_idw()
library(hexbin)
library(lmtest)
library(tibble)
library(dplyr)
library(purrr)
library(geowarp)
library(ggplot2)
library(patchwork)  # make sure it’s loaded
library(scales)   # pretty_breaks()


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

csv  <- "cluster_0_data.csv"

keep_targets <- c("Ag_ppm", "Au_ppb", "Cu_pct", "Ni_pct")
targets <- keep_targets

ranges_list <- setNames(lapply(keep_targets, get_ranges_from_csv), keep_targets)


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



# ── helper: undo min–max normalisation ────────────────────────────────

fits_dir   <- "model_fits"           # sub-directory holding the *.rds files

model_files <- list.files(
  path       = fits_dir,
  pattern    = "\\.rds$",         # anything ending in .rds
  full.names = TRUE
)

model_files <- model_files[!grepl("_cv", basename(model_files), ignore.case = TRUE)]

get_target  <- function(fn)
  sub(".*_([^_]+_[^_]+)\\.rds$", "\\1", basename(fn))        # e.g. Ni_pct

keep_idx    <- grepl(paste(keep_targets, collapse = "|"),
                     get_target(model_files), ignore.case = TRUE)

model_files <- model_files[keep_idx]
if (length(model_files) == 0) stop("No matching model files found.")






get_target <- function(fn)
  sub(".*_([^_]+_[^_]+)\\.rds$", "\\1", basename(fn))

profiles_all <- map_dfr(model_files, function(fn) {
  
  model_nm <- sub("\\.rds$", "", basename(fn))
  if (grepl("_cv", model_nm, ignore.case = TRUE))
    return(NULL)                        # skip any stray CV files
  
  tgt <- get_target(fn)
  
  ## always start from the test set for this target --------------------------
  ts      <- test_list[[tgt]]
  rng_z   <- ranges_list[[tgt]]$z
  rng_tg  <- ranges_list[[tgt]]$target
  depth_vec <- -unnorm_iso_coord(ts$z, ranges_list[[tgt]], "z")
  
  if (grepl("^gw_fit", model_nm, ignore.case = TRUE)) {
    ## ── GeoWarp: predict on *ts* only ───────────────────────────────────────
    fit       <- readRDS(fn)
    mean_vec  <- unnorm_vec(mean_profile(fit, df = ts), rng_tg)
    sd_vec    <- sqrt(marginal_variance_profile(fit, df = ts)) *
      (rng_tg$max - rng_tg$min)
    
  } else {
    ## ── Classical gstat fits (ok / ok_m / sk / idw) ─────────────────────────
    fit_df    <- as.data.frame(readRDS(fn))    # predictions for exactly ts
    mean_vec  <- unnorm_vec(fit_df$var1.pred,  rng_tg)
    sd_vec    <- sqrt(fit_df$var1.var) *
      (rng_tg$max - rng_tg$min)
  }
  
  tibble(
    model  = model_nm,
    target = tgt,
    depth  = depth_vec,
    mean   = mean_vec,
    sd     = sd_vec
  )
})

profiles_all <- profiles_all %>% 
  mutate(
    fit = sub("(_[^_]+){2}$", "", model)   # gw_fit_Ag_ppm → gw_fit
  )





# ── 1.   Small helper: one grid for one fitting method ───────────────────────
make_one_grid <- function(df_fit,
                          col         = "steelblue",   # ribbon + line colour
                          point_col   = "black",       # test-point colour
                          point_alpha = 0.1,
                          point_size  = 0.3) {
  
  # ------------------------------------------------------------------ #
  # 1.  Build a tidy data-frame of test points for every target facet  #
  # ------------------------------------------------------------------ #
  df_points <- purrr::map_dfr(unique(df_fit$target), function(tgt) {
    
    ts      <- test_list[[tgt]]
    rng_z   <- ranges_list[[tgt]]$z
    rng_tg  <- ranges_list[[tgt]]$target
    
    tibble(
      target = tgt,
      depth  = -unnorm_iso_coord(ts$z, ranges_list[[tgt]], "z"),  # <- depth (m)
      obs    = unnorm_vec(ts$target, rng_tg)                     # <- grade units
    )
  })
  
  # ------------------------------------------------------------- #
  # 2.  Plot: ribbon + mean line + scattered test observations    #
  # ------------------------------------------------------------- #
  ggplot(df_fit, aes(x = mean, y = depth)) +
    
    # shaded ±1 SD band
    geom_ribbon(aes(xmin = mean - 2*sd, xmax = mean + 2*sd),
                fill = col, alpha = 0.20, colour = NA) +
    
    # mean profile line
    geom_path(colour = col, linewidth = 0.6) +
    
    # raw test samples
    geom_point(data = df_points,
               aes(x = obs, y = depth),
               colour = point_col, alpha = point_alpha, size = point_size,
               inherit.aes = FALSE) +
    
    scale_y_reverse("Depth (m)", breaks = scales::pretty_breaks(10)) +
    facet_wrap(~ target, scales = "free_x") +
    labs(
      title = df_fit$fit[1],       # “gw_fit”, “ok_m”, …
      x     = "Grade (original units)"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x     = element_text(angle = 45, hjust = 1)
    )
}




# ── 2.   Master function: choose fits → list of grids ────────────────────────
make_grids_by_fit <- function(data, fits = NULL) {
  if (!is.null(fits)) data <- filter(data, fit %in% fits)
  if (nrow(data) == 0)
    stop("Nothing left after filtering for fits = ",
         paste(fits, collapse = ", "))
  
  split(data, data$fit) |>
    lapply(make_one_grid)                  # named list of ggplot objects
}

# ── 3.   How to use it --------------------------------------------------------

grids <- make_grids_by_fit(profiles_all,     # <- your five-column data
                           fits = c("gw_fit", "gw_fit_noWarp"))

# View the GeoWarp grid only
print(grids$gw_fit_noWarp)









# ── 0. Build combined data-frame --------------------------------------------
preds_long <- profiles_all %>% 
  mutate(
    fit = sub(paste0("_", target, "$"), "", model)   # gw_fit_Ag_ppm → gw_fit
  ) %>% 
  select(target, fit, value = mean)

test_long <- map_dfr(names(test_list), function(tg) {
  tibble(
    target = tg,
    fit    = "test points",
    value  = unnorm_vec(test_list[[tg]]$target, ranges_list[[tg]]$target)
  )
})

hist_df <- bind_rows(preds_long, test_long)

# ── 1. Drop SK if requested, rename GeoWarp variants, keep NAs already "test points" ─
drop_sk <- TRUE
hist_df2 <- hist_df %>% 
  { if (drop_sk) filter(., fit != "sk") else . } %>% 
  mutate(
    fit = recode(
      fit,
      "gw_fit"        = "gw_dd",
      "gw_fit_noWarp" = "gw_lw"
    )
  )

# ── 2. Compute 2 %–98 % limits per target for test points --------------------
trim_p <- 0.02
tp_limits <- hist_df2 %>%
  filter(fit == "test points") %>%
  group_by(target) %>%
  summarise(
    lo = quantile(value, trim_p,       na.rm = TRUE),
    hi = quantile(value, 1 - trim_p,   na.rm = TRUE),
    .groups = "drop"
  )

# ── 3. Trim only the test-point rows outside those limits --------------------
hist_df2 <- hist_df2 %>%
  left_join(tp_limits, by = "target") %>%
  filter(
    fit != "test points" |
      (value >= lo & value <= hi)
  ) %>%
  select(-lo, -hi)

# ── 4. Palette: grey for test points, hues for everything else --------------
other_fits <- setdiff(unique(hist_df2$fit), "test points")
pal        <- c(setNames(hue_pal()(length(other_fits)), other_fits),
                "test points" = "grey60")

# ── 5. Plot ------------------------------------------------------------------
n_bins <- 80

ggplot(hist_df2, aes(value, fill = fit, colour = fit)) +
  geom_histogram(position = "identity", alpha = 0.35, bins = n_bins) +
  facet_wrap(~ target, scales = "free") +
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal) +
  labs(
    x      = "Grade (log units)",
    y      = "Count",
    fill   = "Model / Source",
    colour = "Model / Source"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")








###################### WARPINGS #############


# helper (only once in your script)
unnorm_vec <- function(x_norm, rng) x_norm * (rng$max - rng$min) + rng$min

library(patchwork)   # for wrap_plots()

warp_plots <- purrr::map(model_files, function(mf) {
  fit    <- readRDS(mf)
  target <- get_target(mf)
  
  # 1 ── choose a depth grid (test set if available)
  df_norm <- if (exists("test_list") && !is.null(test_list[[target]])) {
    test_list[[target]]
  } else {
    fit$observed_df
  }
  
  # 2 ── ranges for un-normalising depth
  rng_z <- attr(df_norm, "ranges")$z
  if (is.null(rng_z)) rng_z <- list(min = min(df_norm$z), max = max(df_norm$z))  # fallback
  
  # 3 ── build dataframe & compute warped depths
  warping_df <- data.frame(x = 0, y = 0, z = df_norm$z)
  warping_df$z_warped <- geowarp::warped_coordinates(fit, df = warping_df)[, 3]
  
  # 4 ── back-transform to original units
  warping_df$z_raw        <- unnorm_vec(warping_df$z,        rng_z)

  # 5 ── plot in real depth units
  ggplot(warping_df, aes(z_raw, z_warped)) +
    geom_line() +
    scale_x_reverse() +   # depth increases downward
    coord_flip() +        # put depth on vertical axis
    labs(
      x     = "Depth (m)",
      y     = "Warped depth (m)",
      title = paste("Warping function –", target, "\n", basename(mf))
    ) +
    theme_minimal(base_size = 11)
})

wrap_plots(warp_plots, ncol = 2)   # display all warp curves

