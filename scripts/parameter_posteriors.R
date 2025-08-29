.libPaths("C:/Users/23478671/R_libs")
library(dplyr)
library(purrr)
library(geowarp)
library(tidyr)
library(ggplot2)
library(patchwork)
library(stringr)

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
) %>% purrr::map(readRDS)

test_list  <- set_names(
  test_files,
  sub("^test_data_(.*)\\.rds$", "\\1", basename(test_files))
) %>% purrr::map(readRDS)
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


unnorm_z <- function(z_norm, rng) rng$min + z_norm * (rng$max - rng$min)

get_target <- function(fn)
  sub(".*_([^_]+_[^_]+)\\.rds$", "\\1", basename(fn))

profiles_all <- map_dfr(model_files, function(fn) {
  
  model_nm <- sub("\\.rds$", "", basename(fn))
  if (grepl("_cv", model_nm, ignore.case = TRUE))
    return(NULL)                       # skip stray CV files
  
  tgt          <- get_target(fn)
  ts           <- test_list[[tgt]]
  rng_z        <- ranges_list[[tgt]]$z
  rng_tg       <- ranges_list[[tgt]]$target
  depth_vec    <- -unnorm_iso_coord(ts$z, ranges_list[[tgt]], "z")
  depth_vec_norm    <- ts$z
  
  if (grepl("^gw_fit", model_nm, ignore.case = TRUE)) {
    ## ── GeoWarp fit ──────────────────────────────────────────────────────────
    fit          <- readRDS(fn)
    mean_vec     <- unnorm_vec(mean_profile(fit, df = ts), rng_tg)
    sd_vec       <- sqrt(marginal_variance_profile(fit, df = ts)) *
      (rng_tg$max - rng_tg$min)
    # warped_depth <- -unnorm_iso_coord_2(
    #   warped_coordinates(fit, df = ts)[ , 3],
    #   ranges_list[[tgt]], "z") 
    warped_depth <- -unnorm_iso_coord(ts$z, ranges_list[[tgt]], "z")*(warped_coordinates(fit, df = ts)[ , 3]/ts$z)
    warped_depth_norm <- warped_coordinates(fit, df = ts)[ , 3]
    
  } else {
    ## ── gstat / idw / ok / sk fit ────────────────────────────────────────────
    fit_df       <- as.data.frame(readRDS(fn))      # predictions on ts
    mean_vec     <- unnorm_vec(fit_df$var1.pred, rng_tg)
    sd_vec       <- sqrt(fit_df$var1.var) * (rng_tg$max - rng_tg$min)
    warped_depth <- NA_real_                          # not defined for these fits
    warped_depth_norm <- NA_real_
  }
  
  tibble(
    model       = model_nm,
    target      = tgt,
    depth       = depth_vec,
    mean        = mean_vec,
    sd          = sd_vec,
    warp_depth  = warped_depth,
    depth_norm = depth_vec_norm,
    warped_norm = warped_depth_norm
  )
})

profiles_all <- profiles_all %>% 
  mutate(
    fit = sub("(_[^_]+){2}$", "", model)   # gw_fit_Ag_ppm → gw_fit
  )


gw_df <- profiles_all %>%
  # keep only the GeoWarp fits
  filter(str_detect(model, "^gw_")) %>%
  # add a warping flag: “linear” if noWarp in the name, else “nonlinear”
  mutate(
    warping = if_else(
      str_detect(model, "noWarp"), 
      "linear", 
      "nonlinear"
    )
  )

# inspect
distinct(gw_df, model, warping)


  

# --- 1. pivot your gw_df (or profile_df) into long form if you haven’t already:
plot_df <- gw_df %>%
  pivot_longer(
    cols      = c(mean, sd),
    names_to  = "stat",
    values_to = "value"
  ) %>%
  # sort so each stat/target/warping group is ordered by depth
  arrange(stat, target, warping, depth)

# --- 2. define a common theme
my_theme <- theme_minimal() +
  theme(
    text              = element_text(face = "plain"),    # everything plain
    axis.line         = element_line(color = "black", size = 0.8),
    axis.ticks        = element_line(color = "black", size = 0.6),
    axis.ticks.length = unit(4, "pt"),
    axis.text         = element_text(size = 10, color = "black", face = "plain"),
    axis.title        = element_text(size = 12, face = "plain"),
    panel.grid.major  = element_line(color = "grey80", size = 0.3),
    panel.grid.minor  = element_blank(),
    strip.text        = element_text(size = 11, face = "plain"),
    panel.background  = element_rect(fill = "white", color = NA)
  )

# --- 3. build the “mean” row (one panel per metal, free x&y, grouped by warping)
# --- 3. “mean” row ----------------------------------------------------------
# --- label map -------------------------------------------------------------
target_labs <- c(
  Ag_ppm = "log Ag (ppm)",
  Au_ppb = "log Au (ppb)",
  Cu_pct = "log Cu (%)",
  Ni_pct = "log Ni (%)"
)

# --- 3. “mean” row ----------------------------------------------------------
# 0. define a nice palette
warp_cols <- c(
  linear    = "#1b9e77",   # green
  nonlinear = "#d95f02"    # orange
)

# --- 3. “mean” row ----------------------------------------------------------
# --- “mean” row -------------------------------------------------------------
p_mean <- plot_df %>%
  filter(stat == "mean") %>%
  ggplot(aes(x = value, y = depth,
             colour = warping, group = warping)) +
  geom_line(size = 0.6) +
  scale_colour_manual(values = warp_cols,
                      labels  = c("LinearWarp", "DepthDeform"),
                      name    = "Warping") +
  scale_y_reverse() +
  facet_wrap(~ target, nrow = 1, scales = "free",
             strip.position = "top",
             labeller = as_labeller(target_labs)) +
  labs(
    #title = expression(mu(h)),        # ←  μ(h)
    y     = "Depth (m)",
    x     = NULL
  ) +
  my_theme +
  theme(
    legend.position = "none",
    plot.title      = element_text(hjust = 0.5, face = "bold")
  )

# --- “sd” row ---------------------------------------------------------------
p_sd <- plot_df %>%
  filter(stat == "sd") %>%
  ggplot(aes(x = value, y = depth,
             colour = warping, group = warping)) +
  geom_path(size = 0.6) +
  scale_colour_manual(values = warp_cols,
                      labels  = c("LinearWarp", "DepthDeform"),
                      name    = "Warping") +
  scale_y_reverse() +
  facet_wrap(~ target, nrow = 1, scales = "free",
             strip.position = "bottom",
             labeller = as_labeller(target_labs)) +
  labs(
    #title = expression(sigma[delta](h)),   # ←  σ<sub>Δ</sub>(h)
    y     = "Depth (m)",
    x     = NULL
  ) +
  my_theme +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )


# # --- stack & draw ----------------------------------------------------------
# final_plot <- (p_mean / p_sd) +
#   plot_layout(ncol = 1, heights = c(1, 1), guides = "collect") &
#   theme(
#     legend.position  = "bottom",
#     strip.placement  = "outside"
#   )
# 
# print(final_plot)
# 
# 


# 1. common grid of 200 depths
depth_grid <- seq(min(gw_df$depth_norm),
                  max(gw_df$depth_norm),
                  length.out = 200)

# 2. build splines for each fit
splines <- gw_df %>%
  filter(grepl("^gw_fit", model)) %>%
  group_by(model, warping, target) %>%
  arrange(depth_norm, .by_group = TRUE) %>%
  summarise(
    spl = list(smooth.spline(depth_norm, warped_norm, spar = 0.6)),
    .groups = "drop"
  )

# ─────────────────────────────────────────────────────────────────────────────
# A) Warp curves
# ─────────────────────────────────────────────────────────────────────────────
preds_warp <- splines %>%
  mutate(
    warp_df = map(spl, ~{
      pr <- predict(.x, depth_grid)
      tibble(depth_norm = pr$x, warped_norm = pr$y)
    })
  ) %>%
  select(-spl) %>%
  unnest(warp_df)

warp_summary <- preds_warp %>%
  group_by(warping, target, depth_norm) %>%
  summarise(
    med_warp = median(warped_norm),
    lo25     = quantile(warped_norm, .25),
    hi75     = quantile(warped_norm, .75),
    .groups  = "drop"
  )


# ─────────────────────────────────────────────────────────────────────────────
# B) Derivative curves
# ─────────────────────────────────────────────────────────────────────────────
preds_deriv <- splines %>%
  mutate(
    deriv_df = map(spl, ~{
      pr <- predict(.x, depth_grid, deriv = 1)
      tibble(depth_mid = pr$x, deriv = pr$y)
    })
  ) %>%
  select(-spl) %>%
  unnest(deriv_df)




deriv_summary <- preds_deriv %>%
  group_by(warping, target, depth_mid) %>%
  summarise(
    med_der = median(deriv),
    lo25    = quantile(deriv, .25),
    hi75    = quantile(deriv, .75),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 1. convert depth_mid → real depth (m)
# ------------------------------------------------------------

deriv_summary <- deriv_summary %>%             # has depth_mid in [0,1]
  rowwise() %>%                                # because ranges differ by target
  mutate(
    depth_m = unnorm_iso_coord(               # minus sign if “down = positive”
      depth_mid,
      ranges_list[[as.character(target)]],     # ↩ target-specific ranges
      coord = "z"
    )
  ) %>%
  ungroup()

deriv_summary <- deriv_summary %>%          # already has depth_m
  arrange(target, warping, depth_m) %>%          # make sure depth ascends
  group_by(target, warping) %>%
  # find the last index where the derivative still changes
  mutate(idx        = row_number(),
         is_change  = med_der != lag(med_der, default = first(med_der)),
         last_var_i = max(idx[is_change], na.rm = TRUE)) %>%
  # keep rows up to that last-changing index
  filter(idx <= last_var_i) %>%
  select(-idx, -is_change, -last_var_i) %>%      # cleanup
  ungroup()

deriv_summary <- deriv_summary %>%               # already has depth_m
  group_by(target) %>%                           # do the calc per metal
  mutate(
    depth_m = depth_m - (max(depth_m, na.rm = TRUE) +
                           min(depth_m, na.rm = TRUE))
  ) %>%
  ungroup()


# ------------------------------------------------------------
# 2. plot against depth_m instead of depth_mid
# ------------------------------------------------------------
p_deriv <- ggplot(
  deriv_summary,
  aes(med_der, depth_m,
      colour = warping, fill = warping, group = warping)
) +
  geom_ribbon(aes(xmin = lo25, xmax = hi75),
              colour = NA, alpha = 0.25) +
  geom_path(size = 0.6) +
  scale_colour_manual(values = c(linear = "#1b9e77",
                                 nonlinear = "#d95f02"),
                      labels = c("LinearWarp", "DepthDeform"),
                      name   = "Warping") +
  scale_fill_manual(values = c(linear = "#1b9e77",
                               nonlinear = "#d95f02"),
                    labels = c("LinearWarp", "DepthDeform"),
                    name   = "Warping") +
  scale_y_reverse() +
  facet_wrap(
    ~ target,
    nrow   = 1,            # single row
    scales = "free_y"      # ← free Y (depth) for each facet
  ) +
  labs(
    #title = expression(df[3](h)/dh),
    y     = "Depth (m)",
    x     = NULL
  ) +
  my_theme +
  theme(
    legend.position = "bottom",
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    strip.placement = "outside"
  )

#print(p_deriv)





p_mean  <- p_mean  + labs(x = expression(mu(h)))         # row-1 label
p_sd    <- p_sd    + labs(x = expression(sigma[delta](h)))# row-2 label
p_deriv <- p_deriv + labs(x = expression(df[3](h)/dh))   # row-3 label


# leave p_mean as-is (strip.position = "bottom" so labels appear below panels)
# hide them in rows 2 & 3
p_sd    <- p_sd    + theme(strip.text.x = element_blank())
p_deriv <- p_deriv + theme(strip.text.x = element_blank())



warp_cols <- c(linear = "#1b9e77", nonlinear = "#d95f02")

scale_colour <- scale_colour_manual(
  values = warp_cols,
  labels = c("LinearWarp", "DepthDeform"),
  name   = "Warping"          # ← keep this legend
)

scale_fill <- scale_fill_manual(
  values = warp_cols,
  labels = c("LinearWarp", "DepthDeform"),
  name   = "Warping",
  guide  = "none"             # ← suppress the duplicate
)

# use these two scales in every plot ---------------------------
p_mean  <- p_mean  + scale_colour
p_sd    <- p_sd    + scale_colour
p_deriv <- p_deriv + scale_colour + scale_fill   # derivative needs both


# ── 2. remove bold from plot titles on all three rows ───────────────────────
p_mean  <- p_mean  + theme(plot.title = element_text(hjust = 0.5, face = "plain"))
p_sd    <- p_sd    + theme(plot.title = element_text(hjust = 0.5, face = "plain"))
p_deriv <- p_deriv + theme(plot.title = element_text(hjust = 0.5, face = "plain"))

# (everything else in your pipeline stays the same)
library(patchwork)


final_plot <- (p_mean / p_sd / p_deriv) +
  plot_layout(heights = c(1, 1, 1), guides = "collect") &
  theme(
    legend.position  = "bottom",
    strip.placement  = "outside",
    strip.background = element_blank(),
    text = element_text(face = "plain")   # ← forces plain for all text
  )

print(final_plot)












# Pick a color palette with enough hues for 8 lines:
library(RColorBrewer)
line_cols <- brewer.pal(8, "Dark2")

# Create a "metal+warping" label for legend clarity
deriv_summary$curve <- interaction(deriv_summary$target, deriv_summary$warping, sep = " - ")

p_deriv_all <- ggplot(
  deriv_summary,
  aes(x = med_der, y = depth_m, colour = curve, fill = curve, group = curve)
) +
  geom_ribbon(aes(xmin = lo25, xmax = hi75), alpha = 0.15, colour = NA) +
  geom_path(size = 0.8) +
  scale_colour_manual(values = line_cols, name = "Metal & Warping") +
  scale_fill_manual(values = line_cols, name = "Metal & Warping") +
  scale_y_reverse() +
  labs(
    #title = expression(df[3](h)/dh),
    y     = "Depth (m)",
    x     = expression(df[3](h)/dh)
  ) +
  my_theme +
  theme(
    legend.position = "right",
    plot.title      = element_text(hjust = 0.5)
  )

print(p_deriv_all)





# ── keep just the four metals we care about ──────────────────────────
metals <- c("Ag_ppm", "Au_ppb", "Cu_pct", "Ni_pct")

# ── 1. TEST points → original units  ─────────────────────────────────
test_points_df <- imap_dfr(test_list[metals], function(ts, tgt) {
  tibble(
    target = tgt,
    depth  = -unnorm_iso_coord(ts$z,      ranges_list[[tgt]], "z"),
    value  =  unnorm_vec      (ts$target, ranges_list[[tgt]]$target)
  )
})

# ── 2. prediction table (gw_df ALREADY un-normalised) ───────────────
plot_wide <- gw_df %>%                       # columns: target depth mean sd warping
  filter(target %in% metals) %>%             # only four metals
  arrange(target, warping, depth)            # ensure monotone depth per curve

# ── 3. common theme ─────────────────────────────────────────────────
my_theme <- theme_minimal() +
  theme(
    axis.line  = element_line(colour = "black", size = 0.8),
    axis.ticks = element_line(colour = "black", size = 0.6),
    axis.text  = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    panel.grid.major = element_line(colour = "grey80", size = 0.3),
    panel.grid.minor = element_blank(),
    strip.text       = element_text(size = 11, face = "bold")
  )

# ── 4. helper that draws one row (linear OR nonlinear) ──────────────
make_row <- function(df_row, ttl) {
  ggplot(df_row, aes(y = depth)) +
    geom_point(
      data = test_points_df,
      aes(x = value, y = depth),
      inherit.aes = FALSE,
      colour = "steelblue", alpha = .5, size = 0.4
    ) +
    geom_line(aes(x = mean),               colour = "black",  size = 0.6) +
    geom_line(aes(x = mean + 2*sd), linetype = "dashed", colour = "grey40") +
    geom_line(aes(x = mean - 2*sd), linetype = "dashed", colour = "grey40") +
    scale_y_reverse() +
    facet_wrap(~ target, nrow = 1, scales = "free") +
    labs(title = ttl, y = "Depth", x = NULL) +
    my_theme +
    theme(plot.title = element_text(hjust = 0.5))
}

# ── 5. build & stack rows ───────────────────────────────────────────
p_linear    <- make_row(filter(plot_wide, warping == "linear"),    "Linear warping")
p_nonlinear <- make_row(filter(plot_wide, warping == "nonlinear"), "Non-linear warping")

final_plot <- p_linear / p_nonlinear + plot_layout(ncol = 1)

print(final_plot)


