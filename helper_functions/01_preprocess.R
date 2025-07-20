# 01_preprocess.R -----------------------------------------------------------

library(dplyr)
source("helper_functions/00_utils_io.R")

## 1. depth_filter()  -------------------------------------------------------
depth_filter <- function(df, bin = 50, min_cnt = 400, min_samples = 15) {
  
  cnts <- df %>%
    mutate(db = floor(-z / bin) * bin) %>%
    count(db)
  
  cutoff <- cnts %>%
    arrange(desc(db)) %>%
    filter(purrr::map_lgl(db,
                          ~ all(cnts$n[cnts$db < .x] > min_cnt))) %>%
    slice(1) %>% pull(db)
  
  df %>%
    filter(-z <= cutoff) %>%                  # shallower than cutoff
    group_by(name) %>% filter(n() > min_samples) %>% ungroup()
}

## 2. normalise_xyz_target()  ----------------------------------------------
normalise_xyz_target <- function(df) {
  # per-axis min / max
  rx <- range(df$x); ry <- range(df$y); rz <- range(df$z)
  rt <- range(df$target)
  
  # common span = largest axis span
  span <- max(rx[2] - rx[1],
              ry[2] - ry[1],
              rz[2] - rz[1])
  
  # build a data.frame scaled by that common span
  out <- df %>%
    mutate(
      x      = (x - rx[1]) / span,
      y      = (y - ry[1]) / span,
      z      = (z - rz[1]) / span,
      target = (target - rt[1]) / (rt[2] - rt[1])   # target still 0–1
    )
  
  # fake “max” = min + span so the old unscale() keeps working unchanged
  adj_rx <- list(min = rx[1], max = rx[1] + span)
  adj_ry <- list(min = ry[1], max = ry[1] + span)
  adj_rz <- list(min = rz[1], max = rz[1] + span)
  
  attr(out, "ranges") <- list(x = adj_rx, y = adj_ry,
                              z = adj_rz, target = list(min = rt[1], max = rt[2]))
  out
}


## 3. unnormalise_xyz_target()  --------------------------------------------
unnormalise_xyz_target <- function(norm_df) {
  # pull back the per‐axis ranges you stashed in the attr
  rng <- attr(norm_df, "ranges")
  if (is.null(rng)) stop("No 'ranges' attribute found on this data frame.")
  
  # each of rng$x, rng$y, rng$z is a list(min=..., max=...)
  # rng$target likewise has the original target range
  x_min <- rng$x$min; x_max <- rng$x$max
  y_min <- rng$y$min; y_max <- rng$y$max
  z_min <- rng$z$min; z_max <- rng$z$max
  t_min <- rng$target$min; t_max <- rng$target$max
  
  norm_df %>%
    mutate(
      # invert x,y,z all using their common span
      x = x * (x_max - x_min) + x_min,
      y = y * (y_max - y_min) + y_min,
      z = z * (z_max - z_min) + z_min,
      # invert target back to its original [t_min, t_max] interval
      target = target * (t_max - t_min) + t_min
    )
}


unnorm_iso_coord <- function(coord_norm, ranges, coord = "z") {
  
  if (!coord %in% c("x", "y", "z"))
    stop("'coord' must be one of 'x', 'y', 'z'")
  
  # 1.  Pull mins and maxs ----------------------------------------------------
  needed <- c("x", "y", "z")
  if (!all(needed %in% names(ranges)))
    stop("ranges must contain $x, $y, and $z")
  
  get_delta <- function(ax) ranges[[ax]]$max - ranges[[ax]]$min
  scale <- max(vapply(needed, get_delta, numeric(1)))
  
  coord_min <- ranges[[coord]]$min
  
  # 2.  If the vector already looks un-scaled (any value outside [0, 1]),
  #     leave it untouched; otherwise un-scale it.
  if (all(coord_norm >= 0 & coord_norm <= 1, na.rm = TRUE)) {
    coord_min + coord_norm * scale          # un-scaled metres
  } else {
    coord_norm                              # already metres
  }
}

unnorm_iso_coord_2 <- function(coord_norm, ranges, coord = "z") {
  
  if (!coord %in% c("x", "y", "z"))
    stop("'coord' must be one of 'x', 'y', 'z'")
  
  # 1.  Pull mins and maxs ----------------------------------------------------
  needed <- c("x", "y", "z")
  if (!all(needed %in% names(ranges)))
    stop("ranges must contain $x, $y, and $z")
  
  get_delta <- function(ax) ranges[[ax]]$max - ranges[[ax]]$min
  scale <- max(vapply(needed, get_delta, numeric(1)))
  
  coord_min <- ranges[[coord]]$min
  
  # 2.  If the vector already looks un-scaled (any value outside [0, 1]),
  #     leave it untouched; otherwise un-scale it.
  if (all(coord_norm >= 0 & coord_norm <= 1, na.rm = TRUE)) {
    coord_min + coord_norm * scale          # un-scaled metres
  } else {
    coord_min + coord_norm * scale                             # already metres
  }
}








unnorm_vec <- function(x_norm, rng) x_norm * (rng$max - rng$min) + rng$min

