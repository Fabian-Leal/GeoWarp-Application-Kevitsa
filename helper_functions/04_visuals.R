# 04_visuals.R -------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(plotly)
library(tidyr)          # simplest

## ------------------------------------------------------------------------
## 1. depth_profile_plot()
## ------------------------------------------------------------------------
# df        : data.frame with columns name, x, y, z, target  (target is log-grade)
# bin_size  : vertical bin width in metres
depth_profile_plot <- function(df,
                               bin_size     = 50,
                               var_name     = NULL,
                               show_counts  = FALSE) {
  # ── 0. variable label --------------------------------------------------
  if (is.null(var_name))
    var_name <- attr(df, "var") %||% "grade"
  
  # ── 1. depth bins & statistics ----------------------------------------
  binned <- df |>
    dplyr::mutate(depth = -z,
                  bin   = floor(depth / bin_size) * bin_size) |>
    dplyr::group_by(bin) |>
    dplyr::summarise(
      count  = dplyr::n(),
      holes  = dplyr::n_distinct(name),
      avg    = mean(target, na.rm = TRUE),
      sd_val = sd(target,   na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      upper = avg + sd_val,
      lower = avg - sd_val
    ) |>
    dplyr::arrange(bin)                            # shallow → deep
  
  # ── 2. sample-count bars (optional) ------------------------------------
  if (show_counts) {
    rng   <- range(c(binned$lower, binned$upper), na.rm = TRUE)
    bar_w <- diff(rng) * 0.20
    binned <- binned |>
      dplyr::mutate(
        bar_start = rng[1] - bar_w,
        bar_end   = bar_start + (count / max(count)) * bar_w,
        y_lab_up  = bin + 0.20 * bin_size,
        y_lab_dn  = bin - 0.20 * bin_size
      )
  }
  
  # ── 3. depth limits and “nice” 10-tick grid ---------------------------
  y_top <- max(binned$bin) + 0.20 * bin_size
  y_bot <- min(binned$bin) - 0.20 * bin_size
  range_m <- y_top - y_bot                           # profile thickness
  
  # choose spacing: the smallest multiple of 10 that gives ≤ 9 intervals
  step <- ceiling(range_m / 9 / 10) * 10
  step <- max(step, 10)                              # safety
  
  # round top to next multiple of 10, then build 10 ticks downward
  first_tick <- ceiling(y_top / 10) * 10
  y_breaks   <- seq(first_tick, by = -step, length.out = 10)
  
  # ── 4. plot ------------------------------------------------------------
  p <- ggplot2::ggplot(binned, ggplot2::aes(y = bin)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(xmin = lower, xmax = upper),
      fill  = "#b24000", alpha = 0.18
    ) +
    ggplot2::geom_path(
      ggplot2::aes(x = avg, group = 1),
      colour = "#b24000", linewidth = 0.6
    ) +
    ggplot2::scale_y_reverse(
      limits = c(y_top, y_bot),
      breaks = y_breaks,
#      name   = ""
      name = "Depth (m)"
    ) +
    ggplot2::labs(
      x = paste0("log(", var_name, ")"),
      title = NULL
    ) +
    ggplot2::theme_minimal()
  
  # ── 5. add count bars if selected -------------------------------------
  if (show_counts) {
    p <- p +
      ggplot2::geom_segment(
        ggplot2::aes(x = bar_start, xend = bar_end,
                     y = bin,      yend = bin),
        linewidth = 3, colour = "grey80"
      ) +
      ggplot2::geom_text(
        ggplot2::aes(x = bar_end, y = y_lab_up, label = count),
        hjust = 0, size = 3
      ) +
      ggplot2::geom_text(
        ggplot2::aes(x = bar_end, y = y_lab_dn, label = holes),
        hjust = 0, size = 3, colour = "grey40"
      )
  }
  
  return(p)
}


## ------------------------------------------------------------------------
## 2. train_test_scatter3D()
## ------------------------------------------------------------------------
# df          : full (unnormalised) dataset with a 'split' column ("Train"/"Test")
# colour_train, colour_test : colours for the two sets
train_test_scatter3D <- function(df,
                                 ranges = NULL,          # attr from normalise_xyz_target()
                                 colour_train = "blue",
                                 colour_test  = "red") {
  
  # If ranges were supplied → back-transform to metres
  if (!is.null(ranges)) {
    unscale <- function(v, rng) v * (rng$max - rng$min) + rng$min
    df <- df %>%
      mutate(
        x = unscale(x, ranges$x),
        y = unscale(y, ranges$y),
        z = unscale(z, ranges$z)
      )
  }
  
  cols <- c(Train = colour_train, Test = colour_test)
  
  plot_ly(
    df,
    x = ~x, y = ~y, z = ~z,
    color  = ~split, colors = cols,
    type   = "scatter3d", mode = "markers",
    marker = list(size = 3, opacity = 0.8)
  ) %>%
    layout(
      title  = "Train vs Test holes (real metres)",
      legend = list(title = list(text = "Split")),
      scene  = list(xaxis = list(title = "Easting"),
                    yaxis = list(title = "Northing"),
                    zaxis = list(title = "Elevation / depth"))
    )
}


# -------------------------------------------------------------------------
# top_hole_prediction_plot()
# -------------------------------------------------------------------------
# test_df      : test set with (name,x,y,z) *and* un-scaled columns
#                z_raw, target_raw, predicted_mean_raw, krig_pred_raw, idw_pred_raw
# avg_distances: tibble with hole1, hole2, avg_distance (can be NULL)
# n_holes      : number of richest holes to display
# -------------------------------------------------------------------------
top_hole_prediction_plot <- function(test_df,
                                     avg_distances = NULL,
                                     n_holes       = 20,
                                     model_cols    = NULL,
                                     model_names   = NULL) {
  
  # ---- sanity checks ------------------------------------------------------
  needed <- c("z_raw", "target_raw")
  stopifnot(all(needed %in% names(test_df)))
  
  # ---- discover / validate model columns ---------------------------------
  if (is.null(model_cols)) {
    model_cols <- grep("_pred_raw$", names(test_df), value = TRUE)
  }
  if (length(model_cols) == 0)
    stop("No model prediction columns found.")
  
  if (is.null(model_names)) {
    # crude but useful replacements
    model_names <- model_cols %>%
      sub("_pred_raw$", "", .) %>%             # drop suffix
      sub("^predicted_mean",  "GeoWarp", .) %>%
      sub("^krig_pred",       "Kriging",  .) %>%
      sub("^idw_pred",        "IDW",      .)
  }
  
  if (length(model_cols) != length(model_names))
    stop("model_cols and model_names must have the same length")
  
  # rename the chosen columns to friendly names in a copy
  plot_df <- test_df
  names(plot_df)[match(model_cols, names(plot_df))] <- model_names
  
  # ---- choose the top-N holes by #samples --------------------------------
  top_holes <- plot_df %>%
    dplyr::count(name, name = "n_obs") %>%
    dplyr::arrange(dplyr::desc(n_obs)) %>%
    dplyr::slice_head(n = n_holes) %>%
    dplyr::pull(name)
  
  # ---- long format for predictions ---------------------------------------
  pred_long <- plot_df %>%
    dplyr::filter(name %in% top_holes) %>%
    dplyr::select(name, z_raw, target_raw, all_of(model_names)) %>%
    tidyr::pivot_longer(cols = all_of(model_names),
                        names_to = "method",
                        values_to = "prediction") %>%
    dplyr::mutate(depth = -z_raw)          # depth positive downward
  
  obs_long <- plot_df %>%
    dplyr::filter(name %in% top_holes) %>%
    dplyr::transmute(name, depth = -z_raw, target = target_raw)
  
  # ---- facet labels with optional NN info --------------------------------
  if (!is.null(avg_distances)) {
    nn_info <- dplyr::bind_rows(
      dplyr::rename(avg_distances, hole = hole1, neighbour = hole2),
      dplyr::rename(avg_distances, hole = hole2, neighbour = hole1)
    ) %>%
      dplyr::group_by(hole) %>% dplyr::slice_min(avg_distance, n = 1) %>%
      dplyr::ungroup()
    
    label_map <- nn_info %>%
      dplyr::filter(hole %in% top_holes) %>%
      dplyr::mutate(label = paste0(
        hole, "\nNN: ", neighbour, " (", round(avg_distance, 1), " m)"
      )) %>%
      {setNames(.$label, .$hole)}
  } else {
    label_map <- setNames(top_holes, top_holes)
  }
  
  # nice distinct colours
  pal <- scales::hue_pal()(length(model_names))
  names(pal) <- model_names
  
  # ---- final ggplot -------------------------------------------------------
  ggplot() +
    geom_point(data = pred_long,
               aes(prediction, depth, colour = method),
               shape = 16, size = 1.2) +
    geom_point(data = obs_long,
               aes(target, depth, shape = "Observed"),
               colour = "black", size = 0.8) +
    scale_colour_manual(name = "Method", values = pal) +
    scale_shape_manual(name = "", values = c(Observed = 4)) +
    guides(
      colour = guide_legend(order = 1),
      shape  = guide_legend(order = 2,
                            override.aes = list(colour = "black"))
    ) +
    facet_wrap(~ name,
               scales   = "free_x",
               ncol     = 5,
               labeller = as_labeller(label_map)) +
    scale_y_reverse(name = "Depth (m)") +
    labs(title = "Predictions vs Observed Grade (Top Holes)",
         x     = "log(grade)") +
    theme_minimal(base_size = 8) +
    theme(strip.text    = element_text(size = 7),
          panel.spacing = unit(0.5, "lines"))
}




# ------------------------------------------------------------------------
# profile_and_warp_plots()
# ------------------------------------------------------------------------
# fit        – a geowarp fit object (from geowarp_optimise())
# test_df    – dataframe in *scaled* space containing at least x,y,z,target
#              (scaled 0–1, exactly what you used for mean_profile())
# x0, y0     – where to slice the vertical profile.  If NULL take the mid-point
# bin_size   – only used for the axis label / spacing in warping plot
#
# returns a patchwork object:  | Mean profile | SD profile | Warping fn |
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# profile_and_warp_plots()  –  now un-scales depth AND grade
# ------------------------------------------------------------------------
# fit      : geowarp fit object
# test_df  : data.frame in scaled space (x,y,z,target,name)  +  attr(., "ranges")
# x0,y0    : profile slice; NULL → mid-point of test set
# bin_size : only used inside the warping plot axis tick spacing
# ------------------------------------------------------------------------
profile_and_warp_plots <- function(fit, test_df,
                                   x0 = NULL, y0 = NULL,
                                   bin_size = 0.05) {
  
  # ---- safety ------------------------------------------------------------
  rng <- attr(test_df, "ranges")
  stopifnot(!is.null(rng),
            all(c("x","y","z","target") %in% names(test_df)))
  
  unscale_vec <- function(v, r) v * (r$max - r$min) + r$min
  
  # ---- choose slice centre ----------------------------------------------
  if (is.null(x0)) x0 <- mean(test_df$x, na.rm = TRUE)
  if (is.null(y0)) y0 <- mean(test_df$y, na.rm = TRUE)
  
  # ---- build profile (still in scaled space) -----------------------------
  profile_df <- data.frame(
    x = x0,
    y = y0,
    z = sort(unique(test_df$z))
  ) %>%
    dplyr::mutate(
      gwarp_mean  = mean_profile(fit, df = .),
      gwarp_stdev = sqrt(marginal_variance_profile(fit, df = .))
    )
  
  # ---- empirical mean / sd by depth -------------------------------------
  stats_by_z <- test_df %>%
    dplyr::group_by(z) %>%
    dplyr::summarise(
      obs_mean = mean(target, na.rm = TRUE),
      obs_sd   = sd(target,   na.rm = TRUE),
      .groups  = "drop"
    )
  
  profile_df <- dplyr::left_join(profile_df, stats_by_z, by = "z")
  
  # ---- un-scale & un-log -------------------------------------------------
  profile_df <- profile_df %>%
    dplyr::mutate(
      depth_m        = unscale_vec(z,      rng$z),           # metres
      gwarp_mean_raw = exp(unscale_vec(gwarp_mean, rng$target)),
      obs_mean_raw   = exp(unscale_vec(obs_mean,   rng$target)),
      # ±2σ band (convert limits, not the SD itself)
      upper_raw      = exp(unscale_vec(gwarp_mean + 2 * gwarp_stdev,
                                       rng$target)),
      lower_raw      = exp(unscale_vec(gwarp_mean - 2 * gwarp_stdev,
                                       rng$target))
    )
  
  # ---- PLOT 1  –  mean ---------------------------------------------------
  p_mean <- ggplot(profile_df) +
    geom_line(aes(depth_m, gwarp_mean_raw), colour = "blue") +
    # geom_line(aes(depth_m, obs_mean_raw),  colour = "red") +  # optional
    scale_x_reverse() + coord_flip() +
    labs(x = "Depth (m)", y = "Grade",
         title = "Mean profile (original units)")
  
  # ---- PLOT 2  –  ±2 σ band ---------------------------------------------
  p_sd <- ggplot(profile_df) +
    geom_line(aes(depth_m, gwarp_stdev), colour = "blue") +
      scale_x_reverse() + coord_flip() +
    labs(x = "Depth (m)", y = "Grade",
         title = "Sqrt(σ) ")
  
  # ---- PLOT 3  –  warping function --------------------------------------
  warping_df <- data.frame(x = 0, y = 0, z = profile_df$z) %>%
    dplyr::mutate(z_warped = warped_coordinates(fit, df = .)[, 3],
                  depth_m  = unscale_vec(z, rng$z))
  
  p_warp <- ggplot(warping_df, aes(depth_m, z_warped)) +
    geom_line() +
    scale_x_reverse() + coord_flip() +
    labs(x = "Depth (m)", y = "Warped depth",
         title = "Warping function")
  
  # ---- return combined panel --------------------------------------------
  patchwork::wrap_plots(p_mean, p_sd, p_warp, nrow = 1)
}

# ------------------------------------------------------------------------
# target_profile_bandplot()
# ------------------------------------------------------------------------
# fit        – a geowarp fit (from geowarp_optimise())
# test_df    – data.frame in scaled space with x,y,z,target,name
# x0, y0     – horizontal location for the vertical profile (NULL = mid‐point)
# colour_mean, alpha_obs – aesthetic tweaks
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# target_profile_bandplot()  –  now plots in original units
# ------------------------------------------------------------------------
# fit        – geowarp fit object
# test_df    – data.frame in scaled space with x,y,z,target,name
#              and attr(test_df, "ranges") from normalise_xyz_target()
# x0,y0      – slice coordinates (NULL → mid-point of test set)
# colour_mean, alpha_obs – aesthetics
# ------------------------------------------------------------------------
target_profile_bandplot <- function(fit, test_df,
                                    x0 = NULL, y0 = NULL,
                                    colour_mean = "blue",
                                    alpha_obs   = 0.10) {
  
  rng <- attr(test_df, "ranges")
  stopifnot(!is.null(rng),
            all(c("x","y","z","target","name") %in% names(test_df)))
  
  unscale_vec <- function(v, r) v * (r$max - r$min) + r$min
  
  # ── 1. choose slice centre ─────────────────────────────────────────────
  if (is.null(x0)) x0 <- mean(test_df$x, na.rm = TRUE)
  if (is.null(y0)) y0 <- mean(test_df$y, na.rm = TRUE)
  
  # ── 2. build profile in scaled space, then back-transform ──────────────
  profile_df <- data.frame(
    x = x0, y = y0, z = sort(unique(test_df$z))
  ) %>%
    dplyr::mutate(
      gwarp_mean  = mean_profile(fit, df = .),
      gwarp_stdev = sqrt(marginal_variance_profile(fit, df = .)),
      depth_m     = unscale_vec(z,      rng$z),
      mean_raw    = exp(unscale_vec(gwarp_mean, rng$target)),
      upper_raw   = exp(unscale_vec(gwarp_mean + 2 * gwarp_stdev, rng$target)),
      lower_raw   = exp(unscale_vec(gwarp_mean - 2 * gwarp_stdev, rng$target))
    )
  
  # ── 3. back-transform observed samples for spaghetti ───────────────────
  obs_df <- test_df %>%
    dplyr::transmute(
      name,
      depth_m = unscale_vec(z,      rng$z),
      target  = exp(unscale_vec(target, rng$target))
    )
  
  # ── 4. final ggplot ────────────────────────────────────────────────────
  ggplot() +
    # spaghetti – every observed sample
    geom_line(data = obs_df,
              aes(depth_m, target, group = name),
              alpha = alpha_obs) +
    
    # GeoWarp mean and ±2σ band
    geom_line(data = profile_df,
              aes(depth_m, mean_raw), colour = colour_mean) +
    geom_ribbon(data = profile_df,
                aes(depth_m, ymin = lower_raw, ymax = upper_raw),
                fill = colour_mean, alpha = 0.15) +
    
    scale_x_reverse() + coord_flip() +
    labs(x = "Depth (m)", y = "Grade",
         title = "Target profile with GeoWarp ±2σ band") +
    theme_minimal()
}


# ------------------------------------------------------------------------
# collar_xy_plot() – 2D map of drill‐hole collars, coloured by split
# ------------------------------------------------------------------------
# plot_df    : data.frame with columns `name`, `x`, `y`, and `split` ("Train"/"Test")
# point_size : size of the points
# alpha      : transparency level
# ------------------------------------------------------------------------
collar_xy_plot <- function(plot_df, point_size = 2, alpha = 0.8) {
  # lazy‐load dependencies
  library(dplyr)
  library(ggplot2)
  
  collars <- plot_df %>%
    group_by(name, split) %>%
    summarise(
      x = first(x),
      y = first(y),
      .groups = "drop"
    )
  
  ggplot(collars, aes(x = x, y = y, color = split)) +
    geom_point(size = point_size, alpha = alpha) +
    scale_color_manual(values = c(Train = "blue", Test = "red")) +
    labs(
      title = "Collar Locations (Train vs Test)",
      x     = "Easting (m)",
      y     = "Northing (m)",
      color = "Split"
    ) +
    theme_minimal()
}



# ------------------------------------------------------------------------
# example_and_mean_profile_plots()
# ------------------------------------------------------------------------
# df        – data.frame with columns: name (hole ID), z (depth), target (log-grade)
# hole_ids  – optional vector of hole IDs to show as “examples”
#               default = first three unique names
# var       – column name of the log‐grade (as string)
# var_label – axis label for the log‐grade
# ------------------------------------------------------------------------
example_and_mean_profile_plots <- function(df,
                                           hole_ids   = NULL,
                                           var        = "target",
                                           var_label  = "log(grade)") {
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  
  # 1) choose example holes
  if (is.null(hole_ids)) {
    hole_ids <- df %>% distinct(name) %>% slice_head(n = 3) %>% pull(name)
  }
  df_ex <- df %>% filter(name %in% hole_ids)
  
  # 2) compute mean & sd by depth
  stats <- df %>%
    group_by(z) %>%
    summarise(
      avg = mean(.data[[var]], na.rm = TRUE),
      sd  = sd  (.data[[var]], na.rm = TRUE),
      .groups = "drop"
    )
  
  # 3) panel (a): example holes
  p1 <- ggplot(df_ex, aes(x = .data[[var]], y = -z, colour = name)) +
    geom_line() +
    scale_y_continuous(trans = "reverse") +
    labs(x = var_label, y = "Depth [m]", title = "Example Holes") +
    theme_minimal() +
    theme(legend.title = element_blank())
  
  # 4) panel (b): mean ± 2 sd
  p2 <- ggplot(stats, aes(y = -z)) +
    geom_ribbon(aes(xmin = avg - 2*sd, xmax = avg + 2*sd),
                fill = "grey80", alpha = 0.5) +
    geom_line(aes(x = avg), colour = "black") +
    scale_y_continuous(trans = "reverse") +
    labs(x = var_label, y = NULL, title = "Mean ± 2 sd") +
    theme_minimal()
  
  # 5) combine
  p1 + p2 + plot_layout(ncol = 2)
}
