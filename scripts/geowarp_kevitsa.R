.libPaths("C:/Users/23478671/R_libs")
library(dplyr)
library(ggplot2)
library(geowarp)
library(patchwork)
library(caret)
library(sp)
library(gstat)
library(gridExtra)
library(tidyr)
library(plotly)
library(purrr)

set.seed(1)  

data_file <- "C:/Users/23478671/Github/Bayesian-Mixtures-Domaining/data/raw/cluster_0_data.csv"

data <- read.csv(data_file, stringsAsFactors = FALSE)

# ----- 1. Read & prepare raw data -----
raw <- read.csv(data_file, stringsAsFactors = FALSE) %>%
  rename(
    name          = Name,
    x             = X,
    y             = Y,
    z             = Z,
    target        = Cu_pct,
    #rock_quality  = RQD_Pct
  ) %>%
  filter(
    !is.na(x),
    !is.na(y),
    !is.na(z),
    !is.na(target),
    #!is.na(rock_quality),
    target >= 0.0001
  ) %>%
  select(name, x, y, z, target)%>%
  # Transform to log (log(Ni_pct))
  mutate(target = log(target)) %>%
  sample_frac(0.3)







plot_log_target_profile <- function(raw, bin_size = 50) {
  # raw: data frame with columns “name”, “x”, “y”, “z”, and “target” = log(Ni_pct)
  # bin_size: vertical interval in meters (default = 50)
  # ─────────────────────────────────────────────────────────────────────────────
  # Function 1: plot_log_target_profile()
  #   • Bins “raw” by depth (bin_size meters)
  #   • Computes per‐bin:
  #       – count (number of samples)
  #       – holes (number of unique drill‐hole IDs)
  #       – target_avg (mean log(Ni_pct))
  #       – target_sd  (sd   log(Ni_pct))
  #   • Draws a profile with:
  #       – horizontal count bars (grey)
  #       – ±1 sd error bars (steelblue)
  #       – bin‐mean points (steelblue)
  #       – sample‐count labels (above)
  #       – hole‐count labels   (below)
  #   • Inverts the y‐axis so 0 m is at top, deeper bins downwards
  # ─────────────────────────────────────────────────────────────────────────────
  
  # 1. Bin & summarize
  raw_binned <- raw %>%
    mutate(
      depth     = -z,
      depth_bin = floor(depth / bin_size) * bin_size
    ) %>%
    group_by(depth_bin) %>%
    summarise(
      count      = n(),
      holes      = n_distinct(name),
      target_avg = mean(target, na.rm = TRUE),
      target_sd  = sd(target,   na.rm = TRUE)
    ) %>%
    ungroup()
  
  # 2. Find plotting limits
  shallowest_bin <- max(raw_binned$depth_bin, na.rm = TRUE)
  deepest_bin    <- min(raw_binned$depth_bin, na.rm = TRUE)
  
  # 3. Prepare data.frame for plotting (compute bar extents & label positions)
  df_plot <- raw_binned %>%
    mutate(
      x_min = target_avg - target_sd,
      x_max = target_avg + target_sd
    ) %>%
    {
      data_range       <- max(.$x_max, na.rm = TRUE) - min(.$x_min, na.rm = TRUE)
      bar_region_width <- data_range * 0.20
      bar_start        <- min(.$x_min, na.rm = TRUE) - bar_region_width
      
      mutate(
        .,
        bar_start       = bar_start,
        bar_end         = bar_start + (count / max(count, na.rm = TRUE)) * bar_region_width,
        label_x_count   = bar_end + 0.02 * data_range,
        label_y_count   = depth_bin + (bin_size * 0.20),   # 20% of bin_size above
        label_y_holes   = depth_bin - (bin_size * 0.20)    # 20% of bin_size below
      )
    }
  
  # 4. Plot
  ggplot(df_plot, aes(x = target_avg, y = depth_bin)) +
    # a) Horizontal count bars
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
    # b) ±1 sd error bars (horizontal)
    geom_errorbarh(
      aes(
        xmin = x_min,
        xmax = x_max
      ),
      height = bin_size * 0.4,  # 40% of bin_size
      color  = "steelblue",
      alpha  = 0.5
    ) +
    # c) Bin‐mean points
    geom_point(size = 2, color = "steelblue") +
    # d) Sample‐count labels (above)
    geom_text(
      aes(
        x     = label_x_count,
        y     = label_y_count,
        label = count
      ),
      color = "grey20",
      size  = 3,
      hjust = 0
    ) +
    # e) Hole‐count labels (below)
    geom_text(
      aes(
        x     = label_x_count,
        y     = label_y_holes,
        label = holes
      ),
      color = "grey40",
      size  = 3,
      hjust = 0
    ) +
    # f) Invert y‐axis with a small padding so labels aren’t cut off
    scale_y_continuous(
      trans  = "reverse",
      limits = c(shallowest_bin + (bin_size * 0.2),
                 deepest_bin    - (bin_size * 0.2)),
      breaks = seq(shallowest_bin, deepest_bin, by = -bin_size),
      name   = "Depth (m)"
    ) +
    labs(
      title = "log(Ni_pct) Profile vs. Depth\n(with Sample‐Count above & Hole‐Count below)",
      x     = "log(Ni_pct) (bin mean)"
    ) +
    theme_minimal()
}

filter_holes_by_count <- function(raw, bin_size = 50, threshold_count = 400) {
  # raw: data frame with columns “name” and “z”
  # bin_size: vertical interval in meters (default = 50)
  # threshold_count: numeric threshold for bin counts
  # ─────────────────────────────────────────────────────────────────────────────
  # Function 2: filter_holes_by_count()
  #   • Bins “raw” by depth (bin_size meters) and counts samples per bin
  #   • Finds the shallowest bin (cutoff_bin) such that all deeper bins
  #     satisfy count > threshold_count (or adjust logic as needed)
  #   • Returns a list containing:
  #       – cutoff_bin     (numeric)
  #       – kept_bins      (vector of depth_bin ≥ cutoff_bin)
  #       – holes_to_keep  (vector of hole IDs)
  #       – filtered_holes (data.frame of all rows from those holes)
  # ─────────────────────────────────────────────────────────────────────────────
  
  # 1. Bin & count
  
  raw_binned <- raw %>%
    mutate(
      depth     = -z,
      depth_bin = floor(depth / 50) * 50
    ) %>%
    group_by(depth_bin) %>%
    summarise(
      count = n()
    ) %>%
    ungroup()
  
  # 2. Identify the “cutoff” depth_bin: the shallowest bin (largest depth_bin) 
  #    for which all deeper bins (depth_bin < current) have count < 90
  
  # First, sort bins from shallowest (0) to deepest (most negative)
  bins_sorted <- raw_binned %>% arrange(desc(depth_bin))
  
  # For each bin, check if all bins deeper (depth_bin < this) have count < 90
  cutoff_bin <- bins_sorted %>%
    filter(
      # for each row, test:
      sapply(
        depth_bin, 
        function(d) {
          all(raw_binned$count[raw_binned$depth_bin < d] > 400)
        }
      )
    ) %>%
    slice(1) %>% 
    pull(depth_bin)
  
  cutoff_bin
  #> [1] -150   # example result: –150 m
  
  # 3. Filter the binned data to only include bins from cutoff_bin upward (i.e. depth_bin >= cutoff_bin)
  #    because below cutoff (more negative) all counts < 90
  filtered_bins <- raw_binned %>%
    filter(depth_bin >= cutoff_bin) %>%
    arrange(desc(depth_bin))
  kept_bins <- raw_binned %>%
    filter(depth_bin >= cutoff_bin) %>%
    pull(depth_bin)
  
  # 3. Identify hole names that have at least one sample in those kept_bins
  holes_to_keep <- raw %>%
    # compute the same 50 m bin for each sample
    mutate(depth_bin = floor((-z) / 50) * 50) %>%
    # keep only samples in the “kept” bins
    filter(depth_bin %in% kept_bins) %>%
    pull(name) %>%
    unique()
  
  
  filtered_holes <- raw %>%
    # Compute “depth” = –z, then only retain samples where depth >= cutoff_bin
    filter((-z) <= cutoff_bin) %>%
    group_by(name) %>%
    filter(n() > 15) %>%
    ungroup()

  # Return a named list
  list(
    cutoff_bin     = cutoff_bin,
    kept_bins      = kept_bins,
    holes_to_keep  = holes_to_keep,
    filtered_holes = filtered_holes
  )
}


plot_log_target_profile(raw)

res <- filter_holes_by_count(raw, bin_size = 50, threshold_count = 400)

res$cutoff_bin       # the cutoff depth_bin
res$kept_bins       
res$holes_to_keep    
nrow(res$filtered_holes) 


filtered_holes <- res$filtered_holes
all_holes <- unique(filtered_holes$name)

# Sample X% of those hole IDs
set.seed(123) 
sampled_holes <- sample(
  all_holes,
  size = ceiling(length(all_holes) * 0.50)
)

# Keep only rows belonging to those  holes
filtered_holes <- filtered_holes %>%
  filter(name %in% sampled_holes)

length(unique(filtered_holes$name))  
nrow(filtered_holes)   


compute_pairwise_avg_distances <- function(filtered_holes) {
  nested <- filtered_holes %>%
    select(name, x, y) %>%
    group_by(name) %>%
    nest(coords = c(x, y)) %>%
    ungroup()
  
  hole_names <- nested$name
  pairs      <- combn(hole_names, 2, simplify = FALSE)
  
  map_dfr(pairs, function(p) {
    coords_i <- nested %>% filter(name == p[1]) %>% pull(coords) %>% .[[1]]
    coords_j <- nested %>% filter(name == p[2]) %>% pull(coords) %>% .[[1]]
    d2       <- outer(coords_i$x, coords_j$x, "-")^2 +
      outer(coords_i$y, coords_j$y, "-")^2
    tibble(
      hole1        = p[1],
      hole2        = p[2],
      avg_distance = mean(sqrt(d2))
    )
  })
}

avg_distances = compute_pairwise_avg_distances(filtered_holes)





x_min <- min(filtered_holes$x)
x_max <- max(filtered_holes$x)
y_min <- min(filtered_holes$y)
y_max <- max(filtered_holes$y)
z_min <- min(filtered_holes$z)
z_max <- max(filtered_holes$z)

target_min <- min(filtered_holes$target)
target_max <- max(filtered_holes$target)
#rock_quality_min = min(filtered_holes$rock_quality)
#rock_quality_max = max(filtered_holes$rock_quality)


normalised <- filtered_holes %>%
  mutate(
    x_norm      = (x - x_min) / (x_max - x_min),
    y_norm      = (y - y_min) / (y_max - y_min),
    z_norm      = (z - z_min) / (z_max - z_min),
    target_norm = (target - target_min) / (target_max - target_min),
    #rock_quality_norm = (rock_quality - rock_quality_min) / (rock_quality_max - rock_quality_min)
  ) %>%
  select(name, x = x_norm, y = y_norm, z = z_norm, target = target_norm)




unscale <- function(v, v_min, v_max) v * (v_max - v_min) + v_min

# 0) Build a true one‐row‐per‐hole table in real metres
hole_coords <- filtered_holes %>% 
  group_by(name) %>% 
  summarise(
    x_m = unscale(first(x), x_min, x_max),
    y_m = unscale(first(y), y_min, y_max)
  ) %>% 
  ungroup()

# 1) Run k-means on those hole‐level coords with k = 2
set.seed(1)
km <- kmeans(hole_coords[, c("x_m","y_m")], centers = 2)

hole_coords <- hole_coords %>% 
  mutate(cluster = km$cluster)

# 2) Count holes per cluster and fraction
cluster_sizes <- hole_coords %>% 
  count(cluster) %>% 
  mutate(frac = n / sum(n))

# 3) Pick as TEST the cluster whose size is closest to 30%
test_cluster <- cluster_sizes %>% 
  slice_min(abs(frac - 0.30)) %>% 
  pull(cluster)

test_holes  <- hole_coords %>% 
  filter(cluster == test_cluster) %>% 
  pull(name)
train_holes <- hole_coords %>% 
  filter(cluster != test_cluster) %>% 
  pull(name)

# # 4) (Optional) if you want *exactly* 30% of holes in TEST, subsample:
# N      <- nrow(hole_coords)
# n_test <- round(0.30 * N)
# if (length(test_holes) > n_test) {
#   set.seed(1)
#   test_holes  <- sample(test_holes, n_test)
#   train_holes <- setdiff(hole_coords$name, test_holes)
# }

# Plug into your workflow
train_data <- normalised %>% filter(name %in% train_holes)
test_data  <- normalised %>% filter(name %in% test_holes)

cat("#train:", length(train_holes), "\n#test:", length(test_holes), "\n")




# 1) Tag each sample with its split
plot_df <- filtered_holes %>%
  mutate(split = ifelse(name %in% train_holes, "Train", "Test"))

# 2) Define two contrasting colours
split_cols <- c(Train = "blue", Test = "red")

# 3) 3D scatter, coloured by split
plot_ly(
  data  = plot_df,
  x     = ~x,            # normalised X
  y     = ~y,            # normalised Y
  z     = ~z,            # normalised Z
  color = ~split,        # “Train” vs “Test”
  colors = split_cols,
  type   = "scatter3d",
  mode   = "markers",
  marker = list(size = 3, opacity = 0.8)
) %>%
  layout(
    title = "Train vs Test Holes (normalised x,y,z)",
    legend = list(title = list(text = "Split")),
    scene = list(
      xaxis = list(title = "X"),
      yaxis = list(title = "Y"),
      zaxis = list(title = "Z (depth)")
    )
  )






# train_idx  <- createDataPartition(normalised$target, p = 0.7, list = FALSE)
# train_data  <- normalised[train_idx, ]
# test_data   <- normalised[-train_idx, ]



# 
# # HOLE‐BASED
# all_holes   <- unique(normalised$name)
# set.seed(1)
# train_holes <- sample(all_holes, ceiling(length(all_holes) * 0.7))
# test_holes  <- setdiff(sampled_holes, train_holes)
# 
# train_data   <- normalised %>% filter(name %in% train_holes)
# test_data   <- normalised %>% filter((name %in% test_holes))
# 

# --- Deviation: short-range exp + flexible nugget -------------------
dev_mod <- geowarp_vertical_only_deviation_model(
  covariance_function = "matern15",
  axial_warping_unit  = geowarp_bernstein_awu(order = 30), # 
  variance_model      = geowarp_variance_model(
    vertical_basis_functions      = TRUE,
    vertical_basis_function_delta = 10/(z_max-z_min)   # meters
  )
)

horizontal_scaling_x <- x_max - x_min
horizontal_scaling_y <- y_max - y_min




mean_mod <- geowarp_mean_model(
  fixed_formula                = ~ z, #consider not including cov if want to interp
  vertical_basis_functions     = TRUE,
  vertical_basis_function_delta = 1000/(z_max-z_min)   #(this is in meters)
)



dev_mod_full3D <- geowarp_deviation_model(
  covariance_function = "matern15",
  
  axial_warping_units = list(
    geowarp_linear_awu(
      scaling = 1,
      prior = list(
        type  = "inv_uniform",
        shape = 0, 
        rate  = 0,
        lower = 1,
        upper = 90
      )
    ),
    geowarp_linear_awu(
      scaling = 1,
      prior = list(
        type  = "inv_uniform",
        shape = 0, 
        rate  = 0,
        lower = 1,
        upper = 90
      )
    ),
    
    geowarp_bernstein_awu(order = 30)
  ),
  
  geometric_warping_unit = geowarp_geometric_warping_unit(),
  
  variance_model = geowarp_variance_model(
    vertical_basis_functions      = FALSE,
    vertical_basis_function_delta = 1000/(z_max-z_min)   
  )
)



model <- geowarp_model(
  variable               = "target",
  horizontal_coordinates = c("x", "y"),
  horizontal_domains     = list(range(train_data$x), range(train_data$y)),
  vertical_coordinate    = "z",
  #vertical_domain        = range(train_data$z),
  vertical_domain        = c(0, 1), 
  mean_model             = mean_mod,
  deviation_model        = dev_mod_full3D
  #deviation_model        = dev_mod
)

set.seed(1)  
fit <- geowarp_optimise(
  train_data, model,
  n_parents = 32,   
  best_of   = 2,    # x random restarts
  trace     = 3,
  threads   = 8
)




test_data <- test_data %>%
  mutate(
    predicted_mean   = mean_profile(fit, df = .),
    predicted_stdev  = sqrt(marginal_variance_profile(fit, df = .))
  )


errors_gw <- test_data %>%
  mutate(
    residual       = target - predicted_mean,
    abs_error      = abs(residual),
    squared_error  = residual^2,
    within_1sd     = abs(residual) <= predicted_stdev,
    within_2sd     = abs(residual) <= (2 * predicted_stdev)
  )

mae_gw         <- mean(errors_gw$abs_error, na.rm = TRUE)
rmse_gw        <- sqrt(mean(errors_gw$squared_error, na.rm = TRUE))
median_ae_gw   <- median(errors_gw$abs_error, na.rm = TRUE)
r_squared_gw   <- 1 - (sum(errors_gw$squared_error, na.rm = TRUE) /
                         sum((test_data$target - mean(test_data$target, na.rm = TRUE))^2, na.rm = TRUE))
pct_1sd_gw     <- mean(errors_gw$within_1sd, na.rm = TRUE) * 100
pct_2sd_gw     <- mean(errors_gw$within_2sd, na.rm = TRUE) * 100







#Kriging

train_sp <- train_data
coordinates(train_sp) <- ~ x + y + z
proj4string(train_sp) <- CRS("+proj=longlat +datum=WGS84")  # dummy CRS; only relative distances matter

test_sp <- test_data
coordinates(test_sp) <- ~ x + y + z
proj4string(test_sp) <- CRS("+proj=longlat +datum=WGS84")


vgm_cloud <- variogram(target ~ 1, data = train_sp) 
vgm_fit   <- fit.variogram(vgm_cloud, vgm("Sph"), fit.method = 2)


krig_res <- krige(
  formula   = target ~ 1,
  locations = train_sp,
  newdata   = test_sp,
  model     = vgm_fit,
  nmax      = 30,           # <- only use the 30 nearest samples
  debug.level = 0           # <- suppress extra console chatter
)

test_data$krig_pred  <- krig_res$var1.pred
test_data$krig_var   <- krig_res$var1.var
test_data$krig_sd    <- sqrt(test_data$krig_var)


errors_kr <- test_data %>%
  mutate(
    residual_kr       = target - krig_pred,
    abs_error_kr      = abs(residual_kr),
    squared_error_kr  = residual_kr^2,
    within_1sd_kr     = abs(residual_kr) <= krig_sd,
    within_2sd_kr     = abs(residual_kr) <= (2 * krig_sd)
  )

mae_kr       <- mean(errors_kr$abs_error_kr, na.rm = TRUE)
rmse_kr      <- sqrt(mean(errors_kr$squared_error_kr, na.rm = TRUE))
median_ae_kr <- median(errors_kr$abs_error_kr, na.rm = TRUE)
r_squared_kr <- 1 - (sum(errors_kr$squared_error_kr, na.rm = TRUE) /
                       sum((test_data$target - mean(test_data$target, na.rm = TRUE))^2, na.rm = TRUE))
pct_1sd_kr   <- mean(errors_kr$within_1sd_kr, na.rm = TRUE) * 100
pct_2sd_kr   <- mean(errors_kr$within_2sd_kr, na.rm = TRUE) * 100

#IDW

idw_res <- idw(
  formula   = target ~ 1,
  locations = train_sp,    
  newdata   = test_sp,
  idp       = 3            
)

test_data$idw_pred <- idw_res$var1.pred
test_data$idw_sd   <- sqrt(idw_res$var1.var)   

errors_idw <- test_data %>%
  mutate(
    residual_idw       = target - idw_pred,
    abs_error_idw      = abs(residual_idw),
    squared_error_idw  = residual_idw^2,
    within_1sd_idw     = abs(residual_idw) <= idw_sd,
    within_2sd_idw     = abs(residual_idw) <= (2 * idw_sd)
  )

mae_idw       <- mean(errors_idw$abs_error_idw, na.rm = TRUE)
rmse_idw      <- sqrt(mean(errors_idw$squared_error_idw, na.rm = TRUE))
median_ae_idw <- median(errors_idw$abs_error_idw, na.rm = TRUE)
r_squared_idw <- 1 - (sum(errors_idw$squared_error_idw, na.rm = TRUE) /
                        sum((test_data$target - mean(test_data$target, na.rm = TRUE))^2, na.rm = TRUE))
pct_1sd_idw   <- mean(errors_idw$within_1sd_idw, na.rm = TRUE) * 100
pct_2sd_idw   <- mean(errors_idw$within_2sd_idw, na.rm = TRUE) * 100

### ERROR SUMMARY

validation_summary <- data.frame(
  metric          = c("MAE","RMSE","Median_AE","R_squared",
                      "Pct_within_1SD","Pct_within_2SD"),
  GeoWarp         = c(mae_gw, rmse_gw, median_ae_gw, r_squared_gw,
                      pct_1sd_gw, pct_2sd_gw),
  Kriging         = c(mae_kr, rmse_kr, median_ae_kr, r_squared_kr,
                      pct_1sd_kr, pct_2sd_kr),
  IDW             = c(mae_idw, rmse_idw, median_ae_idw, r_squared_idw,
                      pct_1sd_idw, pct_2sd_idw)
)

print(validation_summary)













unscale <- function(v, v_min, v_max) v * (v_max - v_min) + v_min


test_data <- test_data %>%
  mutate(
    # depth
    z_raw = unscale(z, z_min, z_max),
    
    # observed target
    target_raw = unscale(target, target_min, target_max),
    
    # model predictions
    predicted_mean_raw = unscale(predicted_mean, target_min, target_max),
    krig_pred_raw      = unscale(krig_pred,      target_min, target_max),
    idw_pred_raw       = unscale(idw_pred,       target_min, target_max)
  )


hole_counts <- test_data %>%
  count(name, name = "n_obs") %>%
  arrange(desc(n_obs))

top_holes <- hole_counts %>%
  slice(1:20) %>%                 
  pull(name)

pred_long_top <- test_data %>%
  filter(name %in% top_holes) %>%
  select(name, z_raw, target_raw,
         predicted_mean_raw, krig_pred_raw, idw_pred_raw) %>%
  pivot_longer(
    cols      = c(predicted_mean_raw, krig_pred_raw, idw_pred_raw),
    names_to  = "method",
    values_to = "prediction"
  ) %>%
  mutate(
    method = recode(
      method,
      predicted_mean_raw = "GeoWarp",
      krig_pred_raw      = "Kriging",
      idw_pred_raw       = "IDW"
    ),
    depth = -z_raw        # flip so deeper is lower on the plot
  )

obs_top <- test_data %>%
  filter(name %in% top_holes) %>%
  transmute(
    name,
    depth  = -z_raw,
    target = target_raw
  )


nn_info <- bind_rows(
  avg_distances %>% rename(hole = hole1, neighbour = hole2),
  avg_distances %>% rename(hole = hole2, neighbour = hole1)
) %>%
  group_by(hole) %>%
  slice_min(avg_distance, n = 1) %>%
  ungroup()

label_map <- nn_info %>%
  filter(hole %in% top_holes) %>%
  mutate(label = paste0(
    hole,
    "\nNN: ", neighbour,
    " (", round(avg_distance, 1), " m)"
  )) %>%
  {setNames(.$label, .$hole)}


## ----------------  Plot on the original (un-scaled) units  ----------------
ggplot() +
  ## a) GeoWarp / Kriging / IDW predictions – colour legend
  geom_point(
    data  = pred_long_top,
    aes(x = prediction, y = depth, colour = method),
    shape = 16, size = 1.3
  ) +
  
  ## b) Observed data – shape legend
  geom_point(
    data  = obs_top,
    aes(x = target,    y = depth, shape = "Observed"),
    colour = "black", size = 0.8
  ) +
  
  ## colour scale (methods)
  scale_colour_manual(
    name   = "Method",
    values = c(GeoWarp = "#F8766D",
               Kriging = "#00BA38",
               IDW     = "#619CFF")
  ) +
  
  ## shape scale (observed)
  scale_shape_manual(
    name   = "",
    values = c(Observed = 4),
    labels = c(Observed = "Observed data")
  ) +
  
  guides(
    colour = guide_legend(order = 1),
    shape  = guide_legend(order = 2, override.aes = list(colour = "black"))
  ) +
  
  facet_wrap(~ name,
             scales   = "free_x",
             ncol     = 5,
             nrow     = 5,
             labeller = as_labeller(label_map)) +
  
  scale_y_continuous(trans = "reverse", name = "Depth (m)") +
  labs(
    title = "Predictions vs. Observed Ni_pct by Depth (Top 20 Boreholes)",
    x     = "log(Ni_pct)"
  ) +
  theme_minimal() +
  theme(
    strip.text      = element_text(size = 8),
    axis.title      = element_text(size = 10),
    axis.text       = element_text(size = 7),
    legend.position = "bottom",
    panel.spacing   = unit(0.5, "lines")
  )









# ─────────────────────────────────────────────────────────────────────────────
# # 1. define bin width to match your profile step
# bin_width <- 0.005
# 
# # 2. compute observed mean & sd in each bin of test_data
# stats_binned <- test_data %>%
#   # assign each test_data row to the lower edge of its bin
#   mutate(z_bin = floor(z / bin_width) * bin_width) %>%
#   group_by(z_bin) %>%
#   summarise(
#     obs_mean = mean(target, na.rm = TRUE),
#     obs_sd   = sd(target,   na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# # 3. join those binned stats onto profile_df
# profile_df <- profile_df %>%
#   # make sure profile_df has the same z_bin column
#   mutate(z_bin = floor(z / bin_width) * bin_width) %>%
#   left_join(stats_binned, by = "z_bin") %>%
#   select(-z_bin)



## 1.  Choose an (x,y) location in the scaled space
x0 <- mean(train_data$x)           # ≈ 0.5 – the centre of your scaled X
y0 <- mean(train_data$y)           # ≈ 0.5 – the centre of your scaled Y

## 2.  Create a fine Z sequence (still in the scaled space)
profile_df <- data.frame(
  x = x0,
  y = y0,
  z = unique(test_data$z)         # every ~1 % of the depth range
) %>%
  mutate(
    gwarp_mean  = mean_profile(fit, df = .),
    gwarp_stdev = sqrt(marginal_variance_profile(fit, df = .))
)

stats_by_z <- test_data %>%
  group_by(z) %>%
  summarise(
    obs_mean = mean(target, na.rm = TRUE),
    obs_sd   = sd(target,   na.rm = TRUE),
    .groups = "drop"
  )


profile_df <- profile_df %>%
  left_join(stats_by_z, by = "z")
  
wrap_plots(
  ggplot(profile_df) +
    #geom_line(aes(z, obs_mean), col = 'red') +
    geom_line(aes(z, gwarp_mean), col = 'blue') +
    scale_x_reverse() +
    coord_flip() +
    labs(x = 'Depth', y = 'Mean', title = 'Mean profile'),
  ggplot(profile_df) +
    geom_line(aes(z, obs_sd), col = 'red') +
    geom_line(aes(z, gwarp_stdev), col = 'blue') +
    scale_x_reverse() +
    coord_flip() +
    labs(x = 'Depth', y = 'Standard deviation', title = 'Standard deviation'),
  nrow = 1
)




warping_df <- data.frame(
  x = 0,
  y = 0,
  z = test_data$z
) %>%
  mutate(
    z_warped = warped_coordinates(fit, df = .)[, 3]
)

# This should look roughly linear, since the data are actually stationary.
# The way this works is that the data are stationary with a length scale of
# one in the warped space
ggplot(warping_df, aes(z, z_warped)) +
  geom_line() +
  scale_x_reverse() +
  coord_flip() +
  labs(x = 'Depth', y = 'Warped depth', title = 'Warping function')














ggplot() +
  geom_line(
    data = test_data,
    mapping = aes(z, target, group = name),
    alpha = 0.1
  ) +
  geom_line(
    data = profile_df,
    mapping = aes(z, gwarp_mean),
    col = 'blue'
  ) +
  geom_line(
    data = profile_df,
    mapping = aes(z, gwarp_mean - 2 * gwarp_stdev),
    col = 'blue',
    lty = 2
  ) +
  geom_line(
    data = profile_df,
    mapping = aes(z, gwarp_mean + 2 * gwarp_stdev),
    col = 'blue',
    lty = 2
  ) +
  scale_x_reverse() +
  coord_flip() +
  labs(x = 'Depth', y = 'Target', title = 'Target profile 3D deviation')
















library(dplyr)
library(geowarp)
library(progress)

# 1) The worker: fit one model, record its run‐time, return a 1-row tibble
run_one <- function(lx_min, lx_max, ly_min, ly_max,
                    scaling_x, scaling_y,
                    train_data, test_data)
{
  t0 <- Sys.time()
  
  dev_mod <- geowarp_deviation_model(
    covariance_function = "matern15",
    axial_warping_units = list(
      geowarp_linear_awu(
        scaling = 1,
        prior = list(
          type  = "inv_uniform", shape = 0, rate = 0,
          lower = scaling_x/(lx_max),
          upper = scaling_x/(lx_min)
        )
      ),
      geowarp_linear_awu(
        scaling = 1,
        prior = list(
          type  = "inv_uniform", shape = 0, rate = 0,
          lower = scaling_y/(ly_max),
          upper = scaling_y/(ly_min)
        )
      ),
      geowarp_bernstein_awu(order = 30)
    ),
    geometric_warping_unit = NULL,
    variance_model = geowarp_variance_model(
      vertical_basis_functions      = TRUE,
      vertical_basis_function_delta = 2
    )
  )
  
  model <- geowarp_model(
    variable               = "target",
    horizontal_coordinates = c("x","y"),
    horizontal_domains     = list(range(train_data$x),
                                  range(train_data$y)),
    vertical_coordinate    = "z",
    vertical_domain        = c(0,1),
    mean_model             = mean_mod,       # assume you defined mean_mod earlier
    deviation_model        = dev_mod
  )
  
  fit <- geowarp_optimise(train_data, model,
                          n_parents = 32, best_of = 2,
                          trace     = 0, threads = 4)
  
  preds <- test_data %>%
    mutate(
      mu  = mean_profile(fit, df = .),
      sd  = sqrt(marginal_variance_profile(fit, df = .)),
      res = target - mu
    )
  
  # now assign out BEFORE touching it:
  out <- with(preds, data.frame(
    lx_min, lx_max, ly_min, ly_max,
    RMSE    = sqrt(mean(res^2)),
    MAE     = mean(abs(res))
  ))
  
  # append run time in seconds
  out$run_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  out
}

# 2) A grid‐fitting loop that shows progress + ETA
fit_grid <- function(grid_df, scaling_x, scaling_y,
                     train_data, test_data,
                     stage = "coarse")
{
  n  <- nrow(grid_df)
  pb <- progress_bar$new(
    format = sprintf("%s :current/:total (:percent) | elapsed: :elapsedfull  eta: :eta",
                     toupper(stage)),
    total = n, clear = FALSE, width = 80
  )
  
  res <- vector("list", n)
  for(i in seq_len(n)) {
    res[[i]] <- with(grid_df[i, ],
                     run_one(lx_min, lx_max, ly_min, ly_max,
                             scaling_x, scaling_y,
                             train_data, test_data))
    pb$tick()
  }
  bind_rows(res)
}

# 3) Example usage:
scaling_x <- diff(c(x_min,x_max))
scaling_y <- diff(c(y_min,y_max))


min_vals <- c(1, 7)          # will map to lx_min / ly_min
max_vals <- c(100, 500) 

# coarse grid:
grid1 <- expand.grid(
  lx_max = max_vals,
  lx_min = min_vals,
  ly_max = max_vals,
  ly_min = min_vals
) %>% filter(lx_min < lx_max, ly_min < ly_max)

results1 <- fit_grid(grid1, scaling_x, scaling_y,
                     train_data, test_data, stage = "coarse")

# inspect
head(arrange(results1, RMSE), 10)

best1 <- results1 %>% arrange(RMSE) %>% slice(1:10) 
print(best1)










## ---------- 2)  fine search around each top candidate ------------------- ##
## build a bespoke grid per candidate so you zoom in *relative* to each best1 row
make_fine_grid <- function(row, d_min = 0.05, d_max = 0.2, n_min = 5, n_max = 4) {
  expand.grid(
    lx_min = seq(pmax(row$lx_min - d_min, 0.01),
                 row$lx_min + d_min,
                 length.out = n_min),
    lx_max = seq(row$lx_max - d_max,
                 row$lx_max + d_max,
                 length.out = n_max),
    ly_min = seq(pmax(row$ly_min - d_min, 0.01),
                 row$ly_min + d_min,
                 length.out = n_min),
    ly_max = seq(row$ly_max - d_max,
                 row$ly_max + d_max,
                 length.out = n_max)
  ) %>% filter(lx_min < lx_max, ly_min < ly_max)
}

grid2 <- bind_rows(lapply(seq_len(nrow(best1)), function(i) make_fine_grid(best1[i, ]))) %>%
  distinct()

results2 <- do.call(
  rbind,
  lapply(seq_len(nrow(grid2)), function(i) {
    with(grid2[i, ],
         run_one(lx_min, lx_max, ly_min, ly_max,
                 scaling_x, scaling_y,
                 train_data, test_data))
  })
)

best2 <- results2 %>% arrange(RMSE) %>% slice(1:10)
print(best2)

## choose final best across both stages
final_best <- bind_rows(best1, best2) %>%
  arrange(RMSE) %>% slice(1)
print(final_best)


























hole_coords <- normalised %>%
  group_by(name) %>%
  summarise(
    collar_x = first(x),          # already 0–1
    collar_y = first(y),
    collar_z = first(z),
    Dip      = if ("Dip" %in% names(raw)) first(raw$Dip[raw$name == cur_group()$name]) else NA_real_,
    .groups  = "drop"
  )

# -------------------------------------------------------------------------
# 2.  Collar-density score via k-NN in 3-D (normalised space) ------------
# -------------------------------------------------------------------------
k_neigh <- 10
coords  <- as.matrix(hole_coords[, c("collar_x", "collar_y", "collar_z")])
knn_d   <- FNN::knn.dist(coords, k = k_neigh)[, k_neigh]   # distance to kth neighbour

hole_coords <- hole_coords %>%
  mutate(
    density_score   = 1 / knn_d,             # higher = denser drilling
    density_stratum = ntile(density_score, 4)   # quartiles
  )

# -------------------------------------------------------------------------
# 3.  Dip classes (if Dip exists) -----------------------------------------
# -------------------------------------------------------------------------
if (!all(is.na(hole_coords$Dip))) {
  hole_coords <- hole_coords %>%
    mutate(
      dip_class = cut(
        abs(Dip),
        breaks = c(-Inf, 30, 60, Inf),       # shallow / moderate / steep
        labels = c("shallow", "moderate", "steep"),
        right  = FALSE
      )
    )
} else {
  hole_coords <- hole_coords %>% mutate(dip_class = factor("unknown"))
}

# -------------------------------------------------------------------------
# 4.  Stratum ID and 20 % hole-level Test split ---------------------------
# -------------------------------------------------------------------------
hole_coords <- hole_coords %>%
  mutate(stratum_id = interaction(density_stratum, dip_class, drop = TRUE))

test_frac <- 0.20
set.seed(1)

test_holes <- hole_coords %>%
  group_by(stratum_id) %>%
  sample_frac(test_frac) %>%
  pull(name)

train_holes <- setdiff(hole_coords$name, test_holes)

train_data <- normalised %>% filter(name %in% train_holes)
test_data  <- normalised %>% filter(name %in% test_holes)

cat("# holes  (train):", length(train_holes),
    "\n# holes  (test) :", length(test_holes),
    "\n# samples(train):", nrow(train_data),
    "\n# samples(test) :", nrow(test_data), "\n")

# -------------------------------------------------------------------------
# 5.  3-D scatter coloured by Train / Test (unchanged styling) ------------
# -------------------------------------------------------------------------
plot_df <- bind_rows(
  train_data %>% mutate(split = "Train"),
  test_data  %>% mutate(split = "Test")
)

split_cols <- c(Train = "blue", Test = "red")

plot_ly(
  data   = plot_df,
  x      = ~x,
  y      = ~y,
  z      = ~z,
  color  = ~split,
  colors = split_cols,
  type   = "scatter3d",
  mode   = "markers",
  marker = list(size = 3, opacity = 0.8)
) %>%
  layout(
    title  = "Train vs Test Holes (normalised X, Y, Z) – density × dip strata",
    legend = list(title = list(text = "Split")),
    scene  = list(
      xaxis = list(title = "X"),
      yaxis = list(title = "Y"),
      zaxis = list(title = "Z (depth)")
    )
  )








