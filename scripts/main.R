# main.R -------------------------------------------------------------------
setwd("C:/Users/23478671/Github/GeoWarp-Application-Kevitsa")
source("helper_functions/00_utils_io.R")
source("helper_functions/01_preprocess.R")
source("helper_functions/02_split.R")
source("helper_functions/03_models.R")
source("helper_functions/04_visuals.R")



set.seed(1)

csv  <- "cluster_0_data.csv"


# 1. load & pre-process
set.seed(1)
#raw   <- load_raw(csv, sample_frac = 1)
raw <- load_raw(csv, var = "Cu_pct", sample_frac = 1)

depth_profile_plot(raw)

clean <- depth_filter(raw, min_cnt = 800)
normd <- normalise_xyz_target(clean)
z_rng <- attr(normd, "ranges")$z

# 2. choose a split (swap gap_split for strata_split)
set.seed(1)
spl <- km_split(normd, frac = 0.30, centres = 2)
#spl <- strata_split(normd)

train <- spl$train
test  <- spl$test

ranges <- attr(normd, "ranges")    # saved by normalise_xyz_target()

plot_df <- bind_rows(
  train %>% mutate(split = "Train"),
  test  %>% mutate(split = "Test")
)
train_test_scatter3D(plot_df, ranges = ranges)

#train_test_scatter3D(plot_df)
# p_collars_2d <- collar_xy_plot(plot_df)
# print(p_collars_2d)



# 3. fit models
# (a) Full warp, no vertical bases
gw_fit_cv        <- fit_geowarp(train, variable_variance = FALSE)
# (b) Full warp, with vertical bases
#gw_vert_fit        <- fit_geowarp(train, variable_variance = FALSE)

# (c) No warp, no vertical bases
gw_fit_nw_cv <- fit_geowarp_noWarp(train, variable_variance = FALSE)
# (d) No warp, with vertical bases
#gw_noWarp_vert_fit <- fit_geowarp_noWarp(train, variable_variance = FALSE)

kr_pred  <- fit_kriging(train, test)
idw_pred <- fit_idw(train, test)

# 4. attach predictions & evaluate
# test <- test %>% mutate(
#   gw_base_pred        = mean_profile(gw_base_fit,        df = .),
#   gw_base_sd          = sqrt(marginal_variance_profile(gw_base_fit,        df = .)),
#   gw_vert_pred        = mean_profile(gw_vert_fit,        df = .),
#   gw_vert_sd          = sqrt(marginal_variance_profile(gw_vert_fit,        df = .)),
#   gw_noWarp_base_pred = mean_profile(gw_noWarp_base_fit, df = .),
#   gw_noWarp_base_sd   = sqrt(marginal_variance_profile(gw_noWarp_base_fit, df = .)),
#   gw_noWarp_vert_pred = mean_profile(gw_noWarp_vert_fit, df = .),
#   gw_noWarp_vert_sd   = sqrt(marginal_variance_profile(gw_noWarp_vert_fit, df = .)),
#   kr_mean             = kr_pred$var1.pred,
#   kr_sd               = sqrt(kr_pred$var1.var),
#   idw_mean            = idw_pred$var1.pred,
#   idw_sd              = sqrt(idw_pred$var1.var)
# )

test <- test %>% mutate(
  gw_cv_pred        = mean_profile(gw_fit_cv,        df = .),
  gw_cv_sd          = sqrt(marginal_variance_profile(gw_fit_cv,        df = .)),
  gw_noWarp_cv_pred = mean_profile(gw_fit_nw_cv, df = .),
  gw_noWarp_cv_sd   = sqrt(marginal_variance_profile(gw_fit_nw_cv, df = .)),
  kr_mean             = kr_pred$var1.pred,
  kr_sd               = sqrt(kr_pred$var1.var),
  idw_mean            = idw_pred$var1.pred,
  idw_sd              = sqrt(idw_pred$var1.var)
)



eval_metric <- function(obs, pred, pred_sd = NULL) {
  
  # ---- core point-prediction metrics (4 decimals) -------------------------
  mae  <- round(mean(abs(obs - pred), na.rm = TRUE), 4)
  rmse <- round(sqrt(mean((obs - pred)^2, na.rm = TRUE)), 4)
  
  stats <- c(MAE = mae, RMSE = rmse)
  
  # ---- calibration metrics (2-decimal percentages) ------------------------
  if (!is.null(pred_sd)) {
    pct_1sd <- round(
      mean(abs(obs - pred) <= pred_sd,      na.rm = TRUE) * 100, 2)
    pct_2sd <- round(
      mean(abs(obs - pred) <= 2 * pred_sd,  na.rm = TRUE) * 100, 2)
    
    stats <- c(stats,
               Pct_within_1SD = pct_1sd,
               Pct_within_2SD = pct_2sd)
  }
  
  stats
}


validation <- rbind(
  GW_CV        = eval_metric(test$target, test$gw_cv_pred,        test$gw_cv_sd),
#  GW_Vert        = eval_metric(test$target, test$gw_vert_pred,        test$gw_vert_sd),
  GW_NoWarp_CV = eval_metric(test$target, test$gw_noWarp_cv_pred, test$gw_noWarp_cv_sd),
#  GW_NoWarp_Vert = eval_metric(test$target, test$gw_noWarp_vert_pred, test$gw_noWarp_vert_sd),
  Kriging        = eval_metric(test$target, test$kr_mean,             test$kr_sd),
  IDW            = eval_metric(test$target, test$idw_mean,            test$idw_sd)
)

# grab the target name
target_name <- attr(raw, "target_name")

# collect your fits in a named list
model_list <- list(
  gw_fit_cv        = gw_fit_cv,
  gw_fit_noWarp_cv = gw_fit_nw_cv
)

# loop over them and save each
for (nm in names(model_list)) {
  fit_obj <- model_list[[nm]]
  # build a file name like "gw_base_Cu_pct.rds"
  fname   <- sprintf("%s_%s.rds", nm, target_name)
  saveRDS(fit_obj, file = fname)
  message("Saved ", nm, " to ", fname)
}


print(validation)


unscale_vec <- function(v, rng) v * (rng$max - rng$min) + rng$min
# test_data <- test %>% mutate(
#   z_raw              = unscale_vec(z,  attr(normd, "ranges")$z),
#   target_raw         = unscale_vec(target, attr(normd, "ranges")$target),
#   predicted_mean_raw = unscale_vec(gw_base_pred, attr(normd, "ranges")$target),
#   gw_noWarp_pred_raw = unscale_vec(gw_noWarp_base_pred, attr(normd, "ranges")$target),
#   krig_pred_raw      = unscale_vec(kr_mean,       attr(normd, "ranges")$target),
#   idw_pred_raw       = unscale_vec(idw_mean,      attr(normd, "ranges")$target)
# )

test_data <- test %>% mutate(
  z_raw              = unscale_vec(z,  attr(normd, "ranges")$z),
  target_raw         = unscale_vec(target, attr(normd, "ranges")$target),
  predicted_mean_raw = unscale_vec(gw_cv_pred, attr(normd, "ranges")$target),
  krig_pred_raw      = unscale_vec(kr_mean,       attr(normd, "ranges")$target),
  idw_pred_raw       = unscale_vec(idw_mean,      attr(normd, "ranges")$target)
)

top_hole_prediction_plot(
  test_data,
  avg_distances   = NULL,     # or NULL
  n_holes         = 20,
  model_cols      = c("predicted_mean_raw",
                      "krig_pred_raw",
                      "idw_pred_raw"),
  model_names     = c("GeoWarp", "Kriging", "IDW")
)


# static 3-panel plot
profile_and_warp_plots(gw_fit_cv, test)    # uses centre of train-data by default


band_plot <- target_profile_bandplot(gw_fit_cv, test)
print(band_plot)                    # or ggsave("band_plot.pdf", p_band, width = 5, height = 6)






































# -------------------------------------------------------------------------
# master driver – run the full workflow for several target variables
# -------------------------------------------------------------------------
targets <- c("Cu_pct", "Au_ppb", "Ni_pct", "Ag_ppm")

csv_path <- "cluster_0_data.csv"

all_validation <- list()   # collect the per-target metrics

for (tvar in targets) {
  cat("\n==================  modelling ", tvar, "  ==================\n")
  
  ## ----------------------------------------------------------------------
  ## 1. load  &  depth-filter
  ## ----------------------------------------------------------------------
  set.seed(1)                                        # reproducible subsampling
  raw   <- load_raw(csv_path, var = tvar, sample_frac = 1)
  
  clean <- depth_filter(raw, min_cnt = 800)
  if (nrow(clean) == 0) {
    warning("No rows left after depth_filter() for ", tvar, "; skipping.")
    next
  }
  
  normd <- normalise_xyz_target(clean)
  z_rng <- attr(normd, "ranges")$z
  
  ## ----------------------------------------------------------------------
  ## 2. split into train / test
  ## ----------------------------------------------------------------------
  set.seed(1)
  spl <- km_split(normd, frac = 0.30, centres = 2)   # or strata_split(normd)
  
  train <- spl$train
  test  <- spl$test
  
  ## ----------------------------------------------------------------------
  ## 3. fit the two GeoWarp variants + Kriging + IDW
  ## ----------------------------------------------------------------------
  gw_fit      <- fit_geowarp(train,        variable_variance = TRUE)
  gw_fit_nw   <- fit_geowarp_noWarp(train, variable_variance = TRUE)
  
  kr_pred  <- fit_kriging(train, test)
  idw_pred <- fit_idw(train,   test)
  
  ## ----------------------------------------------------------------------
  ## 4. attach predictions and compute metrics
  ## ----------------------------------------------------------------------
  test <- test %>% mutate(
    gw_pred        = mean_profile(gw_fit,        df = .),
    gw_sd          = sqrt(marginal_variance_profile(gw_fit,        df = .)),
    gw_noWarp_pred = mean_profile(gw_fit_nw,     df = .),
    gw_noWarp_sd   = sqrt(marginal_variance_profile(gw_fit_nw,     df = .)),
    kr_mean           = kr_pred$var1.pred,
    kr_sd             = sqrt(kr_pred$var1.var),
    idw_mean          = idw_pred$var1.pred,
    idw_sd            = sqrt(idw_pred$var1.var)
  )
  
  ## ---- helper to compute MAE / RMSE / calibration -----------------------
  do_metrics <- function(obs, pred, sd = NULL) {
    mae  <- round(mean(abs(obs - pred), na.rm = TRUE), 5)
    rmse <- round(sqrt(mean((obs - pred)^2, na.rm = TRUE)), 5)
    out  <- c(MAE = mae, RMSE = rmse)
    if (!is.null(sd)) {
      out <- c(out,
               Pct_within_1SD = round(mean(abs(obs - pred) <= 1*sd, na.rm = TRUE)*100, 2),
               Pct_within_2SD = round(mean(abs(obs - pred) <= 2*sd, na.rm = TRUE)*100, 2)
      )
    }
    out
  }
  
  validation <- rbind(
    GW_CV        = do_metrics(test$target, test$gw_pred,        test$gw_sd),
    GW_NoWarp_CV = do_metrics(test$target, test$gw_noWarp_pred, test$gw_noWarp_sd),
    Kriging      = do_metrics(test$target, test$kr_mean,           test$kr_sd),
    IDW          = do_metrics(test$target, test$idw_mean,          test$idw_sd)
  )
  
  print(validation)


  all_validation[[tvar]] <- validation         # store for later
  
  ## ----------------------------------------------------------------------
  ## 5. save test set + fitted models
  ## ----------------------------------------------------------------------
  # (a) the test set
  test_fname <- sprintf("test_data_%s.csv", tvar)
  write.csv(test, test_fname, row.names = FALSE)
  cat("  ↳ wrote test set  ➜ ", test_fname, "\n")
  
  # (b) fitted GeoWarp models
  saveRDS(gw_fit,      file = sprintf("gw_fit_%s.rds",      tvar))
  saveRDS(gw_fit_nw,   file = sprintf("gw_fit_noWarp_%s.rds", tvar))
  cat("  ↳ saved model RDS files for ", tvar, "\n")
}

## ------------------------------------------------------------------------
## 6. write the big validation summary
## ------------------------------------------------------------------------
val_df <- dplyr::bind_rows(lapply(names(all_validation), function(v)
  cbind(Target = v,
        Method = rownames(all_validation[[v]]),
        all_validation[[v]])))
write.csv(val_df, "validation_summary_all_targets.csv", row.names = FALSE)
cat("\nFinished.  All metrics in validation_summary_all_targets.csv\n")














# ────────────────────────────────────────────────────────────────
# 0.  Saving train/test RDS
# ────────────────────────────────────────────────────────────────

# ── 0.  House–keeping ───────────────────────────────────────────────────
setwd("C:/Users/23478671/Github/Geowarp-Resource-Estimation/GeoWarp_resource_estimation/application_kevitsa")

# library(dplyr) ; library(purrr)
# 
# source("comparison_functions/00_utils_io.R")   # load_raw()
# source("comparison_functions/01_preprocess.R") # depth_filter(), normalise_xyz_target()
# source("comparison_functions/02_split.R")      # km_split()
# 
# # helper to pull the commodity name out of any gw_fit_*.rds file
# get_target <- function(fn)
#   sub(".*_([^_]+_[^_]+)\\.rds$", "\\1", basename(fn))
# 
# # discover which targets you have models for
# model_files <- list.files(pattern = "^gw_fit.*\\.rds$", full.names = TRUE)
# targets     <- unique(purrr::map_chr(model_files, get_target))
# print(targets)   # e.g.  "Cu_pct" "Au_ppb" "Ag_ppm" "Ni_pct" "Density_gcm3"
# 
# # ── 1.  Function that rebuilds train / test for one target ──────────────
# recreate_split <- function(target,
#                            csv_file   = "cluster_0_data.csv",
#                            min_cnt    = 800,
#                            test_frac  = 0.30,
#                            km_centres = 2) {
#   
#   raw    <- load_raw(csv_file, var = target, sample_frac = 1)
#   clean  <- depth_filter(raw, min_cnt = min_cnt)
#   normed <- normalise_xyz_target(clean)
#   
#   set.seed(1)                                       # reproducible!
#   split  <- km_split(normed, frac = test_frac, centres = km_centres)
#   
#   list(train = split$train, test = split$test)      # return both
# }
# 
# # ── 2.  Loop over every target, save the frames ─────────────────────────
# summary_tbl <- map_dfr(targets, function(tg) {
#   
#   cat("→ rebuilding data for", tg, "...\n")
#   ds <- recreate_split(tg)
#   
#   saveRDS(ds$train, file = paste0("train_data_", tg, ".rds"))
#   saveRDS(ds$test , file = paste0("test_data_" , tg, ".rds"))
#   
#   tibble(target = tg,
#          n_train = nrow(ds$train),
#          n_test  = nrow(ds$test))
# })
# 
# cat("\nSaved:\n")
# print(summary_tbl, width = Inf)


