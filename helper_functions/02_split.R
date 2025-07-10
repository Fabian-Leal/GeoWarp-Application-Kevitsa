# 02_split.R ---------------------------------------------------------------

library(dplyr)
library(FNN)          # k-NN for density

## 1. gap_split() – simple x-cutoff ----------------------------------------
gap_split <- function(df, x_cut = 0.5) {
  list(train = df %>% filter(x <= x_cut),
       test  = df %>% filter(x >  x_cut))
}







## 2. strata_split() – density × dip, 20% per stratum ----------------------
strata_split <- function(df, k = 10, test_frac = 0.3, dips = NULL) {
  
  holes <- df %>% group_by(name) %>%
    summarise(
      cx = first(x), cy = first(y), cz = first(z),
      Dip = if (!is.null(dips)) first(dips[name == cur_group()$name]) else NA_real_,
      .groups = "drop")
  
  dist_k <- FNN::knn.dist(as.matrix(holes[, c("cx","cy","cz")]), k = k)[, k]
  holes$density_q <- ntile(1/dist_k, 4)
  
  holes$dip_cls <- if (!all(is.na(holes$Dip))) {
    cut(abs(holes$Dip), c(-Inf, 30, 60, Inf),
        labels = c("shallow","moderate","steep"), right = FALSE)
  } else factor("unknown")
  
  holes$stratum <- interaction(holes$density_q, holes$dip_cls, drop = TRUE)
  
  set.seed(1)
  test_holes <- holes %>%
    group_by(stratum) %>% sample_frac(test_frac) %>% pull(name)
  
  list(train = df %>% filter(!name %in% test_holes),
       test  = df %>% filter( name %in% test_holes))
}



## 3. km_split() – pick the k-means cluster closest to `frac` ---------
km_split <- function(df, frac = 0.30, centres = 2) {
  rng <- attributes(df)$ranges          # saved by normalise_xyz_target()
  if (is.null(rng)) stop("normalised ranges not found")
  
  unscale <- function(v, r) v * (r$max - r$min) + r$min
  
  # one row per hole, back-transform to metres
  holes <- df %>%
    group_by(name) %>%
    summarise(
      x_m = unscale(first(x), rng$x),
      y_m = unscale(first(y), rng$y),
      .groups = "drop"
    )
  
  set.seed(1)
  km <- kmeans(holes[, c("x_m", "y_m")], centers = centres)
  holes$cluster <- km$cluster
  
  # choose the cluster whose size is nearest to `frac`
  cl_sizes <- holes %>% count(cluster) %>%
    mutate(frac_cl = n / sum(n))
  test_cl <- cl_sizes %>% slice_min(abs(frac_cl - frac)) %>% pull(cluster)
  
  test_holes  <- holes %>% filter(cluster == test_cl) %>% pull(name)
  list(
    train = df %>% filter(!name %in% test_holes),
    test  = df %>% filter( name %in% test_holes)
  )
}

