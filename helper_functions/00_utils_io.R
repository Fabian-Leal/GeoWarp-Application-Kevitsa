# 00_utils_io.R -------------------------------------------------------------

library(dplyr)

## load_raw()  --------------------------------------------------------------
# now takes an extra arg `var` giving the column you want to model
load_raw <- function(path, var = "Pb_ppm", sample_frac = 1) {
  df0 <- read.csv(path, stringsAsFactors = FALSE)
  if (! var %in% names(df0)) {
    stop("load_raw(): column '", var, "' not found in the CSV.")
  }
  
  df <- df0 %>%
    rename(
      name   = Name,
      x      = X,
      y      = Y,
      z      = Z,
      target = all_of(var)
    ) %>%
    filter(
      !is.na(x), !is.na(y), !is.na(z), !is.na(target),
      target >= 1e-4
    ) %>%
    # 3) log-transform, subsample, and keep only the columns we want
    mutate(target = log(target)) %>%
    sample_frac(sample_frac) %>%
    select(name, x, y, z, target)
  
  # 4) store the original var name for later use
  attr(df, "target_name") <- var
  df
}


## range_list() helper  ------------------------------------------------------
range_list <- function(x) list(min = min(x), max = max(x))
