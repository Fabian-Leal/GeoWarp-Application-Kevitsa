.libPaths("C:/Users/23478671/R_libs")
setwd("C:/Users/23478671/Github/GeoWarp-Application-Kevitsa")
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
library(tidyr)
library(stringr)

fits_dir   <- "model_fits"
keep_vars  <- c("Ag_ppm", "Au_ppb", "Cu_pct", "Ni_pct")   # <- keep only these

# helper to pull "Au_ppb" from  "gw_fit_noWarp_Au_ppb.rds"
get_target <- function(path)
  sub(".*_([^_]+_[^_]+)\\.rds$", "\\1", basename(path))

# keep GeoWarp fits that match one of the four metals
gw_files <- list.files(
  fits_dir, pattern = "^gw_.*\\.rds$", full.names = TRUE
) %>%
  keep(~ get_target(.x) %in% keep_vars)









# ── helper: split scalar vs vector parameters ---------------------------------
split_params <- function(fit) {
  prm <- fit$parameters
  list(
    scalars = prm[vapply(prm, function(x) length(x) == 1, logical(1))],
    vectors = prm[vapply(prm, function(x) length(x)  > 1, logical(1))]
  )
}

# ── 1)  scalar parameters -----------------------------------------------------
gw_scalars <- map_dfr(gw_files, function(fn) {
  p <- split_params(readRDS(fn))$scalars
  tibble(
    fit_file        = basename(fn),
    tau_mu_rand     = p$tau_squared_mean_random,
    sigma2_nugget   = p$sigma_squared_nugget,
    eta_deviation   = p$eta_deviation,
    ell_dev_random  = p$ell_deviation_random,
    tau2_dev_random = p$tau_squared_deviation_random
  )
})

# ── 2)  vector / matrix parameters in long form ------------------------------
gw_vectors <- map_dfr(gw_files, function(fn) {
  p <- split_params(readRDS(fn))$vectors
  fname <- basename(fn)
  
  bind_rows(
    tibble(
      fit_file = fname,
      param    = "zeta_deviation",
      index    = seq_along(p$zeta_deviation),
      value    = p$zeta_deviation
    ),
    tibble(
      fit_file = fname,
      param    = "gamma_deviation_vertical",
      index    = seq_along(p$gamma_deviation_vertical),
      value    = p$gamma_deviation_vertical
    ),
    tibble(
      fit_file = fname,
      param    = "L_deviation",
      index    = seq_along(as.vector(p$L_deviation)),
      value    = as.vector(p$L_deviation)
    ),  
    tibble(
      fit_file = fname,
      param    = "gamma_deviation_horizontal",
      index    = seq_along(p$gamma_deviation_horizontal),
      value    = as.vector(p$gamma_deviation_horizontal)
    )
  )
})



library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

# 1) Point at all your GeoWarp fits
gw_files <- list.files("model_fits",
                       pattern = "^gw_fit.*\\.rds$",
                       full.names = TRUE)

# 2) Read each fit and unroll all the vector parameters into one long table
long_gw <- map_dfr(gw_files, function(fn) {
  fit     <- readRDS(fn)
  name    <- basename(fn)
  variant <- case_when(
    str_detect(name, "gw_fit_noWarp_cv") ~ "GeoWarp LinearWarp (CV)",
    str_detect(name, "gw_fit_noWarp")    ~ "GeoWarp LinearWarp",
    str_detect(name, "gw_fit_cv")        ~ "GeoWarp DepthDeform (CV)",
    TRUE                                 ~ "GeoWarp DepthDeform"
  )
  metal   <- str_extract(name, "Ag_ppm|Au_ppb|Cu_pct|Ni_pct")
  
  # grab the 3×3 matrix, the two horizontal gammas,
  # the vertical gammas (may be length 0 for LinearWarp),
  # and the zeta spline
  p <- fit$parameters
  tib_L  <- tibble(param = "L_deviation",
                   index = seq_along(as.vector(p$L_deviation)),
                   value = as.vector(p$L_deviation))
  tib_gh <- tibble(param = "gamma_deviation_horizontal",
                   index = seq_along(p$gamma_deviation_horizontal),
                   value = p$gamma_deviation_horizontal)
  tib_gv <- tibble(param = "gamma_deviation_vertical",
                   index = seq_along(p$gamma_deviation_vertical),
                   value = p$gamma_deviation_vertical)
  tib_z  <- tibble(param = "zeta_deviation",
                   index = seq_along(p$zeta_deviation),
                   value = p$zeta_deviation)
  
  bind_rows(tib_L, tib_gh, tib_gv, tib_z) %>%
    mutate(variant = variant, metal = metal)
})

# 3) Pivot it so each metal is its own column
wide_gw2 <- long_gw %>%
  filter(!is.na(metal)) %>%   # drop any stray non-metals
  pivot_wider(
    id_cols     = c(variant, param, index),
    names_from  = metal,
    values_from = value
  ) %>%
  arrange(variant, param, index)

# 4) Inspect
print(wide_gw2, n = 30, digits = 3)



library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)   # pretty_breaks()
library(viridis)  # nice perceptual palette

# ── 1.  Long format ────────────────────────────────────────────────────────
gw_long <- wide_gw2 %>%                                    # 121 × 7
  pivot_longer(cols = c(Ag_ppm, Au_ppb, Cu_pct, Ni_pct),   # -> 484 × 6
               names_to  = "metal",
               values_to = "value")

# Optional: order factors the way you want them to appear
gw_long <- gw_long %>% 
  mutate(
    variant = factor(variant,
                     levels = c("GeoWarp DepthDeform",
                                "GeoWarp DepthDeform (CV)",
                                "GeoWarp LinearWarp",
                                "GeoWarp LinearWarp (CV)")),
    param   = factor(param,
                     levels = c("L_deviation",
                                "gamma_deviation_horizontal",
                                "gamma_deviation_vertical",
                                "zeta_deviation"))
  )

# ── 2.  Heat-map ───────────────────────────────────────────────────────────
p <- ggplot(gw_long,
            aes(x = metal,
                y = factor(index),          # treat index as discrete for tidy spacing
                fill = value)) +
  geom_tile() +
  facet_grid(variant ~ param,               # rows = fit variant, cols = parameter block
             scales = "free_y",
             space  = "free_y") +
  scale_fill_viridis_c(option = "plasma",
                       trans   = "log10",   # comment this line if you prefer linear scale
                       na.value = "grey90",
                       breaks  = pretty_breaks(5)) +
  labs(x    = NULL,
       y    = "index",
       fill = "value") +
  theme_minimal(base_size = 11) +
  theme(
    strip.text.y = element_text(angle = 0),  # keep facet labels horizontal
    axis.text.y  = element_blank(),          # hide long list of indices
    axis.ticks.y = element_blank(),
    panel.spacing = unit(0.6, "lines")
  )

print(p)

# ── 3.  Save to file (optional) ────────────────────────────────────────────
ggsave("wide_gw2_heatmap.pdf", p, width = 9, height = 6, dpi = 300)

















# 1) install xtable if you don’t have it
if (!requireNamespace("xtable", quietly = TRUE)) {
  install.packages("xtable")
}
library(xtable)

# 2) prepare the table: combine param+index, replace NA with “–”
tbl <- wide_gw2 %>%
  mutate(
    Parameter = paste0(param, "_", index),
    across(c(Ag_ppm, Au_ppb, Cu_pct, Ni_pct),
           ~ ifelse(is.na(.), "–", formatC(., format="e", digits=3)))
  ) %>%
  select(Variant = variant, Parameter, Ag_ppm, Au_ppb, Cu_pct, Ni_pct)

# 3) print as a LaTeX longtable
print(
  xtable(
    tbl,
    caption = "GeoWarp Vector Parameters by Variant and Metal",
    label   = "tab:gw_params"
  ),
  include.rownames     = FALSE,
  floating             = FALSE,
  table.placement      = "!htbp",
  tabular.environment  = "longtable",
  booktabs             = TRUE,
  sanitize.text.function = identity
)





