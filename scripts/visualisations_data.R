# main.R -------------------------------------------------------------------
setwd("C:/Users/23478671/Github/GeoWarp-Application-Kevitsa")
source("helper_functions/00_utils_io.R")
source("helper_functions/01_preprocess.R")
source("helper_functions/02_split.R")
source("helper_functions/03_models.R")
source("helper_functions/04_visuals.R")
set.seed(1)
library(purrr)      # map()
library(patchwork)  # wrap_plots()
library(ggplot2)


csv  <- "cluster_0_data.csv"
grade_vars <- c("Ag_ppm", "Au_ppb", "Cu_pct", "Ni_pct")


nice_label <- function(v) {
 switch(v,
        "Ag_ppm" = "log Ag (ppm)",
        "Au_ppb" = "log Au (ppb)",
        "Cu_pct" = "log Cu (%)",
        "Ni_pct" = "log Ni (%)",
        paste0("log ", v))          # fallback
}

plots <- imap(grade_vars, function(v, idx) {
 df <- load_raw(csv, var = v, sample_frac = 1) |>
   depth_filter(min_cnt = 800)
 
 depth_profile_plot(df, bin_size = 2,
                    var_name   = v,
                    show_counts = FALSE) +
   labs(x = nice_label(v))          # replace the x-axis title
})

label_plot <- ggplot() +
 annotate("text", x = 0.5, y = 0.5,
          label = "Depth (m)", angle = 90, size = 4.5) +
 theme_void()

(label_plot | wrap_plots(plots, ncol = 2)) +
 plot_layout(widths = c(0.06, 1))


#outfile <- "depth_profiles.pdf"  # change path/name if desired

# Option 1 ── use ggsave (simplest)
# ggplot2::ggsave(
#   filename = outfile,
#   plot     = label_plot,
#   device   = cairo_pdf,          # crisp text in most viewers
#   width    = 8,                  # inches
#   height   = 10,                 # inches
#   units    = "in"
# )

summary_stats <- function(csv_path, vars, sample_frac = 1, min_cnt = 800) {
  
  purrr::imap_dfr(vars, function(v, i) {
    
    df <- load_raw(csv_path, var = v, sample_frac = sample_frac) |>
      depth_filter(min_cnt = min_cnt)
    
    tibble(
      variable = v,
      N        = nrow(df),
      mean     = mean(df$target, na.rm = TRUE),
      sd       = sd(df$target,   na.rm = TRUE),
      min      = min(df$target, na.rm = TRUE),
      Q1       = quantile(df$target, 0.25, na.rm = TRUE),
      median   = median(df$target,    na.rm = TRUE),
      Q3       = quantile(df$target, 0.75, na.rm = TRUE),
      max      = max(df$target,  na.rm = TRUE)
    )
  })
}

# ---- run it ---------------------------------------------------------------
stats_tbl <- summary_stats(csv, grade_vars)

print(stats_tbl, digits = 3)


