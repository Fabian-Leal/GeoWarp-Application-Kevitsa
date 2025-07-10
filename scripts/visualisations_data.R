# main.R -------------------------------------------------------------------
setwd("C:/Users/23478671/Github/Geowarp-Resource-Estimation/GeoWarp_resource_estimation/application_kevitsa")
source("helper_functions/00_utils_io.R")
source("helper_functions/01_preprocess.R")
source("helper_functions/02_split.R")
source("helper_functions/03_models.R")
source("helper_functions/04_visuals.R")

set.seed(1)
library(purrr)      # map()
library(patchwork)  # wrap_plots()
csv  <- "cluster_0_data.csv"
grade_vars <- c("Ag_ppm", "Au_ppb", "Cu_pct", "Ni_pct")

plots <- map(grade_vars, function(v) {
df <- load_raw(csv, var = v, sample_frac = 1)            # ① load
df <- depth_filter(df, min_cnt = 800)                    # ② filter
depth_profile_plot(df, bin_size = 20, var_name = v, show_counts = FALSE)
})
 
 wrap_plots(plots, ncol = 2)  
 
 
 
 

library(purrr)
library(patchwork)
library(ggplot2)

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


outfile <- "depth_profiles.pdf"  # change path/name if desired

# Option 1 ── use ggsave (simplest)
ggplot2::ggsave(
  filename = outfile,
  plot     = label_plot,
  device   = cairo_pdf,          # crisp text in most viewers
  width    = 8,                  # inches
  height   = 10,                 # inches
  units    = "in"
)
