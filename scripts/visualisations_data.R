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





# ── Histograms for the four grade variables ────────────────────────────────

hist_plots <- imap(grade_vars, function(v, i) {
  load_raw(csv, var = v, sample_frac = 1) |>
    depth_filter(min_cnt = 800) |>
    ggplot(aes(x = target)) +
    geom_histogram(bins = 40, colour = "white", fill = "steelblue", alpha = 0.8) +
    labs(x = nice_label(v), y = "count") +
    theme_minimal(base_size = 10)
})

(histograms <- wrap_plots(hist_plots, ncol = 2))

# Optional: save to PDF
#ggsave("histograms_grades.pdf", histograms,
#       width = 8, height = 6, dpi = 300)




# ── Libraries ──────────────────────────────────────────────────────────────
# install.packages("plotly")   # run once if needed
plot_grade_3d_zero_origin <- function(
    csv_path,
    vars,
    angle_theta = 45,
    angle_phi   = 20,
    sample_frac = 1,
    min_cnt     = 800
) {
  if (!requireNamespace("plot3D", quietly = TRUE)) {
    stop("Please install the 'plot3D' package: install.packages('plot3D')")
  }
  if (!requireNamespace("viridis", quietly = TRUE)) {
    stop("Please install the 'viridis' package: install.packages('viridis')")
  }
  library(plot3D)
  library(viridis)
  
  ##— store & tweak global graphics settings ——
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))                      # ← restore on exit, even on error
  par(
    mfrow = c(2, 2),
    mar   = c(1, 5, 1, 5),               # extra room for titles
#    oma   = c(3, 3, 3, 3),
    cex.axis = 0.8,
    mgp   = c(3.2, 0, 0)                   # move axis titles away from ticks
  )
  
  for (v in vars) {
    ## 1) load & filter
    df <- load_raw(csv_path, var = v, sample_frac = sample_frac) |>
      depth_filter(min_cnt = min_cnt)
    
    ## 2) zero-origin coordinates
    x0 <- df$x - min(df$x, na.rm = TRUE)
    y0 <- df$y - min(df$y, na.rm = TRUE)
    z  <- df$z
    
    ## 3) axis limits
    xlim <- c(0, max(x0, na.rm = TRUE))
    ylim <- c(0, max(y0, na.rm = TRUE))
    
    ## 4) 3D scatter
    scatter3D(
      x      = x0,
      y      = y0,
      z      = z,
      colvar = df$target,
      pch    = 20,
      cex    = 0.6,
      cex.main = 1.2,
      theta  = angle_theta,
      phi    = angle_phi,
      xlab   = "\n \n  Northing (m)",
      ylab   = "\n \n  Easting (m)",
      zlab   = "\n \n  Depth (m)",
      main   = nice_label(v),
      xlim   = xlim,
      ylim   = ylim,
      col    = magma(100),
      colkey = list(length = 0.5, width = 0.6,
                    title = nice_label(v), dist = 0),
      ticktype = "detailed",
      tick.ratio = 0,
      bty      = "g"
      
      
    )
  }
}


# Call it:
plot_grade_3d_zero_origin(
  csv_path    = csv,
  vars        = grade_vars,
  angle_theta = 45,
  angle_phi   = 20
)





# Your list of metals
grade_vars <- c("Au_ppb","Ag_ppm","Cu_pct","Ni_pct")

# Compute min & max depth for each
depth_ranges <- map_dfr(grade_vars, function(v) {
  df <- load_raw(csv, var = v, sample_frac = 1) %>%
    depth_filter(min_cnt = 800)
  
  df %>%
    summarise(
      min_depth = min(z, na.rm = TRUE),
      max_depth = max(z, na.rm = TRUE)
    ) %>%
    mutate(variable = v, .before = min_depth)
})

depth_ranges$magnitude <- depth_ranges$max_depth - depth_ranges$min_depth

print(depth_ranges)





