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















# ------------------------------------------------------------
# Complete two‐panel map with a true geographic square inset
# ------------------------------------------------------------

# Install these once if needed:
# install.packages(c("ggmap","ggplot2","patchwork","maps"))

library(ggmap)
library(ggplot2)
library(patchwork)
library(maps)
library(ggspatial)

# 0. Register your Stadia Maps key ---------------------------------------
register_stadiamaps("af7dc6f7-6cae-4b55-b12b-3f8474e43d04")

# 1. Kevitsa coordinates -------------------------------------------------
kev_lon <- 26.92937428901149
kev_lat <- 67.69053660061874

# Compute half-spans so that the inset is a true square (~55 km × 55 km)
lat_rad      <- kev_lat * pi/180
half_lat_deg <- 0.5                  # 0.5° ≈ 55 km north–south
half_lon_deg <- half_lat_deg / cos(lat_rad)

# 2. Panel A: Finland overview (watercolor) ----------------------------
fin_bbox <- c(left = 15, bottom = 57, right = 35, top = 72)
fin_map  <- get_stadiamap(
  bbox    = fin_bbox,
  zoom    = 4,
  source  = "stamen",
  maptype = "stamen_watercolor"
)

cities <- data.frame(
  name = c("Helsinki", "Oulu", "Rovaniemi", "Sodankylä"),
  lon  = c(24.9384,    25.4651,  25.7482,     26.5895),
  lat  = c(60.1699,    65.0121,  66.5039,     67.4167)
)


# Get Finland boundary from map_data
finland_map <- map_data("world", region = "Finland")

# Red square showing exactly the same geographic square used in Panel B
rectA <- data.frame(
  lon = c(kev_lon - half_lon_deg, kev_lon + half_lon_deg,
          kev_lon + half_lon_deg, kev_lon - half_lon_deg,
          kev_lon - half_lon_deg),
  lat = c(kev_lat - half_lat_deg, kev_lat - half_lat_deg,
          kev_lat + half_lat_deg, kev_lat + half_lat_deg,
          kev_lat - half_lat_deg)
)

pA <- ggmap(fin_map) +
  geom_path(
    data = finland_map,
    aes(x = long, y = lat, group = group),
    colour = "black", size = 0.7
  ) +
  geom_path(
    data = rectA, aes(x = lon, y = lat),
    colour = "red", linewidth = 0.7
  ) +
  annotate(
    "text",
    x     = mean(rectA$lon),
    y     = max(rectA$lat) - 2.05,
    label = "Finland",
    fontface = "bold",
    colour   = "grey20",
    size     = 5
  ) +
  geom_point(data = cities,
             aes(x = lon, y = lat),
             shape = 21, fill = "white",
             colour = "black", size = 2, stroke = 0.6) +
  geom_text(data = cities,
            aes(x = lon, y = lat, label = name),
            hjust = 0.355, vjust = -1.2,
            size = 2.5, colour = "black") +
  # graticules & degree labels
  scale_x_continuous(
    name   = NULL, 
    breaks = seq(15, 35, by = 5),
    labels = function(x) paste0(x, "°E"),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name   = NULL, 
    breaks = seq(57, 72, by = 5),
    labels = function(y) paste0(y, "°N"),
    expand = c(0, 0)
  ) 

# 3. Panel B: Local zoom (watercolor) ----------------------------------
local_bbox <- c(
  left   = kev_lon - half_lon_deg,
  bottom = kev_lat - half_lat_deg,
  right  = kev_lon + half_lon_deg,
  top    = kev_lat + half_lat_deg
)
zoom_map <- get_stadiamap(
  bbox    = local_bbox,
  zoom    = 8,
  source  = "stamen",
  maptype = "stamen_watercolor"
)

# Small red rectangle (same square) around the mine
# … (all the code up to building zoom_map is unchanged) …

# Small red rectangle (same square) around the mine
rectB <- rectA

#  Add a little data.frame of locations to mark:
# … everything up to constructing 'towns' is the same …

# towns data.frame:
towns <- data.frame(
  name  = c("Kevitsa Mine", "Sodankylä"),
  lon   = c(kev_lon,      26.5895),
  lat   = c(kev_lat,      67.4167),
  size  = c(5,            3),    # circle sizes
  shape = c(22,           21)    # 23 = diamond, 21 = circle
)

towns$fill <- c("grey", "white")  # Kevitsa = black, Sodankylä = grey

# 2) map fill in aes(), and use scale_*_identity()
library(ggspatial)

# 1. compute the display limits again
local_bbox <- c(
  left   = kev_lon - half_lon_deg,
  bottom = kev_lat - half_lat_deg,
  right  = kev_lon + half_lon_deg,
  top    = kev_lat + half_lat_deg
)

# 2. rebuild pB with axes, grid, north arrow & scale bar
pB <- ggmap(zoom_map, extent = "panel") +
  # red bounding square
  geom_path(data = rectB, aes(x = lon, y = lat),
            colour = "red", linewidth = 0.3) +
  # points for Kevitsa + Sodankylä
  geom_point(data = towns,
             aes(x = lon, y = lat, shape = shape, size = size, fill = fill),
             colour = "black", stroke = 1.2) +
  scale_shape_identity() +
  scale_size_identity() +
  scale_fill_identity() +
  # labels
  geom_text(data = towns,
            aes(x = lon + 0.03, y = lat + 0.03, label = name),
            hjust = 0, vjust = 0,
            colour = "black", fontface = "bold", size = 3) +
  
  # bring back lon/lat axes
  scale_x_continuous(
    name   = NULL,
    breaks = seq(ceiling(local_bbox["left"]),
                 floor(local_bbox["right"]), by = 0.5),
    labels = function(x) paste0(x, "°E"),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name   = NULL,
    breaks = seq(
      from = floor(local_bbox["bottom"]),
      to   = ceiling(local_bbox["top"]),
      by   = 0.2
    ),
    labels = function(y) paste0(y, "°N"),
    expand = c(0, 0)
  ) +
  
  # add a subtle grid
  theme_minimal(base_family = "sans") +
  theme(
    panel.grid.major = element_line(colour = "grey80", linetype = "dotted"),
    panel.background = element_rect(fill = NA, colour = NA),
    axis.text        = element_text(size = 8, colour = "black"),
    axis.ticks       = element_line(size = 0.3)
  ) +
  
  # north arrow and scale bar
  annotation_north_arrow(
    location    = "tl",
    which_north = "true",
    style       = north_arrow_fancy_orienteering()
  ) +
  coord_fixed(
    ratio  = 1 / cos(lat_rad),
    xlim   = c(local_bbox["left"],  local_bbox["right"]),
    ylim   = c(local_bbox["bottom"], local_bbox["top"]),
    expand = FALSE
  )






final <- pB +
  inset_element(
    kevitsa,
    left   = 0.65,  # start 65% from left
    bottom = 0.65,  # start 65% from bottom
    right  = 1,  # end at 95% from left → width = 0.30
    top    = 1   # end at 95% from bottom → height = 0.30
  )
print(final)




library(png); library(ggplot2)

img <- png::readPNG("figures/kevitsa_mine_2.png")

kevitsa <- ggplot() +
  # draw the raster exactly from 0→1 in both directions
  annotation_raster(img,
                    xmin = 0, xmax = 1,
                    ymin = 0, ymax = 1) +
  # force the scales to exactly [0,1] with no expansion
  scale_x_continuous(expand = c(0,0), limits = c(0,1)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) +
  # draw a 1pt black border around the panel
  theme_void() +
  theme(
    panel.border = element_rect(colour = "black",
                                size   = 1,
                                fill   = NA),
    plot.margin  = unit(c(0,0,0,0), "pt")
  )

print(kevitsa)






# 4. Combine: place pB as an inset on top of pA --------------------------
final_plot_swapped <- pA +
  inset_element(
    pB,
    left   = 0.65,   # 65% in from the left edge of pA
    bottom = -0.3,   #  5% in from the bottom edge of pA
    right  = 1.00,   # full width of pA
    top    = 0.50    # up to 40% of pA’s height
  )

print(final_plot_swapped)





