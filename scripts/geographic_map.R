# ------------------------------------------------------------
# Complete two‐panel map with a true geographic square inset
# ------------------------------------------------------------

# Install these once if needed:
# install.packages(c("ggmap","ggplot2","patchwork","maps"))
library(png); library(ggplot2)
library(ggmap)
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



final <- pB +
  inset_element(
    kevitsa,
    left   = 0.65,  # start 65% from left
    bottom = 0.65,  # start 65% from bottom
    right  = 1,  # end at 95% from left → width = 0.30
    top    = 1   # end at 95% from bottom → height = 0.30
  )
print(final)



print(pA|final)
print(pA|pB/kevitsa)

