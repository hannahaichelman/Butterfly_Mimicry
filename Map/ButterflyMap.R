## Author: Hannah E. Aichelman
## This script creates a map of the western United States where Limenitis butterflies of interest were collected for this study.

# helpful info at this tutorial: https://cran.r-project.org/web/packages/usmap/vignettes/usmap3.html

#### Set-Up ####
#install.packages('usmap')

library(usmap)
library(ggplot2)

#### Create Basic State Map ####
# example of creating a map of the entire US:
plot_usmap(regions = "counties") + 
  labs(title = "US Counties",
       subtitle = "This is a blank map of the counties of the United States.") + 
  theme(panel.background = element_rect(color = "black", fill = "lightblue"))


# plot only western states of interest for butterfly ranges:
map = plot_usmap(include = c("CA", "ID", "NV", "OR", "WA"), fill = c("#74c476", "#fed976", "#fb6a4a", "#fed976", "#fed976"), alpha = 0.7) +
  #labs(title = "Western US States") +
  theme_bw()
map
ggsave(map, file = "/Users/hannahaichelman/Dropbox/BU_Postdoc/ButterflyGenomes/Github/Butterfly_Mimicry/Map/RangeMap.pdf", width=4, height=5, units=c("in"), useDingbats=FALSE)

#### Topographic Map ####
# more info here on how I found this code: https://github.com/milos-agathon/create-crisp-topographic-maps-with-r/

# libraries we need
libs <- c("elevatr", "terra", "tidyverse","sf", "giscoR", "osmdata", "marmap", "maps")

# install missing libraries
#installed_libs <- libs %in% rownames(installed.packages())
#if (any(installed_libs == F)) {
#  install.packages(libs[!installed_libs])
#}

# load libraries
invisible(lapply(libs, library, character.only = T))

# 1. GET SHAPE FILE FOR COUNTRY MAP
#------------------
crsLONGLAT <- "+proj=longlat +datum=WGS84 +no_defs"

# create function to pull shape file of the whole US without state borders
get_country_sf <- function() {
  country_sf <- giscoR::gisco_get_countries(
    year = "2020",
    epsg = "4326",
    resolution = "10",
    country = "USA"
  ) |>
    sf::st_crop(c(xmin = -125, xmax = -60, ymin = 20, ymax = 60)) |>
    sf::st_transform(crs = crsLONGLAT)
  
  return(country_sf)
}

country_sf <- get_country_sf()


# other version of shape file for USA with state borders:
states_sf <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
head(states_sf)

# 2. GET COUNTRY ELEVATION DATA
#------------------------------
get_elevation_data <- function() {
  country_elevation <- elevatr::get_elev_raster(
    locations = states_sf,
    #locations = country_sf,
    z = 7,
    clip = "locations"
  )
  
  return(country_elevation)
}

country_elevation <- get_elevation_data()

# plot elevation data
terra::plot(country_elevation)

# 3. GET ELEVATION DATA FOR BINDING BOX (BBOX)
#------------------------------
get_elevation_data_bbox <- function() {
  country_elevation <- elevatr::get_elev_raster(
    locations = states_sf,
    #locations = country_sf,
    z = 7,
    clip = "bbox"
  )
  
  return(country_elevation)
}

country_elevation <- get_elevation_data_bbox() |>
  terra::rast()

# 4. PLOT
#---------
country_elevation |>
  as.data.frame(xy = T) |>
  ggplot() +
  geom_tile(
    aes(x = x, y = y, fill = file1767ac86bd71)
  ) +
  geom_sf(
    #data = country_sf,
    data = states_sf,
    fill = "transparent", color = "yellow", size = .25
  ) +
  theme_void()

# 5. CROP AREA
#--------------------
get_area_bbox <- function() {
  xmin <- -125
  xmax <- -108
  ymin <- 30
  ymax <- 50
  
  bbox <- sf::st_sfc(
    sf::st_polygon(list(cbind(
      c(xmin, xmax, xmax, xmin, xmin),
      c(ymin, ymin, ymax, ymax, ymin)
    ))),
    crs = crsLONGLAT
  )
  
  return(bbox)
}

bbox <- get_area_bbox()

crop_area_with_polygon <- function() {
  bbox_vect <- terra::vect(bbox)
  bbox_raster <- terra::crop(country_elevation, bbox_vect)
  bbox_raster_final <- terra::mask(
    bbox_raster, bbox_vect
  )
  return(bbox_raster_final)
}

bbox_raster_final <- crop_area_with_polygon()

bbox_raster_final |>
  as.data.frame(xy = T) |>
  ggplot() +
  geom_tile(
    aes(x = x, y = y, fill = file1767ac86bd71)
  ) +
  geom_sf(
    #data = country_sf,
    data = states_sf,
    fill = "transparent", color = "black", size = .25
  ) +
  theme_void()

# # 6. GET REGION LINES
# #--------------------
# region <- "California"
# # define longlat projection
# 
# westcoast_sf <- osmdata::getbb(
#   region,
#   format_out = "sf_polygon"
# )
# 
# westcoast_sf
# 
# ggplot() +
#   geom_sf(
#     data = westcoast_sf$multipolygon,
#     color = "red", fill = "grey80", size = .5
#   ) +
#   theme_void()
# 
# crop_region_with_polygon <- function() {
#   region_vect <- terra::vect(sardinia_sf$multipolygon)
#   region_raster <- terra::crop(country_elevation, region_vect)
#   region_raster_final <- terra::mask(
#     region_raster, region_vect
#   )
#   return(region_raster_final)
# }
# 
# region_raster_final <- crop_region_with_polygon()
# 
# region_raster_final |>
#   as.data.frame(xy = T) |>
#   ggplot() +
#   geom_tile(
#     aes(x = x, y = y, fill = file514862c13e19)
#   ) +
#   geom_sf(
#     data = country_sf,
#     fill = "transparent", color = "black", size = .25
#   ) +
#   theme_void()

# 7. FINAL MAP
#-------------
get_elevation_map <- function() {
  #country_elevation_df <- country_elevation |>
  country_elevation_df <- bbox_raster_final |>
    as.data.frame(xy = T) |>
    na.omit()
  
  names(country_elevation_df)[3] <- "elevation"
  
  country_map <-
    ggplot(data = country_elevation_df) +
    geom_raster(
      aes(x = x, y = y, fill = elevation),
      alpha = 1
    ) +
    marmap::scale_fill_etopo() +
    coord_sf(crs = crsLONGLAT) +
    geom_sf(
      #data = country_sf,
      data = states_sf,
      fill = "transparent", color = "black", size = .25
    ) +
    labs(
      x = "",
      y = "",
      title = "",
      subtitle = "",
      caption = ""
    ) +
    theme_minimal() +
    theme(
      #axis.line = element_blank(),
      #axis.text.x = element_blank(),
      #axis.text.y = element_blank(),
      #axis.ticks = element_blank(),
      #legend.position = "none",
      #panel.grid.major = element_blank(),
      #panel.grid.minor = element_blank(),
      #plot.margin = unit(c(t = 0, r = 0, b = 0, l = 0), "cm"),
      #plot.background = element_blank(),
      #panel.background = element_blank(),
      #panel.border = element_blank()
    )
  
  return(country_map)
}

country_map <- get_elevation_map()

ggsave(
  filename = "/Users/hannahaichelman/Dropbox/BU_Postdoc/ButterflyGenomes/Github/Butterfly_Mimicry/Map/westcoast_topo_map_statelines.png", width = 7, 
  height = 8.5, dpi = 600, device = "png", 
  country_map, bg = "white"
)
