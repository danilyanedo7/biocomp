# Ensure required packages are installed and loaded
if (!require("pacman")) {
  install.packages("pacman")
}
pacman::p_load(pRecipe, giscoR, terra, tidyverse, rayshader, sf, classInt)

# Load country boundary data and subset to the desired region
regions_sf1 <- geodata::gadm(country = "ID",
                             level = 2,
                             path = getwd()) |>
  sf::st_as_sf() #region
region_sf <- subset(regions_sf1, NAME_1 == "Jawa Timur") #sub region
region_sf2 <- subset(regions_sf1, NAME_2 == "Malang") #sub sub region (city)
malang_bbox <- st_bbox(region_sf2)
print(malang_bbox)

# Download and prepare precipitation data
pRecipe::download_data(
  dataset = "mswep",
  path = getwd(),
  domain = "raw",
  timestep = "yearly"
)
list.files() #copy name of the file onto below code
mswep_data <- terra::rast("mswep_tp_mm_global_197902_202301_025_yearly.nc") |>
  terra::crop(region_sf)
terra::plot(mswep_data[[1]])
region_sf2_utm <- st_transform(region_sf, crs = 32749)

# Plot the reprojected geometry
plot(sf::st_geometry(region_sf2_utm), add = TRUE)
plot(sf::st_geometry(region_sf), add = TRUE)

# Prepare data for panel visualization
names(mswep_data) <- 1979:2023
mswep_df <- mswep_data |>
  as.data.frame(xy = TRUE) |>
  tidyr::pivot_longer(!c("x", "y"), names_to = "year", values_to = "precipitation") |>
  dplyr::filter(year >= 2020) |>
  dplyr::filter(year != 2023)
head(mswep_df)

# Define theme and color breaks
theme_for_the_win <- function() {
  # Function to customize plot theme
  theme_minimal() +
    theme(
      axis.line = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 11, color = "grey10"),
      legend.text = element_text(size = 10, color = "grey10"),
      panel.grid.major = element_line(color = NA),
      panel.grid.minor = element_line(color = NA),
      plot.background = element_rect(fill = NA, color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(fill = NA, color = NA),
      plot.margin = unit(c(
        t = 0,
        r = 0,
        b = 0,
        l = 0
      ), "lines")
    )
}

breaks <- classInt::classIntervals(mswep_df$precipitation, n = 5, style = "equal")$brks
colors <- hcl.colors(n = length(breaks),
                     palette = "Temps",
                     rev = TRUE)

# Create panel map for precipitation
map1 <- ggplot(data = mswep_df) +
  geom_raster(aes(x = x, y = y, fill = precipitation)) +
  geom_sf(
    data = region_sf,
    fill = "transparent",
    color = "grey10",
    size = .5
  ) +
  geom_sf(data = region_sf[region_sf$NAME_2 == "Malang", ],
          color = "black",
          fill = "transparent") +
  scale_fill_gradientn(
    name = "mm",
    colors = colors,
    breaks = breaks,
    labels = round(breaks, 0),
    limits = c(min(mswep_df$precipitation), max(mswep_df$precipitation))
  ) +
  facet_wrap( ~ year) + #adding panel to compare each year
  guides(
    fill = guide_colourbar(
      direction = "vertical",
      barheight = unit(50, "mm"),
      barwidth = unit(5, "mm"),
      title.position = "top",
      label.position = "right",
      title.hjust = .5,
      label.hjust = .5,
      ncol = 1,
      byrow = FALSE
    )
  ) +
  geom_magnify(
    from = c(112.366325, 112.776237, -7.409174, -6.711704),
    to = c(115, 116, -7, -5.5),
    colour = "red",
    linewidth = 0.3,
    proj.fill = alpha("red", 0.2)
  ) +
  theme_bw()
map1

# Save the map
ggsave(
  "./map1_presipitasimalang.jpeg",
  map1,
  width = 12,
  height = 4,
  units = "in",
  dpi = 450
)


# Calculate and visualize average precipitation
mswep_average_df <- mswep_df |>
  dplyr::group_by(x, y, .drop = FALSE) |>
  dplyr::summarise(mean = mean(precipitation))
head(mswep_average_df)

breaks <- classInt::classIntervals(mswep_average_df$mean, n = 5, style = "equal")$brks
colors <- hcl.colors(n = length(breaks),
                     palette = "Temps",
                     rev = TRUE)

map2 <- ggplot(data = mswep_average_df) +
  geom_raster(aes(x = x, y = y, fill = mean)) +
  geom_sf(
    data = region_sf,
    fill = "transparent",
    color = "grey10",
    size = .5
  ) +
  scale_fill_gradientn(
    name = "mm",
    colors = colors,
    breaks = breaks,
    labels = round(breaks, 0),
    limits = c(min(mswep_average_df$mean), max(mswep_average_df$mean))
  ) +
  guides(
    fill = guide_colourbar(
      direction = "vertical",
      barheight = unit(50, "mm"),
      barwidth = unit(5, "mm"),
      title.position = "top",
      label.position = "right",
      title.hjust = .5,
      label.hjust = .5,
      ncol = 1,
      byrow = FALSE
    )
  ) +
  geom_magnify(
    from = c(112.366325, 112.776237, -7.409174, -6.711704),
    to = c(115, 116, -7, -5.5),
    colour = "red",
    linewidth = 0.3,
    proj.fill = alpha("red", 0.2)
  ) +
  theme_bw()
map2
