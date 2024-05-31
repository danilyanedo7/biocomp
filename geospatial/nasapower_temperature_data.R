# Ensure required packages are installed and loaded
if (!require("pacman")) {
  install.packages("pacman")
}

pacman::p_load(nasapower, giscoR, terra, tidyverse, rayshader, sf, classInt)

# SHP Study Sites Extent
regions_sf1 <- geodata::gadm(country = "ID",
                             level = 2,
                             path = getwd()) |>
  sf::st_as_sf()

region_sf <- subset(regions_sf1, NAME_1 == "Jawa Timur") #sub region
region_sf2 <- subset(regions_sf1, NAME_2 == "Gresik") #sub sub region (city)
gresik_bbox <- st_bbox(region_sf2) #finding bounding box of the desired region
print(gresik_bbox)

# Calculate yearly average temperature
yearly_avg <- regional_t2m %>%
  rowwise() %>%
  mutate(avg_temp = mean(c_across(JAN:DEC), na.rm = TRUE)) %>%
  ungroup() %>%
  select(YEAR, LAT, LON, avg_temp)

# Convert data for ggplot2
yearly_avg_ggplot <- yearly_avg %>%
  rename(year = YEAR, lat = LAT, lon = LON) %>%
  mutate(year = as.factor(year))  # Convert year to factor for facet_wrap

region_sf <- st_transform(region_sf, crs = 4326)

# Plot with ggplot2
breaks <- classInt::classIntervals(yearly_avg_ggplot$avg_temp, n = 5, style = "equal")$brks

colors <- hcl.colors(n = length(breaks),
                     palette = "Temps",
                     rev = TRUE)
map1 <- ggplot(data = yearly_avg_ggplot) +
  geom_raster(aes(x = lon, y = lat, fill = avg_temp)) +
  geom_sf(
    data = region_sf,
    fill = "transparent",
    color = "grey10",
    size = .5
  ) +
  geom_sf(data = region_sf[region_sf$NAME_2=="Gresik",],color = "white", fill = "transparent",)+
  facet_wrap( ~ year) +
  scale_fill_viridis_c(name = "Average\nTemperature\n(°C)", limits = c(
    min(yearly_avg_ggplot$avg_temp),
    max(yearly_avg_ggplot$avg_temp)
  )) +
  theme(legend.position = "right",
        legend.key.height = unit(1, "inch")) +
  geom_magnify(
    from = c(112.366325, 112.776237, -7.409174, -6.711704),
    to = c(115, 116, -7, -5.5),
    colour = "red",
    linewidth = 0.3,
    proj.fill = alpha("red", 0.2)
  ) +
  theme_bw()
map1

#saving high quality plot
ggsave(
  "./nameoftheplot.jpeg",
  map1,
  width = 12,
  height = 4,
  units = "in",
  dpi = 450
)
