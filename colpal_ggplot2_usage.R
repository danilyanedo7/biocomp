# Step 1: Installation
# Install the RColorBrewer package and other relevant package
install.packages("tidyverse")
install.packages("RColorBrewer")

# Step 2: Loading Required Libraries
# Load the RColorBrewer library
library(RColorBrewer)

# Step 3: Using Pre-defined Color Palettes
# Display available pre-defined palettes
display.brewer.all()

# Choose a pre-defined palette
my_palette <- brewer.pal(8, "Set1")

# Display the colors in the palette
print(my_palette)

# Display color blind friendly palette only
display.brewer.all(colorblindFriendly = TRUE)

# Step 4: Customizing Color Palettes
# Define a custom color palette
custom_palette <- c("blue", "green", "red", "orange") 
# You can also use hexcode 
custom_palette2 <- c("#e76254", "#ef8a47", "#f7aa58", "#ffd06f", "#ffe6b7", "#aadce0", "#72bcd5", "#528fad", "#376795") 

# Display the colors in the custom palette
print(custom_palette)

# Step 5: Using Color Palettes in Plots
# Generate sample data
data <- data.frame(
  x = 1:10,
  y = rnorm(10),
  group = sample(letters[1:4], 10, replace = TRUE)
)

# Create a scatter plot with custom colors
plot(data$x, data$y, col = custom_palette[data$group], pch = 16)
legend("topright", legend = unique(data$group), col = custom_palette, pch = 16)



# Beautiful color palette inspired from arts
install.packages("MetBrewer")
install.packages("devtools")
devtools::install_github("BlakeRMills/MetBrewer")
library(MetBrewer)
met.brewer('Hiroshige', n = 100)


# If you are familiar with ggplot2, you can create your own function layer for ggplot2
library(tidyverse)
# In this part, i named this palette as tropical six palette, i used to have special palette edit for my photography project

tropical_six <- list(
  "prussian blue"    = "#103851",
  "hunter green"     = "#316A4F",
  "satin sheen gold" = "#BAA147",
  "caramel"          = "#BE7E54",
  "cadet gray"       = "#8499AA",
  "caput mortuum"    = "#592822"
)


tropical_pal <- function(
    primary = "prussian blue",
    other = "cadet gray",
    direction = 1
) {
  stopifnot(primary %in% names(tropical_six))
  
  function(n) {
    if (n > 6) warning("Tropical Color Palette only has 6 colors cok.")
    
    if (n == 2) {
      other <- if (!other %in% names(tropical_six)) {
        other
      } else {
        tropical_six[other]
      }
      color_list <- c(other, tropical_six[primary])
    } else {
      color_list <- tropical_six[1:n]
    }
    
    color_list <- unname(unlist(color_list))
    if (direction >= 0) color_list else rev(color_list)
  }
}

scale_colour_tropical <- function(
    primary = "prussian blue",
    other = "cadet gray",
    direction = 1,
    ...
) {
  ggplot2::discrete_scale(
    "colour", "tropical",
    tropical_pal(primary, other, direction),
    ...
  )
}

scale_color_tropical <- scale_colour_tropical

# example of usage
ggplot(mtcars, 
       aes(x = wt, 
           y = mpg, 
           color = factor(carb))) + 
  geom_line(size = 2)+
  scale_colour_tropical()+
  theme_minimal()


# Primary colour on my personal website, i use hugo template and here is the structure
primary = "#26867c" # changes colors of links
"#084c61"
"#3c6e71"
# Menu
menu_primary = "#122140"
menu_text = "rgba(255,255,255,0.6)" # "#34495e" # main text in menu
menu_text_active = "#e58500" # social media icons
menu_title = "#fff" # text in the menu bar

# Home sections
background = "rgb(255, 255, 255)" # page background
home_section_odd = "rgb(255, 255, 255)" # white
home_section_even = "#fafafb" # light gray

# End of script
