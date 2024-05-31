# Install and load required packages
install.packages("vegan")  # Package for ecological analysis
library(vegan)             # Load vegan package for community ecology functions
library(ggplot2)           # Load ggplot2 for plotting capabilities

# Create sample data, manually input in R, or you can import from excel using 'writexl' package

data <- data.frame(
  Site = c("A", "B", "C"),
  Species1 = c(10, 5, 8),
  Species2 = c(7, 3, 5),
  Species3 = c(3, 2, 1),
  Species4 = c(6, 4, 2),
  Species5 = c(8, 1, 0)
)

# Extract site names
site_names <- data$Site

# Extract species data
species_data <- data[, -1]

# Calculate Shannon Diversity Index
shannon_index <- diversity(species_data, index = "shannon")
cat("Shannon Diversity Index:", shannon_index, "\n")

# Calculate Simpson Diversity Index
simpson_index <- diversity(species_data, index = "simpson")
cat("Simpson Diversity Index:", simpson_index, "\n")

# Calculate Species Richness
richness <- specnumber(species_data)
cat("Species Richness:", richness, "\n")

# Create a data frame to store diversity indices
indices_df <- data.frame(
  Site = rep(site_names, 3),
  Index = rep(c("Shannon", "Simpson", "Richness"), each = length(site_names)),
  Value = c(shannon_index, simpson_index, richness)
)

# Print the indices data frame
print(indices_df)


# Using biosampleR to calculate indices such as shannon, simpson, chao1, chao_diff
library(biosampleR)
data <- data.frame(
  Site = c(1, 2, 3), #sites numbering
  Species_A = c(10, 5, 8),
  Species_B = c(15, 9, 6),
  Species_C = c(7, 11, 14),
  Species_D = c(3, 6, 4)
)

biodiv_index <- calc_diversity_indices(data)
biodiv_index
