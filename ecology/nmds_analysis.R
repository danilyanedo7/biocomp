# Load required libraries
library(vegan)
library(tidyverse)

# Sample data
sample_names <- c("Sample1", "Sample2", "Sample3", "Sample4")
sample_type <- c("Type1", "Type2", "Type1", "Type2")
treatment <- c("TreatmentA", "TreatmentB", "TreatmentA", "TreatmentB")
otu_data <- matrix(
  data = sample(1:100, 16, replace = TRUE),
  nrow = 4,
  ncol = 4
)

# Creating the dataframe
df <- data.frame(sample_names, sample_type, treatment, otu_data)

# Selecting only the abundance columns
abundance_columns <- df[, 4:ncol(df)]

# Creating the OTU table
otu_table <- as.matrix(abundance_columns)

# Adding row names (OTU IDs)
rownames(otu_table) <- df$sample_names

# Displaying the OTU table
print(otu_table)

# NMDS Analysis
# NMDS (Non-metric Multidimensional Scaling) is a method commonly used in ecology
# to visualize similarities or dissimilarities between samples based on community
# composition data.
nmds <- metaMDS(otu_table, distance = "bray")

# Plotting the NMDS ordination
# NMDS plots provide insights into the clustering or dispersion of samples, indicating
# similarities or differences in community composition.
plot(nmds)

# Extracting NMDS values: NMDS scores (x and y coordinates)
# The positioning of samples in the NMDS plot can be interpreted in relation to
# environmental variables or treatments to identify patterns or trends.
nmds.scores <- as.data.frame(scores(nmds))

# Adding metadata to the NMDS scores
# This includes information about sample names, treatments, and sample types,
# which can help interpret the NMDS plot.
nmds.scores$sample_names <- df$sample_names
nmds.scores$treatment <- df$treatment
nmds.scores$sample_type <- df$sample_type

# Plotting the NMDS ordination using ggplot2
# By incorporating NMDS analysis into your ecological research workflow, you can gain
# valuable insights into community composition patterns and the factors driving those
# patterns.
ggplot(nmds.scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(size = 2, aes(shape = sample_type, colour = treatment)) +
  theme_minimal() +
  labs(
    x = "NMDS1",
    colour = "Treatment",
    y = "NMDS2",
    shape = "Type"
  ) +
  scale_colour_manual(values = c("#1696d2", "#fdbf11"))
