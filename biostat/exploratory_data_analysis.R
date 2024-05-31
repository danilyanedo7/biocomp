# Load necessary libraries
library(tidyverse)
library(viridis)

# Load diamonds dataset
data(diamonds)

# Display the first few rows of the dataset
head(diamonds)

# Summary statistics of the dataset
summary(diamonds)

# Dimensions of the dataset (number of rows and columns)
dim(diamonds)

# Check for missing values in each column
sapply(diamonds, function(x) sum(is.na(x)))

# Histogram of price values
ggplot(data = diamonds, aes(x = price)) +
  geom_histogram(fill = "#12719e", color = "#12719e", binwidth = 2) +
  ggtitle("Histogram of Price Values")+
  theme_minimal()

# Scatter plot of carat vs. price, colored by cut
ggplot(data = diamonds, aes(x = carat, y = price, color = cut)) + 
  geom_point()+
  scale_color_viridis_d()+
  theme_minimal()

# Boxplot of price by cut
ggplot(data = diamonds, aes(x = cut, y = price)) + 
  geom_boxplot(fill = "#12719e")+
  theme_minimal()

# Scatter plot of carat vs. depth
ggplot(data = diamonds, aes(x = carat, y = depth)) + 
  geom_point()+
  theme_minimal()

# Faceted scatter plot of price vs. carat, faceted by cut
ggplot(data = diamonds, aes(x = carat, y = price)) + 
  geom_point() + 
  facet_wrap(~cut)+
  theme_minimal()

# Correlation matrix plot
correlation_matrix <- cor(diamonds[, c("carat", "depth", "table", "price")])
corrplot::corrplot(correlation_matrix, method = "color")

# Pairwise scatter plots of selected numeric variables
pairs(diamonds[, c("carat", "depth", "table", "price")])
