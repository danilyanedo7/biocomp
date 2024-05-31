## Phylogenetic analysis from numerical excel file
## By: Edo Danilyan
## Website: edodanilyan.com
# ggtree book resource: https://yulab-smu.top/treedata-book/

# Setup: install packages below by removing the "#" symbol
# Setup : install packages di bawah ini dengan menghapus tanda "#"

# install.packages("tidyverse")
# install.packages("here")
# install.packages("remotes")
# remotes::install_github("GuangchuangYu/treeio")
# install.packages("BiocManager")
# BiocManager::install("ggtree")
# BiocManager::install("DECIPHER")
# install.packages("ape")
# install.packages("plotly")
# install.packages("phytools")

# Load required packages
library(tidyverse, warn.conflicts = FALSE) # for data manipulation and visualization
library(seqinr) # for biological sequence analysis
library(adegenet) # for genetic data analysis
library(ape) # for phylogenetic analysis
library(ggtree) # for visualization of phylogenetic trees
library(DECIPHER) # for sequence alignment and analysis
library(viridis) # for color scales in visualization
library(here) # for file path management
library(plotly) # for interactive plots
library(phytools) # for phylogenetic analysis and visualization
library(phangorn) # for phylogenetic analysis

# Read in example data from an Excel file
contoh_data <- readxl::read_excel("./data/contoh_data.xlsx")
here()
# View the data in a pop-up window
View(contoh_data)

# Convert the first column to row names
datamikro <- contoh_data |> column_to_rownames(var = "0")

# Calculate a binary distance matrix based on the data
dist_mikro <- dist.binary(datamikro, method = 1)

## Construct a neighbor joining tree based on the distance matrix
mikro_nj <- nj(dist_mikro)

# Ladderize the tree to improve visualization
mikro_nj <- ladderize(mikro_nj)

# Visualize the tree using various layouts
ggtree(mikro_nj, layout = "rectangular") # rectangular layout (default)
ggtree(mikro_nj, layout = "roundrect") # round rectangular layout
ggtree(mikro_nj, layout = "slanted") # slanted layout
ggtree(mikro_nj, layout = "ellipse") # ellipse layout
ggtree(mikro_nj, layout = "circular") # circular layout
ggtree(mikro_nj, layout = "fan", open.angle = 40) # fan layout
ggtree(mikro_nj, layout = "equal_angle") # equal angle layout
ggtree(mikro_nj, layout = "daylight") # daylight layout
ggtree(mikro_nj, branch.length = 'none') # tree without branch lengths
ggtree(mikro_nj, layout = "ellipse", branch.length = "none") # ellipse layout without branch lengths

# Construct an UPGMA tree based on the distance matrix
mikro_upgma <- upgma(dist_mikro)

# Add graphical elements to the tree visualization
tree_upgma <- ggtree(mikro_upgma, layout = "rectangular") +
  geom_tiplab(size = 3, align = TRUE) + # add labels to the tips of the tree
  geom_nodepoint(color = "green", alpha = 0.5) + # add nodes as orange points
  theme_tree()
tree_upgma

# Correlation-cophenetic analysis
hc_upgma <- hclust(dist_mikro, method = "average", members = NULL) # "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
d2 <- cophenetic(hc_upgma)
cor(dist_mikro, d2) # The accepted correlation index is ≥0.7 or ≥70%. The lower the correlation value, the less the dendrogram represents the similarity matrix

# To save high quality image, you can use ggsave function as follow
ggsave(
  "./tree_upgma.jpeg",
  tree_upgma,
  width = 7,
  height = 4,
  units = "in",
  dpi = 400
)
