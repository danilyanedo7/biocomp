## Phylogenetic analysis from txt file retrieved from GenBank
## By: Edo Danilyan
## Website: edodanilyan.com


# Setup : install packages below by removing the "#" symbol

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

# Load the necessary libraries
library(seqinr)
library(adegenet)
library(ape)
library(ggtree)
library(DECIPHER)
library(viridis)
library(ggplot2)
library(here)
library(plotly)
library(phytools)
library(phangorn)

# Check the working directory
here()

# Import sequence data from the drive to R
# Change "DNA" to "RNA" or "AA" if needed (skipping this part for now)
seqs <- readDNAStringSet("./data/acetobacter.txt", format = "fasta")

# Preview the sequences
seqs

# Note: Sequences must have the same orientation
# If not, the sequences can be reoriented (optional)
seqs <- OrientNucleotides(seqs)

# To construct a phylogenetic tree, we need to align the sequences first
aligned <- AlignSeqs(seqs)

# View the alignment results using your default web browser
BrowseSeqs(aligned, highlight = 0)

# After aligning the sequences, save the result as a single fasta file
writeXStringSet(aligned, file = "acetobacter_aligned.fasta")

# Read the aligned sequences
dna <- read.alignment("acetobacter_aligned.fasta", format = "fasta")

# Create a distance matrix for the aligned file
D <- dist.alignment(dna, matrix = "similarity")
D

# Create a heatmap and tree
temp <- as.data.frame(as.matrix(D))
temp
table.paint(temp,
            cleg = 0,
            clabel.row = .5,
            clabel.col = .5) +
  scale_color_viridis() # darker shades of gray = larger distance
heatmap(as.matrix(D))

# We can observe a pattern as the data is arranged by year
# However, we can't reach a conclusion yet

# Cintruct the tree using NJ method
tre <- nj(D)
class(tre) # all trees created with {ape} package will be of class phylo

tre <- ladderize(tre)

# To plot the tree, we can use base R or ggtree
## Base R plots

plot(tre, cex = 0.6)
title("Acetobacter")

# or
h_cluster <- hclust(D, method = "average", members = NULL) # method = average is used for UPGMA, members can be NULL or a vector with a length of size D
plot(h_cluster, cex = 0.6)

# Tree Plotting in ggtree

# There are so many layout that you can use

# You can fan it out
ggtree(tre, yscale = "NA") +
  geom_tiplab(hjust = -0.3,
              size = 4,
              align = TRUE) +
  xlim(0, 0.5)

# or try this layout
ggtree(tre, layout = "daylight") +
  geom_tiplab(hjust = -0.3,
              size = 4,
              align = TRUE) +
  xlim(0, 0.5)

# or these
ggtree(tre, layout = "roundrect")
ggtree(tre, layout = "slanted")
ggtree(tre, layout = "ellipse")
ggtree(tre, layout = "circular")
ggtree(tre, layout = "fan", open.angle = 120)
ggtree(tre, layout = "equal_angle")
ggtree(tre, layout = "daylight")
ggtree(tre, branch.length = 'none')
ggtree(tre, layout = "ellipse", branch.length = "none")
ggtree(tre, branch.length = 'none', layout = 'circular')
ggtree(tre, layout = "daylight", branch.length = 'none')
ggtree(tre, layout = "rectangular") # Did anything change? This is the default

# Plot a basic tree
treeplot <- ggtree(tre) +
  geom_tiplab(hjust = -0.3,
              size = 4,
              align = TRUE) +
  xlim(0, 0.5)
treeplot

# Setting branch.length = "none" will create a cladogram
ggtree(tre, layout = "rectangular", branch.length = "none")

# Adding a midpoint root (part of phytools package)
ggtree(midpoint.root(tre))

# View plot interactively, this comes in handy when you need to highlight certain node or branch
ggplotly(treeplot)

# Customize the trees
# Plot using ggtree and highlight clusters
# Change the node values for your own data
ggtree(tre) +
  geom_tiplab(hjust = -0.3,
              size = 4,
              align = TRUE) +
  geom_hilight(node = 12,
               fill = "purple",
               alpha = 0.2) +
  geom_hilight(node = 15,
               fill = "firebrick",
               alpha = 0.2) +
  geom_hilight(node = 13,
               fill = "yellow",
               alpha = 0.2) +
  xlim(0, 0.5)

# Adding tip labels
ggtree(tre) +
  geom_treescale(x = 0, y = 0, # x and y position of the treescale
                 width = 0.01) + # width of the scale
  geom_tiplab(size = 4, geom = 'text') + # displaying tip labels
  coord_cartesian(clip = 'off') + # allows drawing outside the plot
  theme(plot.margin = margin(1, 2, 1, 1, "cm")) # add space around the plot

# Highlight clusters and add a vertical line to group clusters
# Change the node values for your own data
p <- ggtree(tre) +
  geom_tiplab(size = 3, align = TRUE) +
  geom_hilight(node = 15,
               fill = "purple",
               alpha = 0.2) +
  geom_hilight(node = 17,
               fill = "green",
               alpha = 0.2) +
  geom_hilight(node = 18,
               fill = "gold",
               alpha = 0.2) +
  geom_cladelabel(
    node = 15,
    label = "Cluster 1",
    color = "purple",
    offset = .1,
    barsize = 2,
    fontsize = 5,
    align = TRUE,
    alpha = 0.5
  ) +
  geom_cladelabel(
    node = 17,
    label = "Cluster 2",
    color = "green",
    offset = .1,
    barsize = 2,
    fontsize = 5,
    align = TRUE,
    alpha = 0.5
  ) +
  geom_cladelabel(
    node = 20,
    label = "Cluster 3",
    color = "gold",
    offset = .1,
    barsize = 2,
    fontsize = 5,
    align = TRUE,
    alpha = 0.5
  ) +
  geom_nodepoint(color = "orange", alpha = 0.5)
p

# Plot the alignment with the tree side by side

# Let's plot the alignment with the tree. To do this, we first have to
# match the names to the tip labels
# Set our tree into a new name
tre.new <- tre
# Change tip labels to full alignment names
tre.new$tip.label <- aligned@ranges@NAMES

# Plot the alignment
dna_new <- as.DNAbin(dna)
dna_new

finalplot <- msaplot(p, dna_new, offset = 0.09, width = 2) +
  scale_fill_viridis_d(alpha = 0.8)
finalplot

# To save high quality image, you can use ggsave function as follow
ggsave(
  "./acetobacter.jpeg",
  finalplot,
  width = 9,
  height = 5,
  units = "in",
  dpi = 320
)