BiocManager::install("msa")
library(msa)

system.file("tex", "texshade.sty", package="msa")
# Load example sequences
my_sequence_file <- system.file("examples", "exampleAA.fasta", package="msa")
my_sequences <- readAAStringSet(my_sequence_file)
my_sequences

# Perform multiple sequence alignment using default substitution matrix
my_first_alignment <- msa(my_sequences)
my_first_alignment

# Print the complete alignment
print(my_first_alignment, show="complete")

# Pretty print the alignment to a PDF
msaPrettyPrint(my_first_alignment, output="pdf", showNames="none",
               showLogo="none", askForOverwrite=FALSE, verbose=FALSE)

# Pretty print a subset of the alignment to a specific format
msaPrettyPrint(my_first_alignment, output="asis", y=c(164, 213),
               subset=c(1:6), showNames="none", showLogo="none",
               consensusColor="ColdHot", showLegend=FALSE,
               askForOverwrite=FALSE)

# Perform multiple sequence alignment using ClustalW
my_clustalw_alignment <- msa(my_sequences, "ClustalW")
my_clustalw_alignment

# Perform multiple sequence alignment using ClustalOmega
my_clustalomega_alignment <- msa(my_sequences, "ClustalOmega")
my_clustalomega_alignment

# Perform multiple sequence alignment using MUSCLE
my_muscle_alignment <- msa(my_sequences, "Muscle")
my_muscle_alignment

# Read Hemoglobin sequences and perform alignment
hemo_seq <- readAAStringSet(system.file("examples/HemoglobinAA.fasta",
                                        package="msa"))
hemo_aln <- msa(hemo_seq)

# Convert alignment to seqinr format
hemo_aln2 <- msaConvert(hemo_aln, type="seqinr::alignment")
library(seqinr)

# Calculate distance matrix for alignment
dist_matrix <- dist.alignment(hemo_aln2, "identity")
as.matrix(dist_matrix)[2:5, "HBA1_Homo_sapiens", drop=FALSE]

# Create and plot phylogenetic tree
library(ape)
hemo_tree <- nj(dist_matrix)
plot(hemo_tree, main="Phylogenetic Tree of Hemoglobin Alpha Sequences")
