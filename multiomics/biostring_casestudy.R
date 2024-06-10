# To learn more about Biostrings, run code below
browseVignettes(package = "Biostrings")

# Step 1: Install Bioconductor and Biostrings
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pwalign", force = T)

# Step 2: Load the Biostrings Library
library(Biostrings)
library(pwalign)

# Step 3: Create DNA, RNA, and Protein Sequences
dna_seq <- DNAString("ACGTACGTGACG")
rna_seq <- RNAString("ACGUACGUGACG")
protein_seq <- AAString("MTEYKLVVVG")

# Print sequences
print(dna_seq)
print(rna_seq)
print(protein_seq)

# Step 4: Perform Basic Operations
# Calculate the length of a sequence
print(length(dna_seq))

# Find a pattern in the sequence
pattern_match <- matchPattern("ACGT", dna_seq)
print(pattern_match)

# Get the reverse complement of a DNA sequence
rev_comp_seq <- reverseComplement(dna_seq)
print(rev_comp_seq)

# Step 5: More Advanced Functions
# Pairwise alignment
align1 <- pairwiseAlignment(DNAString("ACGT"), DNAString("AGT"), type = "global")
print(align1)

# Finding mismatches
mismatch <- mismatchTable(pairwiseAlignment(DNAString("ACGT"), DNAString("AGT"), type = "local"))
print(mismatch)

# Step 6: Reading and Writing Sequence Files
# Reading a FASTA file
# Assuming you have a FASTA file named "example.fasta" in your working directory, uncomment lines below
# fasta_file <- "example.fasta"
# sequences <- readDNAStringSet(fasta_file)

# Writing sequences to a FASTA file
# writeXStringSet(sequences, "output_sequences.fasta")

# Example Workflow
# Create a DNA sequence, uncomment lines below
example_dna_seq <- DNAString("ACGTACGTGACG")

# Get the reverse complement
example_rev_comp_seq <- reverseComplement(example_dna_seq)

# Find a pattern
example_pattern_match <- matchPattern("ACGT", example_dna_seq)

# Display results
print(example_dna_seq)
print(example_rev_comp_seq)
print(example_pattern_match)


# Case study, you have a forward 16s of your bacteria sample
# Raw forward and reverse 16S rRNA sequence of a bacterial sample
reverse_16s <- DNAString("CCCTTCTGTCCACCTTAGGCGGCTGGCTCCAAAGGTTACCCCACCGACTTCGGGTGTTACAAACTCTCGTGGTGTGACGGGCGGTGTGTACAAGGCCCGGGAACGTATTCACCGCGGCATGCTGATCCGCGATTACTAGCGATTCCGGCTTCATGCAGGCGAGTTGCAGCCTGCAATCCGAACTGAGAATGGTTTTATGGGATTCGCTTAACCTCGCGGTCTCGCAGCCCTTTGTACCATCCATTGTAGCACGTGTGTAGCCCAGGTCATAAGGGGCATGATGATTTGACGTCATCCCCACCTTCCTCCGGTTTGTCACCGGCAGTCACCTTAGAGTGCCCAACTGAATGCTGGCAACTAAGATCAAGGGTTGCGCTCGTTGCGGGACTTAACCCAACATCTCACGACACGAGCTGACGACAACCATGCACCACCTGTCATCCTGTCCCCCGAAGGGGAACGCCCTATCTCTAGGGTTGTCAGGAGATGTCAAGACCTGGTAAGGTTCTTCGCGTTGCTTCGAATTAAACCACATGCTCCACCGCTTGTGCGGGCCCCCGTCAATTCCTTTGAGTTTCAGCCTTGCGGCCGTACTCCCCAGGCGGAGTGCTTAATGCGTTTGCTGCAGCACTAAAGGGCGGAAACCCTCTAACACTTAGCACTCATCGTTTACGGCGTGGACTACCAGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGCGCCTCAGCGTCAGTTACAGACCAAAGAGTCGCCTTCGCCACTGGTGTTCCTCCACATCTCTACGCATTTCACCGCTACACGTGGAATTCCACTCTTCTCTTCTGCACTCAAGTTCCCCAGTTTCCAATGACCCTCCCCGGTTGAGCCGGGGGCTTTCACATCAGACTTAAGGAACCGCCTGCGCGCGCTTTACGCCCAATAATTCCGGACAACGCTTGCCACCTACGTATTACCGCGGCTGCTGGCACGTAGTTAGCCGTGGCTTTCTGGTTAGGTACCGTCAAGGTACCGGCAGTTACTCCGGTACTTGTTCTTCCCTAACAACAGAGTTTTACGATCCGAAAACCTTCATCACTCACGCGGCGTTGCTCCGTCAGATCTTCGTCATTGCGGAAGATCCCTACTGCTGCCTCCGTAGGAATCTGGGCCGGGTCTNAGTCCAGTGTGGCGATCACCTNTCAGGTCGGTACCATCGTCGCTTGGGAACCGTACTTCCAACTAACTAATGCCC")

# Reverse 16S rRNA sequence of a bacterial sample
forward_16s <- DNAString("CCAAGCGCTGCTAATACATGCAAGTCGAGCGGACAGATGGGAGCTTGCTCCCTGAAGTCAGCGGCGGACGGGTGAGTAACACGTGGGCAACCTGCCTGTAAGACTGGGATAACTCCGGGAAACCGGGGCTAATACCGGATAATTCTTTCCCTCACATGAGGGAAAGCTGAAAGATGGTTTCGGCTATCACTTACAGATGGGCCCGCGGCGCATTAGCTAGTTGGTGAGGTAACGGCTCACCAAGGCGACGATGCGTAGCCGACCTGAGAGGGTGATCGGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCGCAATGGACGAAAGTCTGACGGAGCAACGCCGCGTGAGTGATGAAGGTTTTCGGATCGTAAAACTCTGTTGTTAGGGAAGAACAAGTACCGGAGTAACTGCCGGTACCTTGACGGTACCTAACCAGAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGCGCGCGCAGGCGGTTCCTTAAGTCTGATGTGAAAGCCCCCGGCTCAACCGGGGAGGGTCATTGGAAACTGGGGAACTTGAGTGCAGAAGAGAAGAGTGGAATTCCACGTGTAGCGGTGAAATGCGTAGAGATGTGGAGGAACACCAGTGGCGAAGGCGACTCTTTGGTCTGTAACTGACGCTGAGGCGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAGTGCTAAGTGTTAGAGGGTTTCCGCCCTTTAGTGCTGCAGCAAACGCATTAAGCACTCCGCCTGGGGAGTACGGCCGCAAGGCTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCTCCTGACAACCCTAGAGATAGGGCGTTCCCCTTCGGGGGACAGGATGACAGGTGGTGCATGGTTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCACCCTTGATCTTAGTTGCCAGCATTCAGTTGGGCACTCTAAGGTGACTGCCGTGACAACCGGAGGAAGGGGGGATGACGTCAATCATCATGCCCTTTGACTGGGCTACCACGTGCTACATGGATGGTACAAGGCTGCAGACGCCAGGTAACNAATCCCTAAACCTTTNCATTCGGATGCAGG")

# Remove sequences containing 'N'
alignment <- pwalign::pairwiseAlignment(forward_16s, reverse_16s, type = "global")

# Get aligned sequences
aligned_forward <- as.character(pattern(alignment))
aligned_reverse <- as.character(subject(alignment))

# Function to replace 'N' with the complementary base from the other sequence
replace_N_with_complement <- function(aligned_forward, aligned_reverse) {
  f_with_complement <- character(length(aligned_forward))
  r_with_complement <- character(length(aligned_reverse))
  
  for (i in seq_along(aligned_forward)) {
    if (aligned_forward[i] == "N") {
      f_with_complement[i] <- aligned_reverse[i]
    } else {
      f_with_complement[i] <- aligned_forward[i]
    }
    
    if (aligned_reverse[i] == "N") {
      r_with_complement[i] <- aligned_forward[i]
    } else {
      r_with_complement[i] <- aligned_reverse[i]
    }
  }
  
  return(list(
    forward = DNAString(paste(f_with_complement, collapse = "")),
    reverse = DNAString(paste(r_with_complement, collapse = ""))
  ))
}

# Apply the function to replace 'N' with the complementary base
aligned_sequences <- replace_N_with_complement(unlist(strsplit(aligned_forward, "")), unlist(strsplit(aligned_reverse, "")))

# Print cleaned sequences
print(aligned_sequences)
print(aligned_sequences$forward)
print(aligned_sequences$reverse)

# Convert cleaned sequences to XStringSet
cleaned_sequences <- DNAStringSet(c(aligned_sequences$forward, aligned_sequences$reverse_complement))

# Use BrowseSeqs to visualize the sequences
BrowseSeqs(cleaned_sequences)


alignment <- DNAStringSet(c(alignment$forward, alignment$reverse))
BrowseSeqs(alignment)
