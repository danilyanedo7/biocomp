# Main code source comes from DADA2 online vignette
# This tutorial is rewrote as my learning materials
# Visit https://benjjneb.github.io/dada2/tutorial.html for more information
# Context for Each Step:
# 1. Package Installation: Ensures all required packages are installed with dependencies.
# 2. Data Input and Quality Profiles: Reads input data and plots quality profiles.
# 3. Filtering and Trimming: Filters and trims the sequence data.
# 4. Error Learning: Learns error rates from the filtered sequences.
# 5. Denoising and Merging: Performs denoising and merging of paired-end reads.
# 6. Sequence Table Construction and Chimera Removal: Constructs the sequence table and removes chimeric sequences.
# 7. Taxonomic Assignment: Assigns taxonomy to sequence variants using the SILVA reference database.
# 8. Species-Level Assignment: Refines taxonomic assignments to the species level.
# 9. DECIPHER Taxonomic Identification: Uses DECIPHER for additional taxonomic identification.
# 10. Mock Community Evaluation: Evaluates the inferred sequences against a mock community.
# 11. Phyloseq Analysis: Uses phyloseq for downstream analysis and visualization.


# Preparations
# Install and load necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install packages with dependencies
BiocManager::install("dada2", dependencies = TRUE, force = TRUE)
BiocManager::install("DECIPHER", dependencies = TRUE, force = TRUE)
BiocManager::install("phyloseq", dependencies = TRUE, force = TRUE)
BiocManager::install("Biostrings", dependencies = TRUE, force = TRUE)

# Load libraries
library(dada2)
library(DECIPHER)
library(phyloseq)
library(Biostrings)
library(ggplot2)

# Set the path to your data and list the files
path <- "./data/MiSeq_SOP"
list.files(path)

# Define forward and reverse read files
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# Define file paths for filtered reads
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Run filterAndTrim and check the output
out <- filterAndTrim(
  fnFs,
  filtFs,
  fnRs,
  filtRs,
  truncLen = c(240, 160),
  maxN = 0,
  maxEE = c(2, 2),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
)
head(out)

# Learn error rates
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
plotErrors(errF, nominalQ = TRUE)

# Inferential analysis
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

# Merging paired-end reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

# Filtering sequences by length
seqtab2 <- seqtab[, nchar(colnames(seqtab)) %in% 250:256]

# Removing chimera
seqtab.nochim <- removeBimeraDenovo(seqtab,
                                    method = "consensus",
                                    multithread = TRUE,
                                    verbose = TRUE)
dim(seqtab.nochim)

# Calculate proportion of non-chimeric sequences
sum(seqtab.nochim) / sum(seqtab)

# Track reads through the pipeline
getN <- function(x)
  sum(getUniques(x))

track <- cbind(
  out,
  sapply(dadaFs, getN),
  sapply(dadaRs, getN),
  sapply(mergers, getN),
  rowSums(seqtab.nochim)
)

colnames(track) <- c("input",
                     "filtered",
                     "denoisedF",
                     "denoisedR",
                     "merged",
                     "nonchim")
rownames(track) <- sample.names
head(track)

# Taxonomic annotation
taxa <- assignTaxonomy(seqtab.nochim,
                       "./data/tax/silva_nr_v132_train_set.fa",
                       multithread = TRUE)

# Species-level assignment
taxa <- addSpecies(taxa, "./data/tax/silva_species_assignment_v132.fa")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# Using DECIPHER for taxonomic identification
dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
load("./data/tax/IDTaxa/SILVA_SSU_r132_March2018.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
ids <- IdTaxa(
  dna,
  trainingSet,
  strand = "top",
  processors = NULL,
  verbose = FALSE
) # use all processors

# Ranks of interest
ranks <- c("domain",
           "phylum",
           "class",
           "order",
           "family",
           "genus",
           "species")

# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks
rownames(taxid) <- getSequences(seqtab.nochim)

# Mock community evaluation
unqs.mock <- seqtab.nochim["Mock", ]
unqs.mock <- sort(unqs.mock[unqs.mock > 0], decreasing = TRUE) # Drop ASVs absent in the Mock
cat(
  "DADA2 inferred",
  length(unqs.mock),
  "sample sequences present in the Mock community.\n"
)
mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x)
  any(grepl(x, mock.ref))))
cat("Of those,",
    sum(match.ref),
    "were exact matches to the expected reference sequences.\n")

# Phyloseq analysis
theme_set(theme_minimal())

samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject, 1, 1)
subject <- substr(subject, 2, 999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject = subject,
                    Gender = gender,
                    Day = day)
samdf$When <- "Early"
samdf$When[samdf$Day > 100] <- "Late"
rownames(samdf) <- samples.out

ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  sample_data(samdf),
  tax_table(taxa)
)
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

# Plot richness
plot_richness(
  ps,
  x = "Day",
  measures = c("Shannon", "Simpson"),
  color = "When"
)

# Bray NMDS ordination
ps.prop <- transform_sample_counts(ps, function(otu)
  otu / sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method = "NMDS", distance = "bray")
plot_ordination(ps.prop, ord.nmds.bray, color = "When", title = "Bray NMDS")

# Plot top 20 taxa
top20 <- names(sort(taxa_sums(ps), decreasing = TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU)
  OTU / sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x = "Day", fill = "Family") + facet_wrap( ~ When, scales =
                                                               "free_x")


#### If you happened to use these package, dont forget to cite their work. To create citation text, use citation()

citation("phyloseq")
citation("dada2")
citation("DECIPHER")
citation("Biostrings")
citation("ggplot2")

