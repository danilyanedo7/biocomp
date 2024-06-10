if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rbwa")

# Set options for knitr graphics
options("knitr.graphics.auto_pdf" = TRUE, eval = TRUE)

# Load the Rbwa library
library(Rbwa)

# Build BWA index for the fasta file
dir <- tempdir()
fasta <- system.file(package = "Rbwa", "fasta/chr12.fa")
index_prefix <- file.path(dir, "chr12")
bwa_build_index(fasta, index_prefix = index_prefix)
list.files(dir)

# Perform BWA alignment
fastq <- system.file(package = "Rbwa", "fastq/sequences.fastq")
bwa_aln(
  index_prefix = index_prefix,
  fastq_files = fastq,
  sai_files = file.path(dir, "output.sai")
)

# Perform BWA alignment with specific parameters
bwa_aln(
  index_prefix = index_prefix,
  fastq_files = fastq,
  sai_files = file.path(dir, "output.sai"),
  n = 3,
  k = 3,
  l = 13
)

# Convert SAI file to SAM file
bwa_sam(
  index_prefix = index_prefix,
  fastq_files = fastq,
  sai_files = file.path(dir, "output.sai"),
  sam_file = file.path(dir, "output.sam")
)

# Read and trim SAM file for display
strtrim(readLines(file.path(dir, "output.sam")), 65)

# Convert SAM file to multi-SAM file
xa2multi(file.path(dir, "output.sam"),
         file.path(dir, "output.multi.sam"))
strtrim(readLines(file.path(dir, "output.multi.sam")), 65)

# Perform BWA-MEM alignment
fastq <- system.file(package = "Rbwa", "fastq/sequences.fastq")
bwa_mem(
  index_prefix = index_prefix,
  fastq_files = fastq,
  sam_file = file.path(dir, "output.sam")
)

# Read and trim SAM file for display
strtrim(readLines(file.path(dir, "output.sam")), 65)
