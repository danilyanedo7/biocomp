install.packages("gsalib")
BiocManager::install("rtracklayer")

# Load necessary libraries
library(gsalib)
library(GenomicRanges)
library(rtracklayer)

# Run GATK command to generate a report (example)
# Replace with actual paths and parameters
system("gatk HaplotypeCaller -R reference.fasta -I input.bam -O output.vcf -ERC GVCF -G Standard")

# Read the GATK report file
report <- read_gatk_report("output.grp")

# Explore the report
print(report)