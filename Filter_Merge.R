# NGS_Q_Filter_Merge.R
# Filter Illumina NGS sequencing reads for quality and merge paired-end reads

# Load required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

required_packages <- c("ggplot2", "tidyverse", "ggrepel", "ShortRead", "dada2", "Biostrings")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("ShortRead", "dada2", "Biostrings")) {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
  library(pkg, character.only = TRUE)
}

# Set path to raw FASTQ files
path <- "raw_fastq/"

# Create output directories if they don't exist
filtered_dir <- file.path(path, "filtered")
if (!dir.exists(filtered_dir)) dir.create(filtered_dir)

plot_dir <- "plots"
if (!dir.exists(plot_dir)) dir.create(plot_dir)

# Get forward and reverse reads
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz$", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz$", full.names = TRUE))

# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Save quality plots
plot_file_F <- file.path(plot_dir, "Forward_Read_Quality.png")
plot_file_R <- file.path(plot_dir, "Reverse_Read_Quality.png")

png(filename = plot_file_F, width = 1000, height = 600)
plotQualityProfile(fnFs[1:min(2, length(fnFs))])
dev.off()

png(filename = plot_file_R, width = 1000, height = 600)
plotQualityProfile(fnRs[1:min(2, length(fnRs))])
dev.off()

# Define file paths for filtered output
filtFs <- file.path(filtered_dir, paste0(sample.names, "_R1_filtered.fastq.gz"))
filtRs <- file.path(filtered_dir, paste0(sample.names, "_R2_filtered.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filter and trim reads
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     trimLeft = c(10, 10),
                     trimRight = c(10, 10),
                     maxN = 0,
                     maxEE = c(2, 4),
                     truncQ = 5,
                     rm.phix = TRUE,
                     compress = TRUE,
                     multithread = TRUE,
                     minQ = 5)

print("Filtering Summary:")
print(head(out))
