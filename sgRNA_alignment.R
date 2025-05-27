# sgRNA_alignment.R

# 1. Load or install required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

required <- c("Biostrings", "Rsubread", "GenomicAlignments")
for (pkg in required) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("Biostrings", "GenomicAlignments")) {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
  library(pkg, character.only = TRUE)
}

# 2. Load and format sgRNA library
lib_df <- read.csv("X.csv", stringsAsFactors = FALSE)

# Extract 19nt sgRNAs (starting from 2nd base)
sgRNA_seqs <- substr(lib_df$Target.sequence, 2, 20)  # 2 to 20 for 19nt
sgRNA_dna <- DNAStringSet(sgRNA_seqs)
names(sgRNA_dna) <- lib_df$Target.sequence

# Save as FASTA
writeXStringSet(sgRNA_dna, filepath = "sgRNA.fa")

# 3. Build sgRNA index for Rsubread alignment
buildindex(basename = "sgRNA", reference = "sgRNA.fa", indexSplit = FALSE, gappedIndex = TRUE)

# 4. Align filtered paired-end reads
input_R1 <- "raw_fastq/filtered/S6-T1-A6-AG-In-2S6-T1-A6-AG-In-2_R1_001_F_filt.fastq.gz"
input_R2 <- "raw_fastq/filtered/S6-T1-A6-AG-In-2S6-T1-A6-AG-In-2_R2_001_R_filt.fastq.gz"
output_bam <- "aligned_reads.bam"

align(index = "sgRNA",
      readfile1 = input_R1,
      readfile2 = input_R2,
      output_file = output_bam,
      nthreads = 4,
      unique = TRUE,
      nBestLocations = 1,
      type = "DNA",
      maxMismatches = 0,
      indels = 0,
      TH1 = 1,
      TH2 = 1,
      minFragLength = 17,
      maxFragLength = 50,
      nTrim5 = 6,
      nTrim3 = 6,
      nsubreads = 16)

# 5. Read alignment and generate sgRNA counts
galn <- readGAlignmentPairs(output_bam)
counts_df <- data.frame(table(seqnames(galn)), row.names = "Var1")
colnames(counts_df) <- "Freq"

# 6. Normalize sgRNA counts (Doench-style)
counts_df$Freqovertotal <- counts_df$Freq / sum(counts_df$Freq)
counts_df$Freqovertotalx1M <- counts_df$Freqovertotal * 1e6
counts_df$Freqlog2 <- log2(1 + counts_df$Freqovertotalx1M)  # avoid log(0)

# 7. Plot histograms
hist(counts_df$Freq, main = "Raw Read Counts per sgRNA", xlab = "Count")
hist(counts_df$Freqlog2, main = "log2(Normalized sgRNA Counts)", xlab = "log2(count + 1)")

# 8. Save output
write.csv(counts_df, "sgRNA_counts_normalized.csv")
