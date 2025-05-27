# NGS_to_sgRNA_R

R scripts to:
- Quality filter Illumina NGS sequencing reads and trim / merge paired-end reads.
- Map NGS sequencing reads to an sgRNA library.


This script helps to:

- Visualize read quality.
- Remove poor-quality sequences.
- Generate cleaned FASTQ files for further processing.
- Map NGS sequencing reads to specific sgRNAs.

## Dependencies

This script uses the following R packages:

- `ggplot2`
- `tidyverse`
- `ggrepel`
- `ShortRead` *(Bioconductor)*
- `dada2` *(Bioconductor)*
- `Biostrings` *(Bioconductor)*

To install missing packages, the script attempts automatic installation.

```r
# For Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("ShortRead", "dada2", "Biostrings"))
