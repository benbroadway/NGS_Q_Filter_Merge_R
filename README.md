# NGS_Q_Filter_Merge_R
R script to quality filter Illumina NGS sequencing reads and then trim and merge paired-end reads.

This script helps to:

- Visualize read quality.
- Remove poor-quality sequences.
- Generate cleaned FASTQ files for further processing.

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
