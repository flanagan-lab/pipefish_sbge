# _Stigmatopora nigra_ transcriptome analysis

This directory documents the analysis of gene expression and tissue specificity in sexually mature male and female _Stigmatopora nigra_.

## Navigating this directory

The `bash_scripts` directory contains the scripts used to QC and assemble the _S. nigra_ transcriptome. The order in which to run them, and details about their usage, are outlined in `docs/SNI_Analysing_RNAseq_data.Rmd`.

The `data` directory contains data that was used in the analyses outlined in `docs/SMA_DE_analysis.Rmd`. This data includes:
  * `txi.salmon_SN.RDS`: This is the gene expression matrix that is generated at the end of the RNAseq pipeline and used for differential expression analysis with DESeq2. This file was generated from the salmon `quant.sf` files with the package tximport.
  * `expression_files`: This directory contains all of the `quant.sf` files that were generated by salmon for each sample.

The `docs` directory contains all Rmarkdown documents used in the various analyses. These documents include:
  * `SNI_Analysing_RNAseq_data.Rmd`: In this .Rmd, we work through the entire RNA-sequencing pipeline starting with raw reads, moving through pre-assembly quality control, _de novo_ assembly generation, post-assembly quality control and thinning, and alignment and abundance estimations. The programs used for each step and how each program was used is highlighted in detail within the document,along with the results of the various steps. The end file generated from this document is a .RDS gene expression matrix that was then used for the differential expression analysis outline in `SNI_DE_analysis.Rmd`.
  * `SNI_DE_analysis.Rmd`: In this .Rmd, we read in the .RDS file, perform some exploratory analyses, including a single-factor analysis that looked at male/female gene expression across all organs, and then a detailed multi-factor analysis that explored male/female expression levels within each organ.
  * `SNI_tissue_specificity.Rmd`: In this .Rmd, we calculate the tissue specificity index tau for all genes. This file includes the filtering of genes and function that was used to calculate tau.

The `figs` directory contains all of the figures that were generated from `SNI_DE_analysis.Rmd` and `SNI_tissue_specificity.Rmd`.
