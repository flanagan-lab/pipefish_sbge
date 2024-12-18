# Analysis of sex-biased gene expression in _Syngnathus fuscus_

Within this repository you will find code related to the transcriptome analysis of the Northern pipefish, _Syngnathus fuscus_. Sequencing data was generated for male and female pipefish from three different organs including the gills, the liver, and the gonads with the Illumina NovaSeq600. As there is no genome for this species, a _de novo_ transcriptome assembly was constructed and used for the subsequent analyses. We set out to answer several questions including:
  1. Are patterns of sex-biased and sex-specific gene expression different for reproductive versus somatic tissues?
  2. Do the functions of sex-biased or sex-specific genes differ between the two sexes?
  3. Does pleiotropy appear to constrain the decoupling of male- and female-specific gene expression? 

## Navigating this directory
### Docs
All Rmarkdown documents used for the various analyses are located in the directory docs/. For each `.Rmd` file there is a rendered `.md` file. They do the following things:

 - `Analyzing_fuscus_RNAseq_data_from_MSU.Rmd`: Works through the entire RNA-sequencing pipeline starting with the raw reads, moving through pre-assembly quality control, _de novo_ assembly generation, post-assembly quality control and thinning, and alignment and abundance estimations. The programs used for each step and how each program was used is highlighted in detail within the document, along with the results of the various steps. The end file generated from this document is a .RDS gene expression matrix that was then used for differential expression analysis.

  - `fuscus_diff_expr_analysis.Rmd`: Read in the .RDS file, perform some exploratory analysis, including a single-factor analysis that looked at MvF gene expression across all organs, and then a detailed multi-factor analysis that explored MvF expression levels within each organ. This file generates all of the datasets that were then used to create figures 1, 2, and 3 for the manuscript as well as results for the Gene ontology analysis.

  - `fuscus_tissue_specificity.Rmd`: Read in the quant.sf files generated from salmon to calculate the tissue specificity index tau ($tau$) for all genes. This file includes the filtering of genes and function that was used to calculate tau. Additionally, this file contains the comparisson analysis of tau vs sex-bias and generates the data that was used to create figure 4 in the manuscript.

### R
The directory R/ contains the script that was used to generate the BUSCO plot.

### Bash
The directory bash/ contains all of the bash scripts referenced to in the Rmarkdown documents. The order in which to run them, and details about their usage, are outlined in `docs/Analyzing_fuscus_RNAseq_data_from_MSU.Rmd` and `docs/fuscus_diff_expr_analysis.Rmd`.

### Data Availability
The raw RNA-sequencing reads and sample meta-data that corresponds to the analyses highlighted in this directory have been deposited in the NCBI Short Read Archive (BioProject PRJNA1193634).
