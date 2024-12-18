Tissue Specificity in *Syngnathus fuscus*
================
Coley Tosto
2024-12-19



- [Differential Expression Analysis](#differential-expression-analysis)
- [Calculating Tissue Specificity](#calculating-tissue-specificity)
  - [Pulling out TPM values](#pulling-out-tpm-values)
  - [Filtering the TPM dataset](#filtering-the-tpm-dataset)
  - [Generating $\tau$](#generating-tau)
  - [Validating the $\tau$
    calculations](#validating-the-tau-calculations)
- [Sex Bias and Tissue Specificity](#sex-bias-and-tissue-specificity)
  - [$\tau$ and logFC in sex-biased genes for each
    organ](#tau-and-logfc-in-sex-biased-genes-for-each-organ)
  - [$\tau$ and logFC for all the signficantly expressed genes in each
    organ](#tau-and-logfc-for-all-the-signficantly-expressed-genes-in-each-organ)
  - [Investigating genes with high tau and high
    logFC\`](#investigating-genes-with-high-tau-and-high-logfc)
    - [Gonads](#gonads)
    - [Liver](#liver)
    - [Gills](#gills)
  - [Categories of sex bias vs $\tau$](#categories-of-sex-bias-vs-tau)
    - [Within the sex-biased genes](#within-the-sex-biased-genes)
    - [Across all significantly expressed
      genes](#across-all-significantly-expressed-genes)

``` r
#This is a cohesive list of all the libraries used in this document
library(DESeq2)
library(spfTools)
library(ggplot2)
library(EnvStats)
library(dplyr)

library(magick)
library(patchwork)
library(pdftools)
library(tidyverse)
```

``` r
#The abundance matrix generated via salmon and tximport to be used for the DE analysis
txi.salmon <- readRDS("data/txi.salmon_FU.RDS")

#The samples file generated for tximport
samples <- read.table("FU_samples.txt", header = TRUE)

samples$group <- factor(paste0(samples$Sex, samples$Organ))

#Changing "Gonad" to be more specific to testis or ovaries
samples$Organ <- ifelse(samples$Sex == "F" & samples$Organ =="Gonad",
                        paste0("Ovaries"),
                        ifelse(samples$Sex == "M" & samples$Organ == "Gonad",
                               paste0("Testis"),
                               paste0(samples$Organ))
                        )

#Make sure the conditions are in the samples file as a factor
samples$Sex <- as.factor(samples$Sex)
samples$Organ <- as.factor(samples$Organ)

#Format colData to be used in the tau function
colData <- as.data.frame(samples)
rownames(colData) <- samples$ID
```

# Differential Expression Analysis

The main goal of this analysis is to relate measure of pleiotropic
constraints (tissue specificity) to sex bias. In order to do that I need
to combine my sex-bias data that I generated in a separate .Rmd
(`fuscus_diff_expr_analysis.Rmd`) with the tau values I will be
generating in this document.

Here I am reading in the sex-bias datasets that correspond to each of
the organ’s single-factor analysis. These datasets include all of the
genes where the p-value was not equal to “NA”.

``` r
liver_res <- read.csv("data/liver_res.csv", 
                      row.names = 1, 
                      header = TRUE)
gill_res <- read.csv("data/gill_res.csv", 
                      row.names = 1, 
                      header = TRUE)
gonad_res <- read.csv("data/gonad_res.csv", 
                      row.names = 1, 
                      header = TRUE)
```

# Calculating Tissue Specificity

To estimate tissue specificity, the TPM estimates are needed which is
the number of transcripts from a particular gene normalized first by
gene length, and then by sequencing depth (in millions) in the sample.
The output quant.sf files from salmon have the following format:

    Name | Length | EffectiveLength | TPM | NumReads

## Pulling out TPM values

From the salmon outputs I pulled out the TPM values for each sample. All
of the necessary .sf files are stored in a folder called
`fuscus_expression_files`.

The function below is pulling out the TPMs for each of the samples into
a column and then binding them together into one dataset. I then remove
the three samples that were not included in the differential expression
analysis.

``` r
#Get the list of file names/paths for all of the quant.sf files
files <- list.files(pattern = ".sf", path = "data/fuscus_expression_files", 
                    full.names = TRUE)

#For each quant.sf file pull out the TPM column
tpms <- do.call(cbind, lapply(files, function(file){

  dat <- read.delim(file, row.names = "Name")
  tpm <- dat["TPM"]
  colnames(tpm) <- gsub("data/fuscus_expression_files/(.*)_quant.sf","\\1",file)
  
  return(tpm)
}))

#Remove the samples that were removed in the diff. expression analysis
tpms <- tpms[, !(colnames(tpms) %in% c("FUT11M4",
                                       "FUG11M4",
                                       "FUL11M4"))]
```

## Filtering the TPM dataset

I know want to filter the TPM dataset to remove a few things. These
include similar things that I filtered out of the DESeq2 datasets prior
to performing the differential expression analysis. This includes:

1.  Keeping only the rows that weren’t filtered out in the DESeq2
    dataset due to low counts (rowSum $\le$ 10).

2.  Removing the rows that corresponded to genes that were “outliers” in
    the DESeq2 analysis.

3.  Removing the samples related to individual 11M4 that were removed in
    the differential expression analysis (see above).

This will be imperfect filtering in this case because I did an
individual single-factor analysis within each organ type. To try and
make life easier I am going to run the multi-factor analysis on the
DESeq dataset and then apply the filtering from that onto the TPM
dataset.

It should also be noted that the DESeq2 filtering was done based off of
the count data and not the TPMs, but as TPMs are just normalized counts
they should be correlated and something that was had low gene counts
should also have a low TPM and anything that was considered an outlier
in the counts could also be an outlier in terms of TPM, but it may not
be perfect.

``` r
#Running the MF diff. expr. analysis
#Create the DESeq dataset
dds_FU <- DESeqDataSetFromTximport(txi.salmon,
                                   colData = samples,
                                   design = ~ group)
```

    ## using counts and average transcript lengths from tximport

``` r
##Remove all 11M4 organs from the dataset
dds_FU <- dds_FU[, !(dds_FU$ID %in% c("FUT11M4", "FUG11M4", "FUL11M4"))]

##Filter the dataset, only keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds_FU)) >= 10 
dds_FU <- dds_FU[keep, ]

#Generate the expression values
dds_FU_exp <- DESeq(dds_FU)
```

    ## estimating size factors

    ## using 'avgTxLength' from assays(dds), correcting for library size

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

``` r
#Filtering the tpm dataset
#Only keeping the rows that weren't filtered out due to low counts
tpms <- tpms[rownames(tpms) %in% rownames(dds_FU), ]

#Pulling out the geneIDs for genes that were categorized as "outliers" by DESeq2
#Calculating the Cooks threshold that would have been used
np <- length(resultsNames(dds_FU_exp))
nsamp <- ncol(dds_FU_exp)
cooks_thresh <- qf(0.99, df1 = np, df2 = nsamp-np)

out_ids <- names(mcols(dds_FU_exp)$maxCooks[mcols(dds_FU_exp)$maxCooks > cooks_thresh])

#Removing the rows in the tpm dataset that were deemed "outliers" by DESeq2
tpms <- tpms[!(rownames(tpms) %in% out_ids), ]
```

## Generating $\tau$

These filtered TPM estimates can then be used to estimate $\tau$, a
tissue specificity estimator that can range from 0 to 1. $\tau$ is
calculated for each tissue, $i$, as follows:

$$
\tau=\frac{\sum_i{[1-ln(TPM_i)/ln(TPM_{max})]}}{N-1}
$$

where $TPM_{max}$ is the maximum **average** TPM for a given tissue
type, and $TPM_i$ is the **average** TPM for tissue $i$. If $\tau = 0$,
that gene is evenly expressed across tissues; if $\tau=1$, the gene is
expressed in an entirely tissue-specific fashion. Because TPM values
approach 0 are impacted by sampling stochasticity, any genes that had an
expression approaching 0 were set to $TPM=2$.

``` r
#Function for estimating tau given the TPM matrix and metadata file
est_tau <- function(geneDat,colDat){
  
  #For each row in the TPM matrix cbind it with the metadata file,
  #this attaches organ type information to the TPM values
  tissue_dat<-data.frame(cbind(colDat,
                               geneDat))
  
  #For the TPM values approaching 0, set them to 2
  tissue_dat$geneDat[tissue_dat$geneDat < 2] <- 2
  
  
  #Get the average TPM for each tissue type (TPMi)
  tissue_avgs <- tapply(tissue_dat$geneDat,tissue_dat$Organ,mean)
  
  #Get the maximum value from the average TPMS (TPMmax)
  tpmMax <- max(tissue_avgs, na.rm=TRUE)
    
  #IF running tau on JUST males of JUST females, this accounts for the
  #fact that ovary or testis will return an NA in the averaging
  if(length(unique(tissue_dat$Organ)) == 3){
    
    #calcualte tau
    tau <- sum(1-(log(tissue_avgs[unique(tissue_dat$Organ)])/log(tpmMax)))/
      (length(unique(tissue_dat$Organ))-1)
    
    #pull out appropriate organ
    if(tau == 0){
    
    organ <- paste0(NA)
    
    }else{
    
    organ <- names(tissue_avgs[tissue_avgs == tpmMax])
    
    }
  
  #Combine tau and the organ information together
  tau_organ <- data.frame(tau = tau, organ = organ)
  
  return(tau_organ)
    
  }
  
  #IF using the WHOLE dataset, calculate tau
  tau <- sum(1-(log(tissue_avgs)/log(tpmMax)))/(length(unique(tissue_dat$Organ))-1)
  
  #Pull out the information about WHICH tissue the expression is specific in
  ##If tau == 0 (not tissue specific) then return NA, but if above zero
  ##pull out the name of the organ with the highest expression
  if(tau == 0){
    
    organ <- paste0(NA)
    
  }else{
    
    organ <- names(tissue_avgs[tissue_avgs == tpmMax])
    
  }
  
  #Combine tau and the organ information together
  tau_organ <- data.frame(tau = tau, organ = organ)
  
  return(tau_organ)
}
```

I then applied that function across each row in the TPMs matrix with my
metadata stored in the object `colData`. The metadata file includes the
sample ID, Sex, and Organ type for every column present in the TPM
matrix.

I ran the $\tau$ function with the whole dataset, and then for just
males and just females.

``` r
colData <- colData[!(colData$ID %in% c("FUT11M4", "FUG11M4", "FUL11M4")),]

tau <- do.call(rbind, apply(tpms, 1, est_tau, colDat=colData))


tau_fem <- do.call(rbind, apply(tpms[,which(colData$Sex=="F")], 1, est_tau,
                                colData[which(colData$Sex=="F"),]))

tau_mal <- do.call(rbind, apply(tpms[,which(colData$Sex=="M")], 1, est_tau,
                                colData[which(colData$Sex=="M"),]))
```

<figure>
<img src="fuscus_tissue_specificity_files/figure-gfm/tauHist-1.png"
alt="Distribution of tissue specificity estimates for only the male samples (left) versus only the female samples (middle) versus the male samples and the female samples (right)." />
<figcaption aria-hidden="true">Distribution of tissue specificity
estimates for only the male samples (left) versus only the female
samples (middle) versus the male samples and the female samples
(right).</figcaption>
</figure>

We can see that the distribution of $\tau$ doesn’t vary drastically when
looking at only male/female samples versus when we included all of the
samples in the calculation. Because of this, I can be confident in using
the $\tau$ with all the samples moving forward.

## Validating the $\tau$ calculations

To make sure $\tau$ is being calculated in a way that makes sense, I am
checking some of the TPM values and plotting the corresponding counts
for the genes with a high tissue specificity index and a low tissue
specificity index.

I also want to make sure that the organ classification I implemented in
the tissue specificity function is accurate.

    ##                              tau   organ
    ## TRINITY_DN5507_c4_g1   0.9320775    Gill
    ## TRINITY_DN41134_c0_g1  0.9317676 Ovaries
    ## TRINITY_DN267909_c0_g1 0.9277302    Gill
    ## TRINITY_DN6930_c1_g1   0.9258133 Ovaries
    ## TRINITY_DN192690_c0_g1 0.9224122    Gill
    ## TRINITY_DN1186_c2_g1   0.9222358 Ovaries
    ## TRINITY_DN31372_c1_g1  0.9220840 Ovaries
    ## TRINITY_DN5203_c1_g1   0.9216494 Ovaries
    ## TRINITY_DN5905_c1_g2   0.9201969   Liver

    ##                       FUG10M2  FUG11F1  FUG11M2  FUG12M1  FUG13F1  FUG13F4
    ## TRINITY_DN5507_c4_g1 19777.27 25944.81 26820.53 30563.74 37748.11 37098.78
    ##                       FUG15M5   FUG2F2  FUG3F2  FUL10M2  FUL11F1  FUL11M2
    ## TRINITY_DN5507_c4_g1 17677.93 31498.43 16204.4 0.408263 0.462261 0.068661
    ##                       FUL12M1 FUL13F1 FUL13F4 FUL15M5 FUL2F2 FUL3F2 FUO11F1
    ## TRINITY_DN5507_c4_g1 0.192583       0       0       0      0      0       0
    ##                      FUO13F1  FUO13F4 FUO2F2 FUO3F2  FUT10M2  FUT11M2  FUT12M1
    ## TRINITY_DN5507_c4_g1       0 0.099184      0      0 1.014572 1.463632 1.859761
    ##                       FUT15M5
    ## TRINITY_DN5507_c4_g1 0.788966

    ##                       FUG10M2 FUG11F1 FUG11M2 FUG12M1 FUG13F1 FUG13F4 FUG15M5
    ## TRINITY_DN41134_c0_g1       0       0       0       0       0       0       0
    ##                       FUG2F2   FUG3F2 FUL10M2 FUL11F1 FUL11M2 FUL12M1 FUL13F1
    ## TRINITY_DN41134_c0_g1      0 0.117978       0       0       0       0       0
    ##                       FUL13F4 FUL15M5   FUL2F2  FUL3F2 FUO11F1  FUO13F1
    ## TRINITY_DN41134_c0_g1       0       0 1.439704 0.17896 32838.6 25016.87
    ##                        FUO13F4   FUO2F2   FUO3F2  FUT10M2  FUT11M2 FUT12M1
    ## TRINITY_DN41134_c0_g1 26196.31 20893.05 24119.07 0.065602 0.203414       0
    ##                        FUT15M5
    ## TRINITY_DN41134_c0_g1 0.826125

    ##                         FUG10M2  FUG11F1  FUG11M2  FUG12M1  FUG13F1  FUG13F4
    ## TRINITY_DN267909_c0_g1 15861.43 18979.82 17814.45 20319.39 27249.73 25673.48
    ##                         FUG15M5   FUG2F2   FUG3F2 FUL10M2 FUL11F1 FUL11M2
    ## TRINITY_DN267909_c0_g1 16809.72 24157.54 18606.81       0       0       0
    ##                        FUL12M1 FUL13F1 FUL13F4 FUL15M5 FUL2F2 FUL3F2 FUO11F1
    ## TRINITY_DN267909_c0_g1       0       0       0       0      0      0       0
    ##                        FUO13F1 FUO13F4 FUO2F2 FUO3F2 FUT10M2 FUT11M2  FUT12M1
    ## TRINITY_DN267909_c0_g1       0       0      0      0       0 1.32893 2.616365
    ##                         FUT15M5
    ## TRINITY_DN267909_c0_g1 1.416825

    ##                      FUG10M2 FUG11F1 FUG11M2 FUG12M1 FUG13F1 FUG13F4 FUG15M5
    ## TRINITY_DN6930_c1_g1       0       0       0       0       0       0       0
    ##                      FUG2F2 FUG3F2 FUL10M2 FUL11F1 FUL11M2 FUL12M1 FUL13F1
    ## TRINITY_DN6930_c1_g1      0      0       0       0       0       0       0
    ##                      FUL13F4 FUL15M5 FUL2F2 FUL3F2  FUO11F1  FUO13F1  FUO13F4
    ## TRINITY_DN6930_c1_g1       0       0      0      0 17354.72 8823.398 11856.55
    ##                       FUO2F2   FUO3F2 FUT10M2 FUT11M2 FUT12M1 FUT15M5
    ## TRINITY_DN6930_c1_g1 18735.4 339.0367       0       0       0       0

    ##                         FUG10M2  FUG11F1  FUG11M2  FUG12M1  FUG13F1  FUG13F4
    ## TRINITY_DN192690_c0_g1 23487.56 24211.85 894.8261 1008.716 1293.544 1019.839
    ##                         FUG15M5   FUG2F2   FUG3F2  FUL10M2 FUL11F1 FUL11M2
    ## TRINITY_DN192690_c0_g1 13587.53 1028.793 39327.73 0.317268       0       0
    ##                        FUL12M1 FUL13F1 FUL13F4 FUL15M5 FUL2F2 FUL3F2 FUO11F1
    ## TRINITY_DN192690_c0_g1       0       0       0       0      0      0       0
    ##                        FUO13F1 FUO13F4 FUO2F2 FUO3F2  FUT10M2  FUT11M2 FUT12M1
    ## TRINITY_DN192690_c0_g1       0       0      0      0 2.226906 0.296346       0
    ##                         FUT15M5
    ## TRINITY_DN192690_c0_g1 2.633663

    ##                       FUG10M2 FUG11F1 FUG11M2 FUG12M1 FUG13F1  FUG13F4 FUG15M5
    ## TRINITY_DN1186_c2_g1 0.219755       0       0       0       0 0.229783       0
    ##                      FUG2F2  FUG3F2 FUL10M2 FUL11F1 FUL11M2 FUL12M1 FUL13F1
    ## TRINITY_DN1186_c2_g1      0 0.12584       0       0       0       0       0
    ##                      FUL13F4 FUL15M5   FUL2F2 FUL3F2  FUO11F1 FUO13F1  FUO13F4
    ## TRINITY_DN1186_c2_g1       0       0 0.152474      0 13155.28 7254.93 5955.645
    ##                        FUO2F2   FUO3F2  FUT10M2  FUT11M2  FUT12M1  FUT15M5
    ## TRINITY_DN1186_c2_g1 4321.204 6469.144 0.109197 0.223647 0.103766 0.526084

    ##                       FUG10M2 FUG11F1 FUG11M2 FUG12M1 FUG13F1 FUG13F4 FUG15M5
    ## TRINITY_DN31372_c1_g1       0       0       0       0       0       0       0
    ##                       FUG2F2 FUG3F2 FUL10M2 FUL11F1 FUL11M2 FUL12M1 FUL13F1
    ## TRINITY_DN31372_c1_g1      0      0       0       0       0       0       0
    ##                       FUL13F4 FUL15M5 FUL2F2 FUL3F2  FUO11F1  FUO13F1  FUO13F4
    ## TRINITY_DN31372_c1_g1       0       0      0      0 8751.549 15466.85 5997.736
    ##                         FUO2F2   FUO3F2 FUT10M2 FUT11M2 FUT12M1 FUT15M5
    ## TRINITY_DN31372_c1_g1 5659.576 640.8449       0       0       0       0

    ##                       FUG10M2  FUG11F1 FUG11M2 FUG12M1  FUG13F1  FUG13F4
    ## TRINITY_DN5203_c1_g1 0.061455 0.045342       0       0 0.192052 0.022836
    ##                      FUG15M5   FUG2F2   FUG3F2 FUL10M2 FUL11F1  FUL11M2
    ## TRINITY_DN5203_c1_g1       0 0.103556 0.080138       0       0 0.681082
    ##                       FUL12M1  FUL13F1  FUL13F4  FUL15M5   FUL2F2 FUL3F2
    ## TRINITY_DN5203_c1_g1 0.018602 0.063119 0.110536 0.049911 0.230039      0
    ##                       FUO11F1  FUO13F1  FUO13F4   FUO2F2   FUO3F2 FUT10M2
    ## TRINITY_DN5203_c1_g1 9576.269 8729.695 5411.466 4883.741 6157.025       0
    ##                      FUT11M2 FUT12M1 FUT15M5
    ## TRINITY_DN5203_c1_g1       0       0       0

    ##                      FUG10M2 FUG11F1 FUG11M2  FUG12M1 FUG13F1 FUG13F4  FUG15M5
    ## TRINITY_DN5905_c1_g2       0       0       0 0.076948       0       0 0.041919
    ##                      FUG2F2 FUG3F2  FUL10M2 FUL11F1 FUL11M2  FUL12M1  FUL13F1
    ## TRINITY_DN5905_c1_g2      0      0 0.053721 0.07994       0 0.064351 7.515366
    ##                      FUL13F4 FUL15M5   FUL2F2   FUL3F2 FUO11F1 FUO13F1  FUO13F4
    ## TRINITY_DN5905_c1_g2 5348.33       0 47872.44 22.06734 0.06111 0.01846 0.014991
    ##                        FUO2F2  FUO3F2  FUT10M2  FUT11M2  FUT12M1  FUT15M5
    ## TRINITY_DN5905_c1_g2 0.090014 0.07987 0.239217 0.437166 0.370539 0.392304

<figure>
<img
src="fuscus_tissue_specificity_files/figure-gfm/ts-plot-counts-high-1.png"
alt="Plots of the counts from the DESeq2 Dataset for genes with the highest tau" />
<figcaption aria-hidden="true">Plots of the counts from the DESeq2
Dataset for genes with the highest tau</figcaption>
</figure>

    ##                       tau organ
    ## TRINITY_DN71114_c0_g1   0    NA
    ## TRINITY_DN71150_c0_g1   0    NA
    ## TRINITY_DN71111_c0_g1   0    NA
    ## TRINITY_DN71157_c0_g2   0    NA
    ## TRINITY_DN71125_c0_g1   0    NA
    ## TRINITY_DN71128_c0_g1   0    NA
    ## TRINITY_DN71103_c0_g1   0    NA
    ## TRINITY_DN71165_c0_g1   0    NA
    ## TRINITY_DN71099_c0_g1   0    NA

    ##                       FUG10M2 FUG11F1 FUG11M2 FUG12M1 FUG13F1 FUG13F4 FUG15M5
    ## TRINITY_DN71114_c0_g1       0       0       0       0       0       0       0
    ##                       FUG2F2 FUG3F2 FUL10M2 FUL11F1  FUL11M2  FUL12M1 FUL13F1
    ## TRINITY_DN71114_c0_g1      0      0       0       0 0.800517 0.839772       0
    ##                       FUL13F4 FUL15M5 FUL2F2 FUL3F2 FUO11F1 FUO13F1 FUO13F4
    ## TRINITY_DN71114_c0_g1       0       0      0      0       0       0       0
    ##                       FUO2F2 FUO3F2 FUT10M2 FUT11M2 FUT12M1 FUT15M5
    ## TRINITY_DN71114_c0_g1      0      0       0       0       0       0

    ##                       FUG10M2 FUG11F1 FUG11M2  FUG12M1 FUG13F1 FUG13F4 FUG15M5
    ## TRINITY_DN71150_c0_g1       0       0       0 0.467926       0 1.03467       0
    ##                       FUG2F2 FUG3F2 FUL10M2 FUL11F1  FUL11M2 FUL12M1 FUL13F1
    ## TRINITY_DN71150_c0_g1      0      0       0       0 1.004014       0       0
    ##                       FUL13F4 FUL15M5 FUL2F2 FUL3F2 FUO11F1 FUO13F1 FUO13F4
    ## TRINITY_DN71150_c0_g1       0       0      0      0       0       0       0
    ##                       FUO2F2 FUO3F2 FUT10M2 FUT11M2  FUT12M1 FUT15M5
    ## TRINITY_DN71150_c0_g1      0      0       0       0 1.503289       0

    ##                       FUG10M2  FUG11F1 FUG11M2 FUG12M1 FUG13F1 FUG13F4  FUG15M5
    ## TRINITY_DN71111_c0_g1       0 0.799107       0       0       0       0 0.927938
    ##                       FUG2F2 FUG3F2 FUL10M2 FUL11F1 FUL11M2 FUL12M1 FUL13F1
    ## TRINITY_DN71111_c0_g1      0      0       0 0.23221       0       0       0
    ##                       FUL13F4 FUL15M5 FUL2F2 FUL3F2 FUO11F1 FUO13F1  FUO13F4
    ## TRINITY_DN71111_c0_g1       0       0      0      0       0       0 0.349856
    ##                       FUO2F2 FUO3F2 FUT10M2  FUT11M2 FUT12M1 FUT15M5
    ## TRINITY_DN71111_c0_g1      0      0       0 1.472693       0       0

    ##                       FUG10M2 FUG11F1  FUG11M2  FUG12M1 FUG13F1  FUG13F4
    ## TRINITY_DN71157_c0_g2       0       0 0.122256 0.824265       0 0.317157
    ##                       FUG15M5 FUG2F2 FUG3F2 FUL10M2 FUL11F1 FUL11M2 FUL12M1
    ## TRINITY_DN71157_c0_g2       0      0      0       0       0       0       0
    ##                       FUL13F1 FUL13F4 FUL15M5 FUL2F2 FUL3F2 FUO11F1 FUO13F1
    ## TRINITY_DN71157_c0_g2       0       0       0      0      0       0       0
    ##                       FUO13F4 FUO2F2 FUO3F2 FUT10M2 FUT11M2 FUT12M1 FUT15M5
    ## TRINITY_DN71157_c0_g2       0      0      0       0       0       0       0

    ##                       FUG10M2 FUG11F1  FUG11M2  FUG12M1  FUG13F1 FUG13F4
    ## TRINITY_DN71125_c0_g1       0       0 0.603372 0.680796 0.481815       0
    ##                       FUG15M5 FUG2F2   FUG3F2  FUL10M2 FUL11F1 FUL11M2 FUL12M1
    ## TRINITY_DN71125_c0_g1       0      0 0.313659 0.536639       0       0       0
    ##                        FUL13F1  FUL13F4 FUL15M5 FUL2F2   FUL3F2 FUO11F1 FUO13F1
    ## TRINITY_DN71125_c0_g1 0.242474 0.419704       0      0 0.686315       0       0
    ##                       FUO13F4   FUO2F2 FUO3F2 FUT10M2 FUT11M2 FUT12M1 FUT15M5
    ## TRINITY_DN71125_c0_g1       0 1.893492      0       0       0 0.92331       0

    ##                       FUG10M2 FUG11F1 FUG11M2 FUG12M1 FUG13F1 FUG13F4 FUG15M5
    ## TRINITY_DN71128_c0_g1       0       0       0       0       0       0       0
    ##                       FUG2F2 FUG3F2 FUL10M2 FUL11F1 FUL11M2 FUL12M1 FUL13F1
    ## TRINITY_DN71128_c0_g1      0      0       0       0       0       0       0
    ##                       FUL13F4 FUL15M5 FUL2F2 FUL3F2 FUO11F1 FUO13F1 FUO13F4
    ## TRINITY_DN71128_c0_g1       0       0      0      0       0       0       0
    ##                       FUO2F2 FUO3F2 FUT10M2 FUT11M2 FUT12M1  FUT15M5
    ## TRINITY_DN71128_c0_g1      0      0       0       0       0 1.635123

    ##                       FUG10M2 FUG11F1 FUG11M2  FUG12M1  FUG13F1 FUG13F4
    ## TRINITY_DN71103_c0_g1       0       0       0 0.580412 0.068605       0
    ##                        FUG15M5  FUG2F2 FUG3F2 FUL10M2 FUL11F1 FUL11M2 FUL12M1
    ## TRINITY_DN71103_c0_g1 0.509153 0.15601      0       0       0       0       0
    ##                       FUL13F1 FUL13F4  FUL15M5 FUL2F2 FUL3F2 FUO11F1 FUO13F1
    ## TRINITY_DN71103_c0_g1       0       0 0.081366      0      0       0       0
    ##                       FUO13F4 FUO2F2 FUO3F2 FUT10M2 FUT11M2 FUT12M1 FUT15M5
    ## TRINITY_DN71103_c0_g1       0      0      0       0       0       0       0

    ##                        FUG10M2 FUG11F1 FUG11M2 FUG12M1 FUG13F1  FUG13F4
    ## TRINITY_DN71165_c0_g1 1.371174       0       0       0       0 1.717833
    ##                        FUG15M5 FUG2F2 FUG3F2 FUL10M2 FUL11F1 FUL11M2 FUL12M1
    ## TRINITY_DN71165_c0_g1 1.086165      0      0       0       0       0       0
    ##                       FUL13F1 FUL13F4 FUL15M5 FUL2F2 FUL3F2 FUO11F1 FUO13F1
    ## TRINITY_DN71165_c0_g1       0       0       0      0      0       0       0
    ##                       FUO13F4 FUO2F2 FUO3F2 FUT10M2 FUT11M2 FUT12M1 FUT15M5
    ## TRINITY_DN71165_c0_g1       0      0      0 0.84644       0       0       0

    ##                       FUG10M2 FUG11F1 FUG11M2 FUG12M1  FUG13F1 FUG13F4  FUG15M5
    ## TRINITY_DN71099_c0_g1       0       0       0       0 0.786263       0 0.120757
    ##                         FUG2F2   FUG3F2 FUL10M2 FUL11F1 FUL11M2 FUL12M1 FUL13F1
    ## TRINITY_DN71099_c0_g1 0.994992 0.252156       0       0       0       0       0
    ##                       FUL13F4 FUL15M5 FUL2F2 FUL3F2 FUO11F1 FUO13F1 FUO13F4
    ## TRINITY_DN71099_c0_g1       0       0      0      0       0       0       0
    ##                       FUO2F2 FUO3F2 FUT10M2 FUT11M2 FUT12M1  FUT15M5
    ## TRINITY_DN71099_c0_g1      0      0       0       0       0 0.363972

<figure>
<img
src="fuscus_tissue_specificity_files/figure-gfm/ts-plot-counts-low-1.png"
alt="Plots of the counts from the DESeq2 Dataset for genes with the lowest tau" />
<figcaption aria-hidden="true">Plots of the counts from the DESeq2
Dataset for genes with the lowest tau</figcaption>
</figure>

We can see that for the genes that have high tau values the plots of the
normalized counts show higher expression in just that tissue, whereas
for the plots of the lower tau values no clear patterns are showing up
in the plots of the counts. This is what we want to see, increasing
confidence that the the calculations of tau are accurate. We can also
see that the assignment of the tissues from the est_tau function match
the outputs from looking at the raw counts.

# Sex Bias and Tissue Specificity

Now that I have calculated $\tau$ I want to relate it back to
information about sex bias to see what, if any, relationship exists. I
am going to do this separately for the different organs to see if the
relationship changes based on organ type.

## $\tau$ and logFC in sex-biased genes for each organ

A gene was determined to be sex-biased if it had a log-fold change
$\ge |2|$ AND an adjusted p-value of $< 0.05$. Earlier I read in four
datasets that contain information for all of the significantly expressed
genes in each organ (i.e. the genes that were not given a p-value of NA
from DESeq).

From these datasets I can pull out the female- and male-biased genes to
plot logFC against $\tau$ within the genes that are sex-biased.

<figure>
<img src="fuscus_tissue_specificity_files/figure-gfm/tau-v-bias-1.png"
alt="Sex bias (in terms of the square root of the absolute value of the log2FoldChange) versus tissue specificity (tau) for all genes that are sex biased in the gonads, liver, skin, and brain." />
<figcaption aria-hidden="true">Sex bias (in terms of the square root of
the absolute value of the log2FoldChange) versus tissue specificity
(tau) for all genes that are sex biased in the gonads, liver, skin, and
brain.</figcaption>
</figure>

The absolute value of the log2 fold change was square root transformed.
There appears to be a relationship between $\tau$ and sex bias for the
gonads, but potentially not for the liver and the gills. This may be due
to a low number of sex-biased genes. We can look further into this
relationship with a Spearman’s rank correlation test. This was done
between $\tau$ and the $|log2FoldChange|$ for each organ separately.

    ## 
    ##  Spearman's rank correlation rho
    ## 
    ## data:  gonad_bias_tau$tau and sqrt(abs(gonad_bias_tau$log2FoldChange))
    ## S = 8.9942e+11, p-value = 0.01769
    ## alternative hypothesis: true rho is not equal to 0
    ## sample estimates:
    ##      rho 
    ## 0.017856

    ## 
    ##  Spearman's rank correlation rho
    ## 
    ## data:  liver_bias_tau$tau and sqrt(abs(liver_bias_tau$log2FoldChange))
    ## S = 893718, p-value = 0.007974
    ## alternative hypothesis: true rho is not equal to 0
    ## sample estimates:
    ##       rho 
    ## 0.1929684

    ## 
    ##  Spearman's rank correlation rho
    ## 
    ## data:  gill_bias_tau$tau and sqrt(abs(gill_bias_tau$log2FoldChange))
    ## S = 53683, p-value = 0.8749
    ## alternative hypothesis: true rho is not equal to 0
    ## sample estimates:
    ##        rho 
    ## 0.01930204

The liver and gonad both show a significant correlation. This isn’t
super obvious from just looking at the plot.

## $\tau$ and logFC for all the signficantly expressed genes in each organ

Only looking at the logFC for the sex-biased genes may be affecting the
relationship that I am seeing between sex bias and $\tau$, especially in
the organs where the samples size is very small if I am only including
sex-biased genes (like the liver and gills).

I used the same steps as above to see how the relationship changes when
I add in ALL of the genes that are significantly expressed for each
organ. In this datasets I have removed all rows where the adjusted
p-value is NA (outliers and low-counts).

<figure>
<img src="fuscus_tissue_specificity_files/figure-gfm/tau-v-FC-1.png"
alt="Sex bias (in terms of the absolute value of the log2FoldChange) versus tissue specificity (tau) for all genes that are signifcantly expressed in the gonads, liver, and gills." />
<figcaption aria-hidden="true">Sex bias (in terms of the absolute value
of the log2FoldChange) versus tissue specificity (tau) for all genes
that are signifcantly expressed in the gonads, liver, and
gills.</figcaption>
</figure>

The absolute value of the log2 fold change was square root transformed.
We can see similar patterns as when we only plotted the sex biased genes
with the strongest relationship showing up in the gonads. When we add in
all of the points, however, it seems that the best fit line for the
liver and gills now points downward. Let’s see how the Spearman’s rank
correlation test have changed. This was done between $\tau$ and the
$|log2FoldChange|$ for each organ separately.

    ## 
    ##  Spearman's rank correlation rho
    ## 
    ## data:  gonad_bias_tau_all$tau and abs(gonad_bias_tau_all$log2FoldChange)
    ## S = 3.3767e+13, p-value < 2.2e-16
    ## alternative hypothesis: true rho is not equal to 0
    ## sample estimates:
    ##       rho 
    ## 0.1126625

    ## 
    ##  Spearman's rank correlation rho
    ## 
    ## data:  liver_bias_tau_all$tau and abs(liver_bias_tau_all$log2FoldChange)
    ## S = 3.6241e+12, p-value < 2.2e-16
    ## alternative hypothesis: true rho is not equal to 0
    ## sample estimates:
    ##       rho 
    ## -0.121105

    ## 
    ##  Spearman's rank correlation rho
    ## 
    ## data:  gill_bias_tau_all$tau and abs(gill_bias_tau_all$log2FoldChange)
    ## S = 1.5541e+13, p-value < 2.2e-16
    ## alternative hypothesis: true rho is not equal to 0
    ## sample estimates:
    ##        rho 
    ## -0.1856574

All correlations are now significant, however, for the gills and the
liver there is a significantly negative correlation between $\tau$ and
the degree of sex bias.

## Investigating genes with high tau and high logFC\`

Because we have information about which tissue each gene is specific in,
I was to see if genes that are very sex-biased in a tissue are also
tissue-specific in that same tissue.

### Gonads

<figure>
<img
src="fuscus_tissue_specificity_files/figure-gfm/gonad-tau-logFC-1.png"
alt="Sex bias (in terms of the absolute value of the log2FoldChange) versus tissue specificity (tau) for all genes that are sex-biased in the gonads. Color represents the tissue that it is specific in which may be either the testis or ovaries (yellow), the liver (pink), the gills (blue), or none (black)." />
<figcaption aria-hidden="true"><em>Sex bias (in terms of the absolute
value of the log2FoldChange) versus tissue specificity (tau) for all
genes that are sex-biased in the gonads. Color represents the tissue
that it is specific in which may be either the testis or ovaries
(yellow), the liver (pink), the gills (blue), or none
(black).</em></figcaption>
</figure>

In the gonads it looks like a lot of the genes that are extremely
sex-biased are also tissue specific in the gonads (whether that is
tissue-specific in the testis or in the ovaries).

### Liver

<figure>
<img
src="fuscus_tissue_specificity_files/figure-gfm/liver-tau-logFC-1.png"
alt="Sex bias (in terms of the absolute value of the log2FoldChange) versus tissue specificity (tau) for all genes that are sex-biased in the liver. Color represents the tissue that it is specific in which may be either the testis or ovaries (yellow), the liver (pink), the gills (blue), or none (black)." />
<figcaption aria-hidden="true"><em>Sex bias (in terms of the absolute
value of the log2FoldChange) versus tissue specificity (tau) for all
genes that are sex-biased in the liver. Color represents the tissue that
it is specific in which may be either the testis or ovaries (yellow),
the liver (pink), the gills (blue), or none (black).</em></figcaption>
</figure>

Similar to the gonads, most of the extremely biased genes in the liver
are also tissue specific in the liver, expect for a few cases. It looks
like some of the most biased genes in the liver are showing more
expression in the gonads and there are two cases where tau is 0. Other
than that, many of the more intermediately sex-biased genes show the
highest expression in gonads.

### Gills

<figure>
<img
src="fuscus_tissue_specificity_files/figure-gfm/gills-tau-logFC-1.png"
alt="Sex bias (in terms of the absolute value of the log2FoldChange) versus tissue specificity (tau) for all genes that are sex-biased in the gills. Color represents the tissue that it is specific in which may be either the testis or ovaries (yellow), the liver (pink), the gills (blue), or none (black)." />
<figcaption aria-hidden="true"><em>Sex bias (in terms of the absolute
value of the log2FoldChange) versus tissue specificity (tau) for all
genes that are sex-biased in the gills. Color represents the tissue that
it is specific in which may be either the testis or ovaries (yellow),
the liver (pink), the gills (blue), or none (black).</em></figcaption>
</figure>

There aren’t many biased genes in the gills, however, of the ones that
are there is quite a mix in terms of which tissue is showing the highest
overall expression.

## Categories of sex bias vs $\tau$

To further investigate the relationship between $\tau$ and sex biased
and possibly get a cleaner idea of what may be going on I categorized
the sex biased genes based on a series of logFC thresholds. After that I
can plot $\tau$ against the bias level.

### Within the sex-biased genes

I am first going to combine all of the individual organ datasets
together in a long format dataset which includes information about the
tissue, the logFC, the geneID, tau, and the organ it has the highest
expression in.

After that I am going to categorize all of the sex-biased genes based on
logFC into different types of bias levels, the same way that I did for
the differential expression.

``` r
#Create a long format dataset that has the tissue information, logFC, tau, and geneID 
logFC_long <- data.frame(tissue=c(rep("Gill", nrow(gill_bias_tau)),
                                  rep("Gonad", nrow(gonad_bias_tau)),
                                  rep("Liver", nrow(liver_bias_tau))),
                         logFC=c(gill_bias_tau$log2FoldChange,
                                 gonad_bias_tau$log2FoldChange,
                                 liver_bias_tau$log2FoldChange),
                         geneID=c(gill_bias_tau$Row.names,
                                  gonad_bias_tau$Row.names,
                                  liver_bias_tau$Row.names),
                         tau=c(gill_bias_tau$tau,
                               gonad_bias_tau$tau,
                               liver_bias_tau$tau),
                         organ_exp=c(gill_bias_tau$organ,
                                     gonad_bias_tau$organ,
                                     liver_bias_tau$organ)
                         )

#Add a column to denote female-biased or male biased
logFC_long$bias <- ifelse(logFC_long$logFC < 0, 
                          paste0("FB"), 
                          paste0("MB"))

#Categorize the degree of sex bias in each row
#Make a vector that contains all of the groupings
biased_bins <- c("low", "med", "high", "extreme", "sex-specific")

#Create a new column in the dataset and use ifelse statements to set the category limits
#abs(logFC) was used to account for the fem-biased genes possessing negative values
logFC_long$bias_cat <- ifelse(abs(logFC_long$logFC) >= 2 
                              & abs(logFC_long$logFC) < 3,
                              biased_bins[1],
                                   ifelse(abs(logFC_long$logFC) >= 3 
                                          & abs(logFC_long$logFC) < 5,
                                          biased_bins[2],
                                          ifelse(abs(logFC_long$logFC) >= 5 &
                                                   abs(logFC_long$logFC) < 10,
                                                 biased_bins[3], 
                                                 biased_bins[4])
                                          )
                                   )

#Making sure the ordering stays correct
logFC_long$bias_cat <- factor(logFC_long$bias_cat,
    levels = biased_bins, ordered = TRUE)

#write.csv(logFC_long, "data/logFC_long_taubiasSB.csv", row.names = FALSE)
```

<figure>
<img
src="fuscus_tissue_specificity_files/figure-gfm/tau-v-biascat-plot-1.png"
alt="Tissue specifcity (tau) across the different bias levels. Color denotes female-biased versus male biased and jitters were added to show all of the raw tau values." />
<figcaption aria-hidden="true"><em>Tissue specifcity (tau) across the
different bias levels. Color denotes female-biased versus male biased
and jitters were added to show all of the raw tau
values.</em></figcaption>
</figure>

### Across all significantly expressed genes

I want to see how $\tau$ in unbiased genes compares to these categories
as well. To do that I recreated the long-style dataset to include all of
the sig. expressed genes rather than just the sex-biased genes.

Similar to last time, I first need to combine all of the tissue specific
datasets together in the long format and then I need to categorize all
of the genes again.

``` r
#Create a long format dataset that has the tissue information, logFC, tau, and geneID 
logFC_long_all <- data.frame(tissue=c(rep("Gill", nrow(gill_bias_tau_all)),
                                      rep("Gonad", nrow(gonad_bias_tau_all)),
                                      rep("Liver", nrow(liver_bias_tau_all))
                                      ),
                             logFC=c(gill_bias_tau_all$log2FoldChange,
                                     gonad_bias_tau_all$log2FoldChange,
                                     liver_bias_tau_all$log2FoldChange
                                     ),
                             geneID=c(gill_bias_tau_all$Row.names,
                                      gonad_bias_tau_all$Row.names,
                                      liver_bias_tau_all$Row.names
                                      ),
                             tau=c(gill_bias_tau_all$tau,
                                   gonad_bias_tau_all$tau,
                                   liver_bias_tau_all$tau
                                   ),
                             padj=c(gill_bias_tau_all$padj,
                                    gonad_bias_tau_all$padj,
                                    liver_bias_tau_all$padj),
                             organ_exp=c(gill_bias_tau_all$organ,
                                         gonad_bias_tau_all$organ,
                                         liver_bias_tau_all$organ)
                             )

#Add a column to denote female-biased or male biased
logFC_long_all$bias <- ifelse(logFC_long_all$logFC <= -2 & 
                                logFC_long_all$padj <= 0.05,
                              paste0("FB"),
                              ifelse(logFC_long_all$logFC >= 2 &
                                       logFC_long_all$padj <= 0.05,
                                     paste0("MB"),
                                     paste0("UB")
                                     )
                              )
  
#Categorize the degree of sex bias in each row
#Make a vector that contains all of the groupings
biased_bins <- c("Unbiased", "Low", "Med", "High", "Extreme", "Sex-specific")

#Create a new column in the dataset and use ifelse statements to set the category limits
#abs(logFC) was used to account for the fem-biased genes possessing negative values
logFC_long_all$bias_cat <- ifelse(logFC_long_all$bias == "UB", 
                              biased_bins[1],
                              ifelse(abs(logFC_long_all$logFC) >= 2 &
                                       abs(logFC_long_all$logFC) < 3,
                                     biased_bins[2],
                                     ifelse(abs(logFC_long_all$logFC) >= 3 &
                                              abs(logFC_long_all$logFC) < 5,
                                            biased_bins[3],
                                            ifelse(abs(logFC_long_all$logFC) >= 5 &
                                                     abs(logFC_long_all$logFC) < 10,
                                                   biased_bins[4],
                                                   biased_bins[5])
                                            )
                                     )
                              )
  

#Making sure the ordering stays correct
logFC_long_all$bias_cat <- factor(logFC_long_all$bias_cat,
    levels = biased_bins, ordered = TRUE)
```

As a bit of a side quest, I am also interested in looking at how the
male expression compares to the female expression for all of the
sex-biased genes across the organs.

``` r
#Looking at the relationship between male expression and female expression
organs <- c(unique(logFC_long_all$tissue))

dds_FU_exp$Organ <- ifelse(dds_FU_exp$Organ %in% c("Ovaries", "Testis"),
                           paste0("Gonad"), 
                           paste0(dds_FU_exp$Organ))

bias_colors <- c("FB" = "#7fc97f", "MB" = "#beaed4", "UB" = "#A9A9A975") 

#pdf("docs/figs/Fig_MFexpression.pdf", width = 10, height=4)

par(mfrow=c(1,3), oma=c(4,4,2,8), mar=c(1,1,1,0))

for (organ in organs) {
  
  tmp <- as.data.frame(counts(dds_FU_exp, normalized = TRUE)[, dds_FU_exp$Organ == organ])
  
  tmp$Fmedian <- apply(tmp[, grep("F\\d", colnames(tmp))], 1, median)
  tmp$Mmedian <- apply(tmp[, grep("M", colnames(tmp))], 1, median)
  
  tmp <- tmp[rownames(tmp) %in% logFC_long_all[logFC_long_all$tissue == organ, ]$geneID, ]
  
  tmp <- merge(tmp, logFC_long_all[logFC_long_all$tissue == organ, ],
               by.x = 0, 
               by.y = "geneID", 
               all = TRUE)
  ymax <- max(log(tmp$Mmedian+0.01), na.rm = TRUE) + 1
  xmax <- max(log(tmp$Fmedian+0.01), na.rm = TRUE) + 1
  
  plot(log(tmp$Fmedian+0.01), log(tmp$Mmedian+0.01),
       #ylab = expression('log'[2] * '(Male Expression)'),
       #xlab = expression('log'[2] * '(Female Expression)'),
       axes=FALSE,
       cex.main=2,
       #main = organ,
       col = bias_colors[tmp$bias],
       pch = 19,
       cex = ifelse(organ != "Gonad" & tmp$bias %in% c("MB", "FB"),
                    0.75,
                    0.25),
       ylim = c(0, ymax),
       xlim = c(0, xmax))
  
  
  
  #Make the axis lines longer and thicker
  axis(1, labels=seq(-2, xmax, 2), 
       line=NA, at=seq(-2, xmax, 2), 
       lwd=2,xlim=c(0, xmax),cex.axis=1.5,las=1)
  axis(2, labels=seq(-2, ymax, 2),
       line=NA, at=seq(-2, ymax, 2), 
       lwd=2, ylim=c(0,ymax),cex.axis=1.5,las=1)
  
  mtext(organ, cex=1.5,outer=FALSE, line=-1)
  
}

# add x axis label
mtext(expression('log'[2] * '(Female Expression)'), side = 1, cex=1.25, outer=TRUE, line=2)

# add y axis label
mtext(expression('log'[2] * '(Male Expression)'), side = 2, cex=1.25, outer=TRUE, line=2)

# add legend
spfTools::outer_legend("right",
                       legend=c("Female\nbiased",
                                "Male\nbiased",
                                "Unbiased"),
                       pt.bg=bias_colors,
                       pch=22,
                       bty='n',
                       ncol=1,
                       cex=2,
                       y.intersp = 1.5,
                       pt.cex=2)
```

<figure>
<img
src="fuscus_tissue_specificity_files/figure-gfm/mvf-expression-1.png"
alt="Compare the male median expression to the female median expression for all male- and female-biased genes in the gonads, liver, and gills." />
<figcaption aria-hidden="true"><em>Compare the male median expression to
the female median expression for all male- and female-biased genes in
the gonads, liver, and gills.</em></figcaption>
</figure>

``` r
#dev.off()
```

I also want to include the sex-specific genes in this comparison against
$\tau$, so I need to pull those out again. The list of these genes was
orignially generated during the differential expression analysis.

``` r
#Specify the directory where Sex-specific gene files are located
ss_genes_path <- "data"

#Create a list of the files I want
FU_sex_specific_genes_files <- list.files(ss_genes_path, 
                                          pattern = "specific")

#Create an empty list to store my datasets in
FU_sex_specific_genes <- list()

#Create a loop to read in all of the sex-specific genes
for(file_name in FU_sex_specific_genes_files){
  
  # Read the file and store it in a new object
  file_data <- read.table(file.path(ss_genes_path, file_name))
  
  #Add column names to the dataset
  colnames(file_data) <- c("trin_geneid")
  
  # Create a new object with a name based on the file name
  ss_name <- sub("_TRgenes.txt$", "", file_name) #Removes the file extension
  FU_sex_specific_genes[[ss_name]] <- file_data
}
```

Now that I have the list of sex-specific genes, I can add that
information to the growing dataset that I built before
(`logFC_long_all`).

``` r
#Changing the bias category for the genes we just labeled as sex-specific
organs <- c("Gill", "Gonad", "Liver",
            "Gill", "Gonad", "Liver")
bias <- c("FB", "FB", "FB",  
          "MB", "MB", "MB")

for(i in 1:length(FU_sex_specific_genes)){
  
  tmp <- FU_sex_specific_genes[[i]]
  
  for(j in 1:nrow(tmp)){
    
    one_gene <- tmp[j, ]
    if (one_gene %in% row.names(tau)){
      
    if(one_gene %in%
       logFC_long_all[logFC_long_all$tissue == organs[[i]] &
                      logFC_long_all$bias == bias[[i]],"gene_ID"]){
       
       logFC_long_all[logFC_long_all$geneID == one_gene &
                      logFC_long_all$tissue == organs[[i]] &
                      logFC_long_all$bias == bias[[i]],"bias_cat"
                      ] <- "Sex-specific"
    }else{
      
      one_gene_dat <- data.frame(matrix(ncol= ncol(logFC_long_all),
                                        nrow=1))
      colnames(one_gene_dat) <- colnames(logFC_long_all)
      
      one_gene_dat$tissue <- organs[[i]]
      one_gene_dat$geneID <- one_gene
      one_gene_dat$tau <- tau[row.names(tau) == one_gene,]$tau
      one_gene_dat$bias <- bias[[i]]
      one_gene_dat$bias_cat <- "Sex-specific"
      rownames(one_gene_dat) <- NULL
      
      logFC_long_all <- rbind(one_gene_dat, logFC_long_all)
      
      rownames(logFC_long_all) <- NULL
    }
    }  
  }
}

#Making sure the ordering stays correct
logFC_long_all$bias_cat <- factor(logFC_long_all$bias_cat,
    levels = biased_bins, ordered = TRUE)

#write.csv(logFC_long_all, "data/logFC_long_taubias_SS.csv", row.names = FALSE)
```

Now that the dataset is properly curated, I can plot all of the
information:

<figure>
<img
src="fuscus_tissue_specificity_files/figure-gfm/tau-v-biascat-plot-all-vp-1.png"
alt="Violin plot showing tissue specifcity (tau) across the different bias levels. Color denotes female-biased versus male biased versus unbiased and jitters were added to show all of the raw tau values." />
<figcaption aria-hidden="true">Violin plot showing tissue specifcity
(tau) across the different bias levels. Color denotes female-biased
versus male biased versus unbiased and jitters were added to show all of
the raw tau values.</figcaption>
</figure>

Let’s also look at the plot without all the jitters added:

<figure>
<img
src="fuscus_tissue_specificity_files/figure-gfm/tau-v-biascat-plot-all-vp2-1.png"
alt="Violin plot showing tissue specifcity (tau) across the different bias levels. Color denotes female-biased versus male biased versus unbiased." />
<figcaption aria-hidden="true">Violin plot showing tissue specifcity
(tau) across the different bias levels. Color denotes female-biased
versus male biased versus unbiased.</figcaption>
</figure>

Remove the violin and boxplots from groups that have a sample size of
less than 10

``` r
#pdf("docs/figs/Fig4B_tau_biascat_violin_jitter2.pdf",width = 12, height=3.75)
# Preprocessing to create a new variable indicating sample size greater than 10
logFC_long_all <- logFC_long_all %>%
  group_by(bias_cat, tissue, bias) %>%
  mutate(sample_size = n()) %>%
  ungroup() %>%
  mutate(plot_type = ifelse(sample_size > 10, "box_violin", "jitter"))

# Plotting
ggplot(logFC_long_all, aes(x = bias_cat, y = tau, fill = bias)) +
  geom_violin(data = filter(logFC_long_all, plot_type == "box_violin"), 
              position = position_dodge(), draw_quantiles = c(0.5)) +
  geom_boxplot(data = filter(logFC_long_all, plot_type == "box_violin"), 
               width = 0.1, color = "black", position = position_dodge(width = 0.9)) +
  geom_point(data = filter(logFC_long_all[logFC_long_all$bias == "MB",],
                           plot_type == "jitter"),
             aes(x = as.numeric(bias_cat) + 0.185, y = tau),
             size=1.5, 
             position = position_jitter(width = 0.10),
             bg=paste0(my_colors["MB"],"75"),
             col=my_colors["MB"], pch = 21) +
  geom_point(data = filter(logFC_long_all[logFC_long_all$bias == "FB",],
                           plot_type == "jitter"), 
             aes(x=as.numeric(bias_cat)-0.185, y=tau),
             size=1.5, 
             position = position_jitter(width = 0.10),
             bg=paste0(my_colors["FB"], "75"),
             col=my_colors["FB"],pch=21) +
  scale_x_discrete(labels= bias_labs) +
  scale_fill_manual(values = my_colors) +
  facet_grid(. ~ tissue) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        strip.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=16)) +
  labs(x = "Bias Category", y = expression(tau["TPM"]))  +
  guides(fill = guide_legend(title = "Bias", order = 3)) +
  guides(color = guide_legend(title = "Bias", order = 3)) +
  stat_n_text(data = logFC_long_all[logFC_long_all$bias == "FB",], 
              aes(x = bias_cat, y = tau),
              y.pos = -0.05,
              color = my_colors["FB"]
  ) +
  stat_n_text(data = logFC_long_all[logFC_long_all$bias == "MB",], 
              aes(x = bias_cat, y = tau),
              y.pos = 0.95,
              color = my_colors["MB"]
  ) +
  stat_n_text(data = logFC_long_all[logFC_long_all$bias == "UB",], 
              aes(x = bias_cat, y = tau),
              y.pos = -0.05,
              color = my_colors["UB"]
  )
```

![](fuscus_tissue_specificity_files/figure-gfm/tau-v-biascat-plot-all-vp3-1.png)<!-- -->

``` r
#dev.off()
```

Generate the $\tau$ v. sex-bias plot in correct format for constructing
the final figure.

``` r
logFC_long_all$tissue <- factor(logFC_long_all$tissue, 
                                levels = c("Gonad", "Liver", "Gill"), 
                                ordered = TRUE)

organ_cols <- c("Gill" = "#20B2AA", 
                "Gonad" = "#EEB422", 
                "Liver" = "#EE8262")

organs <- levels(logFC_long_all$tissue)

#pdf("docs/figs/Fig_tau_sexbias.pdf", width = 8, height=5)
par(oma = c(2,2,1,2),
    mar = c(2,2,1,2),
    xpd = FALSE)

plot(logFC_long_all$tau[logFC_long_all$tissue=="Gonad"]~
       sqrt(abs(logFC_long_all$logFC[logFC_long_all$tissue=="Gonad"])),
     xlim = c(0, 6),
     ylim = c(0, 1),
     xlab = "",
     ylab = "",
     bty = "n",
     type = 'n',
     axes = FALSE)

axis(1, pos = 0, lwd = 2, cex = 2, cex.axis = 1.5, las=1)
axis(2, pos = 0, lwd = 2, cex = 2, cex.axis = 1.5, las = 1)

clip(0, 6, 0, 1)

for(organ in organs){
  
  points(logFC_long_all$tau[logFC_long_all$tissue==organ]~
           sqrt(abs(logFC_long_all$logFC[logFC_long_all$tissue==organ])),
         col = paste0(organ_cols[organ],"50"),
         pch = 19)

}

for(organ in organs){
  
  abline(lm(logFC_long_all$tau[logFC_long_all$tissue==organ]~
              sqrt(abs(logFC_long_all$logFC[logFC_long_all$tissue==organ]))),
         col = ifelse(organ == "Gill",
                      "#6E8B3D",
                      organ_cols[organ]),
         lwd = 3,
         lty = which(organs %in% organ),
         xpd = FALSE)
  
}

outer_legend("top",
             c("Gonad", "Liver", "Gill"),
             col = c("#EEB422", "#EE8262", "#6E8B3D"),
             lwd = 3,
             bty = 'n',
             cex = 1.5,
             lty = 1:3,
             ncol = 3)

outer_legend("top",
             c("               ", "            ", "     "),
             col = c("#EEB422", "#EE8262", "#20B2AA"),
             pch = 21,
             pt.bg = c("#EEB42275", "#EE826275", "#20B2AA75"),
             bty = 'n',
             cex = 1.5,
             ncol = 3)

mtext("|log fold change|",1,cex=2, line=2)
mtext(expression(tau["TPM"]),2,cex=2, line=2.5)
```

![](fuscus_tissue_specificity_files/figure-gfm/tau-v-sexbias-fig-1.png)<!-- -->

``` r
#dev.off()
```

Combine the two plots together for the figure:

``` r
figtaua <- image_ggplot(image_read_pdf('docs/figs/Fig_tau_biascat_violins.pdf'),interpolate = TRUE)
figtaub <- image_ggplot(image_read_pdf('docs/figs/Fig_tau_sexbias.pdf'),interpolate = TRUE)


figtau<-wrap_plots(figtaub,
                   figtaua,
                   ncol=2)

figtau <- figtau + plot_annotation(tag_levels = 'A')
figtau
```

![](fuscus_tissue_specificity_files/figure-gfm/figcomb-1.png)<!-- -->

``` r
ggsave("docs/figs/FigTau.pdf",figtau,height = 4, width=16)
ggsave("docs/figs/FigTau.png",figtau,height = 4, width=16)
```

Investigate the impact of male-biased genes on patterns shown:

``` r
logFC_noMB <- logFC_long_all[logFC_long_all$bias != "MB", ]
logFC_noMB$tissue <- factor(logFC_noMB$tissue, 
                            levels = c("Gonad", "Liver", "Gill"), 
                            ordered = TRUE)

#pdf("docs/figs/Fig_tau_sexbias_noMB.pdf", width = 8, height=5)
par(oma = c(2,2,1,2),
    mar = c(2,2,1,2),
    xpd = FALSE)

plot(logFC_noMB$tau[logFC_noMB$tissue=="Gonad"]~
       sqrt(abs(logFC_noMB$logFC[logFC_noMB$tissue=="Gonad"])),
     xlim = c(0, 6),
     ylim = c(0, 1),
     xlab = "",
     ylab = "",
     bty = "n",
     type = 'n',
     axes = FALSE)

axis(1, pos = 0, lwd = 2, cex = 2, cex.axis = 1.5, las=1)
axis(2, pos = 0, lwd = 2, cex = 2, cex.axis = 1.5, las = 1)

clip(0, 6, 0, 1)

for(organ in organs){
  
  points(logFC_noMB$tau[logFC_noMB$tissue==organ]~
           sqrt(abs(logFC_noMB$logFC[logFC_noMB$tissue==organ])),
         col = paste0(organ_cols[organ],"50"),
         pch = 19)

}

for(organ in organs){
  
  abline(lm(logFC_noMB$tau[logFC_noMB$tissue==organ]~
              sqrt(abs(logFC_noMB$logFC[logFC_noMB$tissue==organ]))),
         col = ifelse(organ == "Gill",
                      "#6E8B3D",
                      organ_cols[organ]),
         lwd = 3,
         lty = which(organs %in% organ),
         xpd = FALSE)
  
}

outer_legend("top",
             c("Gonad", "Liver", "Gill"),
             col = c("#EEB422", "#EE8262", "#6E8B3D"),
             lwd = 3,
             bty = 'n',
             cex = 1.5,
             lty = 1:3,
             ncol = 3)

outer_legend("top",
             c("               ", "            ", ""),
             col = c("#EEB422", "#EE8262", "#20B2AA"),
             pch = 21,
             pt.bg = c("#EEB42275", "#EE826275", "#20B2AA75"),
             bty = 'n',
             cex = 1.5,
             ncol = 3)

mtext("|log fold change|",1,cex=2, line=2)
mtext(expression(tau["TPM"]),2,cex=2, line=2.5)
```

![](fuscus_tissue_specificity_files/figure-gfm/tau-v-sexbas-noMB-1.png)<!-- -->

``` r
#dev.off()
```
