---
title: "Analysing RNAseq Data from MSU"
author: "Coley Tosto"
date: "`r Sys.Date()`"
output: 
    html_document:
        code_folding: show
        toc: yes
        toc_float: yes
        number_sections: yes
---

# Retrieving RNAseq Data from MSU

## Using the command line program wget

In order to access the **RTSF Genomics Core** FTP (File Transfer Protocol) server outside of MSU it requires secure (FTPS) connections. The data will remain available on this serve for 60 days. More details on the RTSF genomics core FTP are available on their [website.](https://rtsf.natsci.msu.edu/genomics/data-retrieval/)


The program **wget** is capable of doing this on command line. However, only versions of wget 1.17 or later are capable of supporting FTPS. If you attempt to run wget and you end up with an error message which stating: 'Unsupported scheme "ftps"' it means that the version of wget you are using does not support FTPS.

1.  Check the version of wget that is running use: `wget -V`

2.  After confirming the version of wget will work execute the command:
    `wget -r -np -nH --ask-password ftps://<username>@<hostname>/<directory_name>`
    \
    `wget -r -np -nH --ask-password ftps://flanagans@titan.bch.msu.edu/20220902_mRNASeq_PE150`

3.  If it run successfully you should be prompted to enter your password.

4.  Wait as all of the files are downloaded.

The username, hostname, password, and directory_name can all be found in the email notifying you that the data is ready. The directory_name is the name of the subdirectory in your FTP account which contains the data for your current run that you want to download.

**NOTE**: This command will create a new sub directory named the same as your directory_name at your **current working directory** containing all the FastQ files for that run. **Make sure there is enough space available wherever these files will be downloaded to**

## Mounting external drives to WSL

1.  Create a directory in the **mnt** folder `mkdir /mnt/d`

2.  Mount the external drive to the directory you just created in the mnt folder (assume the drive shows in Windows as "D:") `mount -t drvfs d: /mnt/d`

**If you need administration privileges use sudo in front of both commands**

# Interpreting the FastQC files

This is a summary of the quality of the samples, for more information about interpreting fastQC files in general see [here](https://rtsf.natsci.msu.edu/sites/_rtsf/assets/File/FastQC_TutorialAndFAQ_080717.pdf) and for RNAseq specific interpretation [here](https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/qc_fastqc_assessment.html).

## What a 'Normal' or 'Good' fastQC file looks like

While there may be an X, !, or Check to indicate the quality of the data, they should not be taken too seriously, particularly in the case of RNA sequencing data. Here are some examples of what deviations can occur without meaning the sample is 'bad'.

| NORMAL READS |           |           |            |            |            |
|:-------------|:----------|:----------|:-----------|:-----------|:-----------|
| FLG2F7_R1    | FLG3M8_R1 | FLO3F1_R2 | FLT3M5_R2  | FUG13F1_R1 | FUO11F1_R2 |
| FLG2F7_R2    | FLG3M8_R2 | FLO3F2_R1 | FLT4M4_R1  | FUG13F4_R1 | FUO13F1_R1 |
| FLG3F2_R1    | FLG4M3_R1 | FLO3F2_R2 | FLT4M4_R2  | FUO2F2_R1  | FUO13F1_R  |
| FLG3M5_R1    | FLG4M3_R2 | FLO3F4_R1 | FLT8M7_R1  | FUO2F2_R2  | FUT10M2_R1 |
| FLG3M5_R2    | FLO2F7_R1 | FLO3F4_R2 | FLT8M7_R2  | FUO3F2_R1  | FUT12M1_R1 |
| FLG3M7_R1    | FLO2F7_R2 | FLT2M3_R1 | FUG3F2_R1  | FUO3F2_R2  |            |
| FLG3M7_R2    | FLO3F1_R1 | FLT2M3_R2 | FUG11F1_R1 | FUO11F1_R1 |            |

### Per base sequence quality and Per sequence quality scores

The **per base sequence quality** alerts us to whether or not there were any problems occurring during sequencing and if these problems mean you need to contact the sequencing facility. The y-axis denotes quality score while the x-axis represents the position in the read.

The **per sequence quality scores** plot gives the average quality score on the x-axis and the number of sequences with that average on the y-axis.

For both reads of ALL 59 samples, there were no issues with either of these modules. It is very common to see a slight decrease in the per base sequence quality plot towards the end of the reads that could be explained by things such as signal decay or phasing.

<p float="left">

<img src="../imgs/per_base_seq_qual.JPG" style="width:350px;"/> <img src="../imgs/per_seq_qual_score.JPG" style="width:350px;"/>

</p>

### Per base sequence content

This plot gives you the percent of bases called for each of the four nucleotides at each position across the read. The proportion of the bases should remain relatively constant across the read.

This plot will **ALWAYS** give a fail for RNA-seq data. This is mainly due to the fact that the first 10-12 bases results from the 'random' hexamer priming that occurs during RNA-seq library preparation. The priming, however, is not as random as one would hope give an enrichment to particular bases for those first few nucleotides (see below for example). This was not taken into account when determining the 'good' samples.

<center>

<img src="../imgs/per_base_seq_content.JPG" style="width:650px;"/>

<center>

### Per sequence GC content

This module plots the number of reads vs. GC% per read. The theoretical distribution is shown by the blue line. If the observed distribution deviates too far from the theoretical it will be called a fail. In RNA-seq there may be greater or lesser distribution of GC content among transcripts causing the plots to be wider or narrower than the normal distribution. XX reads showed a WARNING flag on this module.

<center>

<img src="../imgs/gc_content.JPG" style="width:650px;"/>

<center>

### Per base N content

This module shows the percent of bases at each position with no base call (i.e. 'N'). You don't want to see a point where the curve rises noticeably above zero. None of the reads failed this.

<center>

<img src="../imgs/n_content.JPG" style="width:650px;"/>

<center>

### Sequence duplication levels

This module explores the number of duplicated sequences in the library. When sequencing RNA there will be some very highly abundant transcripts and some lowly abundant. Any that are highly abundant may be observed as duplicate reads. For the majority of the reads the largest peak was seen at \>10 sequence duplication level. A read was considered 'normal' if it showed a similar patter to the plot below (even though it is marked as fail).

<center>

<img src="../imgs/seq_dup_level.JPG" style="width:650px;"/>

<center>

### Overrepresented sequences

This table displays the sequences (at least 20bp) that occur in more that 0.1% of the total number of sequences. It can help aid in identifying contamination such as vector adapter sequences. Each overrepresented sequence is compared to a list of common contaminants to try to identify it. When using RNA-seq data it is possible that some transcripts may be so abundant that they become registered as overrepresented. 

A read was considered 'normal' if it possessed **no overrepresented sequences**.

### Adapter content

A cumulative plot of the fraction of reads where the sequence library adapter sequence is identified at the indicated base position. When using long read lengths it is possible that some of the library inserts are shorter than the read length resulting in read-through to the adapter at the 3' end of the read.

This is more likely to occur for RNA-seq libraries due to the variability of the library insert sizes. This was not taken into account when denoting a 'normal' read.

<center>

<img src="../imgs/adapter_content.JPG" style="width:650px;"/>

<center>

## Reads that raised some concerns

Overall the reads look good (especially in terms of things such as the quality score), however, some weird patterns are present that should be noted.

| GC CONTENT | SEQUENCE DUPLICATION | OVERREPRESENTED SEQ. |
|:-----------|:---------------------|:---------------------|
|            |                      | FLG3F1_R1            |
|            |                      | FLG3F1_R2            |
| FLG3F2_R1  |                      |                      |
|            |                      | FLG4M4_R1            |
| FLG4M4_R2  |                      | FLG4M4_R2            |
|            |                      | FLG8F3_R1            |
|            |                      | FLG8F3_R2            |
| FLL2F7_R1  | FLL2F7_R1            | FLL2F7_R1            |
| FLL2F7_R2  | FLL2F7_R2            | FLL2F7_R2            |
| FLL3F1_R1  | FLL3F1_R1            | FLL3F1_R1            |
|            | FLL3F1_R2            | FLL3F1_R2            |
| FLL3F2_R1  | FLL3F2_R1            | FLL3F2_R1            |
| FLL3F2_R2  | FLL3F2_R2            | FLL3F2_R2            |
| FLL3F4_R1  | FLL3F4_R1            | FLL3F4_R1            |
|            | FLL3F4_R2            | FLL3F4_R2            |
| FLL3M5_R1  | FLL3M5_R1            | FLL3M5_R1            |
| FLL3M5_R2  | FLL3M5_R2            | FLL3M5_R2            |
| FLL3M7_R1  | FLL3M7_R1            | FLL3M7_R1            |
| FLL3M7_R2  | FLL3M7_R2            | FLL3M7_R2            |
| FLL3M8_R1  | FLL3M8_R1            | FLL3M8_R1            |
| FLL3M8_R2  | FLL3M8_R2            | FLL3M8_R2            |
|            | FLL4M3_R1            | FLL4M3_R1            |
|            | FLL4M3_R2            | FLL4M3_R2            |
| FLL4M4_R1  | FLL4M4_R1            | FLL4M4_R1            |
| FLL4M4_R2  | FLL4M4_R2            | FLL4M4_R2            |
| FLL8F3_R1  | FLL8F3_R1            | FLL8F3_R1            |
| FLL8F3_R2  | FLL8F3_R2            | FLL8F3_R2            |
|            |                      | FLO8F3_R1            |
|            |                      | FLO8F3_R2            |
| FLT3M5_R1  |                      |                      |
|            |                      | FLT5M3_R1            |
|            |                      | FLT5M3_R2            |
| FUG2F2_R1  |                      |                      |
| FUG2F2_R2  |                      | FUG2F2_R2            |
| FUG3F2_R2  |                      |                      |
|            |                      | FUG10M2_R1           |
| FUG10M2_R2 |                      | FUG10M2_R2           |
| FUG11F1_R2 |                      | FUG11F1_R2           |
| FUG11M2_R1 |                      |                      |
| FUG11M2_R2 |                      | FUG11M2_R2           |
|            |                      | FUG11M4_R1           |
|            |                      | FUG11M4_R2           |
| FUG12M1_R1 |                      | FUG12M1_R1           |
| FUG12M1_R2 |                      | FUG12M1_R2           |
| FUG13F1_R2 |                      | FUG13F1_R2           |
| FUG13F4_R2 |                      | FUG13F4_R2           |
| FUG15M5_R1 |                      |                      |
| FUG15M5_R2 |                      |                      |
| FUL2F2_R1  | FUL2F2_R1            | FUL2F2_R1            |
| FUL2F2_R2  | FUL2F2_R2            | FUL2F2_R2            |
| FUL3F2_R1  | FUL3F2_R1            | FUL3F2_R1            |
| FUL3F2_R2  | FUL3F2_R2            | FUL3F2_R2            |
| FUL10M2_R1 | FUL10M2_R1           | FUL10M2_R1           |
| FUL10M2_R2 | FUL10M2_R2           | FUL10M2_R2           |
|            | FUL11F1_R1           | FUL11F1_R1           |
|            | FUL11F1_R2           | FUL11F1_R2           |
| FUL11M2_R1 | FUL11M2_R1           | FUL11M2_R1           |
| FUL11M2_R2 | FUL11M2_R2           | FUL11M2_R2           |
| FUL11M4_R1 | FUL11M4_R1           | FUL11M4_R1           |
| FUL11M4_R2 | FUL11M4_R2           | FUL11M4_R2           |
| FUL12M1_R1 | FUL12M1_R1           | FUL12M1_R1           |
|            | FUL12M1_R2           | FUL12M1_R2           |
|            | FUL13F1_R1           | FUL13F1_R1           |
|            | FUL13F1_R2           | FUL13F1_R2           |
|            | FUL13F4_R1           | FUL13F4_R1           |
|            | FUL13F4_R2           | FUL13F4_R2           |
| FUL15M5_R1 | FUL15M5_R1           | FUL15M5_R1           |
| FUL15M5_R2 | FUL15M5_R2           | FUL15M5_R2           |
|            |                      | FUO13F4_R1           |
|            |                      | FUO13F4_R2           |
| FUT10M2_R1 |                      |                      |
| FUT11M2_R1 |                      | FUT11M2_R1           |
|            |                      | FUT11M2_R2           |
|            |                      | FUT11M4_R1           |
|            |                      | FUT11M4_R2           |
| FUT12M1_R1 |                      |                      |
| FUT15M5_R1 |                      | FUT15M5_R1           |
|            |                      | FUT15M5_R2           |

### Violations in the per seqeunce GC content

None of the reads received a [**fail** ]{style="color:indianred"} for this section, they were all just [**warnings** ]{style="color:gold"}. The deviations from normal don't look too extreme for any of the reads, below are some examples of reads that received a [**warning** ]{style="color:gold"} for this module. Reads in which the only thing wrong with them was the GC content can likely be clumped with the "normal" reads as I am not sure there is anything extra we need to do with them. Those reads are [**blue** ]{style="color:mediumblue"} in the table below.

| GC VIOLATIONS                              |           |                                             |                                             |            |                                             |
|:-----------|:-----------|:-----------|:-----------|:-----------|:-----------|
| **FLG3F2_R1** | FLL3M5_R2 | **FLT3M5_R1**  | FUG12M1_R2 | FUL10M2_R1 | FUL15M5_R2 |
| FLG4M4_R2 | FLL3M7_R1 | **FUG2F2_R1**  | FUG13F1_R2 | FUL10M2_R2 | **FUT10M2_R1** |
| FLL2F7_R1 | FLL3M7_R2 | FUG2F2_R2 | FUG13F4_R2 | FUL11M2_R1 | FUT11M2_R1 |
| FLL2F7_R2 | FLL3M8_R1 | **FUG3F2_R2** | **FUG15M5_R1** | FUL11M2_R2 | FUT12M1_R1 |
| FLL3F1_R1 | FLL3M8_R2 | FUG10M2_R2 | **FUG15M5_R2** | FUL11M4_R1 | FUT15M5_R1 |
| FLL3F2_R1 | FLL4M4_R1 | FUG11F1_R2 | FUL2F2_R1 | FUL11M4_R2 |                     |
| FLL3F2_R2 | FLL4M4_R2 | **FUG11M2_R1** | FUL2F2_R2 | FUL12M1_R1 |                 |
| FLL3F4_R1 | FLL8F3_R1 | FUG11M2_R2 | FUL3F2_R1 | FUL10M2_R1 |                     |
| FLL3M5_R1 | FLL8F3_R2 | FUG12M1_R1 | FUL3F2_R2 | FUL15M5_R1 |                     |

<p float="left">

<img src="../imgs/Ex1_gc.JPG" style="width:350px;"/> <img src="../imgs/Ex2_gc.JPG" style="width:350px;"/>

</p>

<p float="left">

<img src="../imgs/Ex3_gc.JPG" style="width:350px;"/> <img src="../imgs/Ex4_gc.JPG" style="width:350px;"/>

</p>

<p float="left">

<img src="../imgs/Ex5_gc.JPG" style="width:350px;"/> <img src="../imgs/Ex6_gc.JPG" style="width:350px;"/>

</p>

### Violations in the seqeunce duplication levels

While all of the reads technically failed this module on the fastQC, there was a small group of reads that showed a completely different pattern than the rest. Interestingly, all of these reads happen to be from the **liver** and it is across BOTH species.

<p float="left">

<img src="../imgs/seq_dup_level.JPG" title="Example of a Normal Read" style="width:350px;"/> <img src="../imgs/seq_Ex.JPG" title="Example of a Liver Read" style="width:350px;"/>

</p>

You can see that for all of the other reads there is the largest peak at the \>10 duplication level then \>100 and then a small peak at \>1K. However for every single liver read the largest peak falls at the \>1K duplication level while the smallest is seen at the \>10. Unsure of why this may be. Possible that there are more sequences that are highly abundant in the liver?

| SEQ. DUP. VIOLATIONS |           |           |            |            |
|:---------------------|:----------|:----------|:-----------|:-----------|
| FLL2F7_R1            | FLL3M5_R1 | FLL4M4_R1 | FUL10M2_R1 | FUL12M1_R1 |
| FLL2F7_R2            | FLL3M5_R2 | FLL4M4_R2 | FUL10M2_R2 | FUL12M1_R2 |
| FLL3F1_R1            | FLL3M7_R1 | FLL8F3_R1 | FUL11F1_R1 | FUL13F1_R1 |
| FLL3F1_R2            | FLL3M7_R2 | FLL8F3_R2 | FUL11F1_R2 | FUL13F1_R2 |
| FLL3F2_R1            | FLL3M8_R1 | FUL2F2_R1 | FUL11M2_R1 | FUL13F4_R1 |
| FLL3F2_R2            | FLL3M8_R2 | FUL2F2_R2 | FUL11M2_R2 | FUL13F4_R2 |
| FLL3F4_R1            | FLL4M3_R1 | FUL3F2_R1 | FUL11M4_R1 | FUL15M5_R1 |
| FLL3F4_R2            | FLL4M3_R2 | FUL3F2_R2 | FUL11M4_R2 | FUL15M5_R2 |

### Violations in the overrepresented sequences

Certain reads possessed sequences that were marked as overrepresented. None of the reads failed this section, there were only warnings. For a lot of the reads, there was contamination from the Illumina Nextera Transposase sequence that was registered as overrepresented. This will need to get removed. The table below outlines all of the reads that received a **warning** for this module and what sequences were marked as overrepresented.

| READ       | Illumina Nextera Transposase Read1 | Illumina Nextera Transposase Read2 | No Hit |
|:------------------|:-----------------|:-----------------|:-----------------|
| FLG3F1_R1  |                                    | **X**                              |        |
| FLG3F1_R2  | **X**                              |                                    |        |
| FLG4M4_R1  |                                    | **X**                              |        |
| FLG4M4_R2  | **X**                              |                                    |        |
| FLG8F3_R1  |                                    | **X**                              |        |
| FLG8F3_R2  | **X**                              |                                    |        |
| FLL2F7_R1  |                                    |                                    | 5      |
| FLL2F7_R2  |                                    |                                    | 5      |
| FLL3F1_R1  |                                    | **X**                              | 2      |
| FLL3F1_R2  | **X**                              |                                    | 2      |
| FLL3F2_R1  |                                    | **X**                              | 4      |
| FLL3F2_R2  | **X**                              |                                    | 2      |
| FLL3F4_R1  |                                    | **X**                              | 1      |
| FLL3F4_R2  | **X**                              |                                    | 3      |
| FLL3M5_R1  |                                    | **X**                              | 4      |
| FLL3M5_R2  | **X**                              |                                    | 9      |
| FLL3M7_R1  |                                    |                                    | 4      |
| FLL3M7_R2  |                                    |                                    | 9      |
| FLL3M8_R1  |                                    |                                    | 3      |
| FLL3M8_R2  |                                    |                                    | 7      |
| FLL4M3_R1  |                                    | **X**                              | 2      |
| FLL4M3_R2  | **X**                              |                                    | 2      |
| FLL4M4_R1  |                                    | **X**                              | 3      |
| FLL4M4_R2  | **X**                              |                                    | 6      |
| FLL8F3_R1  |                                    | **X**                              | 3      |
| FLL8F3_R2  | **X**                              |                                    | 4      |
| FLO8F3_R1  |                                    | **X**                              |        |
| FLO8F3_R2  | **X**                              |                                    |        |
| FLT5M3_R1  |                                    | **X**                              |        |
| FLT5M3_R2  | **X**                              |                                    |        |
| FUG2F2_R2  |                                    |                                    | 2      |
| FUG10M2_R1 |                                    | **X**                              |        |
| FUG10M2_R2 | **X**                              |                                    | 1      |
| FUG11F1_R2 | **X**                              |                                    | 1      |
| FUG11M2_R2 | **X**                              |                                    |        |
| FUG11M4_R1 |                                    | **X**                              |        |
| FUG11M4_R2 | **X**                              |                                    |        |
| FUG12M1_R1 |                                    | **X**                              |        |
| FUG12M1_R2 | **X**                              |                                    | 1      |
| FUG13F1_R2 |                                    |                                    | 3      |
| FUG13F4_R2 |                                    |                                    | 4      |
| FUL2F2_R1  |                                    | **X**                              | 3      |
| FUL2F2_R2  | **X**                              |                                    | 1      |
| FUL3F2_R1  |                                    | **X**                              | 1      |
| FUL3F2_R2  | **X**                              |                                    | 3      |
| FUL10M2_R1 |                                    |                                    | 5      |
| FUL10M2_R2 |                                    |                                    | 8      |
| FUL11F1_R1 |                                    | **X**                              | 3      |
| FUL11F1_R2 | **X**                              |                                    | 8      |
| FUL11M2_R1 |                                    |                                    | 5      |
| FUL11M2_R2 |                                    |                                    | 11     |
| FUL11M4_R1 |                                    |                                    | 5      |
| FUL11M4_R2 |                                    |                                    | 7      |
| FUL12M1_R1 |                                    |                                    | 4      |
| FUL12M1_R2 |                                    |                                    | 12     |
| FUL13F1_R1 |                                    | **X**                              | 3      |
| FUL13F1_R2 | **X**                              |                                    | 5      |
| FUL13F4_R1 |                                    | **X**                              | 1      |
| FUL13F4_R2 | **X**                              |                                    | 4      |
| FUL15M5_R1 |                                    | **X**                              | 2      |
| FUL15M5_R2 | **X**                              |                                    | 12     |
| FUO13F4_R1 |                                    | **X**                              |        |
| FUO13F4_R2 | **X**                              |                                    |        |
| FUT11M2_R1 |                                    | **X**                              |        |
| FUT11M2_R2 | **X**                              |                                    |        |
| FUT11M4_R1 |                                    | **X**                              |        |
| FUT11M4_R2 | **X**                              |                                    |        |
| FUT15M5_R1 |                                    | **X**                              |        |
| FUT15M5_R2 | **X**                              |                                    |        |


Aside from the over representation of the Illumina Nextera Transposase seqeunce, many of the reads showed overrepresented sequences that had no hit. All of these sequences were [blasted](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) through NCBI with a non-redundant database to try and pinpoint what they are from. Many of the "No hit" overrepresented sequences are from reads generated with liver tissue. Because of this, blast results from those sequences were compared to a table of highly expressed genes found in *Syngnathus scovelli* livers ([from Rose et al. 2015](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0139401)). Those results are outlined in the table below. Blast results (for sequences from the liver reads) that did not match the table from Rose et al. 2015 are signified in **bold**.

| READ                                  | GENES                                                                      |                                                                            |                                                                              |                                                                           |                                                                            |                                                     |                                             |                                      |
|:-------|:-------|:-------|:-------|:-------|:-------|:-------|:-------|:-------|
| FLL2F7_R1                             | **S. scovelli small subunit ribosomal RNA gene**                           | **S. louisianae small subunit ribosomal RNA gene**                         | S. acus betaine--homocysteine S-methyltransferase 1                          | S. scovelli betaine--homocysteine S-methyltransferase 1                   | S. scovelli 28S ribosomal RNA                                              |                                                     |                                             |                                      |
| FLL2F7_R2                             | **S. louisianae small subunit ribosomal RNA gene** X2                      | S. acus betaine--homocysteine S-methyltransferase 1 x3                     |                                                                              |                                                                           |                                                                            |                                                     |                                             |                                      |
| FLL3F1_R1                             | **S. scovelli small subunit ribosomal RNA gene**                           | S. scovelli betaine--homocysteine S-methyltransferase 1                    |                                                                              |                                                                           |                                                                            |                                                     |                                             |                                      |
| FLL3F1_R2                             | S. scovelli apolipoprotein A-Ib                                            | S. acus betaine--homocysteine S-methyltransferase 1                        |                                                                              |                                                                           |                                                                            |                                                     |                                             |                                      |
| FLL3F2_R1                             | **S. scovelli small subunit ribosomal RNA gene**                           | **S. louisianae small subunit ribosomal RNA gene**                         | S. scovelli betaine--homocysteine S-methyltransferase 1                      | S. acus betaine--homocysteine S-methyltransferase 1                       |                                                                            |                                                     |                                             |                                      |
| FLL3F2_R2                             | **S. louisianae small subunit ribosomal RNA gene** x2                      |                                                                            |                                                                              |                                                                           |                                                                            |                                                     |                                             |                                      |
| FLL3F4_R1                             | *S. acus apolipoprotein A-Ib*                                              |                                                                            |                                                                              |                                                                           |                                                                            |                                                     |                                             |                                      |
| FLL3F4_R2                             | S. scovelli apolipoprotein A-Ib                                            | **S. scovelli hemopexin a** x2                                             |                                                                              |                                                                           |                                                                            |                                                     |                                             |                                      |
| FLL3M5_R1                             | **S. scovelli leukocyte cell-derived chemotaxin-2** x3                     | *S. acus apolipoprotein A-Ib*                                              |                                                                              |                                                                           |                                                                            |                                                     |                                             |                                      |
| FLL3M5_R2                             | **S. scovelli leukocyte cell-derived chemotaxin-2** x3                     | S. scovelli apolipoprotein A-Ib x2                                         | *S. scovelli leukocyte cell-derived chemotaxin-2*                            | S. scovelli type-4 ice-structuring protein LS-12                          | S. scovelli beta-microseminoprotein J1                                     | **S.scovelli uncharacterized LOC125975224**         |                                             |                                      |
| FLL3M7_R1                             | *S. acus apolipoprotein A-Ib*                                              | **S. scovelli hemopexin a**                                                | S. scovelli betaine--homocysteine S-methyltransferase 1                      | S. acus betaine--homocysteine S-methyltransferase 1                       |                                                                            |                                                     |                                             |                                      |
| FLL3M7_R2                             | **S. scovelli hemopexin a** x3                                             | S. scovelli apolipoprotein A-Ib x2                                         | S. scovelli beta-microseminoprotein J1                                       | S. scovelli type-4 ice-structuring protein LS-12 x2                       | S. acus betaine--homocysteine S-methyltransferase 1                        |                                                     |                                             |                                      |
| FLL3M8_R1                             | *S. acus apolipoprotein A-Ib*                                              | **S. scovelli hemopexin a**                                                | S. scovelli apolipoprotein A-Ib                                              |                                                                           |                                                                            |                                                     |                                             |                                      |
| FLL3M8_R2                             | **S. scovelli hemopexin a** x4                                             | S. scovelli apolipoprotein A-Ib x3                                         |                                                                              |                                                                           |                                                                            |                                                     |                                             |                                      |
| FLL4M3_R1                             | S. scovelli betaine--homocysteine S-methyltransferase 1                    | *S. acus apolipoprotein A-Ib*                                              |                                                                              |                                                                           |                                                                            |                                                     |                                             |                                      |
| FLL4M3_R2                             | S. scovelli apolipoprotein A-Ib                                            | S. acus trypsin-2                                                          |                                                                              |                                                                           |                                                                            |                                                     |                                             |                                      |
| FLL4M4_R1                             | S. scovelli betaine--homocysteine S-methyltransferase 1                    | S. acus betaine--homocysteine S-methyltransferase 1                        | *S. acus apolipoprotein A-Ib*                                                |                                                                           |                                                                            |                                                     |                                             |                                      |
| FLL4M4_R2                             | **S. scovelli hemopexin a** x2                                             | S. scovelli type-4 ice-structuring protein LS-12 x2                        | S. scovelli apolipoprotein A-Ib                                              | S. acus betaine--homocysteine S-methyltransferase 1                       |                                                                            |                                                     |                                             |                                      |
| FLL8F3_R1                             | S. acus betaine--homocysteine S-methyltransferase 1 x2                     | S. scovelli type-4 ice-structuring protein LS-12                           |                                                                              |                                                                           |                                                                            |                                                     |                                             |                                      |
| FLL8F3_R2                             | S. scovelli type-4 ice-structuring protein LS-12 x2                        | S. acus betaine--homocysteine S-methyltransferase 1                        | S. scovelli apolipoprotein A-Ib                                              |                                                                           |                                                                            |                                                     |                                             |                                      |
| [FUG2F2_R2 ]{style="color:orange"}    | S. scovelli galactose-specific lectin nattectin x2                         |                                                                            |                                                                              |                                                                           |                                                                            |                                                     |                                             |                                      |
| [FUG10M2_R2 ]{style="color:orange"}   | S. scovelli galactose-specific lectin nattectin                            |                                                                            |                                                                              |                                                                           |                                                                            |                                                     |                                             |                                      |
| [FUG11F1_R2 ]{style="color:orange"}   | S. scovelli galactose-specific lectin nattectin                            |                                                                            |                                                                              |                                                                           |                                                                            |                                                     |                                             |                                      |
| [ FUG12M1_R2 ]{style="color:orange"}  | S. scovelli galactose-specific lectin nattectin                            |                                                                            |                                                                              |                                                                           |                                                                            |                                                     |                                             |                                      |
| [FUG13F1_R2 ]{style="color:orange"}   | S. scovelli galactose-specific lectin nattectin                            | S. acus lithostathine-1-alpha x2                                           |                                                                              |                                                                           |                                                                            |                                                     |                                             |                                      |
| <span style="color:orange">FUG13F4_R2 | S. scovelli galactose-specific lectin nattectin x3                         | S. acus lithostathine-1-alpha                                              |                                                                              |                                                                           |                                                                            |                                                     |                                             |                                      |
| FUL2F2_R1                             | S. acus vitellogenin-1 x2                                                  | S. scovelli vitellogenin                                                   |                                                                              |                                                                           |                                                                            |                                                     |                                             |                                      |
| FUL2F2_R2                             | S. acus vitellogenin-1                                                     |                                                                            |                                                                              |                                                                           |                                                                            |                                                     |                                             |                                      |
| FUL3F2_R1                             | S. scovelli apolipoprotein A-Ib                                            |                                                                            |                                                                              |                                                                           |                                                                            |                                                     |                                             |                                      |
| FUL3F2_R2                             | S. scovelli apolipoprotein A-Ib                                            | S. scovelli type-4 ice-structuring protein LS-12 x2                        |                                                                              |                                                                           |                                                                            |                                                     |                                             |                                      |
| FUL10M2_R1                            | S. scovelli apolipoprotein A-Ib x3                                         | **S. acus leukocyte cell-derived chemotaxin-2**                            | S. scovelli complement C1q tumor necrosis factor-related protein 3-like      |                                                                           |                                                                            |                                                     |                                             |                                      |
| FUL10M2_R2                            | S. scovelli apolipoprotein A-Ib x2                                         | ***S. acus leukocyte cell-derived chemotaxin-2***                          | S. scovelli complement C1q tumor necrosis factor-related protein 3-like      | S. scovelli type-4 ice-structuring protein LS-12 x2                       | **S. scovelli leukocyte cell-derived chemotaxin-2**                        | **S. scovelli hemopexin a**                         |                                             |                                      |
| FUL11F1_R1                            | S. scovelli apolipoprotein A-Ib x2                                         | S. scovelli complement C1q tumor necrosis factor-related protein 3-like    |                                                                              |                                                                           |                                                                            |                                                     |                                             |                                      |
| FUL11F1_R2                            | S. scovelli apolipoprotein A-Ib x2                                         | S. scovelli complement C1q tumor necrosis factor-related protein 3-like x2 | ***S. scovelli leukocyte cell-derived chemotaxin-2***                        | **S. scovelli hemopexin a** x2                                            | S. scovelli alpha-2-HS-glycoprotein-like                                   |                                                     |                                             |                                      |
| FUL11M2_R1                            | S. scovelli complement C1q tumor necrosis factor-related protein 3-like x3 | S. scovelli apolipoprotein A-Ib                                            | *S. scovelli apolipoprotein A-II*                                            |                                                                           |                                                                            |                                                     |                                             |                                      |
| FUL11M2_R2                            | **S. scovelli hemopexin a** x3                                             | S. scovelli apolipoprotein A-Ib x2                                         | *S. scovelli complement C1q tumor necrosis factor-related protein 3-like* x2 | ***S. acus leukocyte cell-derived chemotaxin-2***                         | S. scovelli complement C1q tumor necrosis factor-related protein 3-like x3 | **S. acus hemopexin a**                             |                                             |                                      |
| FUL11M4_R1                            | S. scovelli complement C1q tumor necrosis factor-related protein 3-like x3 | **S. acus leukocyte cell-derived chemotaxin-2**                            | S. scovelli apolipoprotein A-Ib                                              |                                                                           |                                                                            |                                                     |                                             |                                      |
| FUL11M4_R2                            | S. scovelli apolipoprotein A-Ib x2                                         | **S. scovelli hemopexin a** x2                                             | S. scovelli complement C1q tumor necrosis factor-related protein 3-like      | *S. scovelli complement C1q tumor necrosis factor-related protein 3-like* | **S. scovelli leukocyte cell-derived chemotaxin-2**                        |                                                     |                                             |                                      |
| FUL12M1_R1                            | S. scovelli complement C1q tumor necrosis factor-related protein 3-like x2 | **S. scovelli leukocyte cell-derived chemotaxin-2**                        | S. scovelli apolipoprotein A-Ib                                              |                                                                           |                                                                            |                                                     |                                             |                                      |
| FUL12M1_R2                            | S. scovelli complement C1q tumor necrosis factor-related protein 3-like x3 | **S. scovelli hemopexin a** x3                                             | ***S. scovelli leukocyte cell-derived chemotaxin-2***                        | S. scovelli apolipoprotein A-Ib                                           | **S.scovelli uncharacterized LOC125975224**                                | **S. scovelli leukocyte cell-derived chemotaxin-2** | S. scovelli trypsin-2                       | S. scovelli hatching enzyme 1.2-like |
| FUL13F1_R1                            | S. scovelli apolipoprotein A-Ib x2                                         | S. scovelli alpha-2-HS-glycoprotein-like                                   |                                                                              |                                                                           |                                                                            |                                                     |                                             |                                      |
| FUL13F1_R2                            | S. scovelli apolipoprotein A-Ib                                            | ***S. scovelli leukocyte cell-derived chemotaxin-2***                      | **S. scovelli hemopexin a**                                                  | S. scovelli alpha-2-HS-glycoprotein-like                                  | S. acus betaine--homocysteine S-methyltransferase 1                        |                                                     |                                             |                                      |
| FUL13F4_R1                            | S. scovelli apolipoprotein A-Ib                                            |                                                                            |                                                                              |                                                                           |                                                                            |                                                     |                                             |                                      |
| FUL13F4_R2                            | S. scovelli apolipoprotein A-Ib                                            | S. scovelli beta-microseminoprotein J1                                     | S. scovelli type-4 ice-structuring protein LS-12 x2                          |                                                                           |                                                                            |                                                     |                                             |                                      |
| FUL15M5_R1                            | S. scovelli complement C1q tumor necrosis factor-related protein 3-like x2 |                                                                            |                                                                              |                                                                           |                                                                            |                                                     |                                             |                                      |
| FUL15M5_R2                            | **S. scovelli hemopexin a** x3                                             | S. scovelli complement C1q tumor necrosis factor-related protein 3-like x2 | S. scovelli apolipoprotein A-Ib x2                                           | **S. scovelli leukocyte cell-derived chemotaxin-2**                       | ***S. scovelli leukocyte cell-derived chemotaxin-2***                      | S. scovelli type-4 ice-structuring protein LS-12 x2 | **S. acus genome assembly, chromosome: 16** |                                      |

There were a group of sequences that originally blasted as "No significant similarity found". For all of those sequences, the 'optimize for' setting was changed from **highly similar sequences (megablast)** to **somewhat similar sequences (blastn)** and then those results were recorded. The sequences are *italicized* in the table above and the percent identity for them ranged from 95.83% to 98%.


Overall there were 16 genes that resulted from blasting all of the sequences, 9 of which correspond to previously reported highly abundant genes in *Syngnathus scovelli* livers and 2 that were specific to the gills (NA in table).

| GENE NAME                                                                 | ABUNDANT IN FEMALES? | ABUNDANT IN MALES? |
|:-----------------|:-------------------------:|:-------------------------:|
| *S. scovelli* alpha-2-HS-glycoprotein-like                                |        Y(29)         |       Y(43)        |
| *S. scovelli* apolipoprotein A-Ib                                         |        Y(10)         |        Y(7)        |
| *S. acus/scovelli* betaine--homocysteine S-methyltransferase 1            |         Y(8)         |        Y(9)        |
| *S. scovelli* beta-microseminoprotein J1                                  |        Y(14)         |       Y(16)        |
| *S. scovelli* complement C1q tumor necrosis factor-related protein 3-like |        Y(27)         |       Y(17)        |
| *S. scovelli* galactose-specific lectin nattectin                         |          NA          |         NA         |
| *S. scovelli* hatching enzyme 1.2-like                                    |          N           |         N          |
| *S. scovelli* hemopexin a                                                 |          N           |         N          |
| *S. scovelli* leukocyte cell-derived chemotaxin-2                         |          N           |         N          |
| *S. acus* lithostathine-1-alpha x2                                        |          NA          |         NA         |
| *S. scovelli* 28S ribosomal RNA                                           |         Y(1)         |         N          |
| *S. scovelli/louisianae* small subunit ribosomal RNA gene                 |          N           |         N          |
| *S. acus* trypsin-2                                                       |        Y(21)         |       Y(25)        |
| *S. scovelli* type-4 ice-structuring protein LS-12                        |        Y(33)         |       Y(22)        |
| *S.scovelli* uncharacterized LOC125975224                                 |          N           |         N          |
| *S. acus/scovelli* vitellogenin                                           |        Y(2-3)        |         N          |

    
# Trimming the reads using TRIMMOMATIC

[Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic) is commonly used for Illumina paired-end and single ended data. There are a lot of different trimming steps that can be performed and parameters that correspond to each step. Currently for version **0.40** of trimmomatic these are the different steps:


1.  **ILLUMINACLIP**: Cut adapter and other illumina-specific seuqences from the read.
    `:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>`
    -   fastaWithAdaptersEtc = specifies the path to a fasta file containing all the adapters, PCR sequences etc.
    -   seedMismatches = the max. mismatch count which will still allow a full match to be performed.
    -   palindromeClipThreshold = how accurate the match between the two 'adapter ligated' reads must be for PE palindrome read alignment
    -   simpleClipThreshold = how accurate the match between any adapterect. sequence must be against a read.
2.  **SLIDINGWINDOW**: Perform a sliding window trimming, cutting once the average quality within the window falls below a designated threshold.
    `:<windowSize>:<requiredQuality>`
    -   windowSize = the number of bases to average across
    -   requiredQuality = the average quality required
3.  **LEADING**: Cut bases off the start of a read, if below a designated threshold quality.
    `:<quality>`
    -   quality = the minimum quality required to keep a base
4.  **TRAILING**: Cut bases off the end of a read, if below a designated threshold quality.
    `:<quality>`
    -   quality = the minimum quality required to keep a base
5.  **CROP**: Cut the read to a specified length.
    `:<length>` 
    -   length = the number of bases to keep, from the start of the read
6.  **HEADCROP**: Cut the specified number of bases from the start of a read.
    `:<length>`
    -   length = the number of bases to remove from the start of the read
7.  **MINLEN**: Drop the read if it is below a specified length.
    `:<length>`
    -   length = the minimum length of reads to be kept

For **paired-end** data, two input files are specified and 4 output files. Two for the **paired** output where both of the reads survived the processing, and two for the corresponding **unpaired** output where one read survived but it's partner did not.

**Trimming occurs in the order which the steps are specified on the command line**. Recommended that if the adapter clipping is required it is done as early as possible.

## Examples of using Trimmomatic

[This video](https://www.youtube.com/watch?v=rinGE9TKgKA) works through a very good example of how to format trimmomatic to be used on command line in the form of a bash script. That example is outlined below.

    #!/bin/bash

    FORWARD=$1
    REVERSE=$2

    nohup trimmomatic PE -threads 8 $FORWARD $REVERSE\
        paired_$FORWARD unpaired_$FORWARD\
        paired_$REVERSE unpaired_$REVERSE\
        ILLUMINACLIP:/opt/Timmomatic-0.32/adapters/NexteraPE-PE.fa:2:30:10\
        SLIDINGWINDOW:4:15 MINLEN:36 &


This script sets it up so that when used you must input two arguments with the first argument being the forward set of reads and the second the reverse. To run trimmomatic you would then type: `script_location/script_name.sh XXX_paired-1.fastq.gz XXX_paired-2.fastq.gz`

## Installing Trimmomatic

Trimmomatic was installed on the RCC through conda in a conda environment called **Trim**.

1.  Create the **Trim** environment. Navigate into the location you want to store this environment [anaconda3/envs] and use the code:
    `conda create -n Trim`
2.  After creating the environment you have to activate it:
    `conda activate Trim/`
3.  Now that the **Trim** environment is activated we can install Trimmomatic-0.39.
    `conda install -c bioconda trimmomatic`
4.  To use trimmomatic, you must make sure the **Trim** environment is first activated.
    -   To activate any environment you create, you must be in the location where you first created it [anaconda3/envs].

### Installing conda

1.  Navigate to the [Anaconda](https://www.anaconda.com/products/distribution) website and scroll down to the Anaconda Installers. Copy the link address for the **64-bit (x86) Installer** under the Linux tab.
2.  Paste that installer link into the terminal with **wget**:
    `wget https://repo.anaconda.com/archive/Anaconda3-2022.05-Linux-x86_64.sh`
3.  Run the anaconda shell script:
    `bash Anaconda3-2022.05-Linux-x86_64.sh`
4.  Conda is now installed. Navigate into the anaconda folder that is generated and run `conda activate`. You can now start creating environments and installing packages.
5.  Set up the conda channels as per the [bioconda documentation](https://bioconda.github.io/) so that you will not have to use `-c bioconda` when trying to install packages. You will only need to do this **once**.
```
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
```

## Running Trimmomatic

### One Pair of Reads

The following script was used to trim **one pair of reads** at a time.

```{bash, eval=FALSE}
#!/bin/bash

FORWARD=$1
REVERSE=$2

trimmomatic PE $FORWARD $REVERSE\
        paired_$FORWARD unpaired_$FORWARD\
        paired_$REVERSE unpaired_$REVERSE\
        ILLUMINACLIP:NexteraPE-PE.fa:2:30:10\
        HEADCROP:12 LEADING:3 TRAILING:10\
        SLIDINGWINDOW:4:15 MINLEN:50
```

How the script was used:
`nohup ./trim_script.sh FLG3F2_E1_S52_L004_R1_001.fastq.gz FLG3F2_E1_S52_L004_R2_001.fastq.gz &`

### Multiple Pairs of Reads

In order to process **multiple** paired reads in a given directory, a **for loop** was created:

```{bash, eval=FALSE}
#!/bin/bash

data=/home/rccuser/20220902_mRNASeq_PE150/ ##This is the location of the raw reads

for fq in ${data}*_R1.fastq.gz
        do
        base=$(basename $fq _R1.fastq.gz)
        echo "Running trimmomatic for ${base}..."
        time trimmomatic PE -threads 16 $fq ${data}${base}_R2.fastq.gz \
                ${data}trimmed/${base}_paired_R1.fastq.gz ${data}trimmed/${base}_unpaired_R1.fastq.gz \    
                ${data}trimmed/${base}_paired_R2.fastq.gz ${data}trimmed/${base}_unpaired_R2.fastq.gz \    
                ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 HEADCROP:12 LEADING:3 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:50
done
```

-   Before running the script, change the `trimmomatic` line to `echo "..."` to make sure all of the variables are working correctly.

-   Then remove the `echo ""` and run the script as `nohup bash trim_script.sh >> trim.out &` entering `disown` when prompted.

## Quality checking trimmed reads with FastQC and MultiQC

### Installing and Using the command fastqc
-   fastqc was installed using conda in the **Trim environment** with `conda install -c bioconda fastqc`.

-   There are many different arguments that can be attached to fastqc outlined in the [Ubuntu fastqc manual](https://manpages.ubuntu.com/manpages/trusty/man1/fastqc.1.html)

```{bash, eval=FALSE}
#!/bin/bash

data=/home/rccuser/20220902_mRNASeq_PE150/

for fq in ${data}trimmed_paired/*.gz #The location of the PAIRED trimmed reads
        do
        base=$(basename $fq)
        echo "Running fastqc for ${base} ..."
        time fastqc $fq -t 16 -o ${data}FastQC_trimmed #Where I want the output files to be stored
done
```

 - This script was run with `nohup bash fastqc_script.sh >> fastqc.out &`, entering `disown` when prompted.
     - The output from this should consist of one .html and one .zip file for each read.
 
### Installing and Using multiqc
 - [Multiqc](https://multiqc.info/) is a program that allows you to combine results from different bioinformatic analyses across many samples into one single report. It was installed using conda in the **QC environment** with `conda install multiqc`.
     - Originally tried installing multiqc on just the base environment, ran into errors such as _"Solving environment: failed with initial frozen solve. Retrying with flexible solve"_. Left it installing overnight with nohup and still had not properly finished by the time I came back in the morning. After creating the QC environment and trying to install it there it took just a couple minutes to install properly.
     - Make sure the bioconda channel is set up in conda prior to installing this package (See **3.2.1**).
   
 - Once installed run `multiqc`, followed by a list of directories you want it to search. Multiqc will scan all of the specified directories and produce a report based on details found in any log file that it recognizes. 
 `multiqc FastQC_trimmed/` 
     - IF you are trying to run multiqc on fastqc files you must make sure you have the .zip files and not just the .html files from fastqc.

 - Once MultiQC has finished, you should have a HTML report file called `multiqc_report.html` that contains all of the results!\
 
### Results of the FastQC/MultiQC Analysis of the TRIMMED Reads

The raw reads were trimmed to fix two major issues: 1) Cleaning up the **per base sequence content** to get rid of the first 10-12 bases that result from hexamer priming that occurs during RNA-seq library preparation in EVERY sequence; 2) Remove all leftover **Illumina Nextera Transposase** adapter sequences. 

 - All samples now PASS the per base sequence content module now due to the trimmomatic `HEADCROP:12`
 - Trimming did NOT affect the GC content module so there are still 45 samples with a WARNING flag for it. 
 - Due to the trimming, all samples have a WARNING flag on the sequence length distribution module. This is expected. 
 - Trimming did not affect the sequence duplication level, all samples still FAIL this module. 
 - All samples that just had the Illumina Nextera Transposase sequence in the overrepresented sequences now PASS that module. The samples that also had sequences with 'no hit' still have a WARNING flag for this module (47 samples including ALL livers). 
 - All samples not PASS the adapter content module due to the trimmomatic `ILLUMINACLIP:NexteraPE-PE.fa:2:30:10` 

<center>

<img src="../imgs/fastqc-status-check-heatmap.png" style="width:650px;"/>

<center> 


#### Comparisson of UNTRIMMED vs TRIMMED read 
On the left hand side are images of what modules from untrimmed reads look like vs trimmed in cases where the modules differed.

<p float="left">

<img src="../imgs/untrimmed_base_statistics.JPG" style="width:350px;"/> <img src="../imgs/trimmed_base_statistics.JPG" style="width:350px;"/>

</p>

<p float="left">

<img src="../imgs/per_base_seq_qual.JPG" style="width:350px;"/> <img src="../imgs/trimmed_per_base_seq_qual.JPG" style="width:350px;"/>

</p>

<p float="left">

<img src="../imgs/per_seq_qual_score.JPG" style="width:350px;"/> <img src="../imgs/trimmed_per_seq_qual_score.JPG" style="width:350px;"/>

</p>

<p float="left">

<img src="../imgs/per_base_seq_content.JPG" style="width:350px;"/> <img src="../imgs/trimmed_per_base_seq_content.JPG" style="width:350px;"/>

</p>

<p float="left">

<img src="../imgs/untrimmed_seq_len_dis.JPG" style="width:350px;"/> <img src="../imgs/trimmed_seq_len_dis.JPG" style="width:350px;"/>

</p>

<p float="left">

<img src="../imgs/adapter_content.JPG" style="width:350px;"/> <img src="../imgs/trimmed_adapter_content.JPG" style="width:350px;"/>

</p>

# De novo Transcriptome Assembly using Trinity
Trinity is a commonly used method for de novo reconstruction of transcriptomes from RNA-seq data when there is no genome to map back to. Within Trinity there are three independent software modules that are applied and used sequentially: **Inchworm**, **Chrysalis**, and **Butterfly**.

## Installing Trinity
Trinity requires the installation of several programs, to get around this Trinity was installed via a Docker in the RCC and must be used within one. **The first time that you use Trinity you must pull the latest Docker image for Trinity** like so:
```{bash, eval=FALSE}
docker pull trinityrnaseq/trinityrnaseq

```

## Running Trinity 
Once started, a successful Trinity run will take days to complete. There are a lot of extra components you can add onto a Trinity run (View the [Trinity User Manual](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-Trinity) for more information) but here are the main components you should specify:

1.  **SeqType**: Signifies the type of reads (fa or fq).
    `--seqType <string>`

2.  **Max Memory**: Suggested max memory to use by Trinity where limiting can be enabled.
    `--max_memory <string>`
    -   This should be based on the available memory for your computer. For example `--max_memory 20G`
    
3.  **Read Locations**: Provide exact file locations for the left and right reads. One or more file names separated by commas, NOT spaces.
    `--left <string>`
    `--right <string>` 
    -   The left and right reads must be paired up sequentially, if XXX32_R1 is mentioned first for the left, XXX32_R2 must also be first for the right.
    -   If running across multiple samples `--samples_file <string>` can be used rather than typing out the location for each individual. See tab about running multiple reads for more details.

4.  **CPU**: Tell Trinity the number of CPUs to use, default is 2.
    `--CPU <int>`
    -   This should be based on the number of cores the computer has and it must be an even number.
    
5.  **Output Directory**: Name of the directory for output.This will be created if it doesn't already exist.
    `--output <string>`
    -   Good to include trinity in the name as a safety precaution! Don't forget to include the full file path.

### One Pair of Reads
```{bash, eval=FALSE}
#!/bin/bash

data=/home/rccuser/20220902_mRNASeq_PE150/trimmed_paired/ ##This is the location of the paired trimmed reads

sudo docker run --rm -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq Trinity \
    --seqType fq --max_memory 50G \
    --left ${data}reads_1.fq.gz \
    --right ${data}reads_2.fq.gz \
    --max_memory 20G --CPU 6 --output ${data}trinity_out_dir
    
```

### Multiple Pairs of Reads
```{bash, eval=FALSE}
#!/bin/bash

data=/home/rccuser/20220902_mRNASeq_PE150/ ##This is the location of our data

sudo docker run --rm -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq Trinity --seqType fq --max_memory 188G \
    --samples_file ${data}trimmed_paired/FU_trinity_samples.txt \
    --CPU 12 \
    --output `pwd`/trinity_out_dir_fuscs
```

## Troubleshooting Trinity
The overarching issues I ran into with Trinity surrounded not having enough memory/RAM and space on the RCC. A total of **4.1 TB** of space was added to the disk paired with **192 GB** of RAM in order to get Trinity to work. 

Below is an example of what an error related to RAM might look like:
<p float="left">

<img src="../imgs/trinity_RAM.png" style="width:650px;"/>


</p>

An identifying feature when trying to determine whether or not an error with Trinity is related to memory/RAM is the line `bash: line 1:    219 Killed` along with errors starting with `Thread X terminated abnormally`. 

### Using screen as a substitute for nohup
Because Trinity must be run in the docker, we are unable to use the `nohup` command in order to send the process into the background. Because of this and the nature of the RCC, you may end up with errors such as `error client_loop: send disconnect: Broken pipe.` meaning our RCC session has automatically timed out. To overcome this issue we can use [screen](https://www.gnu.org/software/screen/manual/html_node/index.html), which acts as a sort of RDC for the terminal. Some of the basics for using screen are outlined below:

1. Log in to the RCC and begin a screen session by entering `screen` and then start your process/commands.

2. If you become disconnected at any point from the RCC, the screen session will still go on. **To reconnect to the session** use the command `screen -r`.

3. To disconnect from your screen session (WITHOUT STOPPING IT) use **Crtl+a d**. All commands to screen are started with "Ctrl+a" followed by a letter.

# Differential Expression Analysis using RSEM and EdgeR
## Installing RSEM and Bowtie2
### RSEM
1.  From the [RSEM website](https://deweylab.github.io/RSEM/) download the lastest version of the RSEM source code and move the folder into the RCC (I used `sftp` to do this through my laptop)
    Once in the RCC, unzip the file using the command `gunzip` 
    ```
    gunzip -h
    gunzip -d <filename>
    tar -xf <filename>
    ```

2.  Next you will need to compile RSEM, move into the **unzipped** file and in the command line type `sudo make` to create the executable file. Then to install RSEM run `sudo make install`
    This did not work immediately for me because I didn't have the necessary compiler [GCC](https://linuxconfig.org/how-to-install-g-the-c-compiler-on-ubuntu-20-04-lts-focal-fossa-linux), had to install that first.
    ```
    sudo apt install build-essential
    sudo spt install zlib1g-dev #another error while installing GCC complier, could not find zlib1g-dev
    ```

3.  Add RSEM to your path so that you can run any RSEM programs from any directory. If you have moved, change back into the home directory, run `ls -a` to check that you have a .bashrc script.
    `nano .bashrc`

4.  **Without modifying pre-existing parts of the file**, add the following to the very end:
    ```
    #RSEM
    export PATH=$PATH:location_of_RSEM/RSEM-version.number/
    ```
    
### Bowtie2
Bowtie2 was installed using **conda**. An environment was created with `conda create trinity` and then bowtie2 was installed using `conda install bowtie2`. **To run any further analyses that include bowtie2, the trinity environment MUST be activated** with `conda activate trinity`.

## Running RSEM
RSEM was used for differential expression analysis following the pipeline _rsem-prepare-reference_, _rsem-calculate-expression_, and _rsem-generate-data-matrix_. The data was not pre-processed or filtered after trinity prior to going through RSEM. The data will be filtered  during the statistical analysis portion. To run RSEM about 5-15 GB memory is needed. The _.fasta_ file should be a resulting assembled transcriptome from Trinity, and the read files you use should be **trimmed, paired** reads from Trimmomatic. Each pair of read files will be mapped one at a time against the reference transcriptome. It is easiest if the scripts outlined below are run in a folder that contains your trimmed reads and Trinity output.

### Preparing the Reference 
In the first stage of RSEM you will prepare the reference sequence to run RSEM and build indices for bowtie2 using [rsem-prepare-reference](https://deweylab.biostat.wisc.edu/rsem/rsem-prepare-reference.html). The input for this command is the reference transcriptome from Trinity ( _.fasta_ ) and the output is the RSEM references. You will do this **once per species**.

```{bash, eval=FALSE}
#!/bin/bash

trinity_output_file=$1 ##the .Trinity.fasta file
RSEM_output=$2 ##Ex. TrinityRefFuscus

rsem-prepare-reference --bowtie2 $trinity_output_file $RSEM_output

```
This script was then run using **nohup** as: `nohup bash rsem_ref.sh trinity_out_dir_fuscus.Trinity.fasta TrinityRefFuscus > rsem1_FU.log 2>&1 &`, changing the arguements when necessary.
    By adding the `2>&1` to the command we save BOTH the standard output and the errors into a .log file

### Calculating the expression levels
The expression levels are calculated for **each set of paired reads separately** using [rsem-calculate-expression](https://deweylab.github.io/RSEM/rsem-calculate-expression.html). The input is your RNA-seq reads ( _.fastq_ ) and the output is a genes.results.file, an isoforms.results.file, and .transcript.bam,file and a .log that contains information about the mapping success for every sample.

```{bash, eval=FALSE}
#!/bin/bash

data=/home/rccuser/20220902_mRNASeq_PE150/trimmed_paired/ ##The location on the trimmed paired reads
FU-or_FL=$1 #Which species' reads are you wanting?
rsem_reference_output_file=$2 ##Ex TrinityRefFuscus

for fq in ${data}$1*_R1.fastq.gz
	do
	base=$(basename $fq _R1.fastq.gz)
	echo "Running RSEM for ${base}..."
	time rsem-calculate-expression -p 16 --paired-end --bowtie2 \
		$fq ${data}${base}_R2.fastq.gz \
		$2 ${base}
done
```
This command took several days to run and therefore **nohup** was used: `nohup bash rsem_expression.sh FU TrinityRefFuscus > rsem2_FU.log 2>&1 &`, changing the arguments when necessary.
    There will be multiple outputs from the RSEM reference command with the same prefix (whatever you named it in rsem-prepare-reference) but with multiple extensions (Ex. .idx.fa, .seq, .ti, etc.), when you run rsem-calculate-expression you will just give it the prefix name as it will use all of the file types.

### Generating the matrix 
For this command, rsem-generate-data-matrix, samples sould be listed by group, so if you are interested in differential expression between males and females tissue types then all of the males should be listed first, followed by all of the females. The input files are all of the genes.results.file and the output is a counts.matrix.file. This is the file that will be used in downstream analysis. This will be done **once per species**. 

#### Fuscus 
```{bash, eval=FALSE}
#!/bin/bash

data=/home/rccuser/20220902_mRNASeq_PE150/trimmed_paired/

rsem-generate-data-matrix ${data}FUG10M2_paired.genes.results ${data}FUG11M2_paired.genes.results ${data}FUG11M4_paired.genes.results ${data}FUG12M1_paired.genes.results ${data}FUG15M5_paired.genes.results \
	${data}FUL10M2_paired.genes.results ${data}FUL11M2_paired.genes.results ${data}FUL11M4_paired.genes.results ${data}FUL12M1_paired.genes.results ${data}FUL15M5_paired.genes.results \
	${data}FUT10M2_paired.genes.results ${data}FUT11M2_paired.genes.results ${data}FUT11M4_paired.genes.results ${data}FUT12M1_paired.genes.results ${data}FUT15M5_paired.genes.results \
	${data}FUG11F1_paired.genes.results ${data}FUG13F1_paired.genes.results ${data}FUG13F4_paired.genes.results ${data}FUG2F2_paired.genes.results ${data}FUG3F2_paired.genes.results \
	${data}FUL11F1_paired.genes.results ${data}FUL13F1_paired.genes.results ${data}FUL13F4_paired.genes.results ${data}FUL2F2_paired.genes.results ${data}FUL3F2_paired.genes.results \
	${data}FUO11F1_paired.genes.results ${data}FUO13F1_paired.genes.results ${data}FUO13F4_paired.genes.results ${data}FUO2F2_paired.genes.results ${data}FUO3F2_paired.genes.results \
	> RSEM_fuscus_digi.counts.matrix

```

#### Floridae 
```{bash, eval=FALSE}
#!/bin/bash

data=/home/rccuser/20220902_mRNASeq_PE150/trimmed_paired/

rsem-generate-data-matrix ${data}FLG3M5_paired.genes.results ${data}FLG3M7_paired.genes.results ${data}FLG3M8_paired.genes.results ${data}FLG4M3_paired.genes.results ${data}FLG4M4_paired.genes.results \
        ${data}FLL3M5_paired.genes.results ${data}FLL3M7_paired.genes.results ${data}FLL3M8_paired.genes.results ${data}FLL4M3_paired.genes.results ${data}FLL4M4_paired.genes.results \
        ${data}FLT2M3_paired.genes.results ${data}FLT3M5_paired.genes.results ${data}FLT4M4_paired.genes.results ${data}FLT5M3_paired.genes.results ${data}FLT8M7_paired.genes.results \
        ${data}FLG2F7_paired.genes.results ${data}FLG3F1_paired.genes.results ${data}FLG3F2_paired.genes.results ${data}FLG8F3_paired.genes.results \
        ${data}FLL2F7_paired.genes.results ${data}FLL3F1_paired.genes.results ${data}FLL3F2_paired.genes.results ${data}FLL3F4_paired.genes.results ${data}FLL8F3_paired.genes.results \
        ${data}FLO2F7_paired.genes.results ${data}FLO3F1_paired.genes.results ${data}FLO3F2_paired.genes.results ${data}FLO3F4_paired.genes.results ${data}FLO8F3_paired.genes.results \
        > RSEM_floridae_digi.counts.matrix

```

# Accessing Tools installed within Trinity - Checking the quality of the alignment
Because of the way Trinity was installed on the RCC, all of the additional tools have to be accessed through the docker. There is some info on the [Dockerized Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-in-Docker) webpage about how to do this. However, it should be noted the **path to the tools** in Trinity for this RCC differs from what is on the website. Instead of `/usr/local/bin/trinityrnaseq/util/...` it is `/usr/local/bin/util/...`.

## Trinity Transcriptome Contig Nx Statistics
The Nx statistic is calculated based on the lengths of the assembled transcriptome contigs. This statistic tells us that at least X% of assembled transcript nucleotides are found in contigs that are of Nx length.
 - Traditionally, you compute the N50 where at least **half** of all assembled bases are in transcript contigs of at LEAST the N50 length value.
 - This can be done using the built in TrinityStats.pl script within trinity.
 - The input for this function is the resulting *.Trinity.fasta* file from the Trinity run.

```{bash, eval=FALSE}
#!/bin/bash

data=/home/rccuser/20220902_mRNASeq_PE150/trimmed_paired/
trinity_fasta_file=$1
trinity_stats_output=$2

sudo docker run -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq /usr/local/bin/util/TrinityStats.pl \
        ${data}$trinity_fasta_file \
        >> $trinity_stats_output

```
This script was run as `bash trinity_stats.sh trinity_out_dir_floridae.Trinity.fasta floridae_trinity_stats.txt`

The contig N50 values are often exaggerated due to an assembly program generating too many transcript isoforms. To help mitigate this, we can look at the Nx values based on using only the single longest isoform per 'gene'. These Nx values are often lower than the Nx stats based on all assembled contigs.
### *Syngnathus fuscus* N50 results
################################  
##Counts of transcripts, etc.  
################################  
Total trinity 'genes':	399585  
Total trinity transcripts: 645009  
Percent GC: 46.42  

########################################  
Stats based on ALL transcript contigs:  
########################################  

	Contig N10: 6786
	Contig N20: 4722
	Contig N30: 3520
	Contig N40: 2628
	Contig N50: 1915

	Median contig length: 425
	Average contig: 936.35
	Total assembled bases: 603957249


#####################################################  
##Stats based on ONLY LONGEST ISOFORM per 'GENE':  
#####################################################  

	Contig N10: 5180
	Contig N20: 3180
	Contig N30: 1976
	Contig N40: 1276
	Contig N50: 857

	Median contig length: 338
	Average contig: 614.57
	Total assembled bases: 245572519

For *S. fuscus* we can see that 50% of the assembles bases are found in transcript contigs that are at least 857 bases in length (going off of the stats based on the longest isoform).

### *Syngnathus floridae* N50 results
################################  
##Counts of transcripts, etc.  
################################  
Total trinity 'genes':	347758  
Total trinity transcripts: 606899  
Percent GC: 46.42  

########################################  
Stats based on ALL transcript contigs:  
########################################  

	Contig N10: 6684
	Contig N20: 4682
	Contig N30: 3513
	Contig N40: 2663
	Contig N50: 1977

	Median contig length: 436
	Average contig: 961.91
	Total assembled bases: 583784427


#####################################################  
##Stats based on ONLY LONGEST ISOFORM per 'GENE':  
#####################################################  

	Contig N10: 5450
	Contig N20: 3559
	Contig N30: 2412
	Contig N40: 1592
	Contig N50: 1028

	Median contig length: 339
	Average contig: 655.95
	Total assembled bases: 228111197

For *S. floridae* we can see that 50% of the assembles bases are found in transcript contigs that are at least 1028 bases in length (going off of the stats based on the longest isoform).

## Contig Ex90N50 Statistic and Ex90 Gene Count
This statistic is an alternative to the contig Nx statistic that can be considered more appropriate for transcriptome assembly data. With the test, the N50 statistic is computed as it is above BUT limited to the top most highly expressed genes that represent x% of the total normalized expression data. The gene expression is taken as the sum of the transcript isoform expression and the gene length is computed as the expression-weighted mean of the isoform lengths. Prior to calculating this you must have first performed the transcript abundance estimation.

### Building Transcript and Gene Expression Matrices within Trinity
To be able to calculate the Ex90N50 statistic we need to make sure out gene espression matrix is in the right format. To ensure this we can use the *abundance_estimates_to_matrix.pl* script from [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Transcript-Quantification). This script takes the *.gene.results.file from RSEM* (or any other transcript abundance estimation software) and *.fasta.gene_trans_map* from your original Trinity run as input and gives you a RSEM.isoform and RSEM.gene matrices as outputs.

```{bash, eval=FALSE}
#!/bin/bash

data=/home/rccuser/20220902_mRNASeq_PE150/trimmed_paired/ ##the location of the .genes.results.file

sudo docker run -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq /bin/sh -c "cd /home/rccuser/20220902_mRNASeq_PE150/trimmed_paired/ && /usr/local/bin/util/abundance_estimates_to_matrix.pl \
	--est_method RSEM \
	--gene_trans_map ${data}trinity_out_dir_fuscus.Trinity.fasta.gene_trans_map \
	${data}FUG10M2_paired.genes.results ${data}FUG11M2_paired.genes.results ${data}FUG11M4_paired.genes.results ${data}FUG12M1_paired.genes.results ${data}FUG15M5_paired.genes.results \
	${data}FUL10M2_paired.genes.results ${data}FUL11M2_paired.genes.results ${data}FUL11M4_paired.genes.results ${data}FUL12M1_paired.genes.results ${data}FUL15M5_paired.genes.results \
	${data}FUT10M2_paired.genes.results ${data}FUT11M2_paired.genes.results ${data}FUT11M4_paired.genes.results ${data}FUT12M1_paired.genes.results ${data}FUT15M5_paired.genes.results \
	${data}FUG11F1_paired.genes.results ${data}FUG13F1_paired.genes.results ${data}FUG13F4_paired.genes.results ${data}FUG2F2_paired.genes.results ${data}FUG3F2_paired.genes.results \
	${data}FUL11F1_paired.genes.results ${data}FUL13F1_paired.genes.results ${data}FUL13F4_paired.genes.results ${data}FUL2F2_paired.genes.results ${data}FUL3F2_paired.genes.results \
	${data}FUO11F1_paired.genes.results ${data}FUO13F1_paired.genes.results ${data}FUO13F4_paired.genes.results ${data}FUO2F2_paired.genes.results ${data}FUO3F2_paired.genes.results"

```
When the script was originally run no output files could be found. To fix this a shell wrapper (`/bin/sh -c`) was added as the main container process to the command. This issue was likely that the docker was running under its root directory and not the current working directory.

### Calculating the ExN50 Statistic

```{bash, eval=FALSE}
#!/bin/bash

data=/home/rccuser/20220902_mRNASeq_PE150/trimmed_paired/
isoform_expression_matrix=$1
trinity_fasta_output=$2
ExN50_output=$3

sudo docker run -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq /usr/local/bin/util/misc/contig_ExN50_statistic.pl \
        ${data}$isoform_expression_matrix \
        ${data}$trinity_fasta_output | tee $ExN50_output

```
This script was run as the following: `bash exn50.sh RSEM.isoform.TMM.EXPR.matrix trinity_out_dir_fuscus.Trinity.fasta ExN50.stats`. This results in a table that looks like this:

| Ex      | ExN50   | num_genes |
|:--------|:--------|:----------|
|    1    |   304   |     1     |
|    3    |   304   |     2     |
|    4    |   304   |     3     |
|    5    |   304   |     4     |
|   ...   |   ...   |    ...    |
|   99    |   801   |  306213   |
|   100   |   677   |  399585   |

We can then plot the Ex value (first column) against the ExN50 values:

`sudo docker run -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq /bin/sh -c "cd /home/rccuser/20220902_mRNASeq_PE150/trimmed_paired/ && /usr/local/bin/util/misc/plot_ExN50_statistic.Rscript /home/rccuser/20220902_mRNASeq_PE150/trimmed_paired/ExN50.stats" `  
`xpdf ExN50.stats.plot.pdf`