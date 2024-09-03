Analyzing *Syngnathus fuscus* RNAseq Data from MSU
================
Coley Tosto
2024-09-03

- <a href="#pre-assembly-quality-control-and-filtering"
  id="toc-pre-assembly-quality-control-and-filtering">Pre-assembly Quality
  Control and Filtering</a>
  - <a href="#trimming-the-raw-reads-with-trimmomatic"
    id="toc-trimming-the-raw-reads-with-trimmomatic">Trimming the raw reads
    with Trimmomatic</a>
  - <a href="#using-kraken2-to-remove-biological-contamination"
    id="toc-using-kraken2-to-remove-biological-contamination">Using Kraken2
    to remove biological contamination</a>
  - <a href="#using-sortmerna-to-remove-rrna-contamination"
    id="toc-using-sortmerna-to-remove-rrna-contamination">Using SortMeRNA to
    remove rRNA contamination</a>
  - <a href="#doing-a-k-mer-based-correction-with-rcorrector"
    id="toc-doing-a-k-mer-based-correction-with-rcorrector">Doing a k-mer
    based correction with RCorrector</a>
- <a href="#checking-quality-of-trimmed-and-filtered-reads"
  id="toc-checking-quality-of-trimmed-and-filtered-reads">Checking quality
  of trimmed and filtered reads</a>
- <a href="#mapping-reads-to-s-scovelli-genome"
  id="toc-mapping-reads-to-s-scovelli-genome">Mapping reads to <em>S.
  scovelli</em> genome</a>
- <a href="#alignment-and-abundance-estimations"
  id="toc-alignment-and-abundance-estimations">Alignment and abundance
  estimations</a>
  - <a href="#generate-the-index" id="toc-generate-the-index">Generate the
    index</a>
  - <a href="#quantify-the-samples" id="toc-quantify-the-samples">Quantify
    the samples</a>
- <a href="#assembly-thinning-and-redundancy-reduction---supertranscripts"
  id="toc-assembly-thinning-and-redundancy-reduction---supertranscripts">Assembly
  thinning and redundancy reduction - SuperTranscripts</a>
  - <a href="#re-evaluating-the-quality-of-the-assembly"
    id="toc-re-evaluating-the-quality-of-the-assembly">Re-evaluating the
    quality of the assembly</a>
    - <a href="#busco" id="toc-busco">BUSCO</a>
    - <a href="#examining-rna-seq-read-representation-via-salmon"
      id="toc-examining-rna-seq-read-representation-via-salmon">Examining
      RNA-Seq read representation (via Salmon)</a>
- <a href="#prepping-for-differential-expression-analysis"
  id="toc-prepping-for-differential-expression-analysis">Prepping for
  differential expression analysis</a>
  - <a href="#importing-transcript-abundance-with-tximport"
    id="toc-importing-transcript-abundance-with-tximport">Importing
    transcript abundance with tximport</a>

# Pre-assembly Quality Control and Filtering

These reads were obtained from the **RTSF Genomics Core** at Michigan
State University. They were downloaded with the File Transfer Protocol
provided on their
[website.](https://rtsf.natsci.msu.edu/genomics/data-retrieval/)

Quality scores were provided by MSU, so FastQC was not run on the raw
reads. The provided scores were used to assist in choosing the settings
for trimming and filtering. The trimming and filtering for these reads
will follow the pipeline that was laid out for *Syngnathus floridae*.
Any changes to scripts or steps will be highlighted in this document.
Analysis was conducted in a Remote Computing Cluster (RCC) at the
University of Canterbury.

For *Syngnathus fuscus* we started with an average of 50 ± 6.9 million
reads per sample for a total of 2.99 billion reads.

## Trimming the raw reads with Trimmomatic

Trimmomatic was install via a conda environment `Trim` on the RCC.
`trimmomatic v0.39` was used for the following script.

``` bash
#!/bin/bash

#Create the arguements
input_dir=$1 #This should be the location of the RAW reads to be trimmed

for fq in $1*_R1.fastq.gz
    do
    base=$(basename $fq _R1.fastq.gz)
    echo "Running trimmomatic for ${base}..."
    time trimmomatic PE -threads 16 $fq $1${base}_R2.fastq.gz \
        $1trimmed/${base}_paired_R1.fastq.gz $1trimmed/${base}_unpaired_R1.fastq.gz \
        $1trimmed/${base}_paired_R2.fastq.gz $1trimmed/${base}_unpaired_R2.fastq.gz \
        ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 HEADCROP:12 LEADING:3 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:50
done
```

- Before running the script, change the `trimmomatic` line to
  `echo "..."` to make sure all of the variables are working correctly.

- Then remove the `echo ""` and run the script as
  `nohup bash trim_script.sh > trim.out 2>&1 &`.  

- The **NexteraPE-PE.fa** file was pulled from the [trimmomatic
  github](https://github.com/usadellab/Trimmomatic/tree/main/adapters)
  using
  `wget https://github.com/usadellab/Trimmomatic/blob/main/adapters/NexteraPE-PE.fa`

## Using Kraken2 to remove biological contamination

Kraken2 was installed in the conda environment `kraken2` on the RCC.
`kraken2 v2.1.2` was used.

``` bash
#!/bin/bash

#Create arguments
ref_fasta=$1 #desired reference database
input_dir=$2 #Input directory with the location on the reads
output_dir=$3 #Desired location for the output

## Loop through all pairs of reads in the input directory
for pair in $2/*_R1.fq.gz
        do

        #Extract the sample name from the file name
        sample=$(basename $pair _R1.fq.gz)

        ##Echo the sample name it is currently running
        echo "Running Kraken2 for ${sample}..."

        #Define the paths to the input and output files for this sample
        input1=$2/${sample}_R1.fq.gz
        input2=$2/${sample}_R2.fq.gz

        #Run Kraken2 on this pair of reads
        time kraken2 --threads 16 --db $1 --paired $input1 $input2 --unclassified-out $3/${sample}#.fq --report $3/${sample}.log

done
```

This script was run as
`nohup bash bash_scripts/kraken2.sh ../kraken2_pluspfp/ fuscus_trimmed fuscus_nobio > kraken2.log 2>&1 &`

The kraken2 database used for this analysis included the standard
database (archea, bacteria, viral) plus plant, fungi, and protozoan
datbases. Only reads that did not map back to these databases were
retained.

## Using SortMeRNA to remove rRNA contamination

SortmeRNA was installed via a conda environment `sortmerna` on the
RCC.`sortmerna v4.3.6` was used

``` bash
#!/bin/bash

#Create arguments
input_dir=$1 #location of the reads
ref_fasta=$2 #desired reference fasta
output_dir_rrna=$3 #desired location for the reads that are rRNA
output_dir_norrna=$4 #desired location for the reads that are NOT rRNA

## Loop through all pairs of reads in the input directory
for pair in $1/*_R1.fastq.gz
        do

        #Extract the sample name from the file name
        sample=$(basename $pair _R1.fastq.gz)

        #Extract Fish ID
        ID=$(basename $pair _paired_R1.fastq.gz)

        ##Echo the sample name it is currently running
        echo "Running SortMeRNA for ${ID}..."

        #Define the paths to the input and output files for this sample
        input1=$1/${sample}_R1.fastq.gz
        input2=$1/${sample}_R2.fastq.gz

        #Run SortMeRNA on this pair of reads
        time sortmerna --threads 16 --ref $2 --reads $input1 --reads $input2 --fastx --aligned $3/${ID} --other $4/${ID} --out2

        #Remove SortMeRNA intermediate files before running again
        rm -r /home/rccuser/sortmerna/run/kvdb

done
```

This script was run as
`nohup bash bash_scripts/sortmerna.sh fuscus_trimmed ../rRNA_databases_v4.3.4/smr_v4.3_fast_db.fasta fuscus_rrna fuscus_norrna > sortmerna.log 2>&1 &`.

The chosed reference FASTA file contains **a subset of sequences** from
the default SortMeRNA database. Here, the number of sequences in each
database is reduced to improve the speed of the analysis.

## Doing a k-mer based correction with RCorrector

The Rcorrector github repo was cloned and Rcorrector was installed in
the /shared folder on the RCC.

``` bash
#!/bin/bash

#Create arguments
rcorrector_path=$1 #Path to the run_rcorrector.pl file
input_dir=$2 #Path to the input directory containing the reads
output_dir=$3 #Path to the desired output location

##Loop through all pairs of reads in the directory
for pair in $2/*_R1.fq
    do

    #Extract sample name from the file name
    sample=$(basename $pair _R1.fq)

    #Echo the sample name that is currently running
    echo "Running rcorrector for ${sample} ..."

    #Run rcorrector on this pair of reads
    time perl $1/run_rcorrector.pl -t 16 -1 $2/${sample}_R1.fq -2 $2/${sample}_R2.fq -od $3

done
```

This script was run as
`nohup bash bash_scripts/rcor.sh ../../rcorrector fuscus_nobio fuscus_kmer_corrected/ > rcor.log 2>&1 &`.
After this step the fasta files were g-zipped.

# Checking quality of trimmed and filtered reads

The quality of the reads once they finished going through the filtering
and trimming pipeline outlined above was assessed with FastQC and the
results were compiled using MultiQC.

``` bash
#!/bin/bash

input_dir=$1 #location of the reads
output_dir=$2 #name of the desired output directory

for fq in $1/*.gz
    do
    base=$(basename $fq)
    echo "Running fastqc for ${base} ..."
    time fastqc $fq -t 16 -o $2
done
```

This script was run as
`nohup bash bash_scripts/fastqc_script.sh fuscus_kmer_corrected fuscus_FastQC > fastqc.log 2>&1 &`.
Once finished MultiQC was run: `multiqc fuscus_FastQC`.

From the General stats reported in the MultiQC report we can see that
following this filtering/trimming process we end up with an average of
42 ± 6.2 million reads per sample for a total of 2.54 billion reads.

# Mapping reads to *S. scovelli* genome

There is currently a genome that exists for a closely related species
within the same genus as *S. fuscus*. As a result of this, rather than
constructing a *de novo* transcriptome assembly, I am going to use that
genome to map the reads back to.

There are several programs that have the capability of doing this. The
two main things that I will be looking to increase are run time (how
much space and time will running these programs take up) and accuracy
(what proportion of my reads are mapping back to the genome). Some of
the alignment programs include:

- `STAR`
- `HISAT2`
- `Bowtie2`
- `BWA`

A previous study comparing these different alignment software ([Musich
et al.,
2021](https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2021.657240/full))
highlighted that `Bowtie2` with local alignment mode and `BWA`had the
highest alignment rate followed by `STAR` and lastly `HISAT2`. However,
`Bowtie2` with the local alignment took the longest to run followed by
`STAR` and `BWA` which were similar and lastly, `HISAT2` which was by
far the fastest. For `HISAT2` and `STAR` it appears that they are able
to better handle reads that are \>1,000bp compared to the other
aligners.

Outside of these components, one other important component to keep in
mind is the fact that one data source is genomic while the other is
transcriptomic. Where this presents potential issues is in the
fragmentation associated with transcriptomic data, which when
unrecognized, translates to large gaps in comparison to the genome. To
mediate this, some of the aligners have been designed to recognize this,
including `STAR` and `HISAT2`. They are able to do so by allowing the
option of inputting a transcriptome alongside the reference.

# Alignment and abundance estimations

Salmon was installed locally on the RCC. It must be run as
`/usr/local/bin/salmon`. `Salmon v 1.10.1` was used.

## Generate the index

`salmon index` was used for this step. The index is built only **once**
per transcriptome.

``` bash
#!/bin/bash

#Create arguments
transcriptome_file=$1 #This is the output .Trinity.fasta file from the txome assembly
desired_index_name=$2
kmer_length=$3 #Recommended 31 on the Salmon user guide

/usr/local/bin/salmon index -t $transcriptome_file \
            -i $desired_index_name \
            -k $kmer_length -p 16
```

The script was run as
`nohup bash bash_scripts/salmon_txome_indicies.sh trinity_out_dir_fuscus_June2023.Trinity.fasta fuscus_salmon_index 31 > salmon_index.log 2>&1 &`.
A kmer length of 31 was used upon recommendation from the Salmon
website.

## Quantify the samples

After building the index this is run on every individual read pair.

``` bash
#!/bin/bash

#Create arguments
input_dir=$1
index_file=$2
output_dir=$3

for fq in $input_dir*_R1.fq.gz
    do

    #Extract sample name from the file
    sample=$(basename $fq _R1.fq.gz)

    #Echo sample name that is currently running
    echo "Processing sample ${sample} ..."

    #Quantify this pair of reads
    /usr/local/bin/salmon quant -i $index_file -l A \
        -1 $input_dir${sample}_R1.fq.gz \
        -2 $input_dir${sample}_R2.fq.gz \
        -p 16 -o $output_dir${sample}_quant

done
```

This script was run as
`nohup bash bash_scripts/salmon_quant.sh fuscus_kmer_corrected/ fuscus_salmon_index/ fuscus_salmon_quant/ > salmon_quant.log 2>&1 &`.

The `quant.sf` and `.log` files were moved out of their nested location
to a shared folder and renamed to include the sample name.

``` bash
#!/bin/bash

#Create arguements
input_dir=$1 #This is the location of all the Salmon quant folders that were generated
out_dir1=$2 #Desired location for the new renamed quant.sf files
out_dir2=$3 #Desired location for the new renamed salmon_quant.log files

for dir in $input_dir*_quant
    do

    echo "$dir"

    #Extract the sample name from the directory
    sample=$(basename $dir _quant)

    #Echo the sample that is being processed
    echo "Rename sample ${sample} ..."

    #Rename the quant.sf file
    mv $dir/quant.sf $out_dir1/${sample}_quant.sf

    #Rename the log file
    mv $dir/logs/salmon_quant.log $out_dir2/${sample}_salmon_quant.log

done
```

The mapping rate was pulled from all of the log files with
`grep "Mapping rate" *log > FUmap.txt` and the results are below.
Generally, we want to see a mapping rate of at least 80% for our
samples.

    FUG10M2_salmon_quant.log:[2023-10-09 23:13:27.146] [jointLog] [info] Mapping rate = 99.8809%
    FUG11F1_salmon_quant.log:[2023-10-09 23:43:18.088] [jointLog] [info] Mapping rate = 99.8791%
    FUG11M2_salmon_quant.log:[2023-10-10 00:14:21.567] [jointLog] [info] Mapping rate = 99.8895%
    FUG11M4_salmon_quant.log:[2023-10-10 00:45:11.435] [jointLog] [info] Mapping rate = 99.5617%
    FUG12M1_salmon_quant.log:[2023-10-10 01:04:04.579] [jointLog] [info] Mapping rate = 99.8942%
    FUG13F1_salmon_quant.log:[2023-10-10 01:43:52.069] [jointLog] [info] Mapping rate = 99.8912%
    FUG13F4_salmon_quant.log:[2023-10-10 02:18:06.249] [jointLog] [info] Mapping rate = 99.9204%
    FUG15M5_salmon_quant.log:[2023-10-10 02:51:42.322] [jointLog] [info] Mapping rate = 99.93%
    FUG2F2_salmon_quant.log:[2023-10-10 03:17:58.246] [jointLog] [info] Mapping rate = 99.8933%
    FUG3F2_salmon_quant.log:[2023-10-10 03:41:01.692] [jointLog] [info] Mapping rate = 99.8957%
    FUL10M2_salmon_quant.log:[2023-10-10 04:12:32.747] [jointLog] [info] Mapping rate = 99.9229%
    FUL11F1_salmon_quant.log:[2023-10-10 04:28:38.157] [jointLog] [info] Mapping rate = 99.935%
    FUL11M2_salmon_quant.log:[2023-10-10 04:47:41.021] [jointLog] [info] Mapping rate = 99.9017%
    FUL11M4_salmon_quant.log:[2023-10-10 05:02:14.331] [jointLog] [info] Mapping rate = 99.9295%
    FUL12M1_salmon_quant.log:[2023-10-10 05:18:34.970] [jointLog] [info] Mapping rate = 99.9099%
    FUL13F1_salmon_quant.log:[2023-10-10 05:39:55.564] [jointLog] [info] Mapping rate = 99.9373%
    FUL13F4_salmon_quant.log:[2023-10-10 05:58:09.029] [jointLog] [info] Mapping rate = 99.9419%
    FUL15M5_salmon_quant.log:[2023-10-10 06:14:57.386] [jointLog] [info] Mapping rate = 99.9355%
    FUL2F2_salmon_quant.log:[2023-10-10 06:31:12.803] [jointLog] [info] Mapping rate = 78.6675%
    FUL3F2_salmon_quant.log:[2023-10-10 06:40:34.545] [jointLog] [info] Mapping rate = 99.94%
    FUO11F1_salmon_quant.log:[2023-10-10 07:02:53.060] [jointLog] [info] Mapping rate = 99.8901%
    FUO13F1_salmon_quant.log:[2023-10-10 07:38:11.441] [jointLog] [info] Mapping rate = 99.905%
    FUO13F4_salmon_quant.log:[2023-10-10 08:11:16.544] [jointLog] [info] Mapping rate = 99.9083%
    FUO2F2_salmon_quant.log:[2023-10-10 08:58:16.701] [jointLog] [info] Mapping rate = 99.9399%
    FUO3F2_salmon_quant.log:[2023-10-10 09:48:39.310] [jointLog] [info] Mapping rate = 99.832%
    FUT10M2_salmon_quant.log:[2023-10-10 10:28:14.802] [jointLog] [info] Mapping rate = 99.9139%
    FUT11M2_salmon_quant.log:[2023-10-10 11:10:17.014] [jointLog] [info] Mapping rate = 99.9051%
    FUT11M4_salmon_quant.log:[2023-10-10 11:44:42.531] [jointLog] [info] Mapping rate = 99.8352%
    FUT12M1_salmon_quant.log:[2023-10-10 12:10:03.402] [jointLog] [info] Mapping rate = 99.872%
    FUT15M5_salmon_quant.log:[2023-10-10 12:38:18.215] [jointLog] [info] Mapping rate = 99.8842%

A very high mapping rate is seen for all but one sample, **FUL2F2**. If
that stays the same after the assembly thinning I will go back and
investigate why that may be. Overall, very happy with these results.

# Assembly thinning and redundancy reduction - SuperTranscripts

SuperTranscripts were generated to reduce redundancy while still
retaining all of the sequence information. These were generated with the
tool installed with Trinity.

``` bash
#!/bin/bash

trinity_out=$1 #The .fasta file generated by the Trinity run


sudo docker run -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq /bin/sh -c "cd /home/rccuser/shared/coley_files/ && /usr/local/bin/Analysis/SuperTranscripts/Trinity_gene_splice_modeler.py \
        --incl_malign \
        --trinity_fasta $trinity_out"
```

This script was run in a screen session as
`bash bash_scripts/trin_super_transcripts.sh trinity_out_dir_fuscus_June2023.Trinity.fasta`.
It took about 11 minutes to run.

## Re-evaluating the quality of the assembly

### BUSCO

#### Short summary

    # BUSCO version is: 5.2.2 
    # The lineage dataset is: actinopterygii_odb10 (Creation date: 2021-02-19, number of genomes: 26, number of BUSCOs: 3640)
    # Summarized benchmarking in BUSCO notation for file /home/rccuser/shared/coley_files/trinity_supertran_fuscus.fasta
    # BUSCO was run in mode: transcriptome

        ***** Results: *****

        C:92.5%[S:90.0%,D:2.5%],F:3.2%,M:4.3%,n:3640       
        3367    Complete BUSCOs (C)            
        3277    Complete and single-copy BUSCOs (S)    
        90  Complete and duplicated BUSCOs (D)     
        117 Fragmented BUSCOs (F)              
        156 Missing BUSCOs (M)             
        3640    Total BUSCO groups searched        

    Dependencies and versions:
        hmmsearch: 3.1
        metaeuk: 6.a5d39d9

#### Graphing results

<p float="center">

<img src="../imgs/busco_fuscus_super.png" style="width:650px;"/>

</p>

### Examining RNA-Seq read representation (via Salmon)

#### Mapping Rate

| Sample  | Trinity Assembly | SuperTranscript Assembly |
|:-------:|:----------------:|:------------------------:|
| FUG10M2 |     99.8809%     |         96.3946%         |
| FUG11F1 |     99.8791%     |         96.3215%         |
| FUG11M2 |     99.8895%     |         96.2211%         |
| FUG11M4 |     99.5617%     |         93.7851%         |
| FUG12M1 |     99.8942%     |         96.3071%         |
| FUG13F1 |     99.8912%     |         96.3934%         |
| FUG13F4 |     99.9204%     |         96.4825%         |
| FUG15M5 |      99.93%      |         96.5262%         |
| FUG2F2  |     99.8933%     |         96.4757%         |
| FUG3F2  |     99.8957%     |         96.4263%         |
| FUL10M2 |     99.9229%     |         96.1121%         |
| FUL11F1 |     99.935%      |         96.3962%         |
| FUL11M2 |     99.9017%     |         96.771%          |
| FUL11M4 |     99.9295%     |         96.8967%         |
| FUL12M1 |     99.9099%     |         96.3373%         |
| FUL13F1 |     99.9373%     |         96.2684%         |
| FUL13F4 |     99.9419%     |         96.8745%         |
| FUL15M5 |     99.9355%     |         96.7522%         |
| FUL2F2  |     78.6675%     |         76.6031%         |
| FUL3F2  |      99.94%      |         96.163%          |
| FUO11F1 |     99.8901%     |         96.0597%         |
| FUO13F1 |     99.905%      |         96.1805%         |
| FUO13F4 |     99.9083%     |         96.1208%         |
| FUO2F2  |     99.9399%     |         96.4821%         |
| FUO3F2  |     99.832%      |         95.8775%         |
| FUT10M2 |     99.9139%     |         95.9261%         |
| FUT11M2 |     99.9051%     |         96.117%          |
| FUT11M4 |     99.8352%     |         95.7566%         |
| FUT12M1 |     99.872%      |         96.0704%         |
| FUT15M5 |     99.8842%     |         96.1064%         |

# Prepping for differential expression analysis

## Importing transcript abundance with tximport

``` r
##Create a file containing information about the samples
ID <- c("FUG10M2", "FUG11F1", "FUG11M2", "FUG11M4", "FUG12M1", "FUG13F1", "FUG13F4", "FUG15M5", "FUG2F2", "FUG3F2",
        "FUL10M2", "FUL11F1", "FUL11M2", "FUL11M4", "FUL12M1", "FUL13F1", "FUL13F4", "FUL15M5", "FUL2F2", "FUL3F2",
        "FUO11F1", "FUO13F1", "FUO13F4", "FUO2F2", "FUO3F2",
        "FUT10M2", "FUT11M2", "FUT11M4", "FUT12M1", "FUT15M5")
Sex <- c("M", "F", "M", "M", "M", "F", "F", "M", "F", "F",
         "M", "F", "M", "M", "M", "F", "F", "M", "F", "F",
         "F", "F", "F", "F", "F",
         "M", "M", "M", "M", "M")
Organ <- rep(c("Gill", "Liver", "Gonad"), times = c(10, 10, 10))
samples <- cbind(ID, Sex, Organ)
samples <- as.data.frame(samples)
write.table(samples, file = "FU_samples.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
```

Generating the transcript abundance was done on the RCC via R command
line:

``` r
#Enter R
conda activate R
R

#Load required package
library(tximport)

#Read in sample file
samples <- read.table("FU_samples.txt", header=TRUE)

#Create list of quant.sf files from salmon
files <- list.files("/home/rccuser/shared/coley_files/fuscus_super_salmon_quant/expression_files", full.names = TRUE)
names(files) <- paste0(samples$ID)
all(file.exists(files))
  ## [1] TRUE
  
#Pull the transcript-gene relationship from the .gtf file generated during the SuperTranscripts step
gtf <- read.table("trinity_genes.gtf", header = FALSE)
tx2gene <- gtf[,c(10, 10)] #Using the gene_id twice since there are no longer isoforms in the SuperTranscript assembly
tx2gene <- unique(tx2gene)
colnames(tx2gene) <- c("gene_id", "transcript_id")

#Building the matrix
txi.salmon.FU <- tximport(files, type = "salmon", tx2gene = tx2gene)
head(txi.salmon.FU$counts)
saveRDS(txi.salmon.FU, "txi.salmon_FU.RDS") #Move the file off the RCC to continue the Differential expression analysis in R
```
