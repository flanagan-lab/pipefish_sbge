Analyzing *Syngnathus fuscus* RNAseq Data from MSU
================
Coley Tosto
2024-07-16

- [Pre-assembly Quality Control and
  Filtering](#pre-assembly-quality-control-and-filtering)
  - [Trimming the raw reads with
    Trimmomatic](#trimming-the-raw-reads-with-trimmomatic)
  - [Using Kraken2 to remove biological
    contamination](#using-kraken2-to-remove-biological-contamination)
  - [Using SortMeRNA to remove rRNA
    contamination](#using-sortmerna-to-remove-rrna-contamination)
  - [Doing a k-mer based correction with
    RCorrector](#doing-a-k-mer-based-correction-with-rcorrector)
- [Checking quality of trimmed and filtered
  reads](#checking-quality-of-trimmed-and-filtered-reads)
- [De novo transcriptome assembly](#de-novo-transcriptome-assembly)
  - [Running Trinity](#running-trinity)
- [Post-assembly Quality Control - On the Trinity
  assembly](#post-assembly-quality-control---on-the-trinity-assembly)
  - [Using BUSCO to assess composition or ‘completeness’ of
    assembly](#using-busco-to-assess-composition-or-completeness-of-assembly)
  - [Trinity transcriptome contig Nx
    Statistics](#trinity-transcriptome-contig-nx-statistics)
  - [Trinity transcriptome contig Ex90N50 Statistic and gene
    count](#trinity-transcriptome-contig-ex90n50-statistic-and-gene-count)
- [Alignment and abundance
  estimations](#alignment-and-abundance-estimations)
  - [Generate the index](#generate-the-index)
  - [Quantify the samples](#quantify-the-samples)
- [Assembly thinning and redundancy reduction -
  SuperTranscripts](#assembly-thinning-and-redundancy-reduction---supertranscripts)
  - [Re-evaluating the quality of the
    assembly](#re-evaluating-the-quality-of-the-assembly)
    - [BUSCO](#busco)
    - [Examining RNA-Seq read representation (via
      Salmon)](#examining-rna-seq-read-representation-via-salmon)
- [Prepping for differential expression
  analysis](#prepping-for-differential-expression-analysis)
  - [Importing transcript abundance with
    tximport](#importing-transcript-abundance-with-tximport)

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

# De novo transcriptome assembly

## Running Trinity

Due to the timeline, and to ensure all three species of *Syngnathus* are
processed in the same way, the assembly was also generated with Trinity
as I did for *Syngnathus floridae*.

``` bash
#!/bin/bash

#Create arguements
samples_file=$1 #File containing sample names and locations
out_dir_name=$2 #Desired name for the output directory

sudo docker run --rm -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq Trinity --seqType fq --max_memory 188G \
    --samples_file $1 \
    --CPU 16 \
    --output `pwd`/$2
```

This script was run in a screen session as
`bash bash_scripts trinity.sh FU_trinity_samples.txt trinity_out_dir_fuscus_June2023`.

# Post-assembly Quality Control - On the Trinity assembly

## Using BUSCO to assess composition or ‘completeness’ of assembly

`BUSCO` was installed in the `busco` conda environment on the RCC.
`BUSCO v 5.2.2` was used.

``` bash
#!/bin/bash

#Create arguments
transcriptome=$1 #Output fasta file from Trinity
lineage=$2 #chosen dataset for assessment
out_dir_name=$3 #Desired name for the output directory

busco -i $1 -l $2 -m tran -o $3 -c 16
```

The lineage chosen for *S fuscus* was `actinopterygii_odb10` since it is
the closest clade provided with BUSCO. This script was run as
`nohup bash bash_scripts/busco_tran.sh trinity_out_dir_fuscus_June2023.Trinity.fasta actinopterygii_odb10 busco_out_fuscus > busco.log 2>&1 &`.

#### Short Summary

    # BUSCO version is: 5.2.2
    # The lineage dataset is: actinopterygii_odb10 (Creation date: 2021-02-19, number of genomes: 26, number of BUSCOs: 3640)
    # Summarized benchmarking in BUSCO notation for file /home/rccuser/shared/coley_files/trinity_out_dir_fuscus_June2023.Trinity.fasta
    # BUSCO was run in mode: transcriptome

            ***** Results: *****

            C:94.7%[S:7.6%,D:87.1%],F:2.4%,M:2.9%,n:3640
            3447    Complete BUSCOs (C)
            277     Complete and single-copy BUSCOs (S)
            3170    Complete and duplicated BUSCOs (D)
            87      Fragmented BUSCOs (F)
            106     Missing BUSCOs (M)
            3640    Total BUSCO groups searched

    Dependencies and versions:
            hmmsearch: 3.1
            metaeuk: 6.a5d39d9

Very similar results to what was found with *S. floridae* prior to any
thinning. Good completeness but high duplication (likely due to the
multiple isoforms that are outputted from Trinity).

#### Plot

<p float="center">

<img src="../imgs/busco_fuscus_trinity.png" style="width:650px;"/>

</p>

## Trinity transcriptome contig Nx Statistics

``` bash
#!/bin/bash

##Create the arguments
trinity_fasta_file=$1 #This is the output .fasta file from your Trinity run
trinity_stats_output=$2 #Desired name for your output results

sudo docker run -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq /usr/local/bin/util/TrinityStats.pl \
        $trinity_fasta_file \
        >> $trinity_stats_output
```

#### Results

    ################################
    ## Counts of transcripts, etc.
    ################################
    Total trinity 'genes':  323048
    Total trinity transcripts:  534363
    Percent GC: 46.61

    ########################################
    Stats based on ALL transcript contigs:
    ########################################

        Contig N10: 7775
        Contig N20: 5458
        Contig N30: 4161
        Contig N40: 3192
        Contig N50: 2414

        Median contig length: 456
        Average contig: 1088.98
        Total assembled bases: 581911276


    #####################################################
    ## Stats based on ONLY LONGEST ISOFORM per 'GENE':
    #####################################################

        Contig N10: 5947
        Contig N20: 3778
        Contig N30: 2391
        Contig N40: 1502
        Contig N50: 979

        Median contig length: 346
        Average contig: 656.63
        Total assembled bases: 212122880

Based on the longest isoform we can see that 50% of the assembled bases
are found in transcript contigs that are at least 979 bases in length.
This is lower than the 1307 N50 from *Syngnathus floridae*.

## Trinity transcriptome contig Ex90N50 Statistic and gene count

#### Create an Abundance Matrix

    #Get a list of the salmon quant.sf files so you don't have to list them individually
    find . -maxdepth 2 -name "*quant.sf" | tee FUsalmon.quant_files.txt

``` bash
#!/bin/bash

trinity_gene_map=$1
salmon_quant_files=$2
out_file_name=$3

sudo docker run -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq /bin/sh -c "cd /home/rccuser/shared/coley_files/ && /usr/local/bin/util/abundance_estimates_to_matrix.pl \
        --est_method salmon \
        --gene_trans_map $trinity_gene_map \
        --quant_files $salmon_quant_files \
        --out_prefix $out_file_name"
```

This script was run as
`nohup bash bash_scripts/trinity_abund_est.sh`pwd`/trinity_out_dir_fuscus_June2023.Trinity.fasta.gene_trans_map`pwd`/fuscus_salmon_quant/FUsalmon.quant_files.txt fuscus_salmon > trin_abund.log 2>&1 &`.
The shell wrapper (`/bin/sh -c`) was added to the script to ensure the
ouput files were sent to the right place.

#### Calculate the Ex90N50 stat

``` bash
#!/bin/bash

isoform_expression_matrix=$1
trinity_fasta_output=$2
ExN50_output=$3

sudo docker run -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq /usr/local/bin/util/misc/contig_ExN50_statistic.pl \
        $isoform_expression_matrix \
        $trinity_fasta_output | tee $ExN50_output
```

This script was run as
`bash bash_scripts/exn50.sh`pwd`/fuscus_salmon.isoform.TMM.EXPR.matrix`pwd`/trinity_out_dir_fuscus_June2023.Trinity.fasta fuscus_exN50.stats`

| Ex  | ExN50 | num_genes |
|:----|:------|:----------|
| 2   | 281   | 1         |
| 4   | 2370  | 2         |
| 6   | 2370  | 3         |
| 7   | 1570  | 4         |
| …   | …     | …         |
| 90  | 2981  | 25899     |
| 99  | 2431  | 492686    |
| 100 | 2414  | 534363    |

We can then plot the Ex value (first column) against the ExN50 values:

`sudo docker run -v`pwd`:`pwd`trinityrnaseq/trinityrnaseq /bin/sh -c "cd /home/rccuser/shared/coley_files/ && /usr/local/bin/util/misc/plot_ExN50_statistic.Rscript /home/rccuser/shared/coley_files/fuscus_exN50.stats"`

the Ex90N50 for *S. fuscus* (2981) is slightly higher but very similar
to the Ex90N50 for *S. floridae* (2918).

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
