Investigating selection on genes expressed in S. fuscus
================
Sarah P. Flanagan
2025-10-16



- [Translating the transcriptome](#translating-the-transcriptome)
- [Identifying orthologs using genome
  sequences](#identifying-orthologs-using-genome-sequences)
  - [single-copy orthologs](#single-copy-orthologs)
- [Aligning single-copy orthologs for selection
  analysis](#aligning-single-copy-orthologs-for-selection-analysis)
- [Selection analysis](#selection-analysis)
  - [Are selected loci pleiotropically
    constrained?](#are-selected-loci-pleiotropically-constrained)
  - [what are the selected loci?](#what-are-the-selected-loci)

``` r
parse_gff<-function(filename){
  gff <- read.delim(filename,
                  comment.char = "#", 
                  header = FALSE)
  colnames(gff) <- c("chromosome","method","type","start","stop","x","strand","y","description")
  return(gff)
}
```

``` r
sex_bias_colors <- c("FB" = "#7fc97f", "MB" = "#beaed4", "UB" = "darkgray")
```

The goal is to better understand the sex specific genes – what’s the
evolutionary rate?

## Translating the transcriptome

I will translate the transcriptome to be able to use it in the
orthofinder analysis with genome sequences

``` bash
#transdecoder is installed on trinity-big (not on trinity-em)
# from shared/coley_files/all_sex_specific_genes/
for fa in *fasta; do ../TransDecoder-TransDecoder-v5.7.1/TransDecoder.LongOrfs -t $fa -m 50 --complete_orfs_only; done
```

``` bash
# to run on the RCH
module load prefix/2025
module load prefix/latest
module load TransDecoder/5.7.1-GCC-12.3.0


../TransDecoder-TransDecoder-v5.7.1/TransDecoder.LongOrfs -t /home/spf50/trinity_out_dir_fuscus_June2023.Trinity.fasta -m 50 --complete_orfs_only
```

Some of them have multiple sequences and they’re all in sub-folders now,
which is kind of annoying

The manual says it works better on a whole transcriptome than on single
genes, so I’ll run it on the transcriptome file, I guess.

``` bash
# now in shared/coley_files
../TransDecoder-TransDecoder-v5.7.1/TransDecoder.LongOrfs -t  trinity_supertran_fuscus.fasta -m 50 --complete_orfs_only
```

I had to lower the length threshold to include at one of the
sex-specific genes (931), if not more (I only checked that one).

## Identifying orthologs using genome sequences

So I will first identify orthologs. I’ve downloaded protein sequences
from GenBank:

Based on the OrthoFinder tutorial on best practices, I’m using the
proteomes where available. I have named these S_scovelli.protein.zip
etc. Unfortunately Syngnathus floridae did not have a protein file
available.

I’ve downloaded two Hippocampus protein files because in case we are
really interested in changes at the Syngnathus fuscus branch (as opposed
to comparing across clades) – according to the [OrthoFinder
tutorial](https://davidemms.github.io/orthofinder_tutorials/orthofinder-best-practices.html):

    If you’re interested in what happened on a particular branch of the species tree, then you should likewise ensure good species sampling—ideally at least two species below the branch, at least two species on the closest branch above and two or more species in the outgroup.

I unzipped the ncbi downloads and pulled all the fasta files out into
one folder. I’m skipping S. floridae for now.

I used the orthofinder tool to identify primary transcripts (not sure it
did anything, since these are not ensembl files)

``` bash
for f in *faa ; do python ~/Programs/OrthoFinder/tools/primary_transcript.py $f ; done
```

This resulted in a folder called primary_transcripts containing the
protein fasta files from the genomes. I also need to put the outputs
from transdecoder into the folder with these primary transcripts

Then, I need to run OrthoFinder:

``` bash
# on C001KR -- can't use dropbox because the "(" breaks things
# without the transdecoder transcriptome
/home/sarah/Programs/OrthoFinder/orthofinder -f primary_transcripts/ #Feb 11 option

# copy all *protein.faa files into orthofinder/ along with longest_orfs.pep
./orthofinder -f /mnt/BigData/genomes/orthofinder/ #Aug 15 option

# with the non-supertranscriptome option, with trinity_out_dir_fuscus_June2023.Trinity.fasta.transdecoder.pep (removed other fas, faa, pep, fa, and fasta files from the dir)
/home/sarah/Programs/OrthoFinder/orthofinder -f /mnt/BigData/genomes/orthofinder/ # Sep 12 and 16
```

(genomes only) The overall statistics report that of 726348 genes,
319862 were placed in orthogroups (44.0370181 %), which is fairly good
considering the messiness of the transcriptome (there were ~112k genes
without the transcriptome). The total number of single-copy orthogroups
is 2278 (down from 10,504 without the transcriptome).

The species tree looks to correctly reflect the known phylogenetic tree
structure (albeit with a long branch to *S. fuscus*):

``` r
library(ape)
ortho_tree <- ape::read.tree(text="((H_comes.protein:0.0300949,H_zosterae.protein:0.0355681)0.947398:0.0607965,((S_scovelli.protein:0.0112549,longest_orfs:0.274165)0.438367:0.0217156,(S_typhle.protein:0.0191163,S_acus.protein:0.0120793)0.490145:0.0126776)0.947398:0.0607965);")
plot(ortho_tree)
```

![](sfu_sex_specific_evolution_files/figure-gfm/orthoTree-1.png)<!-- -->

### single-copy orthologs

``` bash
cd OrthoFinder/Results_Aug15/
ls Single_Copy_Orthologue_Sequences/*fa > Single_Copy_Orthologue_Sequences.txt

cd OrthoFinder/Results_Sep12/
ls Single_Copy_Orthologue_Sequences/*fa > Single_Copy_Orthologue_Sequences.txt

cd OrthoFinder/Results_Sep16/
ls Single_Copy_Orthologue_Sequences/*fa > Single_Copy_Orthologue_Sequences.txt
```

``` r
# these should be intermediate_data/ but are also on Dropbox (not in 08_)

single_copy_groups <- read.table("08_sex_specific_evolution/OrthoFinder/Results_Sep16/Orthogroups/Orthogroups.tsv", sep='\t', header=TRUE)

single_copy_list <- read.table("08_sex_specific_evolution/OrthoFinder/Results_Sep16/Single_Copy_Orthologue_Sequences.txt")
single_copy_list <- gsub("Single_Copy_Orthologue_Sequences\\/(.*)\\.fa","\\1",single_copy_list$V1)
single_copy_groups <- single_copy_groups[single_copy_groups$Orthogroup %in% single_copy_list,]
```

Are any of the single copy orthologs the sex-biased or sex-specific
genes?

``` r
# get the info about which genes were sex specific in which tissues
# now includes ones that did not align
# in FlanaganLab/intermediate_data/SFU_DE_pipeline_analysis_2023
# on linux: setwd("/run/user/1000/gvfs/smb-share:server=file,share=research/FlanaganLab/intermediate_data/SFU_DE_pipeline_analysis_2023/")
files<-c(list.files(path="06_diff_expr_analysis/01_high_exp_trinity_gene_names",
           pattern="specific",full.names = TRUE),
         list.files(path="06_diff_expr_analysis/01_high_exp_trinity_gene_names",
           pattern="biased",full.names = TRUE) )
         
degs<-do.call(rbind,lapply(files,function(file){
  dat<-read.delim(file, header=FALSE)
  colnames(dat)<-"Trinity_gene"
  dat$sex <- gsub("^.*FU_(\\w{3})(\\w+).*$","\\1",file)
  dat$organ <- gsub("^.*FU_(\\w{3})(\\w+)_.*_.*$","\\2",file)
  dat$DEG_cat<-gsub("^.*FU_(\\w{3})(\\w+)_(\\w+)_.*$","\\3",file)
  return(dat)
}))

# prep the single copy groups for merging
#single_copy_groups$S_fuscus_gene<-gsub("(g\\d+)\\.p\\d+","\\1",single_copy_groups$longest_orfs) for Aug one
# single_copy_groups$S_fuscus_gene<-gsub("(g\\d+)\\.p\\d+","\\1",single_copy_groups$trinity_out_dir_fuscus_June2023.Trinity.fasta.transdecoder) # for Sep 12 one - saved in single_copy_DEG_notsuper.csv
single_copy_groups$S_fuscus_gene<-gsub("(g\\d+)\\.p\\d+","\\1",single_copy_groups$trinity_supertran_fuscus.fasta.transdecoder) # for Sep 16 one



deg_sc_orthologs<-merge(single_copy_groups,
                        degs,
                                 by.x="S_fuscus_gene",
                                 by.y = "Trinity_gene",
                        all.x=TRUE) # keep all single copy orthogroups but not all DEGs


# fill in gaps to avoid NAs
deg_sc_orthologs$DEG_cat[is.na(deg_sc_orthologs$DEG_cat)]<-"unbiased"
deg_sc_orthologs$sex[is.na(deg_sc_orthologs$sex)]<-"unbiased"
deg_sc_orthologs$organ[is.na(deg_sc_orthologs$organ)]<-"unbiased"

# remove and condense duplicates
# Some orthogroups have multiple combinations of sex/organ/DEG categories. Most of these are the case of being labelled as both biased and specific in the same sex and organ, which would be better to be condensed into just 'specific'. 
dup_ogs<-deg_sc_orthologs$Orthogroup[duplicated(deg_sc_orthologs$Orthogroup)]
deduped<-do.call(rbind,lapply(unique(deg_sc_orthologs$Orthogroup), function(og){
  this_og<-deg_sc_orthologs[deg_sc_orthologs$Orthogroup==og,]
  # check if it's duplicated
  if(nrow(this_og)>1){
    # check that all other elements are duplicated
    if(nrow(unique(this_og[,1:8]))==1){
      # if it's the same organ and sex then we probably just need to convert to "biased"
      if(length(unique(this_og$sex)) ==1 & length(unique(this_og$organ))==1){
        # if it is the case of biased and specific being the only difference
        if(sum(this_og$DEG_cat %in% c("biased","specific")) == nrow(this_og)){
          this_og$DEG_cat[this_og$DEG_cat=="biased"]<-"specific"
          if(nrow(unique(this_og))==1){ # sanity check that this fixed it
            return(unique(this_og))
          } else{
            browser()
          }
        }
      } else{
        # we'll need to combine info somehow
        this_og$sex<-paste0(this_og$sex, collapse=";")
        this_og$organ<-paste0(this_og$organ, collapse=";")
        if(nrow(unique(this_og))==1){ # sanity check that this fixed it
            return(unique(this_og))
          } else{
            browser()
          }
      }
    } else{
      browser()
    }
  } else{
    return(this_og) # just return it if it's not duplicated
  }
}))

# save to file
write.csv(deduped,"single_copy_DEG.csv",quote=FALSE, row.names = FALSE)
```

``` r
deg_sc_orthologs<-read.csv("single_copy_DEG.csv")
# factorise for easier investigation
deg_sc_orthologs$DEG_cat<-as.factor(deg_sc_orthologs$DEG_cat)
deg_sc_orthologs$sex<-as.factor(deg_sc_orthologs$sex)
deg_sc_orthologs$organ<-as.factor(deg_sc_orthologs$organ)
```

Now we can summarise the information

``` r
summary(deg_sc_orthologs$DEG_cat)
table(deg_sc_orthologs$DEG_cat, deg_sc_orthologs$sex)
```

One interesting (?) development is that some of the trinity genes that
seem to have the same name have different proteins in the single copy
groups. For example:

``` r
deg_sc_orthologs[deg_sc_orthologs$S_fuscus_gene=="TRINITY_DN0_c2_g1",]
```

## Aligning single-copy orthologs for selection analysis

To run Hyphy, I will need the nucleotide alignments. The orthogroups are
amino acids, so I need to get the original sequences from the files
available from NCBI. Once that is done, I will save all of the protein
names to files to then process them.

``` r
# dropbox
write.table(single_copy_groups$S_scovelli.protein, "generating_alignments/S_scovelli_protein_names.txt",quote=FALSE,row.names = FALSE, col.names = FALSE, eol='\n')
write.table(single_copy_groups$H_comes.protein, "generating_alignments/H_comes_protein_names.txt",quote=FALSE,row.names = FALSE, col.names = FALSE, eol='\n')
write.table(single_copy_groups$H_zosterae.protein, "generating_alignments/H_zosterae_protein_names.txt",quote=FALSE,row.names = FALSE, col.names = FALSE, eol='\n')
write.table(single_copy_groups$S_acus.protein, "generating_alignments/S_acus_protein_names.txt",quote=FALSE,row.names = FALSE, col.names = FALSE, eol='\n')
write.table(single_copy_groups$S_typhle.protein, "generating_alignments/S_typhle_protein_names.txt",quote=FALSE,row.names = FALSE, col.names = FALSE, eol='\n')
write.table(single_copy_groups$trinity_supertran_fuscus.fasta.transdecoder,
            "generating_alignments/S_fuscus_protein_names.txt", quote=FALSE, row.names = FALSE, col.names = FALSE, eol='\n')
```

Then I use this script to extract the headers and remove the “\>” symbol
at the start of the line:

``` bash
#!/bin/bash
for fa in *genomic.fna; 
do
  spp=$(echo $fa | perl -pe 's/([A-Z]\_[a-z]+)\_cds_from_genomic.fna/\1/g')
  grep -f ${spp}_protein_names.txt $fa | perl -pe 's/^>//g' >> ${spp}_fna_headers.txt
done
```

Now I can extract the genes from each genome cds using my custom
`extract_genes_fasta` program. For it to work, I had to edit the program
to allow the names to be simplified (remove anything after the first
space in the header). *Note*: this will not find matches if the “\>” is
in the list of headers

``` bash
# in generating_alignments
for fa in *genomic.fna; do spp=$(echo $fa | perl -pe 's/([A-Z]\_[a-z]+)\_cds_from_genomic.fna/\1/g'); ../extract_genes_fasta -f $fa -l ${spp}_fna_headers.txt -o ${spp}_ -s; done
```

I will also need to do this with the $S. fuscus$ TransDecoder output

``` bash
# extract the header info
grep -f S_fuscus_protein_names.txt ../trinity_supertran_fuscus.fasta.transdecoder_dir/longest_orfs.cds | perl -pe 's/^>//g' > S_fuscus_headers.txt
# extract the CDS files from the fuscus transcriptome
../extract_genes_fasta -f ../trinity_supertran_fuscus.fasta.transdecoder_dir/longest_orfs.cds -l S_fuscus_headers.txt -o S_fuscus_ -s
```

Now, I need a set of files with all the fasta names per orthogroup, so
that I can then use those to grep across all samples.

``` r
apply(single_copy_groups, 1, function(this_row){
  write(c(this_row[2:6], paste0(this_row[7],".fasta")), # exclude the Orthogroup (1) and S_fuscus_gene (8) columns
        file=paste0("generating_alignments/",this_row["Orthogroup"],".txt"))
})
```

I’ll use a cheeky grep + cat + sed to combine all six fasta nucleotide
sequences per orthogroup.

``` bash
# in generating_alignments
ls *fasta > all_fastas.txt
# remove windows endlines from the og lists (if OGs were created on windows R)
sed -i 's/\r//g' OG*.txt

# now use grep plus cat plus sed to merge files into one
for og in OG*.txt; do og_id=$(basename $og .txt); grep -f $og all_fastas.txt | xargs cat | sed 's/\w>/\n>/g' > $og_id.fasta; done 
```

Do a quick check to make sure they all have six sequences as expected

``` bash
for file in *fasta; do count=$(grep ">" $file | wc -l); if [ $count -gt '6' ]; then echo $file; fi; done
```

I’m going to move these fasta files into a new folder for better file
organisation

``` bash
mkdir ../OrthoFinder/Results_Sep16/Single_Copy_Orthologue_Seqs_NT
mv OG*fasta ../OrthoFinder/Results_Sep16/Single_Copy_Orthologue_Seqs_NT/
```

I installed prank v.250331 from source (from github) and copied the
executable to /usr/bin on C001KR.

Now I will try aligning the single copy orthologues using PRANK. I will
use the codon-aware option, but I am letting it infer the tree itself
(mainly because the species tree created by orthofinder and the fasta
files have different names). I have to run this in a directory that does
not have weird characters (so I’m running it in
/mnt/BigData/genomes/orthofinder) and I need to use the -shortnames flag
to avoid an allocation error (`std::bad_alloc()`).

``` bash
#!/bin/bash

# paths relative to location of script
FASTA_DIR="/mnt/BigData/genomes/orthofinder/OrthoFinder/Results_Sep16/Single_Copy_Orthologue_Seqs_NT"
OUT_DIR="/mnt/BigData/genomes/orthofinder/OrthoFinder/Results_Sep16/Single_Copy_Orthologue_Alignments"

# move to location of script
cd "${0%/*}" 

# check if the output directory does not exist
if [ ! -d "$OUT_DIR" ]; then
  mkdir $OUT_DIR # make the dir if it doesn't already exist
fi

for fa in $FASTA_DIR/*fasta
do
    ID=$(basename $fa .fasta) # remove the extension
    prank -d=$fa -o=$OUT_DIR/$ID -translate -shortnames
    rm -r tmpdirprank* #remove the temporary directory that's created
done
```

I will likely need to check some of these alignments.

## Selection analysis

I am installing hyphy v 2.5.78 on C001KR from source. My goal is to run
it on all the single copy orthologs and then we can compare estimates of
selection on the genes in different categories.

Before running the analysis program, I need to remove stop codons (using
hyphy’s tool) and change the names to match the species tree (using my
own bash/regex). I will use the `single_copy_DEG.csv` file to get the
names of the genes for each orthogroup.

``` bash
#!/bin/bash

ALIGNMENTS="/mnt/BigData/genomes/orthofinder/OrthoFinder/Results_Sep16/Single_Copy_Orthologue_Alignments"
ref=single_copy_DEG.csv
species=("H_comes" "H_zosterae" "S_acus" "S_scovelli" "S_typhle" "S_fuscus")

for fa in $ALIGNMENTS/*best.nuc.fas
do
  og=$(basename $fa .best.nuc.fas)
  
  # remove stop codons
  hyphy ~/Programs/hyphy-2.5.78/res/TemplateBatchFiles/CleanStopCodons.bf Universal $fa "Disallow stops" $ALIGNMENTS/${og}_aln.fa

  # get the gene names for each species 
  genes=()
  genes+=( $(grep $og $ref | perl -pe 's/^.*\,OG\d+\,(.*)\,(.*)\,(.*)\,(.*)\,(.*)\,(.*)\,(.*)\,(.*),(.*)/\1/g') )
  genes+=( $(grep $og $ref | perl -pe 's/^.*\,OG\d+\,(.*)\,(.*)\,(.*)\,(.*)\,(.*)\,(.*)\,(.*)\,(.*),(.*)/\2/g') )
  genes+=( $(grep $og $ref | perl -pe 's/^.*\,OG\d+\,(.*)\,(.*)\,(.*)\,(.*)\,(.*)\,(.*)\,(.*)\,(.*),(.*)/\3/g') )
  genes+=( $(grep $og $ref | perl -pe 's/^.*\,OG\d+\,(.*)\,(.*)\,(.*)\,(.*)\,(.*)\,(.*)\,(.*)\,(.*),(.*)/\4/g') )
  genes+=( $(grep $og $ref | perl -pe 's/^.*\,OG\d+\,(.*)\,(.*)\,(.*)\,(.*)\,(.*)\,(.*)\,(.*)\,(.*),(.*)/\5/g') )
  genes+=( $(grep $og $ref | perl -pe 's/^.*\,OG\d+\,(.*)\,(.*)\,(.*)\,(.*)\,(.*)\,(.*)\,(.*)\,(.*),(.*)/\6/g') )


  # create a file with the regex to change the file names
  echo "#!/bin/bash" > fix_names_tmp.sh
  # for each species name, match it to ghe correct gene
  for i in "${!genes[@]}"; do
    spp=$( echo "${species[$i]}" )
    g=$( echo  "${genes[$i]}" )
    echo "sed -ie 's/^>.*$g.*$/>$spp/g' $ALIGNMENTS/${og}_aln.fa" >> fix_names_tmp.sh
  done
  # run the file to change the names
  ./fix_names_tmp.sh

done




    
```

I think the tree will also need to be annotated so that the S. fuscus
branch has the the {FG} label.

``` bash
# move tree to syngnathus_fuscus/ on Dropbox
cp /mnt/BigData/genomes/orthofinder/OrthoFinder/Results_Sep16/Species_Tree/SpeciesTree_rooted.txt ./
# remove '.protein' after the downloaded species' info
sed -ie 's/\.protein//g' SpeciesTree_rooted.txt
# change and annotate the S_fuscus label
sed -ie 's/trinity_supertran_fuscus.fasta.transdecoder/S_fuscus\{FG\}/g' SpeciesTree_rooted.txt 
```

``` bash
#!/bin/bash
ALIGNMENTS="/mnt/BigData/genomes/orthofinder/OrthoFinder/Results_Sep16/Single_Copy_Orthologue_Alignments"
OUTDIR="/mnt/BigData/genomes/orthofinder/OrthoFinder/Results_Sep16/Single_Copy_Orthologue_absrel2"
for aln in $ALIGNMENTS/*_aln.fa
do
  og=$(basename $aln _aln.fa)
  hyphy absrel --alignment $aln --tree SpeciesTree_rooted.txt  --code Universal --branches FG --output ${OUTDIR}/${og}.absrel
done
```

``` r
# in \\file\Research\FlanaganLab\intermediate_data\SFU_DE_pipeline_analysis_2023\08_sex_specific_evolution\OrthoFinder\Results_Sep16
library(jsonlite)
absrel_files<-list.files(pattern="absrel",
                         path = "Single_Copy_Orthologue_absrel2",
                         full.names = TRUE) 
# remove empty files (investigate why they're empty)
empty <- absrel_files[file.size(absrel_files) == 0L]
absrel_files<-absrel_files[! absrel_files %in% empty]

test_results<-do.call(rbind,lapply(absrel_files, function(file){
  absrel_output<-read_json(file)

  test_result <- absrel_output$`test results`$`positive test results`
  attributes<-absrel_output$`branch attributes`$`0`
  
  fuscus<-attributes[grep("fuscus", names(attributes))]
  fuscus_omega <- fuscus[[1]]$`Baseline MG94xREV omega ratio`
  og<-gsub("Single_Copy_Orthologue_absrel2/(OG\\d+)\\.absrel","\\1",file)
  return(c(og,test_result, fuscus_omega))

}))
rownames(test_results)<-absrel_files
colnames(test_results)<-c("orthogroup","selection", "omega")
```

I would like to know what type of expression these orthogroups had, so I
will merge them ..

``` r
orthogroups_selection_all<-merge(test_results, deg_sc_orthologs, by.x="orthogroup",by.y="Orthogroup")
write.csv(orthogroups_selection_all,
          "orthogroups_absrel_results.csv",
          quote = FALSE,
          row.names = FALSE)
```

``` r
orthogroups_selection_all<-read.csv("orthogroups_absrel_results.csv")

table(orthogroups_selection_all$sex,
      orthogroups_selection_all$selection)
```

    ##           
    ##               0    1
    ##   fem       238  117
    ##   fem;mal     2    0
    ##   mal       222   94
    ##   mal;mal     1    0
    ##   unbiased 1122  567

``` r
table(paste0(orthogroups_selection_all$sex,"_" ,orthogroups_selection_all$organ),
      orthogroups_selection_all$selection, orthogroups_selection_all$DEG_cat)
```

I would like to also match this up with the expression data for these
genes.

``` r
logFC_tau<-read.csv("./logFC_long_taubias_SS.csv") 

orthogroups_selection_all<-do.call(rbind,apply(orthogroups_selection_all,
                            1,
                            function(gene, expr_info){
                            
                             this_gene<-data.frame(rbind(gene))
                             #print(gene["S_fuscus_gene"])
                             this_info<-expr_info[expr_info$geneID==gene["S_fuscus_gene"],]
                              
                              #this_info<-this_info[!is.na(this_info$logFC),]
                             if(nrow(this_info)>3 | nrow(this_info) < 1){
                               
                               this_gene$tau <- NA
                               this_gene$Gill_logFC<-NA
                               this_gene$Gonad_logFC<-NA
                               this_gene$Liver_logFC<-NA
                               this_gene$sex<-NA
                               this_gene$organ<-NA
                               this_gene$DEG_cat<-NA
                             
                             } else {
                               
                               this_gene$tau <- mean(this_info$tau, na.rm=TRUE)
                               this_gene$Gill_logFC<-ifelse("Gill" %in% this_info$tissue,this_info$logFC[this_info$tissue=="Gill"], NA)
                               this_gene$Gonad_logFC<-ifelse("Gonad" %in% this_info$tissue,this_info$logFC[this_info$tissue=="Gonad"], NA)
                               this_gene$Liver_logFC<-ifelse("Liver" %in% this_info$tissue,this_info$logFC[this_info$tissue=="Liver"], NA)
                             
                              
                             }
                             return(this_gene)
                            }, expr_info=logFC_tau))

# make sure it's numeric!
orthogroups_selection_all$omega<-as.numeric(orthogroups_selection_all$omega)
```

A total of
`nrow(orthogroups_selection_all[is.na(orthogroups_selection_all$sex),])`
Trinity genes assigned to orthogroups had been removed from the
differential expression analysis (either due to low overall expression
or having been identified as outliers), so they are removed from any
analyses linking selection to expression levels.

``` r
table(orthogroups_selection_all$sex,
           orthogroups_selection_all$selection)
```

    ##           
    ##               0    1
    ##   fem       238  117
    ##   fem;mal     2    0
    ##   mal       222   92
    ##   mal;mal     1    0
    ##   unbiased 1036  550

``` r
# for now, remove those with nonsensical omegas (inspect/revisit later?)
orthogroups_selection<-orthogroups_selection_all[which(orthogroups_selection_all$omega<100),]
# remove the 3 observations that are a switch
orthogroups_selection$sex[orthogroups_selection$sex=="mal;mal"]<-"mal"
orthogroups_selection<-orthogroups_selection[orthogroups_selection$sex %in% c("fem","mal","unbiased"),]

# calculate ns
ns<-table(orthogroups_selection$sex,
      orthogroups_selection$selection)
```

We removed a total of 147 that had non-sensical omega estimates.

``` r
kable(ns,"latex",col.names = c("no selection","selection"))
```

We can use a chi-squared test to see if there is a higher proportion of
female-biased genes than unbiased?

``` r
(selchi<-chisq.test(ns))
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  ns
    ## X-squared = 2.3101, df = 2, p-value = 0.315

``` r
#selchi$expected
```

are the omega values different?

The residuals look ok despite the skewed distributions so I’ll just
continue without transformation

``` r
mod<-lm(rank(orthogroups_selection$omega)~orthogroups_selection$sex*orthogroups_selection$selection)
par(mfrow=c(2,2))
plot(mod)
```

``` r
anova(mod)
```

    ## Analysis of Variance Table
    ## 
    ## Response: rank(orthogroups_selection$omega)
    ##                                                             Df    Sum Sq
    ## orthogroups_selection$sex                                    2    339321
    ## orthogroups_selection$selection                              1  80447832
    ## orthogroups_selection$sex:orthogroups_selection$selection    2     46630
    ## Residuals                                                 2112 698756864
    ##                                                            Mean Sq  F value
    ## orthogroups_selection$sex                                   169661   0.5128
    ## orthogroups_selection$selection                           80447832 243.1544
    ## orthogroups_selection$sex:orthogroups_selection$selection    23315   0.0705
    ## Residuals                                                   330851         
    ##                                                           Pr(>F)    
    ## orthogroups_selection$sex                                 0.5989    
    ## orthogroups_selection$selection                           <2e-16 ***
    ## orthogroups_selection$sex:orthogroups_selection$selection 0.9320    
    ## Residuals                                                           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
mod_aov<-aov(rank(orthogroups_selection$omega)~as.factor(orthogroups_selection$sex)*as.factor(orthogroups_selection$selection))
tuk<-TukeyHSD(mod_aov)
interactions<-tuk[[3]]
round(interactions,4)
```

    ##                           diff       lwr      upr  p adj
    ## mal:0-fem:0             5.4999 -156.2132 167.2129 1.0000
    ## unbiased:0-fem:0      -44.9053 -168.1399  78.3292 0.9046
    ## fem:1-fem:0           389.6127  201.5871 577.6382 0.0000
    ## mal:1-fem:0           399.2729  195.2991 603.2467 0.0000
    ## unbiased:1-fem:0      368.1379  236.7633 499.5125 0.0000
    ## unbiased:0-mal:0      -50.4052 -179.4134  78.6030 0.8756
    ## fem:1-mal:0           384.1128  192.2536 575.9721 0.0000
    ## mal:1-mal:0           393.7731  186.2600 601.2861 0.0000
    ## unbiased:1-mal:0      362.6380  225.8329 499.4432 0.0000
    ## fem:1-unbiased:0      434.5180  273.7523 595.2837 0.0000
    ## mal:1-unbiased:0      444.1783  265.0211 623.3355 0.0000
    ## unbiased:1-unbiased:0 413.0432  325.0018 501.0847 0.0000
    ## mal:1-fem:1             9.6603 -218.9533 238.2738 1.0000
    ## unbiased:1-fem:1      -21.4748 -188.5620 145.6125 0.9991
    ## unbiased:1-mal:1      -31.1350 -215.9859 153.7158 0.9968

``` r
anova(mod)
```

    ## Analysis of Variance Table
    ## 
    ## Response: rank(orthogroups_selection$omega)
    ##                                                             Df    Sum Sq
    ## orthogroups_selection$sex                                    2    339321
    ## orthogroups_selection$selection                              1  80447832
    ## orthogroups_selection$sex:orthogroups_selection$selection    2     46630
    ## Residuals                                                 2112 698756864
    ##                                                            Mean Sq  F value
    ## orthogroups_selection$sex                                   169661   0.5128
    ## orthogroups_selection$selection                           80447832 243.1544
    ## orthogroups_selection$sex:orthogroups_selection$selection    23315   0.0705
    ## Residuals                                                   330851         
    ##                                                           Pr(>F)    
    ## orthogroups_selection$sex                                 0.5989    
    ## orthogroups_selection$selection                           <2e-16 ***
    ## orthogroups_selection$sex:orthogroups_selection$selection 0.9320    
    ## Residuals                                                           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
labs<-multcompLetters4(mod_aov,tuk)
```

The key takeaways here are the genes with selection have higher omega
values (unsurprisingly), but that of those that are under selection in
S. fuscus, female-biased genes \> unbiased and female-biased genes \>
male-biased genes (though unbiased==male-biased).

``` r
orthogroups_selection$plot_group<-paste0(orthogroups_selection$selection, "_", orthogroups_selection$sex)
par(mfrow=c(1,1), mar=c(4,3,2,1))
vioplot(orthogroups_selection$omega ~ orthogroups_selection$plot_group,
        xlab="",
        ylab="",
        main="",
        ylim=c(0,10),
        col=rep(sex_bias_colors, 2),
        cex.axis=1,
        cex=1.5,
        frame.plot=FALSE,
        xaxt='n',
        pchMed=21,
        colMed2=rep(sex_bias_colors,2),
        colMed="black",
        las=1)
axis(1,pos=0,at=0:6,labels=c("",
                             paste0("FB\nn=", ns["fem",1]),
                             paste0("MB\nn=", ns["mal",1]),
                             paste0("UB\nn=", ns["unbiased",1]),
                             paste0("FB\nn=", ns["fem",2]),
                             paste0("MB\nn=", ns["mal",2]),
                             paste0("UB\nn=", ns["unbiased",2])),
     tick=FALSE)
# add letters
text(x=1:6,y=rep(9.25,3),
     c(labs[[3]]$Letters[["fem:0"]],
       labs[[3]]$Letters[["mal:0"]],
       labs[[3]]$Letters[["unbiased:0"]],
       labs[[3]]$Letters[["fem:1"]],
       labs[[3]]$Letters[["mal:1"]],
       labs[[3]]$Letters[["unbiased:1"]]), xpd=TRUE)

# add group labels
lines(x=c(1,3),y=c(10,10))
text(x=2,y=10.5,"no selection", xpd=TRUE)
lines(x=c(4,6),y=c(10,10))
text(x=5,y=10.5,"selection", xpd=TRUE)

# add axis labels
mtext("direction of bias",1,outer=TRUE,line=-2)
mtext(expression(italic("S. fuscus")~"branch"~italic("\u03C9")), 2, outer=TRUE, line=-1)
```

![](sfu_sex_specific_evolution_files/figure-gfm/fig5A-1.pdf)<!-- -->

### Are selected loci pleiotropically constrained?

``` r
mod2<-lm(orthogroups_selection$tau[orthogroups_selection$selection==1]~orthogroups_selection$plot_group[orthogroups_selection$selection==1])
anova(mod2)
```

    ## Analysis of Variance Table
    ## 
    ## Response: orthogroups_selection$tau[orthogroups_selection$selection == 
    ## Response:     1]
    ##                                                                         Df
    ## orthogroups_selection$plot_group[orthogroups_selection$selection == 1]   2
    ## Residuals                                                              754
    ##                                                                         Sum Sq
    ## orthogroups_selection$plot_group[orthogroups_selection$selection == 1]  3.4268
    ## Residuals                                                              26.3890
    ##                                                                        Mean Sq
    ## orthogroups_selection$plot_group[orthogroups_selection$selection == 1]  1.7134
    ## Residuals                                                               0.0350
    ##                                                                        F value
    ## orthogroups_selection$plot_group[orthogroups_selection$selection == 1]  48.957
    ## Residuals                                                                     
    ##                                                                           Pr(>F)
    ## orthogroups_selection$plot_group[orthogroups_selection$selection == 1] < 2.2e-16
    ## Residuals                                                                       
    ##                                                                           
    ## orthogroups_selection$plot_group[orthogroups_selection$selection == 1] ***
    ## Residuals                                                                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
TukeyHSD(aov(mod2))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = mod2)
    ## 
    ## $`orthogroups_selection$plot_group[orthogroups_selection$selection == 1]`
    ##                         diff        lwr         upr     p adj
    ## 1_mal-1_fem      -0.04779568 -0.1090130  0.01342163 0.1594029
    ## 1_unbiased-1_fem -0.16893158 -0.2136736 -0.12418956 0.0000000
    ## 1_unbiased-1_mal -0.12113589 -0.1706346 -0.07163720 0.0000000

``` r
mod2b<-lm(orthogroups_selection$tau~orthogroups_selection$sex*orthogroups_selection$selection)
anova(mod2b)
```

    ## Analysis of Variance Table
    ## 
    ## Response: orthogroups_selection$tau
    ##                                                             Df Sum Sq Mean Sq
    ## orthogroups_selection$sex                                    2 12.869  6.4346
    ## orthogroups_selection$selection                              1  0.013  0.0127
    ## orthogroups_selection$sex:orthogroups_selection$selection    2  0.170  0.0848
    ## Residuals                                                 2112 81.527  0.0386
    ##                                                            F value Pr(>F)    
    ## orthogroups_selection$sex                                 166.6903 <2e-16 ***
    ## orthogroups_selection$selection                             0.3301 0.5657    
    ## orthogroups_selection$sex:orthogroups_selection$selection   2.1960 0.1115    
    ## Residuals                                                                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
labsB<-multcompLetters4(aov(mod2b),TukeyHSD(aov(mod2b)))
```

``` r
par(mar=c(4,4,3,1))
vioplot(orthogroups_selection$tau~orthogroups_selection$plot_group,
        xlab="bias category",
        ylab=expression(italic("\u03C4")),
        col=rep(sex_bias_colors,2),
        xaxt='n',
        ylim=c(0,1),
        cex.axis=1,
        cex=1.5,
        frame.plot=FALSE,
        pchMed=21,
        colMed2=rep(sex_bias_colors,2),
        colMed="black",
        las=1)
axis(1,pos=0,at=0:6,labels=c("",
                             paste0("FB\nn=", ns["fem",1]),
                             paste0("MB\nn=", ns["mal",1]),
                             paste0("UB\nn=", ns["unbiased",1]),
                             paste0("FB\nn=", ns["fem",2]),
                             paste0("MB\nn=", ns["mal",2]),
                             paste0("UB\nn=", ns["unbiased",2])),tick = FALSE)

text(rep(labsB$`orthogroups_selection$sex`$Letters,2),x=1:6,y=1)

# add group labels
lines(x=c(1,3),y=c(1.05,1.05),xpd=TRUE)
text(x=2,y=1.1,"no selection", xpd=TRUE)
lines(x=c(4,6),y=c(1.05,1.05),xpd=TRUE)
text(x=5,y=1.1,"selection", xpd=TRUE)
```

``` r
vioplot(orthogroups_selection$tau[orthogroups_selection$selection==1]~orthogroups_selection$plot_group[orthogroups_selection$selection==1],
        xlab="bias category",
        ylab=expression(italic("\u03C4")),
        col=sex_bias_colors,
        xaxt='n',
        ylim=c(0,1),
        cex.axis=1,
        cex=1.5,
        frame.plot=FALSE,
        pchMed=21,
        colMed2=sex_bias_colors,
        colMed="black",
        las=1)
axis(1,pos=0,at=0:3,labels=c("",
                             paste0("FB\nn=", ns["fem",1]),
                             paste0("MB\nn=", ns["mal",1]),
                             paste0("UB\nn=", ns["unbiased",1])),tick = FALSE)

text(c("a","a","b"),x=1:3,y=1)
```

``` r
par(mfrow=c(2,1), mar=c(4,4,1.5,1),oma=c(1,1,1,1))
vioplot(orthogroups_selection$omega ~ orthogroups_selection$plot_group,
        xlab="",
        ylab="",
        main="",
        ylim=c(0,10),
        col=rep(sex_bias_colors, 2),
        cex.axis=1,
        cex=1.5,
        frame.plot=FALSE,
        xaxt='n',
        pchMed=21,
        colMed2=rep(sex_bias_colors,2),
        colMed="black",
        las=1)
axis(1,pos=0,at=0:6,labels=c("",
                             paste0("FB\nn=", ns["fem",1]),
                             paste0("MB\nn=", ns["mal",1]),
                             paste0("UB\nn=", ns["unbiased",1]),
                             paste0("FB\nn=", ns["fem",2]),
                             paste0("MB\nn=", ns["mal",2]),
                             paste0("UB\nn=", ns["unbiased",2])),
     tick=FALSE)
# add letters
text(x=1:6,y=rep(9.25,3),
     c(labs[[3]]$Letters[["fem:0"]],
       labs[[3]]$Letters[["mal:0"]],
       labs[[3]]$Letters[["unbiased:0"]],
       labs[[3]]$Letters[["fem:1"]],
       labs[[3]]$Letters[["mal:1"]],
       labs[[3]]$Letters[["unbiased:1"]]), xpd=TRUE)

# add group labels
lines(x=c(1,3),y=c(10,10))
text(x=2,y=10.5,"no selection", xpd=TRUE)
lines(x=c(4,6),y=c(10,10))
text(x=5,y=10.5,"selection", xpd=TRUE)

# add axis labels

mtext(expression(italic("S. fuscus")~"branch"~italic("\u03C9")), 2, outer=FALSE, line=2)
text(x=-0.1,y=11,"A",cex=2,xpd=TRUE)

# part B
vioplot(orthogroups_selection$tau~orthogroups_selection$plot_group,
        xlab="bias category",
        ylab=expression(italic("\u03C4")),
        col=rep(sex_bias_colors,2),
        xaxt='n',
        ylim=c(0,1),
        cex.axis=1,
        cex=1.5,
        frame.plot=FALSE,
        pchMed=21,
        colMed2=rep(sex_bias_colors,2),
        colMed="black",
        las=1)
axis(1,pos=0,at=0:6,labels=c("",
                             paste0("FB\nn=", ns["fem",1]),
                             paste0("MB\nn=", ns["mal",1]),
                             paste0("UB\nn=", ns["unbiased",1]),
                             paste0("FB\nn=", ns["fem",2]),
                             paste0("MB\nn=", ns["mal",2]),
                             paste0("UB\nn=", ns["unbiased",2])),tick = FALSE)

text(rep(labsB$`orthogroups_selection$sex`$Letters,2),x=1:6,y=1)

# add group labels
lines(x=c(1,3),y=c(1.05,1.05),xpd=TRUE)
text(x=2,y=1.1,"no selection", xpd=TRUE)
lines(x=c(4,6),y=c(1.05,1.05),xpd=TRUE)
text(x=5,y=1.1,"selection", xpd=TRUE)

text(x=-0.1,y=1.1,"B",cex=2,xpd=TRUE)
```

### what are the selected loci?

I used the web-hosted version of eggNOG to predict KEGG annotations from
the S. scovelli genome. This can then be fed into the getEnrich program
(I’m using the online interface). For this to work I also need a list of
background files and a list of enrichment testing files. I will do this
with the single copy orthogroups as the background and then separately
for female-biased, male-biased, and unbiased genes.

``` r
write.table(orthogroups_selection_all$S_scovelli.protein,
            "background_genes_scovelli.txt",
            sep='\t',quote = FALSE,row.names = FALSE,col.names="gene")
write.table(orthogroups_selection$S_scovelli.protein[orthogroups_selection$selection==1 & orthogroups_selection$sex=="unbiased"],
            "unbiased_selection_genes_scovelli.txt",
            sep='\t',quote = FALSE,row.names = FALSE,col.names="gene")
write.table(orthogroups_selection$S_scovelli.protein[orthogroups_selection$selection==1 & orthogroups_selection$sex=="mal"],
            "malebiased_selection_genes_scovelli.txt",
            sep='\t',quote = FALSE,row.names = FALSE,col.names="gene")
write.table(orthogroups_selection$S_scovelli.protein[orthogroups_selection$selection==1 & orthogroups_selection$sex=="fem"],
            "fembiased_selection_genes_scovelli.txt",
            sep='\t',quote = FALSE,row.names = FALSE,col.names="gene")
```

I need to reformat the KEGG IDs too

``` r
keggs<-read.table("MM_at20c6re.emapper.annotations.tsv",
                  sep='\t',
                  comment.char="#")
keggs_for_enrich<-keggs[,c(1,12)]
colnames(keggs_for_enrich)<-c("gene","term")
keggs_for_enrich$term<-gsub("ko:","",keggs_for_enrich$term)
write.table(keggs_for_enrich,
            "S_scovelli_keggs.tsv",
            quote=FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep='\t')
```

``` r
keggs_to_orthogroups<-merge(keggs, orthogroups_selection_all, by.x="V1",by.y="S_scovelli.protein") #1373
```

A total of 1373 orthogroups were annotated for KEGG pathways.

Female-biased enrichment job ID: 82da111a-3324-42ba-86cb-0b1223d31caa
Male-biased enrichment job ID: 5e0c540d-1b47-4cdf-805a-78f9ecc02292
Unbiased enrichment job ID:2b7a5b6a-4bc6-4c26-8741-809f32661785

I would also like to know if the single copy orthologs are enriched for
different categories, so I will need to output the list of all S.
scovelli genes.

``` bash
# in /mnt/BigData/genomes/orthofinder
echo "gene" > all_scovelli_prot_names.txt
grep ">" S_scovelli.protein.faa | sed -e 's/>//g' | sed -e 's/ .*$//g' | uniq - >> all_scovelli_prot_names.txt
```

orthogroups enrichment job ID: 9b85eddb-eb3c-46a5-a2b3-038305675933

``` r
# Let's look up the sex-specific genes under selection
fem_specific_sel<-orthogroups_selection[orthogroups_selection$selection==1 & orthogroups_selection$sex=="fem" & orthogroups_selection$DEG_cat=="specific",]

mal_specific_sel<-orthogroups_selection[orthogroups_selection$selection==1 & orthogroups_selection$sex=="mal" & orthogroups_selection$DEG_cat=="specific",]

keggs[keggs$V1%in% mal_specific_sel$S_scovelli.protein,]
```

``` r
ssc_gff<-parse_gff("genomic.gff")

scov_annotations<-do.call(rbind,lapply(orthogroups_selection_all$S_scovelli.protein, function(prot, gff){
  
  IDs<-grep(prot, gff$description, value=TRUE)
  IDs<-gsub("^.*product=(.*);protein_id.*$","\\1",IDs)
  IDs <- unique(IDs)
  if(length(IDs)==1){
    return(cbind(prot, IDs))
  } else{
    browser()
  }
}, gff=ssc_gff))
write.csv(scov_annotations, "scovelli_ortholog_annotations.csv", quote=FALSE, row.names = FALSE)
```

``` r
scov_annotations<- read.csv("scovelli_ortholog_annotations.csv")
orthogroups_selection_all<-merge(orthogroups_selection_all,
                                 scov_annotations,
                                 by.x="S_scovelli.protein",
                                 by.y="prot")
```

Let’s take a look at those that were sex-specific and under selection –
I can report these in the text.

``` r
orthogroups_selection_all[orthogroups_selection_all$selection==1 & 
                            orthogroups_selection_all$DEG_cat=="specific" &
                            orthogroups_selection_all$sex=="fem",]
```

``` r
orthogroups_selection_all[orthogroups_selection_all$selection==1 & 
                            orthogroups_selection_all$DEG_cat=="specific" &
                            orthogroups_selection_all$sex=="mal",]
```

``` r
# reorder the columns
orthogroups_selection_all<-orthogroups_selection_all[,c(
  "trinity_supertran_fuscus.fasta.transdecoder",
  "tau","Gill_logFC", "Gonad_logFC", "Liver_logFC",
  "DEG_cat", "sex",    "organ",
  "orthogroup",
  "H_comes.protein", "H_zosterae.protein",
   "S_acus.protein", "S_typhle.protein",
  "S_scovelli.protein", 
  "IDs",
  "selection", "omega"
)]
write.csv(orthogroups_selection_all, 
          "orthogroups_selection_annotations.csv",
          quote=FALSE,
          row.names=FALSE)
```

It might also be interesting to have a table of the genes that have
omega \> 1 (and are under selection). Let’s look at that.

``` r
orthogroups_selection_all$omega<-as.numeric(orthogroups_selection_all$omega)
orthogroups_selection_all$selection<-as.numeric(orthogroups_selection_all$selection)

top_omega<-orthogroups_selection_all[orthogroups_selection_all$selection==1 & 
                                       orthogroups_selection_all$omega > 1 &
                                       orthogroups_selection_all$omega <100,]

# sort by omega value
top_omega<- top_omega[order(top_omega$omega, decreasing = TRUE),]

# rename categories
top_omega$DEG_cat <- paste0(top_omega$sex, "-", top_omega$DEG_cat, "_", top_omega$organ)
top_omega$DEG_cat[top_omega$DEG_cat=="unbiased-unbiased_unbiased"]<-"unbiased"

# reorder the columns and remove unnecessary ones
top_omega<-top_omega[,c("trinity_supertran_fuscus.fasta.transdecoder",
                        "DEG_cat",
                        "tau",
                        "Gonad_logFC",
                        "omega",
                        "orthogroup",
                        "IDs"
                        )]

# rename the columns
colnames(top_omega)<-c("S. fuscus Trinity protein ID",
                        "Differential expression category",
                       "tau",
                       "log Fold Change in gonad",
                        "S. fuscus branch omega",
                        "orthogroup ID",
                        "Product name from S. scovelli genome"
                        )
write.csv(top_omega, "genes_with_omega_above_1.csv", row.names=FALSE, quote=FALSE)
```

``` r
kable(top_omega, "latex")
```

``` r
point_cols<-orthogroups_selection$sex
point_cols[point_cols=="mal"]<-"MB"
point_cols[point_cols=="fem"]<-"FB"
point_cols[point_cols=="unbiased"]<-"UB"

shapes<-rep(15,length(point_cols))
shapes[which(point_cols=="MB")]<-16
shapes[which(point_cols=="UB")]<-2
plot(rank(orthogroups_selection$omega)~rank(orthogroups_selection$tau),
     col=sex_bias_colors[point_cols],
     pch=shapes)
```
