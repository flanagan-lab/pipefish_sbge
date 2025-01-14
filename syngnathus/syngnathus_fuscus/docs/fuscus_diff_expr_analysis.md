Differential Expression Analysis in *Syngnathus fuscus*
================
Coley Tosto
2024-12-16



- [Single factor analysis - Comparing Males v
  Females](#single-factor-analysis---comparing-males-v-females)
  - [Visualizing the results](#visualizing-the-results)
    - [MA-plot - MvF Diff.](#ma-plot---mvf-diff)
    - [Heatmap - Looking at patterns of overall
      expression](#heatmap---looking-at-patterns-of-overall-expression)
    - [Sample-dist Heatmap](#sample-dist-heatmap)
    - [PCA plots](#pca-plots)
- [Multifactor design - Comparing M v F across the diff tissue
  types](#multifactor-design---comparing-m-v-f-across-the-diff-tissue-types)
  - [Invesitgate the results of the differential
    expression](#invesitgate-the-results-of-the-differential-expression)
    - [M-F Liver Comparisson](#m-f-liver-comparisson)
    - [M-F Gill Comparisson](#m-f-gill-comparisson)
    - [M-F Gonad Comparisson](#m-f-gonad-comparisson)
    - [MA Plots](#ma-plots)
  - [Variation in FC across sex-bias and tissue
    type](#variation-in-fc-across-sex-bias-and-tissue-type)
- [Running a single factor analysis between males and females WITHIN
  each
  organ](#running-a-single-factor-analysis-between-males-and-females-within-each-organ)
  - [Liver M-F differential
    expression](#liver-m-f-differential-expression)
  - [Gonad M-F differential
    expression](#gonad-m-f-differential-expression)
  - [Gill M-F differential
    expression](#gill-m-f-differential-expression)
  - [Re-looking at variation in logFC across both sex-bias categories
    and tissue
    types](#re-looking-at-variation-in-logfc-across-both-sex-bias-categories-and-tissue-types)
  - [Investigating expression across samples in the sex-biased
    genes](#investigating-expression-across-samples-in-the-sex-biased-genes)
  - [Looking at shared sex-biased genes with an Upset
    plot](#looking-at-shared-sex-biased-genes-with-an-upset-plot)
  - [Categorizing sex-specific genes](#categorizing-sex-specific-genes)
    - [Investigating the results of the sex-specific
      subsetting](#investigating-the-results-of-the-sex-specific-subsetting)
    - [Pulling out the gene IDs for all of the sex-specific
      genes](#pulling-out-the-gene-ids-for-all-of-the-sex-specific-genes)
  - [Creating categories and binning the sex-biased genes based on
    degree of
    logFC](#creating-categories-and-binning-the-sex-biased-genes-based-on-degree-of-logfc)
  - [Create a combined figure](#create-a-combined-figure)
- [Gene Ontologogy Analysis and BLASTing the different
  genes](#gene-ontologogy-analysis-and-blasting-the-different-genes)
  - [BLASTING against the zebrafish
    proteome](#blasting-against-the-zebrafish-proteome)
  - [Read in and filter the BLAST
    results](#read-in-and-filter-the-blast-results)
  - [Collect gene names corresponding to each protein -sex-biased
    genes](#collect-gene-names-corresponding-to-each-protein--sex-biased-genes)
  - [Pull out the gene names and protein ID for sex-specific
    BLAST](#pull-out-the-gene-names-and-protein-id-for-sex-specific-blast)
  - [Reading in the PANTHER results - Sex-biased
    Genes](#reading-in-the-panther-results---sex-biased-genes)
  - [Looking at GO groups for the male- and female-biased
    genes](#looking-at-go-groups-for-the-male--and-female-biased-genes)
  - [Create combined Gene Ontology
    Figure](#create-combined-gene-ontology-figure)

``` r
#The abundance matrix generated via salmon and tximport to be used for the DE analysis
FU_txi.salmon <- readRDS("data/txi.salmon_FU.RDS")

#The samples file generated for tximport
FU_samples <- read.table("FU_samples.txt", header = TRUE)

#Make sure the conditions are in the samples file as a factor
FU_samples$Sex <- as.factor(FU_samples$Sex)
FU_samples$Organ <- as.factor(FU_samples$Organ)
```

The package `DESeq2` was used for the differential expression analysis
outlined below.

# Single factor analysis - Comparing Males v Females

I am starting with a single factor analysis to look at overall
differences between male and female expression across the organs and to
do some exploratory analysis of the data.

To analyze my data with DESeq2 I must first generate the DESeqDataSet.
In order to do this we need the abundance matrix generated with
`tximport` which is stored in the object “FU_txi.salmon” and was
generated at the end of the `Analyzing_fuscus_RNAseq_data_from_MSU.Rmd`
document and a `samples` file which is stored in the object “FU_samples”
that lays out all of the conditions.

``` r
#Create the DESeq dataset
dds_FU <- DESeqDataSetFromTximport(FU_txi.salmon, 
                                   colData = FU_samples,
                                   design = ~ Sex)
```

I am then going to pre-filter the data to remove low gene counts before
running further DESeq2 functions. By doing this we remove rows in which
there are very few reads thus reducing the memory size of the `dds`
object and increasing the speed at which we can use the transformation
and testing functions in DESeq2.

The cutoff used here was to remove reads that had a **total** count of
less than 10 in each row.

``` r
#only keeping rows that have at lead 10 reads total
keep <- rowSums(counts(dds_FU)) >= 10
dds_FU <- dds_FU[keep, ]
```

After filtering we can now perform the standard differential expression
analysis that is wrapped into DESeq2.

``` r
#Generate the expression values
dds_FU_exp <- DESeq(dds_FU)

#Compile the results
res <- results(dds_FU_exp)
res
```

    ## log2 fold change (MLE): Sex M vs F 
    ## Wald test p-value: Sex M vs F 
    ## DataFrame with 193637 rows and 6 columns
    ##                         baseMean log2FoldChange     lfcSE       stat    pvalue
    ##                        <numeric>      <numeric> <numeric>  <numeric> <numeric>
    ## TRINITY_DN0_c0_g1     1224.15627     -0.7105364  0.255547 -2.7804536 0.0054283
    ## TRINITY_DN0_c1_g1     4528.38058     -0.4239128  0.297621 -1.4243368 0.1543490
    ## TRINITY_DN0_c10_g1     552.52273     -0.0126323  0.182706 -0.0691402 0.9448780
    ## TRINITY_DN0_c11_g1      28.55379      0.7130072  0.506622  1.4073758 0.1593160
    ## TRINITY_DN0_c125_g1      2.62197      0.2234249  0.820431  0.2723262 0.7853712
    ## ...                          ...            ...       ...        ...       ...
    ## TRINITY_DN99991_c0_g1   0.644610       1.419240  1.615112   0.878725  0.379550
    ## TRINITY_DN99993_c0_g1   2.586928       0.660824  0.812081   0.813741  0.415793
    ## TRINITY_DN99995_c0_g1   0.361104       0.996595  2.973076   0.335206  0.737469
    ## TRINITY_DN99997_c0_g1   1.841278      -0.990609  1.293794  -0.765662  0.443878
    ## TRINITY_DN99999_c0_g1   0.452691      -2.016813  2.530238  -0.797084  0.425402
    ##                            padj
    ##                       <numeric>
    ## TRINITY_DN0_c0_g1      0.059560
    ## TRINITY_DN0_c1_g1      0.460392
    ## TRINITY_DN0_c10_g1     0.987621
    ## TRINITY_DN0_c11_g1     0.468337
    ## TRINITY_DN0_c125_g1    0.943650
    ## ...                         ...
    ## TRINITY_DN99991_c0_g1        NA
    ## TRINITY_DN99993_c0_g1  0.752957
    ## TRINITY_DN99995_c0_g1        NA
    ## TRINITY_DN99997_c0_g1  0.775208
    ## TRINITY_DN99999_c0_g1        NA

Once that has finished we can now start exploring some of the
single-factor analysis results.

``` r
##Ordering our results based on p-value
resOrdered <- res[order(res$pvalue),]
resOrdered
```

    ## log2 fold change (MLE): Sex M vs F 
    ## Wald test p-value: Sex M vs F 
    ## DataFrame with 193637 rows and 6 columns
    ##                        baseMean log2FoldChange     lfcSE      stat      pvalue
    ##                       <numeric>      <numeric> <numeric> <numeric>   <numeric>
    ## TRINITY_DN22310_c0_g1  804.2311        27.6241   1.89885   14.5477 6.03587e-48
    ## TRINITY_DN13824_c0_g1   92.1347        25.5625   1.92144   13.3039 2.19822e-40
    ## TRINITY_DN35004_c1_g2  129.5764        26.0372   2.10520   12.3680 3.89194e-35
    ## TRINITY_DN34941_c0_g1  260.2356        27.0195   2.24241   12.0493 1.95608e-33
    ## TRINITY_DN9279_c0_g1    23.3464        23.3489   1.95493   11.9436 7.00895e-33
    ## ...                         ...            ...       ...       ...         ...
    ## TRINITY_DN99935_c0_g1         0              0         0         0           1
    ## TRINITY_DN99938_c0_g1         0              0         0         0           1
    ## TRINITY_DN99956_c0_g1         0              0         0         0           1
    ## TRINITY_DN99961_c0_g1         0              0         0         0           1
    ## TRINITY_DN99973_c0_g1         0              0         0         0           1
    ##                              padj
    ##                         <numeric>
    ## TRINITY_DN22310_c0_g1 3.82523e-43
    ## TRINITY_DN13824_c0_g1 6.96562e-36
    ## TRINITY_DN35004_c1_g2 8.22172e-31
    ## TRINITY_DN34941_c0_g1 3.09916e-29
    ## TRINITY_DN9279_c0_g1  8.88384e-29
    ## ...                           ...
    ## TRINITY_DN99935_c0_g1          NA
    ## TRINITY_DN99938_c0_g1          NA
    ## TRINITY_DN99956_c0_g1          NA
    ## TRINITY_DN99961_c0_g1          NA
    ## TRINITY_DN99973_c0_g1          NA

``` r
summary(res)
```

    ## 
    ## out of 174116 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 4003, 2.3%
    ## LFC < 0 (down)     : 3581, 2.1%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 130262, 75%
    ## (mean count < 1)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
#Looking at  an alpha=0.05, the default is 0.1
res05 <- results(dds_FU_exp, alpha = 0.05)
summary(res05)
```

    ## 
    ## out of 174116 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 2626, 1.5%
    ## LFC < 0 (down)     : 2488, 1.4%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 123551, 71%
    ## (mean count < 1)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
sum(res05$padj < 0.05, na.rm = TRUE)
```

    ## [1] 5114

Looking at the summary of our results we can see that there are 2,626
male-biased genes and 2,488 female-biased genes when all of the organs
are organized together. While this is informative, I am more interested
in a multi-factor approach that takes into account both **Sex** and
**Tissue** type. Before the multi-factor analysis, we can visualize some
results from the single-factor analysis and use it for exploratory data
analysis.

Interestingly, the summary is saying that there are 0 outliers, which I
have a hard time believing.

## Visualizing the results

### MA-plot - MvF Diff.

Generate an MA-plot to show the log2 fold changes attributable to sex
over the mean of normalized counts for all of the samples in the `dds`.
Points will be colored if the adjusted p-value is less that 0.05.

<figure>
<img src="fuscus_diff_expr_analysis_files/figure-gfm/MA-plot-1.png"
alt="LogFC versus the mean normalized count for all of the genes. Points that are blue have a p-value less than 0.05. The two plots are the same with the only difference coming from an adjusted y-limit. Points that are triangles represent genes with a fold change higher/lower than the y-limit." />
<figcaption aria-hidden="true">LogFC versus the mean normalized count
for all of the genes. Points that are blue have a p-value less than
0.05. The two plots are the same with the only difference coming from an
adjusted y-limit. Points that are triangles represent genes with a fold
change higher/lower than the y-limit.</figcaption>
</figure>

### Heatmap - Looking at patterns of overall expression

We can also generate a heat map to look ar overall expression levels
across our samples. Note, this is not differentially expressed genes.

<figure>
<img src="fuscus_diff_expr_analysis_files/figure-gfm/heatmap-1.png"
alt="Heatmap showing the expression level across all organs for the top 20 genes with the highest expression levels." />
<figcaption aria-hidden="true">Heatmap showing the expression level
across all organs for the top 20 genes with the highest expression
levels.</figcaption>
</figure>

From the heatmap we can see that for a lot of the top expressed genes
the highest expression are found in the liver. There appears to be a few
genes that are highly expressed across all tissues and a few genes that
are highly expressed in the gill tissues relative to the other tissues.
With this heatmap we can also pull out the names of the Trinity genes
that are showing this high expression and BLAST the corresponding
sequences to see what the genes are.

``` r
#Pull out the corresponding trinity_geneIDS that are plotted in the heatmap
heatmap_TG <- cbind(row.names(assay(vsd)[select,]))

#Export the gene names
write.table(heatmap_TG,
            'FU_heatmap_trinitygenes.txt',
            sep = "",
            quote=FALSE,
            row.names = FALSE,
            col.names = FALSE)
```

### Sample-dist Heatmap

We can then generate a sample-to-sample distances heatmap that will give
an overview of the similarities and differences between our samples.

<figure>
<img src="fuscus_diff_expr_analysis_files/figure-gfm/sample-dist-1.png"
alt="Sample-to-sample distances heatmap showing the similarities and differences between all of the samples. The darker the color, the more similar the samples are. The diagonal line of very dark blue represents the comparissons between the same samples." />
<figcaption aria-hidden="true">Sample-to-sample distances heatmap
showing the similarities and differences between all of the samples. The
darker the color, the more similar the samples are. The diagonal line of
very dark blue represents the comparissons between the same
samples.</figcaption>
</figure>

We can see that the highest similarities (aside from same samples
comparisons) is between samples from the same organ. After that we can
see high similarities between the male gonads and the female gonads.
There are also a couple samples that may be an issue for use as they are
not similar to any of the other samples from the same group. This
includes FUG11M4 and FUT11M4. I may end up removing all organs related
to 11M4.

### PCA plots

Next we can look at PCA plots for our samples to see how our different
samples are grouping.

    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.
    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.
    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.
    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.
    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.
    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.

<figure>
<img
src="fuscus_diff_expr_analysis_files/figure-gfm/pca-pairs-plot-1.png"
alt="Principal components analysis pairsplot showing PCA axis 1 - 4." />
<figcaption aria-hidden="true">Principal components analysis pairsplot
showing PCA axis 1 - 4.</figcaption>
</figure>

From the pairs plot we can see that PC1 appears to be explaining
differences between the liver and the other tissues, PC2 appears to be
explaining differences between gonadal and somatic tissues, PC3 seems to
account for differences between male and female samples, especially in
the gonads and in PC4 we can see the one problematic gill sample from
above sticking out from all the rest. I am going to plot all of PC2 - 4
against PC1 below for a clearer picture.

<figure>
<img src="fuscus_diff_expr_analysis_files/figure-gfm/pca-plot-1.png"
alt="Principal components analysis reflects that most variation corresponds to differences in expression between organs (green: gill; yellow: gonad; pink: liver), rather than variation due to sex (circle = female; triangle = male). The first axis explains 31% of variation in the dataset, and largely explains differences between the liver and the other tissues. The second axis explains 20% of the variation in gene expression, with somatic tissues on the opposite side of the axis from the gonads." />
<figcaption aria-hidden="true">Principal components analysis reflects
that most variation corresponds to differences in expression between
organs (green: gill; yellow: gonad; pink: liver), rather than variation
due to sex (circle = female; triangle = male). The first axis explains
31% of variation in the dataset, and largely explains differences
between the liver and the other tissues. The second axis explains 20% of
the variation in gene expression, with somatic tissues on the opposite
side of the axis from the gonads.</figcaption>
</figure>

<figure>
<img src="fuscus_diff_expr_analysis_files/figure-gfm/pca1v3-1.png"
alt="The same Principal Components Analysis as above but with PC3. This third axis explains 10% of the variation, largely attributed to differences between the sexes, particularily in the gonads." />
<figcaption aria-hidden="true">The same Principal Components Analysis as
above but with PC3. This third axis explains 10% of the variation,
largely attributed to differences between the sexes, particularily in
the gonads.</figcaption>
</figure>

<figure>
<img src="fuscus_diff_expr_analysis_files/figure-gfm/pca1v4-1.png"
alt="The same Principal Components Analysis as above but with PC4. This fourth axis explains 6% of the variation, largely attributed to one of the gill samples." />
<figcaption aria-hidden="true">The same Principal Components Analysis as
above but with PC4. This fourth axis explains 6% of the variation,
largely attributed to one of the gill samples.</figcaption>
</figure>

#### Saving the PCA plots to use in figure creation

I will be including the PCAs describing the relationship between PC1 - 3
in the main text, as they account for the majority of the differences
between our samples. The plot with PC1 and PC4 I still want to export,
and will likely include as a supplemental to highlight why samples
related to 11M4 will be removed from the subsequent analyses.

``` r
#Setting the shapes I want to show up for the different organs
organ_shapes <- c(Gill = 16, Liver = 17, Gonad = 15)

#Setting the colors I want for each sex
sex_cols <- c("F" = "#7fc97f", "M" = "#beaed4")
```

``` r
#Create the blank pdf to store the plot in
pdf("docs/figs/Fig_PCA1v2.pdf",height = 6,width=6)

#Set the plotting parameters
par(mar=c(4,5,4,1), oma=c(2,2,2,2))

#Create the plot for PC1 v PC2
plot(pca_plotting_data$PC1,
     pca_plotting_data$PC2,
     col = paste0(sex_cols[pca_plotting_data$Sex],"75"),
     pch = organ_shapes[pca_plotting_data$Organ],
     cex = 2,
     cex.lab = 2,
     cex.axis = 1.75,
     xlab = paste0("PC1: ",percentVar[1], "% variance"),
     ylab = paste0("PC2: ",percentVar[2], "% variance"),
     bty = 'l',
     xpd = TRUE)

#Add a legend describing the sex
outer_legend("top",
             c("Female","Male"),
             pch = 18,
             bty = 'n',
             col=paste0(sex_cols,"75"),
             cex = 1.75,
             ncol = 2,
             pt.cex = 3)

dev.off()
```

    ## png 
    ##   2

``` r
#Create the blank pdf to store the plot in
pdf("docs/figs/Fig_PCA1v3.pdf",height = 6,width=6)

#Set the plotting parameters
par(mar=c(4,5,4,1), oma=c(2,2,2,2))

#Create the plot for PC1 v PC3
plot(pca_plotting_data$PC1,
     pca_plotting_data$PC3,
     col = paste0(sex_cols[pca_plotting_data$Sex],"75"),
     pch = organ_shapes[pca_plotting_data$Organ],
     cex = 2,
     cex.lab = 2,
     cex.axis = 1.75,
     xlab = paste0("PC1: ",percentVar[1], "% variance"),
     ylab = paste0("PC3: ",percentVar[3], "% variance"),
     bty = 'l',
     xpd = TRUE)

#Add  legend describing the organs
outer_legend("top",
             c("Gill","Gonad","Liver"),
             pch = c(16,15,17),
             bty = 'n',
             col = 'darkgrey',
             cex = 1.75,
             ncol = 3,
             pt.cex = 3)

dev.off()
```

    ## png 
    ##   2

``` r
#Create the blank pdf to store the plot in
pdf("docs/figs/Fig_PCA1v4.pdf",height = 6,width=8)

#Set the plotting parameters
par(mar=c(4,5,1,5),oma=c(2,2,2,2))

#Create the plot for PC1 v PC4
plot(pca_plotting_data$PC1,
     pca_plotting_data$PC4,
     col = paste0(sex_cols[pca_plotting_data$Sex],"75"),
     pch = organ_shapes[pca_plotting_data$Organ],
     cex = 2,
     cex.lab = 2,
     cex.axis = 1.75,
     xlab = paste0("PC1: ",percentVar[1], "% variance"),
     ylab = paste0("PC4: ",percentVar[4], "% variance"),
     bty = 'l',
     xpd = TRUE)

text(pca_plotting_data$PC1[pca_plotting_data$ID == "FUG11M4"] + 50,
     pca_plotting_data$PC4[pca_plotting_data$ID == "FUG11M4"] - 25,
     labels = pca_plotting_data$ID[pca_plotting_data$ID == "FUG11M4"])

#Add  legend describing the organs and sex
outer_legend("right",
             c("Gill","Gonad","Liver"),
             pch = c(16,15,17),
             bty = 'n',
             col = 'darkgrey',
             cex = 1.5,
             ncol = 1,
             pt.cex = 2.25)

outer_legend("topright",
             c("Female","Male"),
             pch = 18,
             bty = 'n',
             col=paste0(sex_cols,"75"),
             cex = 1.5,
             ncol = 1,
             pt.cex = 2.25)
dev.off()
```

    ## png 
    ##   2

#### Creating heatmaps based on the PCA axes

I now want to generate heatmaps based on the genes that have the highest
loadings for each of the PCA axes. I will be using a cut off of 0.02 to
determine which genes will be includd in the heatmaps. Once again I am
doing this for PCA1-4, but will likely only inlcude PCA1-3 in a figure
for the maintext and PC4 will go in the supplemental.

``` r
df <- as.data.frame(FU_samples[,c("Sex", "Organ")])
rownames(df) <- FU_samples$ID
df$Organ <- as.character(df$Organ)

organ_cols<-c("Gill" = "#20B2AA",
              "Gonad" = "#EEB422", 
              "Liver" = "#EE8262")

ann_colors = list(Sex=sex_cols,
                  Organ=organ_cols)

col_order <- c(rownames(df[df$Sex=="F" & df$Organ=="Gill",]),
               rownames(df[df$Sex=="M" & df$Organ=="Gill",]),
               rownames(df[df$Sex=="F" & df$Organ=="Gonad",]),
               rownames(df[df$Sex=="M" & df$Organ=="Gonad",]),
               rownames(df[df$Sex=="F" & df$Organ=="Liver",]),
               rownames(df[df$Sex=="M" & df$Organ=="Liver",]))

pca_rotation <- pca$rotation[, 1:4]
```

``` r
#Function to use to export the heatmaps to pdfs
# from https://stackoverflow.com/questions/43051525/how-to-draw-pheatmap-plot-to-screen-and-also-save-to-file
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
```

``` r
pc1 <- pheatmap(assay(vsd)[which(abs(pca_rotation[,1]) >= 0.02), col_order], 
                cluster_rows = FALSE, 
                show_rownames = FALSE, 
                cluster_cols = FALSE, 
                show_colnames = FALSE,
                annotation_col = df,
                annotation_colors = ann_colors,
                cellwidth = 9,
                fontsize = 16,
                annotation_legend = FALSE,
                main = "Top loading genes on PC1")
```

![](fuscus_diff_expr_analysis_files/figure-gfm/figheatmap-PC1-1.png)<!-- -->

``` r
save_pheatmap_pdf(pc1, "docs/figs/Fig_pc1_heatmap.pdf",
                  width=6,
                  height=6)
```

    ## png 
    ##   2

``` r
pc2 <- pheatmap(assay(vsd)[which(abs(pca_rotation[,2]) >= 0.02), col_order], 
                cluster_rows = FALSE, 
                show_rownames = FALSE, 
                cluster_cols = FALSE, 
                show_colnames = FALSE,
                annotation_col = df,
                annotation_colors=ann_colors,
                cellwidth = 9,
                fontsize = 16,
                annotation_legend = FALSE,
                main = "Top loading genes on PC2")
```

![](fuscus_diff_expr_analysis_files/figure-gfm/figheatmap-PC2-1.png)<!-- -->

``` r
save_pheatmap_pdf(pc2, "docs/figs/Fig_pc2_heatmap.pdf",
                  width=6,
                  height=6)
```

    ## png 
    ##   2

``` r
pc3 <- pheatmap(assay(vsd)[which(abs(pca_rotation[,3]) >= 0.02), col_order], 
                cluster_rows = FALSE, 
                show_rownames = FALSE, 
                cluster_cols = FALSE, 
                show_colnames = FALSE,
                annotation_col = df,
                annotation_colors=ann_colors,
                cellwidth = 9,
                fontsize = 16,
                border_color = NA,
                main = "Top loading genes on PC3")
```

![](fuscus_diff_expr_analysis_files/figure-gfm/figheatmap-PC3-1.png)<!-- -->

``` r
save_pheatmap_pdf(pc3, "docs/figs/Fig_pc3_heatmap.pdf",
                  width=6,
                  height=6)
```

    ## png 
    ##   2

``` r
pc4 <- pheatmap(assay(vsd)[which(abs(pca_rotation[,4]) >= 0.02), col_order], 
                cluster_rows = FALSE, 
                show_rownames = FALSE, 
                cluster_cols = FALSE, 
                show_colnames = FALSE,
                annotation_col = df,
                annotation_colors=ann_colors,
                cellwidth = 9,
                fontsize = 16,
                border_color = NA,
                main = "Top loading genes on PC4")
```

![](fuscus_diff_expr_analysis_files/figure-gfm/figheatmap-PC4-1.png)<!-- -->

``` r
save_pheatmap_pdf(pc4, "docs/figs/Fig_pc4_heatmap.pdf",
                  width=6,
                  height=6)
```

    ## png 
    ##   2

#### Making the combined figure

``` r
figPCAa <- image_ggplot(image_read_pdf('docs/figs/Fig_PCA1v2.pdf'),
                        interpolate = TRUE)
figPCAb <- image_ggplot(image_read_pdf('docs/figs/Fig_PCA1v3.pdf'),
                        interpolate = TRUE)
figPCAc <- image_ggplot(image_read_pdf('docs/figs/Fig_pc1_heatmap.pdf'),
                        interpolate = TRUE)
figPCAd <- image_ggplot(image_read_pdf('docs/figs/Fig_pc2_heatmap.pdf'),
                        interpolate = TRUE)
figPCAe <- image_ggplot(image_read_pdf('docs/figs/Fig_pc3_heatmap.pdf'),
                        interpolate = TRUE)

# make two patchworks
pcas <- figPCAa + figPCAb
hms <- figPCAc + figPCAd + figPCAe

# put the patchworks together
figPCA <-  wrap_plots(pcas,
                      hms,
                      ncol=1) + 
  plot_annotation(tag_levels = 'A') 

ggsave("docs/figs/FigPCA.pdf", figPCA, height=4, width=5)
ggsave("docs/figs/FigPCA.png", figPCA, height=4, width=5) # also save as a png
```

``` r
figPCA4 <- image_ggplot(image_read_pdf('docs/figs/Fig_PCA1v4.pdf'),
                        interpolate = TRUE)
figPCA4_hm <- image_ggplot(image_read_pdf('docs/figs/Fig_pc4_heatmap.pdf'),
                           interpolate = TRUE)

# put the patchworks together
figPCA_supp <-  wrap_plots(figPCA4, figPCA4_hm, ncol=2) + 
  plot_annotation(tag_levels = 'A') 

ggsave("docs/figs/FigPCA_supp.pdf", figPCA_supp, height=4, width=8)
ggsave("docs/figs/FigPCA_supp.png", figPCA_supp, height=4, width=8) # also save as a png
```

``` r
figPCA
```

![](fuscus_diff_expr_analysis_files/figure-gfm/showAssembledFigPCA-1.png)<!-- -->

# Multifactor design - Comparing M v F across the diff tissue types

If we investigate the column data of our DESeq dataset we can see that
each sample has both a sex and organ type attached to it, we will be
using these two factors in our multi-factor analysis. Prior to running
the multi-factor analysis I am going to remove the samples that
correspond to the individual 11M4.

``` r
#Create an additional column the has the two conditions combined(sex and organ type)
FU_samples$group <- factor(paste0(FU_samples$Sex, FU_samples$Organ))

##Create a copy of the DESeq dataset
ddsMF_FU <- DESeqDataSetFromTximport(FU_txi.salmon,
                                     colData = FU_samples,
                                     design = ~ group)

##Remove all 11M4 organs from the dataset
ddsMF_FU <- ddsMF_FU[, !(ddsMF_FU$ID %in% c("FUT11M4", "FUG11M4", "FUL11M4"))]

##Filter the dataset, only keeping rows that have at least 10 reads total
keep <- rowSums(counts(ddsMF_FU)) >= 10
ddsMF_FU <- ddsMF_FU[keep, ]

#Run the differential expression analysis
ddsMF_FU_exp <- DESeq(ddsMF_FU)
resultsNames(ddsMF_FU_exp)
```

    ## [1] "Intercept"             "group_FGonad_vs_FGill" "group_FLiver_vs_FGill"
    ## [4] "group_MGill_vs_FGill"  "group_MGonad_vs_FGill" "group_MLiver_vs_FGill"

``` r
#Looking at the Normalized read counts for each sample
cbind(FU_samples[!(FU_samples$ID %in% 
                     c("FUT11M4", "FUG11M4", "FUL11M4")),],
      normalizedReads=colSums(counts(ddsMF_FU_exp, 
                                     normalized=TRUE)
                              )
      )
```

    ##         ID Sex Organ  group normalizedReads
    ## 1  FUG10M2   M  Gill  MGill        24127581
    ## 2  FUG11F1   F  Gill  FGill        25620273
    ## 3  FUG11M2   M  Gill  MGill        24773886
    ## 5  FUG12M1   M  Gill  MGill        25685953
    ## 6  FUG13F1   F  Gill  FGill        26212779
    ## 7  FUG13F4   F  Gill  FGill        24803053
    ## 8  FUG15M5   M  Gill  MGill        23866915
    ## 9   FUG2F2   F  Gill  FGill        28033979
    ## 10  FUG3F2   F  Gill  FGill        25512447
    ## 11 FUL10M2   M Liver MLiver        94425427
    ## 12 FUL11F1   F Liver FLiver        82595903
    ## 13 FUL11M2   M Liver MLiver       118909536
    ## 15 FUL12M1   M Liver MLiver       114375220
    ## 16 FUL13F1   F Liver FLiver        97875307
    ## 17 FUL13F4   F Liver FLiver        90282028
    ## 18 FUL15M5   M Liver MLiver       117231749
    ## 19  FUL2F2   F Liver FLiver       241462147
    ## 20  FUL3F2   F Liver FLiver        69932839
    ## 21 FUO11F1   F Gonad FGonad        32378203
    ## 22 FUO13F1   F Gonad FGonad        29331475
    ## 23 FUO13F4   F Gonad FGonad        25761811
    ## 24  FUO2F2   F Gonad FGonad        26870920
    ## 25  FUO3F2   F Gonad FGonad        29802282
    ## 26 FUT10M2   M Gonad MGonad        21635189
    ## 27 FUT11M2   M Gonad MGonad        21114528
    ## 29 FUT12M1   M Gonad MGonad        20609565
    ## 30 FUT15M5   M Gonad MGonad        20385206

If we look at the summary of normalized reads across the different
tissue types, we can see that there are considerably more in the liver
relative to the other tissue types. Because during a multi-factor
analysis normalization occurs across all of the tissues, I am opting to
also remove some of the extremely expressed genes to try and offset the
differences.

``` r
##Filter the dataset, only keeping rows that have at least 10 reads total AND less than one million
keep <- rowSums(counts(ddsMF_FU)) >= 10 & rowSums(counts(ddsMF_FU)) < 1e6
ddsMF_FU <- ddsMF_FU[keep, ]

#Re-run the differential expression analysis
ddsMF_FU_exp <- DESeq(ddsMF_FU)

#Re investigate the Normalized read counts for each sample
cbind(FU_samples[!(FU_samples$ID %in% 
                     c("FUT11M4", "FUG11M4", "FUL11M4")),],
      normalizedReads=colSums(counts(ddsMF_FU_exp, 
                                     normalized=TRUE)
                              )
      )
```

    ##         ID Sex Organ  group normalizedReads
    ## 1  FUG10M2   M  Gill  MGill        19883120
    ## 2  FUG11F1   F  Gill  FGill        20865455
    ## 3  FUG11M2   M  Gill  MGill        20543583
    ## 5  FUG12M1   M  Gill  MGill        20858795
    ## 6  FUG13F1   F  Gill  FGill        21316492
    ## 7  FUG13F4   F  Gill  FGill        20229015
    ## 8  FUG15M5   M  Gill  MGill        19923606
    ## 9   FUG2F2   F  Gill  FGill        22070324
    ## 10  FUG3F2   F  Gill  FGill        21008945
    ## 11 FUL10M2   M Liver MLiver        38972714
    ## 12 FUL11F1   F Liver FLiver        35640412
    ## 13 FUL11M2   M Liver MLiver        42259174
    ## 15 FUL12M1   M Liver MLiver        45477421
    ## 16 FUL13F1   F Liver FLiver        40097766
    ## 17 FUL13F4   F Liver FLiver        36257031
    ## 18 FUL15M5   M Liver MLiver        43583313
    ## 19  FUL2F2   F Liver FLiver        67407162
    ## 20  FUL3F2   F Liver FLiver        33148924
    ## 21 FUO11F1   F Gonad FGonad        29957898
    ## 22 FUO13F1   F Gonad FGonad        27326742
    ## 23 FUO13F4   F Gonad FGonad        23429206
    ## 24  FUO2F2   F Gonad FGonad        24666955
    ## 25  FUO3F2   F Gonad FGonad        28109050
    ## 26 FUT10M2   M Gonad MGonad        19290194
    ## 27 FUT11M2   M Gonad MGonad        18967681
    ## 29 FUT12M1   M Gonad MGonad        18757421
    ## 30 FUT15M5   M Gonad MGonad        18429583

With this change we can now see that there is more even levels of
expression across the different tissue types.

## Invesitgate the results of the differential expression

Thanks to the multi-factor analysis we can now explore differential
expression between all of the different combinations of tissues and
sexes:

- Male Liver v. Female Liver
- Male Gill v. Female Gill  
- Male Gonad v. Female Gonad  
- All of the within sex tissue comparisons (e.g. Male Liver v. Male
  Gill, etc.)

### M-F Liver Comparisson

``` r
##Pulling out the liver M-F results with an alpha of 0.05
FU_liver_con_res <- results(ddsMF_FU_exp, 
                            contrast = c("group", "MLiver", "FLiver"), 
                            alpha = 0.05)
FU_liver_con_res$trin_geneid <- rownames(FU_liver_con_res)
head(FU_liver_con_res)
```

    ## log2 fold change (MLE): group MLiver vs FLiver 
    ## Wald test p-value: group MLiver vs FLiver 
    ## DataFrame with 6 rows and 7 columns
    ##                       baseMean log2FoldChange     lfcSE      stat     pvalue
    ##                      <numeric>      <numeric> <numeric> <numeric>  <numeric>
    ## TRINITY_DN0_c0_g1   1259.36897      -0.946137  0.310797 -3.044227 0.00233279
    ## TRINITY_DN0_c1_g1   4676.25499       0.230724  0.265546  0.868866 0.38492048
    ## TRINITY_DN0_c10_g1   563.87777       0.199304  0.339963  0.586252 0.55770625
    ## TRINITY_DN0_c11_g1    29.98589       0.847210  0.747662  1.133146 0.25715309
    ## TRINITY_DN0_c125_g1    2.23427      -1.006939  1.512383 -0.665797 0.50554118
    ## TRINITY_DN0_c13_g1    24.89771       0.986249  0.775161  1.272314 0.20326140
    ##                          padj         trin_geneid
    ##                     <numeric>         <character>
    ## TRINITY_DN0_c0_g1   0.0284562   TRINITY_DN0_c0_g1
    ## TRINITY_DN0_c1_g1   0.7002120   TRINITY_DN0_c1_g1
    ## TRINITY_DN0_c10_g1  0.8367827  TRINITY_DN0_c10_g1
    ## TRINITY_DN0_c11_g1  0.5683710  TRINITY_DN0_c11_g1
    ## TRINITY_DN0_c125_g1        NA TRINITY_DN0_c125_g1
    ## TRINITY_DN0_c13_g1  0.5018550  TRINITY_DN0_c13_g1

``` r
summary(FU_liver_con_res)
```

    ## 
    ## out of 157005 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 2302, 1.5%
    ## LFC < 0 (down)     : 1289, 0.82%
    ## outliers [1]       : 5538, 3.5%
    ## low counts [2]     : 114872, 73%
    ## (mean count < 6)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

I will be classifying sex-biased genes the same way that I did for
*Syngnathus floridae* as **genes with a p-value \< 0.05 AND a logFC
$\ge$ \|2\|**. With that criteria in the liver there are 1805
male-biased genes and 446 female-biased genes.

I have then pulled out all of the male-biased, female-biased, and
non-biased genes based on the criteria outlined above. I determined
non-biased genes as a p-vaule \> 0.05.

``` r
#Removing the rows where padj. is NA in results
FU_liver_con_res_noNA <- FU_liver_con_res[!is.na(FU_liver_con_res$padj),]
summary(FU_liver_con_res_noNA) #We can now see that there are no outliers or low counts since the NAs have been removed
```

    ## 
    ## out of 36595 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 2302, 6.3%
    ## LFC < 0 (down)     : 1289, 3.5%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 0, 0%
    ## (mean count < 6)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
#Creating a vector that contains all of the male-biased and female-biased genes in the liver
FU_liver_mal_biased <- FU_liver_con_res_noNA[which(FU_liver_con_res_noNA$log2FoldChange >= 2 &
                                                     FU_liver_con_res_noNA$padj <= 0.05),]

FU_liver_fem_biased <- FU_liver_con_res_noNA[which(FU_liver_con_res_noNA$log2FoldChange <= -2 &
                                                     FU_liver_con_res_noNA$padj <= 0.05),]

#Creating an object that contains all of the non-biased genes in the liver
FU_liver_non_biased <- FU_liver_con_res_noNA[which(FU_liver_con_res_noNA$padj > 0.05),]
```

I will be generating a table that outlines the genes that are the most
male or female biased in each organ type. To do this I need to pull the
top 50 differentially expressed genes in males or females, get the
Trinity gene IDs and then BLAST the corresponding sequences. Here I am
pulling out the Trinity gene IDs for those genes.

``` r
#Creating a subset of the results where the p-value is less than 0.05
FU_liver_con_res_p05 <- FU_liver_con_res_noNA[FU_liver_con_res_noNA$padj <= 0.05, ]


#Pulling the top 50 diff expressed genes in Males
top50_FU_MLiver <- head(FU_liver_con_res_p05[order(FU_liver_con_res_p05$log2FoldChange, 
                                                   decreasing = TRUE), ], 
                        n = 50)

#Run once to generate the file and then commented out
write.table(cbind(rownames(top50_FU_MLiver)),
            'FU_maleL_top50TRgenes.txt', 
             sep = "", 
            quote=FALSE, 
            row.names = FALSE, 
            col.names = FALSE)


#Pulling the top 50 diff expressed genes in Females
top50_FU_Fliver <- head(FU_liver_con_res_p05[order(FU_liver_con_res_p05$log2FoldChange, 
                                                   decreasing = FALSE), ], 
                        n = 50)

#Run once to generate the file and then commented out
write.table(cbind(rownames(top50_FU_Fliver)),
            'FU_femL_top50TRgenes.txt', 
             sep = "", 
            quote=FALSE, 
            row.names = FALSE, 
            col.names = FALSE)
```

### M-F Gill Comparisson

``` r
##Pulling out the gill M-F results
gill_con_res_FU <- results(ddsMF_FU_exp, 
                           contrast = c("group", "MGill", "FGill"), 
                           alpha = 0.5)
gill_con_res_FU$trin_geneid <- rownames(gill_con_res_FU)
head(gill_con_res_FU)
```

    ## log2 fold change (MLE): group MGill vs FGill 
    ## Wald test p-value: group MGill vs FGill 
    ## DataFrame with 6 rows and 7 columns
    ##                       baseMean log2FoldChange     lfcSE        stat    pvalue
    ##                      <numeric>      <numeric> <numeric>   <numeric> <numeric>
    ## TRINITY_DN0_c0_g1   1259.36897    -0.21026123  0.309346 -0.67969513  0.496698
    ## TRINITY_DN0_c1_g1   4676.25499    -0.00217425  0.265296 -0.00819557  0.993461
    ## TRINITY_DN0_c10_g1   563.87777    -0.04637597  0.333922 -0.13888259  0.889543
    ## TRINITY_DN0_c11_g1    29.98589     0.51752973  0.723823  0.71499451  0.474612
    ## TRINITY_DN0_c125_g1    2.23427    -1.09088174  1.238903 -0.88052246  0.378576
    ## TRINITY_DN0_c13_g1    24.89771     0.26373722  0.701413  0.37600843  0.706911
    ##                          padj         trin_geneid
    ##                     <numeric>         <character>
    ## TRINITY_DN0_c0_g1           1   TRINITY_DN0_c0_g1
    ## TRINITY_DN0_c1_g1           1   TRINITY_DN0_c1_g1
    ## TRINITY_DN0_c10_g1          1  TRINITY_DN0_c10_g1
    ## TRINITY_DN0_c11_g1          1  TRINITY_DN0_c11_g1
    ## TRINITY_DN0_c125_g1        NA TRINITY_DN0_c125_g1
    ## TRINITY_DN0_c13_g1          1  TRINITY_DN0_c13_g1

``` r
summary(gill_con_res_FU)
```

    ## 
    ## out of 157005 with nonzero total read count
    ## adjusted p-value < 0.5
    ## LFC > 0 (up)       : 71, 0.045%
    ## LFC < 0 (down)     : 46, 0.029%
    ## outliers [1]       : 5538, 3.5%
    ## low counts [2]     : 109489, 70%
    ## (mean count < 4)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

In the gills there were 14 genes that we can consider male-biased and 16
female-biased genes based on our criteria for sex-bias.

I have then pulled out all of the male-biased, female-biased, and
non-biased genes to save them into their own objects.

``` r
#Removing the rows where padj. is NA in results
FU_gill_con_res_noNA <- gill_con_res_FU[!is.na(gill_con_res_FU$padj), ]
summary(FU_gill_con_res_noNA) #We can now see that there are no outliers or low counts since the NAs have been removed
```

    ## 
    ## out of 41978 with nonzero total read count
    ## adjusted p-value < 0.5
    ## LFC > 0 (up)       : 71, 0.17%
    ## LFC < 0 (down)     : 46, 0.11%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 0, 0%
    ## (mean count < 4)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
#Creating a vector that contains all of the male-biased and female-biased genes in the gills
FU_gill_mal_biased <- FU_gill_con_res_noNA[which(FU_gill_con_res_noNA$log2FoldChange >= 2 
                                                 & FU_gill_con_res_noNA$padj <= 0.05),]
FU_gill_fem_biased <- FU_gill_con_res_noNA[which(FU_gill_con_res_noNA$log2FoldChange <= -2 
                                                 & FU_gill_con_res_noNA$padj <= 0.05),]

#Creating an object that contains all of the non-biased genes in the gills, p>0.05
FU_gill_non_biased <- FU_gill_con_res_noNA[which(FU_gill_con_res_noNA$padj > 0.05),]
```

Here I am getting the trinity gene IDs that correspond to the top 50
differentially expressed genes for males and females.

``` r
#Creating a subset of results where p-value is less than 0.05
FU_gill_con_res_p05 <- FU_gill_con_res_noNA[FU_gill_con_res_noNA$padj <= 0.05,]

#Pulling the top 50 diff expressed genes in Males
top50_FU_MGill <- head(FU_gill_con_res_p05[order(FU_gill_con_res_p05$log2FoldChange, 
                                                 decreasing = TRUE),], 
                       n=50)

#Run once to generate the file and then commented out
write.table(cbind(rownames(top50_FU_MGill)),
            'FU_maleG_top50TRgenes.txt', 
            sep = "", 
            quote=FALSE, 
            row.names = FALSE, 
            col.names = FALSE)

#Pulling the top 50 diff expressed genes in Females
top50_FU_FGill <- head(FU_gill_con_res_p05[order(FU_gill_con_res_p05$log2FoldChange, 
                                                 decreasing = FALSE),], 
                       n=50)

#Run once to generate the file and then commented out
write.table(cbind(rownames(top50_FU_FGill)),
            'FU_femG_top50TRgenes.txt', 
            sep = "", 
            quote=FALSE, 
            row.names = FALSE, 
            col.names = FALSE)
```

### M-F Gonad Comparisson

``` r
##Pulling out the gonad M-F results
gonad_con_res_FU <- results(ddsMF_FU_exp, 
                            contrast = c("group", "MGonad", "FGonad"), 
                            alpha = 0.05)
gonad_con_res_FU$trin_geneid <- rownames(gonad_con_res_FU)
head(gonad_con_res_FU)
```

    ## log2 fold change (MLE): group MGonad vs FGonad 
    ## Wald test p-value: group MGonad vs FGonad 
    ## DataFrame with 6 rows and 7 columns
    ##                       baseMean log2FoldChange     lfcSE      stat      pvalue
    ##                      <numeric>      <numeric> <numeric> <numeric>   <numeric>
    ## TRINITY_DN0_c0_g1   1259.36897     -1.0427427  0.309652 -3.367466 7.58623e-04
    ## TRINITY_DN0_c1_g1   4676.25499     -1.3428381  0.264994 -5.067420 4.03243e-07
    ## TRINITY_DN0_c10_g1   563.87777     -0.0888395  0.333925 -0.266046 7.90204e-01
    ## TRINITY_DN0_c11_g1    29.98589      1.7127089  0.745281  2.298071 2.15578e-02
    ## TRINITY_DN0_c125_g1    2.23427      3.7538363  1.508099  2.489118 1.28060e-02
    ## TRINITY_DN0_c13_g1    24.89771     -1.9197578  0.710697 -2.701232 6.90831e-03
    ##                            padj         trin_geneid
    ##                       <numeric>         <character>
    ## TRINITY_DN0_c0_g1   3.87893e-03   TRINITY_DN0_c0_g1
    ## TRINITY_DN0_c1_g1   4.52179e-06   TRINITY_DN0_c1_g1
    ## TRINITY_DN0_c10_g1  9.67759e-01  TRINITY_DN0_c10_g1
    ## TRINITY_DN0_c11_g1  6.37398e-02  TRINITY_DN0_c11_g1
    ## TRINITY_DN0_c125_g1 4.19374e-02 TRINITY_DN0_c125_g1
    ## TRINITY_DN0_c13_g1  2.52487e-02  TRINITY_DN0_c13_g1

``` r
summary(gonad_con_res_FU)
```

    ## 
    ## out of 157005 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 12259, 7.8%
    ## LFC < 0 (down)     : 7309, 4.7%
    ## outliers [1]       : 5538, 3.5%
    ## low counts [2]     : 90084, 57%
    ## (mean count < 2)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

Between the ovaries and testis (gonads) there were 8886 genes that we
can consider male-biased and 3884 female-biased genes based on our
criteria for sex-bias.

From the gonad M-F results I have then pulled out all of the
male-biased, female-biased, and non-biased genes and saved them to their
own objects.

``` r
#Removing the rows where padj. is NA in results
FU_gonad_con_res_noNA <- gonad_con_res_FU[!is.na(gonad_con_res_FU$padj), ]
summary(FU_gonad_con_res_noNA) #We can now see that there are no outliers or low counts since the NAs have been removed
```

    ## 
    ## out of 61383 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 12259, 20%
    ## LFC < 0 (down)     : 7309, 12%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 0, 0%
    ## (mean count < 2)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
#Creating an object that contains all of the male-biased and female-biased genes in the gonads
FU_gonad_mal_biased <- FU_gonad_con_res_noNA[which(FU_gonad_con_res_noNA$log2FoldChange >= 2 
                                                   & FU_gonad_con_res_noNA$padj <= 0.05),]
FU_gonad_fem_biased <- FU_gonad_con_res_noNA[which(FU_gonad_con_res_noNA$log2FoldChange <= -2 
                                                   & FU_gonad_con_res_noNA$padj <= 0.05),]

#Creating an object that contains all of the non-biased genes in the gonads, p>0.05
FU_gonad_non_biased <- FU_gonad_con_res_noNA[which(FU_gonad_con_res_noNA$padj > 0.05),]
```

I then pulled out the trinity gene IDs that corresponded to the top 50
differentially expressed genes in males and females so that I could
BLAST the sequences.

``` r
#Creating a subset of the results where the p-value is less than 0.05
FU_gonad_con_res_p05 <- FU_gonad_con_res_noNA[FU_gonad_con_res_noNA$padj <= 0.05,]

#Pulling the top 50 diff expressed genes in Males
top50_FU_MGonad <- head(FU_gonad_con_res_p05[order(FU_gonad_con_res_p05$log2FoldChange, 
                                                   decreasing = TRUE),], 
                        n=50)

#Run once to generate the file and then commented out
write.table(cbind(rownames(top50_FU_MGonad)),
            'FU_maleGon_top50TRgenes.txt', 
            sep = "", 
            quote=FALSE, 
            row.names = FALSE, 
            col.names = FALSE)

#Pulling the top 50 diff expressed genes in Females
top50_FU_FGonad <- head(FU_gonad_con_res_p05[order(FU_gonad_con_res_p05$log2FoldChange, 
                                                   decreasing = FALSE),], 
                        n=50)

#Run once to generate the file and then commented out
write.table(cbind(head(top50_FU_FGonad)),
            'FU_femGon_top50TRgenes.txt', 
            sep = "", 
            quote=FALSE, 
            row.names = FALSE, 
            col.names = FALSE)
```

### MA Plots

<figure>
<img src="fuscus_diff_expr_analysis_files/figure-gfm/MA-plot-MF-1.png"
alt="MA-plots generated for each organ type that compares logFC to the mean expression. Female-biased and Male-biased genes are represented by color (green and purple correspondingly) and are determing with a fc cutoff of 2 and a p-value cutoff of 0.05." />
<figcaption aria-hidden="true">MA-plots generated for each organ type
that compares logFC to the mean expression. Female-biased and
Male-biased genes are represented by color (green and purple
correspondingly) and are determing with a fc cutoff of 2 and a p-value
cutoff of 0.05.</figcaption>
</figure>

## Variation in FC across sex-bias and tissue type

I want to create a plot the highlights the variation we see in
fold-change both across the different biases groups (male-biased,
female-biased, and non-biased) within one tissue type and across all of
the tissue types. To do this I first need to create a “long” dataset in
the following format:

| Tissue_type |   Bias    | Fold-Change |
|:-----------:|:---------:|:-----------:|
|    Gill     | male_bias |    XXXXX    |
|    Gill     | male_bias |    XXXXX    |
|    Gill     | fem_bias  |    XXXXX    |
|    Gill     |  no_bias  |    XXXXX    |
|      …      |    ….     |      …      |
|    Liver    | male_bias |    XXXXX    |
|    Liver    | male_bias |    XXXXX    |
|    Liver    | fem_bias  |    XXXXX    |
|    Liver    |  no_bias  |    XXXXX    |
|      …      |    ….     |      …      |
|    Gonad    | male_bias |    XXXXX    |
|    Gonad    | male_bias |    XXXXX    |
|    Gonad    | fem_bias  |    XXXXX    |
|    Gonad    |  no_bias  |    XXXXX    |
|      …      |    ….     |      …      |

That dataset is created here, with the Trinity gene IDs included as
well:

``` r
FU_logFC_long <- data.frame(
  tissue=c(rep("Gill",nrow(FU_gill_fem_biased)),
           rep("Gill", nrow(FU_gill_mal_biased)),
           rep("Gill", nrow(FU_gill_non_biased)),
           rep("Gonad", nrow(FU_gonad_fem_biased)),
           rep("Gonad", nrow(FU_gonad_mal_biased)),
           rep("Gonad", nrow(FU_gonad_non_biased)),
           rep("Liver", nrow(FU_liver_fem_biased)),
           rep("Liver", nrow(FU_liver_mal_biased)),
           rep("Liver", nrow(FU_liver_non_biased))
         ),
  bias=c(rep("FB",nrow(FU_gill_fem_biased)),
         rep("MB",nrow(FU_gill_mal_biased)),
         rep("NB", nrow(FU_gill_non_biased)),
         rep("FB", nrow(FU_gonad_fem_biased)),
         rep("MB", nrow(FU_gonad_mal_biased)),
         rep("NB", nrow(FU_gonad_non_biased)),
         rep("FB", nrow(FU_liver_fem_biased)),
         rep("MB", nrow(FU_liver_mal_biased)),
         rep("NB", nrow(FU_liver_non_biased))
         ),
  logFC=c(FU_gill_fem_biased$log2FoldChange,
          FU_gill_mal_biased$log2FoldChange,
          FU_gill_non_biased$log2FoldChange,
          FU_gonad_fem_biased$log2FoldChange,
          FU_gonad_mal_biased$log2FoldChange,
          FU_gonad_non_biased$log2FoldChange,
          FU_liver_fem_biased$log2FoldChange,
          FU_liver_mal_biased$log2FoldChange,
          FU_liver_non_biased$log2FoldChange
          ),
  geneID=c(rownames(FU_gill_fem_biased),
           rownames(FU_gill_mal_biased),
           rownames(FU_gill_non_biased),
           rownames(FU_gonad_fem_biased),
           rownames(FU_gonad_mal_biased),
           rownames(FU_gonad_non_biased),
           rownames(FU_liver_fem_biased),
           rownames(FU_liver_mal_biased),
           rownames(FU_liver_non_biased)
           )
  
)
```

With the dataset in the proper format now, we can generate the plot

    ## png 
    ##   2

There are some weird patterns popping up in this analysis. The first one
is the magnitude of differences in male- and female-biased genes in the
liver. This weirdness is then confirmed in the above plot, we can see
this clumping pattern showing um in the male-biased genes in the liver.
There is a cluster of genes with a logFC between 5 -15 and then another
clump at 25 and above, and very little in between. This may be a result
of normalizing across all of the organs. I am going to try re-running
the DE analysis within each organ and see if that helps anything.

# Running a single factor analysis between males and females WITHIN each organ

## Liver M-F differential expression

``` r
#Create a subset of the DESeq dataset that only includes the liver samples
dds_FU_liver <- ddsMF_FU[, ddsMF_FU$Organ == "Liver"]

#Filtering to remove any rows that have less than 10 reads total
keep <- rowSums(counts(dds_FU_liver)) >= 10
dds_FU_liver <- dds_FU_liver[keep, ]

design(dds_FU_liver) <- ~ Sex

#Run the differential expression analysis
dds_FU_liver_exp <- DESeq(dds_FU_liver)

#Compiling the results
res_liver <- results(dds_FU_liver_exp, alpha = 0.05)

summary(res_liver)
```

    ## 
    ## out of 57172 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 139, 0.24%
    ## LFC < 0 (down)     : 65, 0.11%
    ## outliers [1]       : 6294, 11%
    ## low counts [2]     : 23989, 42%
    ## (mean count < 4)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

The single factor analysis from the liver showed 130 male-biased genes
and 58 female-biased genes. I am going to pull out all of the female-
and male-biased genes to store them in their own objects as well as
genes that are non-biased. Non-biased genes will be categorized as genes
with a p-value \> 0.05, regardless of what the logFC may be.

``` r
#Removing the rows where padj. is NA in results
liver_res_noNA <- res_liver[!is.na(res_liver$padj),]
summary(liver_res_noNA) #We can now see that there are no outliers or 
```

    ## 
    ## out of 26889 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 139, 0.52%
    ## LFC < 0 (down)     : 65, 0.24%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 0, 0%
    ## (mean count < 4)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
                            #low counts since the NAs have been removed
#write.csv(as.data.frame(liver_res_noNA), "data/liver_res.csv", row.names = TRUE)

#Creating a vector that contains all of the male-biased and female-biased genes in the brain
liver_mal_biased <- liver_res_noNA[which(liver_res_noNA$log2FoldChange >= 2
                                             & liver_res_noNA$padj <= 0.05),]
liver_fem_biased <- liver_res_noNA[which(liver_res_noNA$log2FoldChange <= -2 
                                             & liver_res_noNA$padj <= 0.05),]

#Creating an object that contains all of the non-biased genes in the brain
liver_non_biased <- liver_res_noNA[which(liver_res_noNA$padj > 0.05),]
```

## Gonad M-F differential expression

``` r
#Create a subset of the DESeq dataset that only includes the gonad samples
dds_FU_gonad <- ddsMF_FU[, ddsMF_FU$Organ == "Gonad"]

#Filtering to remove any rows that have less than 10 reads total
keep <- rowSums(counts(dds_FU_gonad)) >= 10
dds_FU_gonad <- dds_FU_gonad[keep, ]

design(dds_FU_gonad) <- ~ Sex

#Run the differential expression analysis
dds_FU_gonad_exp <- DESeq(dds_FU_gonad)

#Compiling the results
res_gonad <- results(dds_FU_gonad_exp, alpha = 0.05)

summary(res_gonad)
```

    ## 
    ## out of 88396 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 15718, 18%
    ## LFC < 0 (down)     : 10814, 12%
    ## outliers [1]       : 3091, 3.5%
    ## low counts [2]     : 23993, 27%
    ## (mean count < 2)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

The single factor analysis from the gonad showed 11557 male-biased genes
and 6126 female-biased genes. I am going to pull out all of the female-
and male-biased genes to store them in their own objects as well as
genes that are non-biased. Non-biased genes will be categorized as genes
with a p-value \> 0.05, regardless of what the logFC may be.

``` r
#Removing the rows where padj. is NA in results
gonad_res_noNA <- res_gonad[!is.na(res_gonad$padj),]
summary(gonad_res_noNA) #We can now see that there are no outliers or 
```

    ## 
    ## out of 61312 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 15718, 26%
    ## LFC < 0 (down)     : 10814, 18%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 0, 0%
    ## (mean count < 2)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
                            #low counts since the NAs have been removed
#write.csv(as.data.frame(gonad_res_noNA), "data/gonad_res.csv", row.names = TRUE)

#Creating a vector that contains all of the male-biased and female-biased genes in the brain
gonad_mal_biased <- gonad_res_noNA[which(gonad_res_noNA$log2FoldChange >= 2
                                             & gonad_res_noNA$padj <= 0.05),]
gonad_fem_biased <- gonad_res_noNA[which(gonad_res_noNA$log2FoldChange <= -2 
                                             & gonad_res_noNA$padj <= 0.05),]

#Creating an object that contains all of the non-biased genes in the brain
gonad_non_biased <- gonad_res_noNA[which(gonad_res_noNA$padj > 0.05),]
```

## Gill M-F differential expression

``` r
#Create a subset of the DESeq dataset that only includes the gill samples
dds_FU_gill <- ddsMF_FU[, ddsMF_FU$Organ == "Gill"]

#Filtering to remove any rows that have less than 10 reads total
keep <- rowSums(counts(dds_FU_gill)) >= 10
dds_FU_gill <- dds_FU_gill[keep, ]

design(dds_FU_gill) <- ~ Sex

#Run the differential expression analysis
dds_FU_gill_exp <- DESeq(dds_FU_gill)

#Compiling the results
res_gill <- results(dds_FU_gill_exp, alpha = 0.05)

summary(res_gill)
```

    ## 
    ## out of 92188 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 116, 0.13%
    ## LFC < 0 (down)     : 53, 0.057%
    ## outliers [1]       : 1278, 1.4%
    ## low counts [2]     : 47927, 52%
    ## (mean count < 5)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

The single factor analysis from the gill showed 40 male-biased genes and
29 female-biased genes. I am going to pull out all of the female- and
male-biased genes to store them in their own objects as well as genes
that are non-biased. Non-biased genes will be categorized as genes with
a p-value \> 0.05, regardless of what the logFC may be.

``` r
#Removing the rows where padj. is NA in results
gill_res_noNA <- res_gill[!is.na(res_gill$padj),]
summary(gill_res_noNA) #We can now see that there are no outliers or 
```

    ## 
    ## out of 42983 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 116, 0.27%
    ## LFC < 0 (down)     : 53, 0.12%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 0, 0%
    ## (mean count < 5)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
                            #low counts since the NAs have been removed
#write.csv(as.data.frame(gill_res_noNA), "data/gill_res.csv", row.names = TRUE)

#Creating a vector that contains all of the male-biased and female-biased genes in the brain
gill_mal_biased <- gill_res_noNA[which(gill_res_noNA$log2FoldChange >= 2
                                             & gill_res_noNA$padj <= 0.05),]
gill_fem_biased <- gill_res_noNA[which(gill_res_noNA$log2FoldChange <= -2 
                                             & gill_res_noNA$padj <= 0.05),]

#Creating an object that contains all of the non-biased genes in the brain
gill_non_biased <- gill_res_noNA[which(gill_res_noNA$padj > 0.05),]
```

## Re-looking at variation in logFC across both sex-bias categories and tissue types

I want to now re-create the plot the highlights the variation we see in
fold-change both across the different biases groups (male-biased,
female-biased, and non-biased) within one tissue type and across all of
the tissue types to see if the weird patterns are still there. To do
this I first need to re-generate the “long” style dataset with the new
information:

``` r
logFC_long <- data.frame(
  tissue=c(rep("Gill",nrow(gill_fem_biased)),
           rep("Gill", nrow(gill_mal_biased)),
           rep("Gill", nrow(gill_non_biased)),
           rep("Gonad", nrow(gonad_fem_biased)),
           rep("Gonad", nrow(gonad_mal_biased)),
           rep("Gonad", nrow(gonad_non_biased)),
           rep("Liver", nrow(liver_fem_biased)),
           rep("Liver", nrow(liver_mal_biased)),
           rep("Liver", nrow(liver_non_biased))
         ),
  bias=c(rep("FB",nrow(gill_fem_biased)),
         rep("MB",nrow(gill_mal_biased)),
         rep("NB", nrow(gill_non_biased)),
         rep("FB", nrow(gonad_fem_biased)),
         rep("MB", nrow(gonad_mal_biased)),
         rep("NB", nrow(gonad_non_biased)),
         rep("FB", nrow(liver_fem_biased)),
         rep("MB", nrow(liver_mal_biased)),
         rep("NB", nrow(liver_non_biased))
         ),
  logFC=c(gill_fem_biased$log2FoldChange,
          gill_mal_biased$log2FoldChange,
          gill_non_biased$log2FoldChange,
          gonad_fem_biased$log2FoldChange,
          gonad_mal_biased$log2FoldChange,
          gonad_non_biased$log2FoldChange,
          liver_fem_biased$log2FoldChange,
          liver_mal_biased$log2FoldChange,
          liver_non_biased$log2FoldChange
          ),
  geneID=c(rownames(gill_fem_biased),
           rownames(gill_mal_biased),
           rownames(gill_non_biased),
           rownames(gonad_fem_biased),
           rownames(gonad_mal_biased),
           rownames(gonad_non_biased),
           rownames(liver_fem_biased),
           rownames(liver_mal_biased),
           rownames(liver_non_biased)
           )
  
)

#Exporting this data, run once and then commented out
#write.csv(logFC_long, "data/logFC_long_sexbias.csv", row.names = FALSE)
```

After re-creating the long version of the dataset I can re-generate the
plot.

    ## png 
    ##   2

The weird pattern does appear to still be there for the male-biased
genes in the liver. I am going to look more in detail to the expression
level of the male- and female-biased genes within each organ across all
of the tissues using heatmaps. What does change, however, is the degree
of the logFC. When doing the multi-facotr analysis, there were several
genes that had a logFC of 35. Now, the highest logFC we see is around
25.

## Investigating expression across samples in the sex-biased genes

For each organ type I am going to pull out either all of the sex-biased
genes, or just the top 200 if there are more than that, and look at the
expression levels of those genes across all of the samples.

![](fuscus_diff_expr_analysis_files/figure-gfm/heatmap_gonads-1.png)<!-- -->![](fuscus_diff_expr_analysis_files/figure-gfm/heatmap_gonads-2.png)<!-- -->

![](fuscus_diff_expr_analysis_files/figure-gfm/heatmap_liver-1.png)<!-- -->![](fuscus_diff_expr_analysis_files/figure-gfm/heatmap_liver-2.png)<!-- -->

![](fuscus_diff_expr_analysis_files/figure-gfm/heatmap_gill-1.png)<!-- -->![](fuscus_diff_expr_analysis_files/figure-gfm/heatmap_gill-2.png)<!-- -->

In the liver, for both males and females, expression in the top biased
genes appears to be driven almost exclusively by two individuals. I am
going to continue on with the SF analysis over the MF analysis for a
more conservative approach. I am also choosing to leave those
individuals in the analysis as they may represent true variation in
expression of genes across the samples.

Looking back at the logFC plot made from the single-factor analysis,
there appears to be differences in the logFC across both the organs and
the sex-bias groups (male-biased, fem-biased, and non-biased). In the
gill and liver, it looks like the female-biased genes are more biased on
average than the male-biased genes, but this pattern reverses in the
gonads. In every case the non-biased genes have the lowest median logFC
even though there is some variation.

I am going to fit a linear model to the data and explore some of the
summary statistics to see what is happening

``` r
#Looking at some summary statistics for logFC between the different groups
tapply(abs(logFC_long$logFC), list(logFC_long$bias, logFC_long$tissue), median)
```

    ##         Gill    Gonad     Liver
    ## FB 5.6611821 3.860896 7.1918331
    ## MB 5.4242997 4.097896 5.2275218
    ## NB 0.3223442 1.156665 0.5509466

``` r
tapply(abs(logFC_long$logFC), list(logFC_long$bias, logFC_long$tissue), sd)
```

    ##         Gill    Gonad    Liver
    ## FB 4.0337608 2.211555 7.289603
    ## MB 2.1559514 1.605536 9.545186
    ## NB 0.6580531 1.328525 1.293784

``` r
model <- lm(abs(logFC_long$logFC)~ logFC_long$tissue * logFC_long$bias)
summary(model)
```

    ## 
    ## Call:
    ## lm(formula = abs(logFC_long$logFC) ~ logFC_long$tissue * logFC_long$bias)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -9.6366 -0.6204 -0.2737  0.3451 20.8980 
    ## 
    ## Coefficients:
    ##                                          Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)                                5.8766     0.2374  24.753  < 2e-16
    ## logFC_long$tissueGonad                    -1.5441     0.2380  -6.488 8.70e-11
    ## logFC_long$tissueLiver                     3.8900     0.2908  13.378  < 2e-16
    ## logFC_long$biasMB                         -0.5912     0.3118  -1.896    0.058
    ## logFC_long$biasNB                         -5.3421     0.2375 -22.494  < 2e-16
    ## logFC_long$tissueGonad:logFC_long$biasMB   0.4855     0.3125   1.554    0.120
    ## logFC_long$tissueLiver:logFC_long$biasMB   2.4880     0.3715   6.698 2.12e-11
    ## logFC_long$tissueGonad:logFC_long$biasNB   2.5539     0.2381  10.724  < 2e-16
    ## logFC_long$tissueLiver:logFC_long$biasNB  -3.4922     0.2909 -12.003  < 2e-16
    ##                                             
    ## (Intercept)                              ***
    ## logFC_long$tissueGonad                   ***
    ## logFC_long$tissueLiver                   ***
    ## logFC_long$biasMB                        .  
    ## logFC_long$biasNB                        ***
    ## logFC_long$tissueGonad:logFC_long$biasMB    
    ## logFC_long$tissueLiver:logFC_long$biasMB ***
    ## logFC_long$tissueGonad:logFC_long$biasNB ***
    ## logFC_long$tissueLiver:logFC_long$biasNB ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.278 on 122210 degrees of freedom
    ## Multiple R-squared:  0.5029, Adjusted R-squared:  0.5028 
    ## F-statistic: 1.545e+04 on 8 and 122210 DF,  p-value: < 2.2e-16

The summary statistics confirm the assumptions made from the boxplots,
female-biased genes in the liver and gill have a higher logFC compared
to the male-biased genes, but male-biased genes in the gonads have a
higher logFC compared to the female-biased. In all three organs,
unbiased genes possess the lowest logFC. Interestingly, the highest
logFC across tissues can be found in the liver, followed by the gill,
and lastly the gonads.

## Looking at shared sex-biased genes with an Upset plot

I now want to see if there are genes that are female/male-biased in
multiple organs or if there are genes that are female-biased in one
organ and then male-biased in a different organ, etc. In order to do
this I will be creating upset plots first for the sexes individually and
then across all biased genes.

<figure>
<img
src="fuscus_diff_expr_analysis_files/figure-gfm/upset-plot-mal-1.png"
alt="Upset plots to show the number of shared sex-biased genes across the organs in males." />
<figcaption aria-hidden="true">Upset plots to show the number of shared
sex-biased genes across the organs in males.</figcaption>
</figure>

<figure>
<img
src="fuscus_diff_expr_analysis_files/figure-gfm/upset-plot-fem-1.png"
alt="Upset plots to show the number of shared sex-biased genes across the organs in females." />
<figcaption aria-hidden="true">Upset plots to show the number of shared
sex-biased genes across the organs in females.</figcaption>
</figure>

<figure>
<img
src="fuscus_diff_expr_analysis_files/figure-gfm/upset-plot-all-1.png"
alt="Upset plots to show the number of shared sex-biased genes across the organs across both sexes." />
<figcaption aria-hidden="true">Upset plots to show the number of shared
sex-biased genes across the organs across both sexes.</figcaption>
</figure>

In females, the highest number of shared sex-biased genes is found
between the gonads and the liver, but in males, the highest number of
shared sex-biased genes is found between the gills and the gonads. There
are no genes that are shared across the same organ between males and
females (i.e. biased in males and females in the liver) which is good.
There are 3 and 2 genes that are biased in all three organs for males
and females respectively. Overall, however, there is not a high amount
of genes that are sex-biased in more than one organ type.

## Categorizing sex-specific genes

Now that we have investigated what the sex-bias is like across the
different tissue types I want to dive further into genes with
**sex-specific expression**. These are genes expressed only in one sex
or the other. I am classifying a sex-specific genes within each tissue
as ones where the expression of that gene is less than 10 for all of one
sex and there is a median of $\ge 20$ in the other sex.

Prior to running the for loop, any “outliers” that were determined by
DESeq were removed from the pool of genes. Additionally, normalized
counts were used for the classification of sex-specific genes as DESeq2
uses the moralized counts for all of the logFC/sexbias calculations.

``` r
#Store all of the single factor analysis outputs in a list
DE_models <- list(dds_FU_gill_exp, 
                  dds_FU_gonad_exp, 
                  dds_FU_liver_exp)

#Create an empty list to store my sex-specific datasets in
FU_sex_specific_genes <- list()

#For each DESeq dataset in DE_models pull out the Mspecific and Fspecific
#genes based on medians
for (dataset in DE_models) {
  
  #Pull out the geneIDs for genes that were categorized as "outliers" by DESeq2
  ##Calculate the Cooks threshold that would have been used
  np <- length(resultsNames(dataset))
  nsamp <- ncol(dataset)
  cooks_thresh <- qf(0.99, df1 = np, df2 = nsamp - np)
  
  ##Apply threshold calculated above to the cooks values in the dataset
  out_ids <- names(mcols(dataset)$maxCooks[mcols(dataset)$maxCooks >
                                                cooks_thresh])
  
  #Filtering the dds dataset to remove the outliers identified by DESeq
  dataset_filtered <- dataset[!(rownames(dataset) %in% out_ids), ]
  
  
  #Male-specific Genes
  ##Pull out all of the rows where fem count <=10 in every female sample
  fem10_organ_names <- which(rowSums(t(apply(counts(dataset_filtered,
                                                    normalized = TRUE)[, dataset_filtered$Sex == "F"],
                                             1,
                                             function(x) x <= 10)
                                       )
                                     ) == ncol(counts(dataset_filtered, 
                                                      normalized = TRUE)[, dataset_filtered$Sex == "F"])
                             )
  
  fem10_organ <- counts(dataset_filtered,
                        normalized = TRUE)[rownames(counts(dataset_filtered,
                                                          normalized = TRUE)) %in%
                                            names(fem10_organ_names),]
  
  ##Pull out the rows where median of male count >=50  
  mal50_organ <- apply(counts(dataset_filtered,
                              normalized = TRUE)[, dataset_filtered$Sex == "M"],
                       1,
                       function(x) median(x) >= 50)
  
  ##Keep only rows where all female samples <=10 and the male median >= 50
  fem10_mal50_organ <- fem10_organ[rownames(fem10_organ) %in%
                                     names(mal50_organ[mal50_organ == TRUE]),
                                   ]
  
  ##Create a new object with a name based on the organ type
  organ_malsp <- sub("$", "_male_specific", dataset$Organ[1])
  FU_sex_specific_genes[[organ_malsp]] <- fem10_mal50_organ
  
  
  #Female-specific Genes
  ##Pull out ll the rows where male count is <=10 in every male sample
  mal10_organ_names <- which(rowSums(t(apply(counts(dataset_filtered,
                                                    normalized = TRUE)[, dataset_filtered$Sex == "M"],
                                             1,
                                             function(x) x <= 10)
                                       )
                                     ) == ncol(counts(dataset_filtered,
                                                      normalized = TRUE)[, dataset_filtered$Sex == "M"])
                             )
  
  mal10_organ <- counts(dataset_filtered,
                        normalized = TRUE)[rownames(counts(dataset_filtered,
                                                           normalized = TRUE)) %in%
                                             names(mal10_organ_names), ]
  
  ##Pull out rows where median of female count is >=50 
  fem50_organ <- apply(counts(dataset_filtered,
                              normalized = TRUE)[, dataset_filtered$Sex == "F"],
                       1,
                       function(x) median(x) >= 50)
  
  ##keep only the rows where male <=10 and fem median >= 50
  mal10_fem50_organ <- mal10_organ[rownames(mal10_organ) %in%
                                     names(fem50_organ[fem50_organ == TRUE]),
                                   ]
  
  ##Create a new object with a name based on the organ type
  organ_femsp <- sub("$", "_female_specific", dataset$Organ[1])
  FU_sex_specific_genes[[organ_femsp]] <- mal10_fem50_organ
  
}
```

### Investigating the results of the sex-specific subsetting

Let’s now take a look at how many sex-specific genes we have in each
tissue type:

| tissue | sex | num_genes |
|:-------|:----|----------:|
| Gill   | M   |        12 |
| Gonad  | M   |       624 |
| Liver  | M   |        18 |
| Gill   | F   |         9 |
| Gonad  | F   |       447 |
| Liver  | F   |        14 |

To get a better idea what is going on between the sex-specific genes vs
the sex-biased genes I first determined how many overlaps there were
between our sex-biased and sex-specific genes and then used plotCounts
to plot the counts for each gene separately.

In the **Male Gills** there are 2 genes that overlap between sex-biased
and sex-specific, there are 570 overlapping in the **gonads** and 18
overlapping in the **liver**. For females we have 3 overlapping
sex-specific and sex-biased genes in the **gills**, 412 genes in the
**gonads**, and 13 genes in the **liver**.

### Pulling out the gene IDs for all of the sex-specific genes

``` r
#Run once to create the file and then commented out
write.table(cbind(rownames(FU_sex_specific_genes$Gill_male_specific)),
            'data/FU_malG_specific_TRgenes.txt', 
            sep = "", 
            quote=FALSE, 
            row.names = FALSE, 
            col.names = FALSE)
write.table(cbind(rownames(FU_sex_specific_genes$Gill_female_specific)),
            'data/FU_femG_specific_TRgenes.txt', 
            sep = "", 
            quote=FALSE, 
            row.names = FALSE, 
            col.names = FALSE)
write.table(cbind(rownames(FU_sex_specific_genes$Gonad_male_specific)),
            'data/FU_malGon_specific_TRgenes.txt', 
            sep = "", 
            quote=FALSE, 
            row.names = FALSE, 
            col.names = FALSE)
write.table(cbind(rownames(FU_sex_specific_genes$Gonad_female_specific)),
            'data/FU_femGon_specific_TRgenes.txt', 
            sep = "", 
            quote=FALSE, 
            row.names = FALSE, 
            col.names = FALSE)
write.table(cbind(rownames(FU_sex_specific_genes$Liver_male_specific)),
            'data/FU_malL_specific_TRgenes.txt', 
            sep = "", 
            quote=FALSE, 
            row.names = FALSE, 
            col.names = FALSE)
write.table(cbind(rownames(FU_sex_specific_genes$Liver_female_specific)),
            'data/FU_femL_specific_TRgenes.txt', 
            sep = "", 
            quote=FALSE, 
            row.names = FALSE, 
            col.names = FALSE)
```

## Creating categories and binning the sex-biased genes based on degree of logFC

This binned was only done on the sex-biased genes, will not have a
category for the unbiased genes. The cutoffs for the different groups
are as follows:

    1. Low - biased = LFC 2 - 3
    2. Medium - biased = LFC 3 - 5
    3. High - biased = LFC 5 - 10
    4. Extreme sex bias = LFC > 10

Counts of the sex-specific genes will be added on separately as they
were not classified based on fold-change but rather a presence in one
sex and an absence in the other sex. A for loop will then be used to
make sure genes are not counted twice (i.e. as both extremely biased and
sex-specific).

``` r
#Make a vector that contains all of the groupings
biased_bins <- c("Unbiased", "Low", "Med", "High", "Extreme", "Sex-specific")

#Create a new column in the dataset and use ifelse statements to set the category limits
#abs(logFC) was used to account for the fem-biased genes possessing negative values
logFC_long$bias_cat <- ifelse(logFC_long$bias == "NB",
                              biased_bins[1],
                              ifelse(abs(logFC_long$logFC) >= 2 & abs(logFC_long$logFC) < 3,
                                     biased_bins[2],
                                     ifelse(abs(logFC_long$logFC) >= 3 & abs(logFC_long$logFC) < 5,
                                            biased_bins[3],
                                            ifelse(abs(logFC_long$logFC) >= 5 & abs(logFC_long$logFC) < 10,
                                                   biased_bins[4],
                                                   biased_bins[5])
                                            )
                                     )
                              )

#Making sure the genes we categorized as sex-specific are labeled as sex-specific for their 
#bias cat. in the dataset
organs <- c("Gill", "Gill", "Gonad", "Gonad", "Liver", "Liver")
bias <- c("MB", "FB", "MB", "FB", "MB", "FB")

for(i in 1:length(FU_sex_specific_genes)){

  tmp <- FU_sex_specific_genes[[i]]
  tmp <- as.data.frame(tmp)
  tmp$geneID <- rownames(tmp)
  
  for(j in 1:nrow(tmp)){
    
    one_gene <- tmp[j, ]
    
    if(one_gene[["geneID"]] %in%
       logFC_long[logFC_long$tissue == organs[[i]] &
                  logFC_long$bias == bias[[i]],"geneID"]){
       
      logFC_long[logFC_long$geneID == one_gene[["geneID"]] &
                   logFC_long$tissue == organs[[i]] &
                   logFC_long$bias == bias[[i]],
                 "bias_cat"] <- "Sex-specific"
    }else{
      
      one_gene_dat <- data.frame(matrix(ncol= ncol(logFC_long),
                                        nrow=1))
      colnames(one_gene_dat) <- colnames(logFC_long)
      
      one_gene_dat$tissue <- organs[[i]]
      one_gene_dat$geneID <- one_gene[["geneID"]]
      one_gene_dat$bias <- bias[[i]]
      one_gene_dat$bias_cat <- "Sex-specific"
      rownames(one_gene_dat) <- NULL
      
      logFC_long <- rbind(one_gene_dat, logFC_long)
      
      rownames(logFC_long) <- NULL
    }
  }
 }

#Create  subset of our long dataset that does not include the non-biased genes
logFC_long_noNB <- logFC_long[logFC_long$bias_cat != "Unbiased",]

#Make a table to count the number of genes in each category for each organ
bias_cat_table <- table(logFC_long_noNB$bias, logFC_long_noNB$bias_cat, logFC_long_noNB$tissue)

#Seperate out the different organs
bias_cat_gill <- bias_cat_table[,, "Gill"]
#write.table(bias_cat_gill, "data/bias_cat_gills.txt", row.names = TRUE)

bias_cat_gonad <- bias_cat_table[,, "Gonad"]
#write.table(bias_cat_gonad, "data/bias_cat_gonad.txt", row.names = TRUE)

bias_cat_liver <- bias_cat_table[,, "Liver"]
write.table(bias_cat_liver, "data/bias_cat_liver.txt", row.names = TRUE)

#Plot the counts for each tissue type
sex_cols <- c("F" = "#7fc97f", "M" = "#beaed4" )
labs <- c("Low", "Med", "High", "Extreme", "Specific")

pdf("docs/figs/FigSB_biasCat_counts.pdf",width = 10, height=4)

par(mfrow=c(1, 3), oma=c(6,4,2,8), mar=c(1,2.5,1,0), 
    cex.main=2,
    cex.axis=2)

bp <- barplot(bias_cat_gill[,biased_bins[biased_bins!="Unbiased"]], 
        beside = TRUE,
        xaxt='n',
        ylim = c(0, max(bias_cat_gill)+5), 
        col = sex_cols,
        main = "")
mtext("Gill",3,outer = FALSE,cex=1.5,line=-1)
axis(2,lwd=2)
text(cex=2, x=colMeans(bp), y=-0.5, labs, xpd=NA, srt=35, adj = 1)

barplot(bias_cat_gonad[,biased_bins[biased_bins!="Unbiased"]], 
        beside = TRUE, 
        ylim = c(0, max(bias_cat_gonad)+1000), 
        xaxt='n',
        col = sex_cols,
        main = "",
        cex.main=2,
        cex.axis=2)
mtext("Gonad",3,outer = FALSE,cex=1.5,line=-1)
axis(2,lwd=2,labels = NA)
text(cex=2, x=colMeans(bp), y=-100, labs, xpd=NA, srt=35, adj=1)

barplot(bias_cat_liver[,biased_bins[biased_bins!="Unbiased"]], 
        beside = TRUE, 
        ylim = c(0, max(bias_cat_liver)+20), 
        xaxt='n',
        col = sex_cols,
        main = "",
        cex.main=2,
        cex.axis=2)
axis(2,lwd=2, labels = NA)
text(cex=2, x=colMeans(bp), y=-2, labs, xpd=NA, srt=35,adj=1)
mtext("Liver",3,outer = FALSE,cex=1.5,line=-1)

mtext("Number of Genes",2,outer=TRUE, cex=1.5, line=2.25)
mtext("Bias category",1, outer=TRUE, cex=1.5, line=4)

outer_legend("right", 
       legend = c("Female\nbiased", "Male\nbiased"), 
       pt.bg = sex_cols,
       pch=22,
       bty='n',
       ncol=1,
       cex=2,
       y.intersp = 1.5,
       pt.cex=2)

dev.off()
```

    ## png 
    ##   2

There is a weird pattern showing up in the liver where most of the male
biased genes are showing an “Extreme” level of sex bias. With the
heatmaps we created above, we likely guess that these extreme biased
genes are the ones driven by those two male individuals.

## Create a combined figure

``` r
figSBa <- image_ggplot(image_read_pdf('docs/figs/FigSF_logFC_boxplots.pdf'),interpolate = TRUE)
figSBb <- image_ggplot(image_read_pdf('docs/figs/FigSB_biasCat_counts.pdf'),interpolate = TRUE)
figSBc <- image_ggplot(image_read_pdf('docs/figs/FigSB_sexbias_upset.pdf'),interpolate = TRUE)


figSB <- wrap_plots(figSBa,
                    figSBb,
                    nrow = 2)
figSB <- figSB + plot_annotation(tag_levels = 'A')

figSB_all <- wrap_plots(figSB,
                        figSBc,
                        ncol = 2)
figSB_all <- figSB_all + plot_annotation(tag_levels = 'A')
figSB_all
```

![](fuscus_diff_expr_analysis_files/figure-gfm/figSB-1.png)<!-- -->

``` r
ggsave("docs/figs/FigSB.pdf", figSB_all, height=8, width=10)
ggsave("docs/figs/FigSB.png", figSB_all, height=8, width=10)
```

``` r
figHMa <- image_ggplot(image_read_pdf('docs/figs/liverMB_heatmap.pdf'),interpolate = TRUE)
figHMb <- image_ggplot(image_read_pdf('docs/figs/liverFB_heatmap.pdf'),interpolate = TRUE)

figSBHM <- wrap_plots(
  figSBa,
  figHMa,
  figSBb,
  figHMb,
  ncol = 2,
  nrow = 2)

figSBHM <- figSBHM + plot_annotation(tag_levels = 'A')
figSBHM

ggsave("docs/figs/FigSBHM.pdf", figSBHM, height=8, width=10)
ggsave("docs/figs/FigSBHM.png", figSBHM, height=8, width=10)
```

# Gene Ontologogy Analysis and BLASTing the different genes

There are several things that I want to accomplish with the gene
ontology analysis, to start with I want to classify the sex-biased genes
and then the sex-specific genes.

## BLASTING against the zebrafish proteome

The first step to the GO analysis is to BLAST the trinity sequences
against one of the PANTHER organisms. The zebrafish, *Danio reiro*, was
used for this as it was the only fish in the list.

A BLAST database was generated with the *D. reiro* proteome as:
`makeblastdb -in ../ncbi_dataset/data/GCF_000002035.6/protein.faa -out d_rerio_prot -dbtyp prot`

The trinity gene IDs that correspond to ALL of the female- or
male-biased genes were pulled out and saved to `.txt` files.

``` r
write.table(cbind(rownames(gill_fem_biased)),
            'FU_femG_biased_TRgenes.txt', 
            sep = "", 
            quote=FALSE, 
            row.names = FALSE, 
            col.names = FALSE)
write.table(cbind(rownames(gill_mal_biased)),
            'FU_malG_biased_TRgenes.txt', 
            sep = "", 
            quote=FALSE, 
            row.names = FALSE, 
            col.names = FALSE)
write.table(cbind(rownames(gonad_fem_biased)),
            'FU_femGon_biased_TRgenes.txt', 
            sep = "", 
            quote=FALSE, 
            row.names = FALSE, 
            col.names = FALSE)
write.table(cbind(rownames(gonad_mal_biased)),
            'FU_malGon_biased_TRgenes.txt', 
            sep = "", 
            quote=FALSE, 
            row.names = FALSE, 
            col.names = FALSE)
write.table(cbind(rownames(liver_fem_biased)),
            'FU_femL_biased_TRgenes.txt', 
            sep = "", 
            quote=FALSE, 
            row.names = FALSE, 
            col.names = FALSE)
write.table(cbind(rownames(liver_mal_biased)),
            'FU_malL_biased_TRgenes.txt', 
            sep = "", 
            quote=FALSE, 
            row.names = FALSE, 
            col.names = FALSE)
```

The `.txt` files were imported into the RCC and then BLASTed against the
zebrafish proteome using the following script:

``` bash
#!/bin/bash

#Create the arguements
input_dir_TR=$1 #Location of the .txt files that contain the trinity gene IDs
subset_fasta=$2 #Path to the subset_fasta_file script
assembly_file=$3 #Name/location of the de novo assembly
blast_database=$4
output_dir=$5 #Desired output directory for .fasta and BLAST files

## Loop through all the Trinity gene ID .txt files
for file in $input_dir_TR/*TRgenes.txt
    do

    #Extract the sample name from the file name
    sample=$(basename $file .txt)

    #Get the corresponding sequences for the Trinity Gene IDs
    echo "Running subset_fasta_file for ${sample}..."
    $subset_fasta -c -f $assembly_file -l $input_dir_TR/${sample}.txt

    #Rename the automated output from the subset_fasta_file script
    fasta_out=$(basename $assembly_file fasta)
    mv $fasta_out.subset.fasta ${sample}.fasta

    #Blast the sequences
    echo "Running BLAST for ${sample}..."
    blastn -db $blast_database -query ${sample}.fasta -out ${sample}_blast.txt \
        -evalue 0.001 \
        -num_threads 12 \
        -outfmt "6 qseqid qstart qend stitle sstart send evalue bitscore length pident gaps"

    #Move the outputs into desired output directory
    mv ${sample}* $output_dir

done
```

This script was run as
`bash bash_scripts/blast_pipeline.sh fuscus_high_exp_genes/high_expTG_names other_scripts/subset_fasta_file trinity_supertran_fuscus.fasta blastx ../genomes/d_rerio_prot/d_reiro_prot _reiro_blast.txt fuscus_high_exp_genes/high_exp_results/`.

The trinity gene IDs that correspond to all of the sex-specific genes
were pulled out above. Those were BLASTed with the script above using
**blastn** against the *Syngnathus scovelli* genome as:
`bash bash_scripts/blast_pipeline.sh fuscus_high_exp_genes/high_expTG_names other_scripts/subset_fasta_file trinity_supertran_fuscus.fasta blastn ../genomes/s_scov/s_scov_genome _blast.txt fuscus_high_exp_genes/high_exp_results/`.

## Read in and filter the BLAST results

I then exported all of the blast results as .txt files from the RCC and
need to first read them into R so that I can filter the results and make
sure only one hit is kept per gene.

Not that I have all of the BLAST files for SS genes and SB genes in R I
am going to filter them:

``` r
#Use lapply to apply the function to each dataset stored in the list created above
blast_output_filtered <- lapply(FU_blast_list, function(data_frame){
 
  #Pull out the Unique Trinity gene IDs
  uniqueID <- unique(data_frame[1])
  
  #Create an empty dataframe to store intermediate results in
  output <- data.frame(matrix(data = NA, ncol = ncol(data_frame), nrow = 0))
  colnames(output) <- c("trin_geneid", "trin_gene_start", "trin_gene_end", "reiro_prot_info", "prot_id", "reiro_prot_start", "reiro_prot_end", "evalue", "bit_score", "length", "pident", "gaps")
  
  #Generate a for loop that pulls out the lowest e-value + highest % identity for each gene
  for(gene in uniqueID$trin_geneid){
    
    #Subset the dataset for each Trinity gene
    this_trin_gene<- subset(data_frame, trin_geneid==gene)
    #Pull out gene with smallest e-value
    uni_gene_subset_lowe <- this_trin_gene[which(this_trin_gene$evalue == min(this_trin_gene$evalue)),]
    #In case mult. genes have same e-value, pull out gene with highest % identity
    uni_gene_subset_lowe_highpid <- uni_gene_subset_lowe[which(uni_gene_subset_lowe$pident == max(uni_gene_subset_lowe$pident)),]
    
    # keep only one of multiple identical rows
    uni_gene_subset_lowe_highpid <- unique(uni_gene_subset_lowe_highpid)

    # check that there is a single gene ID in scovelli, if not, only first row is kept
    if(length(gsub("^.*(GeneID:\\d+)\\].*$","\\1",uni_gene_subset_lowe_highpid$reiro_prot_info))>1){
    
    
      uni_gene_subset_lowe_highpid<-uni_gene_subset_lowe_highpid[1,]
    }
    
    #Add the final gene into the empty dataframe from above
    output <- rbind(output, uni_gene_subset_lowe_highpid)
  }
  
  return(output)
})
```

``` r
#Use lapply to apply the function to each dataset stored in the list created above
blast_output_filtered_SS <- lapply(FU_blast_list_sex_specific, function(data_frame){
 
  #Pull out the Unique Trinity gene IDs
  uniqueID <- unique(data_frame[1])
  
  #Create an empty dataframe to store intermediate results in
  output <- data.frame(matrix(data = NA, ncol = ncol(data_frame), nrow = 0))
  colnames(output) <- c("trin_geneid", "trin_gene_start", "trin_gene_end",
                        "scov_gene_info", "scov_prot_info",
                        "scov_gene_start", "scov_gene_end", 
                        "evalue", "bit_score", "length", "pident", "gaps")
  
  #Generate a for loop that pulls out the lowest e-value + highest % identity for each gene
  for(gene in uniqueID$trin_geneid){
    
    #Subset the dataset for each Trinity gene
    this_trin_gene<- subset(data_frame, trin_geneid==gene)
    #Pull out gene with smallest e-value
    uni_gene_subset_lowe <- this_trin_gene[which(this_trin_gene$evalue ==
                                                   min(this_trin_gene$evalue)),]
    #In case mult. genes have same e-value, pull out gene with highest % identity
    uni_gene_subset_lowe_highpid <- 
      uni_gene_subset_lowe[which(uni_gene_subset_lowe$pident ==
                                   max(uni_gene_subset_lowe$pident)),]
    
    # keep only one of multiple identical rows
    uni_gene_subset_lowe_highpid <- unique(uni_gene_subset_lowe_highpid)

    # check that there is a single gene ID in scovelli, if not, only first row is kept
    if(length(gsub("^.*(GeneID:\\d+)\\].*$",
                   "\\1",
                   uni_gene_subset_lowe_highpid$scov_gene_info))>1){
      
      uni_gene_subset_lowe_highpid <- uni_gene_subset_lowe_highpid[1,]
    }
    
    #Add the final gene into the empty dataframe from above
    output <- rbind(output, uni_gene_subset_lowe_highpid)
  }
  
  return(output)
})
```

## Collect gene names corresponding to each protein -sex-biased genes

In order to run PANTHER we need IDs that follow one of their supported
formats. From NCBI, the gene name works (e.g. rpl27). That information
is not automatically supplied with the BLAST output, but we can use the
.gff file from NCBI to pull out the gene names that correspond to the
proteinID (e.g. NP_956018.1).

``` r
#Read in the GFF file for zebrafish that has info about the geneID
reiro_gff <- read.delim("data/d_reiro_genomic.gff",
                        header = FALSE,
                        comment.char = "#")

#Keep only the columns I want and rename them
reiro_gff <- reiro_gff[,c(1:5,9)]
colnames(reiro_gff) <- c("seqname", "source", "feature", "start", "end", "gene_info")

#Subset the dataset to only include the rows we're interested in
genes <- reiro_gff[reiro_gff$feature == "CDS",]

#Pull out the gene name that corresponds to each protein ID
genes$gene_name <- gsub("^(.*;)(gene=)(\\w+\\d*);(.*$)", "\\3", genes$gene_info)
genes$prot_id <- gsub("^(.*;)(protein_id=)(.*)$", "\\3", genes$gene_info)

#Merge the gene names pulled out above with the rest of the BLAST datasets based on the proteinID
blast_output_merged <- lapply(blast_output_filtered, function(dataframe){
  
   output <- merge(dataframe, unique(genes[,7:8]), by = "prot_id")
  
   return(output)
})
```

I now need to save only the gene_id for all of the BLASTed genes to then
upload to PANTHER.

``` r
write.table(c(blast_output_merged$FU_femG_biased_TRgenes_reiro_blast$gene_name, 
              blast_output_merged$FU_femGon_biased_TRgenes_reiro_blast$gene_name,
              blast_output_merged$FU_femL_biased_TRgenes_reiro_blast$gene_name, 
              blast_output_merged$FU_malG_biased_TRgenes_reiro_blast$gene_name, 
              blast_output_merged$FU_malGon_biased_TRgenes_reiro_blast$gene_name, 
              blast_output_merged$FU_malL_biased_TRgenes_reiro_blast$gene_name),
           'FU_GOnames_SBG.txt', 
            sep = " ", 
            quote=FALSE, 
            row.names = FALSE, 
            col.names = FALSE)

write.table(c(blast_output_merged$FU_femL_biased_TRgenes_reiro_blast$gene_name, 
              blast_output_merged$FU_malL_biased_TRgenes_reiro_blast$gene_name),
           'FUL_GOnames_SBG.txt', 
            sep = " ", 
            quote=FALSE, 
            row.names = FALSE, 
            col.names = FALSE)
```

## Pull out the gene names and protein ID for sex-specific BLAST

Once all the BLAST results for the sex-specific genes have been filtered
I want to pull out specifically the name of the gene and the protein ID
from the `scov_gene_info` column and save them as their own column.

``` r
for (i in 1:length(blast_output_filtered_SS)) {
  
  #Pulling out the gene name and adding it to a new column
  blast_output_filtered_SS[[i]]$gene <- gsub("^(.*)(gene=)(\\w+\\d*)(.*$)",
                                             "\\3",
                                             blast_output_filtered_SS[[i]]$scov_gene_info)
  
  #Pulling out the protein name and adding it to a new column
  blast_output_filtered_SS[[i]]$protein <- gsub("^(.*)(protein=)(.*)([[:graph:]]\\s[[:graph:]])(protein_id=)(.*$)",
                                                "\\3",
                                                blast_output_filtered_SS[[i]]$scov_gene_info)
  
  
}
```

Now that I have isolated the gene names and protein IDS I want to reform
the BLAST results so that they can all be compiled into one single
dataset.

``` r
#Create a long format dataset that has the tissue information, BLAST gene and protein, 
#geneID, and sex 
blast_long_SS <- data.frame(
  tissue=c(rep("Gill", 
               sum(nrow(blast_output_filtered_SS$FU_femG_specific_TRgenes_blast),
                   nrow(blast_output_filtered_SS$FU_malG_specific_TRgenes_blast))),
           rep("Gonad", 
               sum(nrow(blast_output_filtered_SS$FU_femGon_specific_TRgenes_blast),
                   nrow(blast_output_filtered_SS$FU_malGon_specific_TRgenes_blast))),
           rep("Liver",
               sum(nrow(blast_output_filtered_SS$FU_femL_specific_TRgenes_blast),
                   nrow(blast_output_filtered_SS$FU_malL_specific_TRgenes_blast)))
         ),
  sex=c(rep("Female",
            nrow(blast_output_filtered_SS$FU_femG_specific_TRgenes_blast)),
        rep("Male",
            nrow(blast_output_filtered_SS$FU_malG_specific_TRgenes_blast)),
        rep("Female",
            nrow(blast_output_filtered_SS$FU_femGon_specific_TRgenes_blast)),
        rep("Male",
            nrow(blast_output_filtered_SS$FU_malGon_specific_TRgenes_blast)),
        rep("Female",
            nrow(blast_output_filtered_SS$FU_femL_specific_TRgenes_blast)),
        rep("Male",
            nrow(blast_output_filtered_SS$FU_malL_specific_TRgenes_blast))
        ),
  gene=c(blast_output_filtered_SS$FU_femG_specific_TRgenes_blast$gene,
         blast_output_filtered_SS$FU_malG_specific_TRgenes_blast$gene,
         blast_output_filtered_SS$FU_femGon_specific_TRgenes_blast$gene,
         blast_output_filtered_SS$FU_malGon_specific_TRgenes_blast$gene,
         blast_output_filtered_SS$FU_femL_specific_TRgenes_blast$gene,
         blast_output_filtered_SS$FU_malL_specific_TRgenes_blast$gene
         ),
  protein=c(blast_output_filtered_SS$FU_femG_specific_TRgenes_blast$protein,
            blast_output_filtered_SS$FU_malG_specific_TRgenes_blast$protein,
            blast_output_filtered_SS$FU_femGon_specific_TRgenes_blast$protein,
            blast_output_filtered_SS$FU_malGon_specific_TRgenes_blast$protein,
            blast_output_filtered_SS$FU_femL_specific_TRgenes_blast$protein,
            blast_output_filtered_SS$FU_malL_specific_TRgenes_blast$protein
            ),
  trinityID=c(blast_output_filtered_SS$FU_femG_specific_TRgenes_blast$trin_geneid,
              blast_output_filtered_SS$FU_malG_specific_TRgenes_blast$trin_geneid,
              blast_output_filtered_SS$FU_femGon_specific_TRgenes_blast$trin_geneid,
              blast_output_filtered_SS$FU_malGon_specific_TRgenes_blast$trin_geneid,
              blast_output_filtered_SS$FU_femL_specific_TRgenes_blast$trin_geneid,
              blast_output_filtered_SS$FU_malL_specific_TRgenes_blast$trin_geneid
              )
  )

#write.csv(blast_long_SS, "FU_BLAST_results_sex_specific.csv", row.names = FALSE)
```

## Reading in the PANTHER results - Sex-biased Genes

In PANTHER I uploaded the `FU_GOnames_SBG.txt` file, selected *Danio
rerio* as the organism, and then chose “Functional classification viewed
in gene list” as the analysis type.

I exported those results to a .txt file and will now be reading that
into R and visualizing the results.

``` r
#Read in the .txt file containing PANTHER results
panther <- read.delim("data/pantherGeneList_SBG.txt", 
                      header = FALSE)

#Add the approporiate column names for the PANTHER data
colnames(panther) <- c("GeneID", "MappedID", 
                       "GeneName", "pantherFAM", 
                       "panther_prot_class", "species")

#Add the identified PANTHER protein classes to the BLAST datasets
blast_output_merged$FU_femG_biased_TRgenes_reiro_blast <-
  merge(blast_output_merged$FU_femG_biased_TRgenes_reiro_blast, 
        panther[, c(2,5)], 
        by.x = "gene_name", 
        by.y = "MappedID", 
        all.x = TRUE)
blast_output_merged$FU_femGon_biased_TRgenes_reiro_blast <-
  merge(blast_output_merged$FU_femGon_biased_TRgenes_reiro_blast, 
        panther[, c(2,5)], 
        by.x = "gene_name", 
        by.y = "MappedID", 
        all.x = TRUE)
blast_output_merged$FU_femL_biased_TRgenes_reiro_blast <-
  merge(blast_output_merged$FU_femL_biased_TRgenes_reiro_blast, 
        panther[, c(2,5)], 
        by.x = "gene_name", 
        by.y = "MappedID", 
        all.x = TRUE)
blast_output_merged$FU_malG_biased_TRgenes_reiro_blast <-
  merge(blast_output_merged$FU_malG_biased_TRgenes_reiro_blast, 
        panther[, c(2,5)], 
        by.x = "gene_name", 
        by.y = "MappedID", 
        all.x = TRUE)
blast_output_merged$FU_malGon_biased_TRgenes_reiro_blast <-
  merge(blast_output_merged$FU_malGon_biased_TRgenes_reiro_blast, 
        panther[, c(2,5)], 
        by.x = "gene_name", 
        by.y = "MappedID", 
        all.x = TRUE)
blast_output_merged$FU_malL_biased_TRgenes_reiro_blast <-
  merge(blast_output_merged$FU_malL_biased_TRgenes_reiro_blast, 
        panther[, c(2,5)], 
        by.x = "gene_name", 
        by.y = "MappedID", 
        all.x = TRUE)
```

Now that we have the PANTHER protein classes attached to all of the
successfully BLASTed sex-biased genes within each organ, we can plot the
proportions of genes falling into each category to see what is most
important for each of the organ types.

``` r
#Count the number of genes that are falling into each protein class
go_tables <- lapply(blast_output_merged, function(dat){
  
  #Create a table that counts the number of genes falling within each PANTHER prot. class
  panther_table <- data.frame(table(dat$panther_prot_class))
  
  #Make sure the protein names are in there as a character
  panther_table$Var1 <- as.character(panther_table$Var1)
  
  #Change all of the null results to say Unclassified
  panther_table$Var1[panther_table$Var1==""] <- "Unclassified"

  return(panther_table)
  })

#Merge all of the different organ data into one dataset
#Create an empty dataframe to store all of the values in
all_go_dat <- data.frame(matrix(ncol=4,nrow=0))

#Add column names to that empty dataset
colnames(all_go_dat) <- c(
  colnames(go_tables$FU_femG_biased_TRgenes_reiro_blast),
  "bias_cat", "tissue")

#Create a vector of you tissues, in the same order that they appear in the go_table list
tissues <- c("Gill", "Gonad", "Liver", 
             "Gill", "Gonad", "Liver")

#Loop through all of the go_tables and combine into one
for(i in 1:length(go_tables)){
  #browser()
  tmp <- go_tables[[i]]
  if (nrow(tmp) == 0) {
    print(tmp)
  }else{
  
  tmp$bias_cat <- names(go_tables)[i]
  tmp$tissue <- tissues[[i]]
  all_go_dat <- rbind(all_go_dat, tmp)}
  
}

#Calculate the frequency at which each protein class occurs within the different tissue types
##Create an empty dataset to store the values in
all_go_sums <- data.frame(matrix(ncol = 4,
                                 nrow = 0))

##Add appropriate column names
colnames(all_go_sums) <- c("Freq", "prot_class", 
                           "tissue", "prop")

##Loop through each organ type to calculate frequencies
for(organ in unique(all_go_dat$tissue)){
  
  #subset the dataset based on the organ
  tmp <- all_go_dat[all_go_dat$tissue == organ, ]
  
  #Sum up all of the genes in each category
  go_sums <- as.data.frame(tapply(tmp$Freq, tmp$Var1, sum))
  
  #Add the protein class as a column rather than row names
  go_sums$prot_class <- rownames(go_sums)
  
  #Change the column name to "Frequency"
  colnames(go_sums)[1] <- "Freq"
  
  #Add the organ in a separate column and remove rownames
  go_sums$tissue <- organ
  rownames(go_sums) <- NULL
  
  #Calculate the respective proportion for each protein class
  go_sums$prop <- go_sums$Freq/sum(go_sums$Freq)
  
  #Export the data to the previously created dataframe
  all_go_sums <- rbind(all_go_sums, go_sums)
}

#Change all of the protein classes that have a frequency of less than 0.02 to "Other"
for(class in unique(all_go_sums$prot_class)){
  
  #Pull out the rows that contain that protein class
  match_rows <- which(all_go_sums$prot_class==class)
  
  #If any of the proportions associated with those rows is less than 0.02, classify as other
  if(any(all_go_sums$prop[match_rows]>0.02)==FALSE){
    
    all_go_sums$prot_class[match_rows] <- "other"
    
  }
}


#write.csv(all_go_sums, "data/GO_freq_SBG.csv", row.names = FALSE)

#Generate the GO plot
organ_cols <- c("Gill" = "#20B2AA", "Gonad" = "#EEB422",
                "Liver" = "#EE8262")

#pdf("docs/figs/FigGO_SBGbarplot.pdf",height=8, width=12)

all_go_sums$prot_class <- str_wrap(all_go_sums$prot_class, 
                                       width = 60)

ggplot(all_go_sums[!(all_go_sums$prot_class %in% c("other", "Unclassified")),], 
       aes(fct_rev(prot_class), prop, fill = tissue)) +   
  geom_bar(position = "dodge", stat="identity") +
  scale_fill_manual(values = organ_cols) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=20),
        axis.text.y = element_text(size=15, #hjust = 0.5, 
                                   margin = margin(b = 20),
                                   vjust = 0.5),
        text=element_text(size=20),
        legend.position = "bottom",
        axis.ticks.length.y = unit(.25, "cm")) +
  coord_flip() + 
  labs(y="Proportion of sex-biased genes", x="")
```

![](fuscus_diff_expr_analysis_files/figure-gfm/plot-panther-SBG-1.png)<!-- -->

``` r
#dev.off()
```

## Looking at GO groups for the male- and female-biased genes

``` r
#Calculate the frequency at which each protein class occurs within the different tissue types
##Create an empty dataset to store the values in
all_go_sums_liver <- data.frame(matrix(ncol = 5,
                                      nrow = 0))

##Add appropriate column names
colnames(all_go_sums_liver) <- c(colnames(all_go_dat),
                                "prop")

##Loop through each bias type to calculate frequencies
for (bias in c("FU_femL_biased_TRgenes_reiro_blast",
               "FU_malL_biased_TRgenes_reiro_blast")) {
  
  tmp <- all_go_dat[all_go_dat$bias_cat == bias, ]
  
  #Calculate the respective proportion for each protein class
  tmp$prop <- tmp$Freq/sum(tmp$Freq)
  
  #Export the data to the previously created dataframe
  all_go_sums_liver <- rbind(tmp, all_go_sums_liver)
  
}


#Change all of the protein classes that have a frequency of less than 0.02 to "Other"
for(class in unique(all_go_sums_liver$Var1)){
  
  #Pull out the rows that contain that protein class
  match_rows <- which(all_go_sums_liver$Var1==class)
  
  #If any of the proportions associated with those rows is less than 0.02, classify as other
  if(any(all_go_sums_liver$prop[match_rows]>0.01)==FALSE){
    
    all_go_sums_liver$Var1[match_rows] <- "other"
    
  }
}

#Generate the GO plot
bias_col <- c("FU_femL_biased_TRgenes_reiro_blast" = "#7fc97f", 
               "FU_malL_biased_TRgenes_reiro_blast" = "#beaed4")


#pdf("docs/figs/FigGO_SBGliver.pdf",height=10, width=12)

ggplot(all_go_sums_liver[!(all_go_sums_liver$Var1 %in% c("other", "Unclassified")),], 
       aes(fct_rev(Var1), prop, fill = bias_cat)) +   
  geom_bar(position = "dodge", stat="identity") +
  scale_fill_manual(values = bias_col) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1,size=20),
        axis.text.y = element_text(size=15, hjust = 1,
                                   vjust = 0.3),
        text=element_text(size=20),
        legend.position = "bottom",
        axis.ticks.length.y = unit(.25, "cm")) +
  coord_flip() + 
  labs(y="Proportion of sex-biased genes", x="")
```

![](fuscus_diff_expr_analysis_files/figure-gfm/plot-panther-SBG-liver-1.png)<!-- -->

``` r
#dev.off()
```

``` r
#Calculate the frequency at which each protein class occurs within the different tissue types
##Create an empty dataset to store the values in
all_go_sums_gonad <- data.frame(matrix(ncol = 5,
                                      nrow = 0))

##Add appropriate column names
colnames(all_go_sums_gonad) <- c(colnames(all_go_dat),
                                "prop")

##Loop through each bias type to calculate frequencies
for (bias in c("FU_femGon_biased_TRgenes_reiro_blast",
               "FU_malGon_biased_TRgenes_reiro_blast")) {
  
  tmp <- all_go_dat[all_go_dat$bias_cat == bias, ]
  
  #Calculate the respective proportion for each protein class
  tmp$prop <- tmp$Freq/sum(tmp$Freq)
  
  #Export the data to the previously created dataframe
  all_go_sums_gonad <- rbind(tmp, all_go_sums_gonad)
  
}


#Change all of the protein classes that have a frequency of less than 0.02 to "Other"
for(class in unique(all_go_sums_gonad$Var1)){
  
  #Pull out the rows that contain that protein class
  match_rows <- which(all_go_sums_gonad$Var1==class)
  
  #If any of the proportions associated with those rows is less than 0.02, classify as other
  if(any(all_go_sums_gonad$prop[match_rows]>0.01)==FALSE){
    
    all_go_sums_gonad$Var1[match_rows] <- "other"
    
  }
}

#Generate the GO plot
bias_col <- c("FU_femGon_biased_TRgenes_reiro_blast" = "#7fc97f", 
               "FU_malGon_biased_TRgenes_reiro_blast" = "#beaed4")


#pdf("docs/figs/FigGO_SBGgonad.pdf",height=10, width=12)

ggplot(all_go_sums_gonad[!(all_go_sums_gonad$Var1 %in% c("other", "Unclassified")),], 
       aes(fct_rev(Var1), prop, fill = bias_cat)) +   
  geom_bar(position = "dodge", stat="identity") +
  scale_fill_manual(values = bias_col) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1,size=20),
        axis.text.y = element_text(size=15, hjust = 1,
                                   vjust = 0.3),
        text=element_text(size=20),
        legend.position = "bottom",
        axis.ticks.length.y = unit(.25, "cm")) +
  coord_flip() + 
  labs(y="Proportion of sex-biased genes", x="")
```

![](fuscus_diff_expr_analysis_files/figure-gfm/plot-panther-SBG-gonads-1.png)<!-- -->

``` r
#dev.off()
```

## Create combined Gene Ontology Figure

``` r
#Combine the GO figures
figGOa <- image_ggplot(image_read_pdf('docs/figs/FigGO_SBGbarplot.pdf'),
                        interpolate = TRUE)
figGOb <- image_ggplot(image_read_pdf('docs/figs/FigGO_SBGliver.pdf'),
                        interpolate = TRUE)
figGOc <- image_ggplot(image_read_pdf('docs/figs/FigGO_SBGgonad.pdf'),
                           interpolate = TRUE)


figGO <- wrap_plots(figGOa,
                    figGOb,
                    figGOc,
                    ncol = 3)

figGO <- figGO + plot_annotation(tag_levels = 'A')

ggsave("docs/figs/FigGO_all.pdf", figGO, height=5, width=12)
ggsave("docs/figs/FigGO_all.png", figGO, height=5, width=12)
```

``` r
write.table(blast_output_merged$FU_femGon_biased_TRgenes_reiro_blast$gene_name,
            'FU_GOnames_SBG_ovary.txt',
            sep = " ",
            quote=FALSE,
            row.names = FALSE,
            col.names = FALSE)

write.table(blast_output_merged$FU_malGon_biased_TRgenes_reiro_blast$gene_name,
            'FU_GOnames_SBG_testes.txt',
            sep = " ",
            quote=FALSE,
            row.names = FALSE,
            col.names = FALSE)

write.table(blast_output_merged$FU_femL_biased_TRgenes_reiro_blast$gene_name,
            'FU_GOnames_SBG_fLiver.txt',
            sep = " ",
            quote=FALSE,
            row.names = FALSE,
            col.names = FALSE)

write.table(blast_output_merged$FU_malL_biased_TRgenes_reiro_blast$gene_name,
            'FU_GOnames_SBG_mLiver.txt',
            sep = " ",
            quote=FALSE,
            row.names = FALSE,
            col.names = FALSE)
```
