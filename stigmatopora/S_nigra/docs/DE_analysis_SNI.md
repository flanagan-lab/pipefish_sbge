Differential Expression Analysis in *Stigmatopora nigra*
================



-   [Single factor analysis - Comparing Males v Females across all
    organs](#single-factor-analysis---comparing-males-v-females-across-all-organs)
    -   [Visualizing the results](#visualizing-the-results)
        -   [MA-plot - MvF Differntial
            expression](#ma-plot---mvf-differntial-expression)
        -   [Heatmap - overall
            expression](#heatmap---overall-expression)
        -   [Sample-dist Heatmap](#sample-dist-heatmap)
        -   [PCA plots](#pca-plots)
-   [Multifactor design - Comparing M v F across the diff tissue
    types](#multifactor-design---comparing-m-v-f-across-the-diff-tissue-types)
    -   [Invesitgate the results of the differential
        expression](#invesitgate-the-results-of-the-differential-expression)
        -   [M-F Liver Comparisson](#m-f-liver-comparisson)
        -   [M-F Gill Comparisson](#m-f-gill-comparisson)
        -   [M-F Gonad Comparisson](#m-f-gonad-comparisson)
        -   [MA Plots](#ma-plots)
    -   [Creating an Upset Plot](#creating-an-upset-plot)
    -   [Variation in FC across sex-bias and tissue
        type](#variation-in-fc-across-sex-bias-and-tissue-type)
    -   [Categorizing sex-specific
        genes](#categorizing-sex-specific-genes)
        -   [Investigating the results of the sex-specific
            subsetting](#investigating-the-results-of-the-sex-specific-subsetting)
        -   [Pulling out the gene IDs for all of the sex-specific
            genes](#pulling-out-the-gene-ids-for-all-of-the-sex-specific-genes)
    -   [Creating categories and binning the sex-biased genes based on
        degree of
        logFC](#creating-categories-and-binning-the-sex-biased-genes-based-on-degree-of-logfc)
    -   [Gene Ontology Analysis](#gene-ontology-analysis)
        -   [Using BLAST in command line](#using-blast-in-command-line)
        -   [Filtering Blast results](#filtering-blast-results)
        -   [Merging the BLAST filtered Data frames with the sex-biased
            information](#merging-the-blast-filtered-data-frames-with-the-sex-biased-information)
        -   [Blasting the sex-specific genes against the *C.
            intestinalis*
            genome](#blasting-the-sex-specific-genes-against-the-c-intestinalis-genome)
        -   [Blasting against the Zebrafish
            genome](#blasting-against-the-zebrafish-genome)
    -   [Looking at GO groups for the male- and female-biased
        genes](#looking-at-go-groups-for-the-male--and-female-biased-genes)

``` r
setwd("~/Documents/GitHub/pipefish_sbge/stigmatopora/S_nigra")
# This is a cohesive list of all the libraries used in this
# document
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(PCAtools)
library(spfTools)
library(pdftools)
library(magick)
library(patchwork)
library(ggpubr)
library(cowplot)
library(UpSetR)
library(ggplot2)
library(kableExtra)
library(tidyverse)
library(knitr)
```

``` r
# The abundance matrix generated via salmon and tximport to
# be used for the DE analysis
txi.nigra <- readRDS("data/txi.salmon_SN.RDS")

# The samples file generated for tximport
samples <- read.table("data/SN_samples.txt", header = TRUE)

# Make sure the conditions are in the samples file as a
# factor
samples$sex <- as.factor(samples$sex)
samples$organ <- as.factor(samples$organ)
```

The package `DESeq2` was used for the differential expression analysis
outlined below.

# Single factor analysis - Comparing Males v Females across all organs

To analyze your data with DESeq2 you must first generate the
DESeqDataSet. In order to do this we need the abundance matrix generated
with `tximport` and a `samples` file that lays out all of the
conditions. The model for this single factor analysis was run as counts
\~ Sex.

``` r
# Create the DESeq dataset
dds_SN <- DESeqDataSetFromTximport(txi.nigra, colData = samples,
    design = ~sex)
```

The data is then pre-filtered to remove low gene counts before running
further DESeq2 functions. By doing this we remove rows in which there
are very few reads thus reducing the memory size of the `dds` object and
increasing the speed at which we can use the transformation and testing
functions in DESeq2.

The cutoff here was to remove rows that had counts fewer than 10 across
all samples.

``` r
# only keeping rows that have at lead 10 reads total
keep <- rowSums(counts(dds_SN)) >= 10
dds_SN <- dds_SN[keep, ]
```

After filtering we can now perform the standard differential expression
analysis that is wrapped into DESeq2.

``` r
# Generate the expression values
dds_SN_exp <- DESeq(dds_SN)

# Compile the results
res <- results(dds_SN_exp)
res
```

    ## log2 fold change (MLE): sex M vs F 
    ## Wald test p-value: sex M vs F 
    ## DataFrame with 100796 rows and 6 columns
    ##                         baseMean log2FoldChange     lfcSE       stat    pvalue
    ##                        <numeric>      <numeric> <numeric>  <numeric> <numeric>
    ## TRINITY_DN0_c0_g1     1609.17100      0.0811555  0.151157   0.536896 0.5913396
    ## TRINITY_DN0_c0_g2       33.63765     -0.4878635  0.479273  -1.017925 0.3087137
    ## TRINITY_DN0_c0_g3        2.95197     -2.5194164  1.246570  -2.021078 0.0432717
    ## TRINITY_DN0_c1_g1       74.34900      0.1295842  0.961197   0.134815 0.8927578
    ## TRINITY_DN0_c2_g1        1.80104      0.9685934  0.944291   1.025736 0.3050159
    ## ...                          ...            ...       ...        ...       ...
    ## TRINITY_DN99980_c0_g1   0.722176     -1.0040679  1.369190 -0.7333298  0.463357
    ## TRINITY_DN9999_c0_g1    2.776003      0.0136508  0.972717  0.0140337  0.988803
    ## TRINITY_DN99991_c0_g1   0.152641     -0.7853673  1.943352 -0.4041303  0.686117
    ## TRINITY_DN99993_c0_g1   0.181944      0.2055342  1.480007  0.1388738  0.889550
    ## TRINITY_DN99998_c0_g1   0.185763      1.0252108  2.956045  0.3468184  0.728728
    ##                            padj
    ##                       <numeric>
    ## TRINITY_DN0_c0_g1      0.775985
    ## TRINITY_DN0_c0_g2      0.536273
    ## TRINITY_DN0_c0_g3      0.160673
    ## TRINITY_DN0_c1_g1      0.951185
    ## TRINITY_DN0_c2_g1      0.532517
    ## ...                         ...
    ## TRINITY_DN99980_c0_g1        NA
    ## TRINITY_DN9999_c0_g1   0.995264
    ## TRINITY_DN99991_c0_g1        NA
    ## TRINITY_DN99993_c0_g1        NA
    ## TRINITY_DN99998_c0_g1        NA

Once that has finished we can now start exploring some of the
single-factor analysis results

``` r
## Ordering our results based on p-value
resOrdered <- res[order(res$pvalue), ]
resOrdered
```

    ## log2 fold change (MLE): sex M vs F 
    ## Wald test p-value: sex M vs F 
    ## DataFrame with 100796 rows and 6 columns
    ##                        baseMean log2FoldChange     lfcSE      stat      pvalue
    ##                       <numeric>      <numeric> <numeric> <numeric>   <numeric>
    ## TRINITY_DN16454_c0_g1   66.2318        25.3751   1.54968   16.3744 2.91396e-60
    ## TRINITY_DN2584_c0_g1  1126.5583        29.0515   1.95240   14.8799 4.45084e-50
    ## TRINITY_DN5911_c0_g1    64.7404        25.4467   1.72192   14.7781 2.02887e-49
    ## TRINITY_DN4321_c0_g3    43.4312        24.8885   1.68537   14.7674 2.37589e-49
    ## TRINITY_DN2003_c1_g1    36.1188        24.6296   1.70952   14.4074 4.65183e-47
    ## ...                         ...            ...       ...       ...         ...
    ## TRINITY_DN98855_c0_g1         0              0         0         0           1
    ## TRINITY_DN99164_c0_g1         0              0         0         0           1
    ## TRINITY_DN99189_c0_g1         0              0         0         0           1
    ## TRINITY_DN99610_c0_g1         0              0         0         0           1
    ## TRINITY_DN99658_c0_g1         0              0         0         0           1
    ##                              padj
    ##                         <numeric>
    ## TRINITY_DN16454_c0_g1 1.16130e-55
    ## TRINITY_DN2584_c0_g1  8.86896e-46
    ## TRINITY_DN5911_c0_g1  2.36716e-45
    ## TRINITY_DN4321_c0_g3  2.36716e-45
    ## TRINITY_DN2003_c1_g1  3.70779e-43
    ## ...                           ...
    ## TRINITY_DN98855_c0_g1          NA
    ## TRINITY_DN99164_c0_g1          NA
    ## TRINITY_DN99189_c0_g1          NA
    ## TRINITY_DN99610_c0_g1          NA
    ## TRINITY_DN99658_c0_g1          NA

``` r
summary(res)
```

    ## 
    ## out of 99810 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 4765, 4.8%
    ## LFC < 0 (down)     : 3915, 3.9%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 60943, 61%
    ## (mean count < 1)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
# How many ADJUSTED p-values were less than 0.1?
sum(res$padj < 0.1, na.rm = TRUE)
```

    ## [1] 8680

``` r
# Looking at an alpha=0.05, the default is 0.1
res05 <- results(dds_SN_exp, alpha = 0.05)
summary(res05)
```

    ## 
    ## out of 99810 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 3437, 3.4%
    ## LFC < 0 (down)     : 3154, 3.2%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 64811, 65%
    ## (mean count < 1)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
sum(res05$padj < 0.05, na.rm = TRUE)
```

    ## [1] 6591

out of 99810 with nonzero total read count adjusted p-value \< 0.05 LFC
\> 0 (up) : 3437, 3.4% LFC \< 0 (down) : 3154, 3.2% outliers \[1\] : 0,
0% low counts \[2\] : 64811, 65% (mean count \< 1) \[1\] see
‘cooksCutoff’ argument of ?results \[2\] see ‘independentFiltering’
argument of ?results

sum(res05\$padj \< 0.05, na.rm = TRUE) 6591

## Visualizing the results

### MA-plot - MvF Differntial expression

Generate an MA-plot to show the log2 fold changes attributable to sex
over the mean of normalized counts for all of the samples in the `dds`.
Points will be colored if the adjusted p-value is less that 0.1.

![LogFC versus the mean normalized count for all of the genes. Points
that are blue have a p-value less than 0.1. The two plots are the same
with the only difference coming from an adjusted y-limit. Points that
are triangles represent genes with a fold change higher/lower than the
y-limit.](DE_analysis_SNI_files/figure-gfm/MA-plot-1.png)

### Heatmap - overall expression

We can also generate a heat map to look at overall expression levels
across our samples. Note, this is not differentially expressed genes.

![Heatmap showing the expression level across all organs for the top 20
genes with the highest expression
levels.](DE_analysis_SNI_files/figure-gfm/heatmap-1.png)

From the heatmap we can see that for all of these top 20 expressed genes
the highest expression in found in the liver. S34 (female liver) is
showing expression patterns that are different from the other livers.
There are also a few genes that appear to be highly expressed across all
of the tissue types. With this heatmap we can also pull out the names of
the Trinity genes that are showing this high expression and BLAST the
corresponding sequences to see what the genes are.

``` r
# Pull out the corresponding trinity_geneIDS that are
# plotted in the heatmap
heatmap_TG <- cbind(row.names(assay(vsd)[select, ]))
# Run once to create the file and then commented out
# write.table(heatmap_TG,
#'heatmap_trinitygenes.txt',
# sep = '', quote=FALSE, row.names = FALSE, col.names =
# FALSE)
```

### Sample-dist Heatmap

We can then generate a sample-to-sample distances heatmap that will give
an overview of the similarities and differences between our samples.

![Sample-to-sample distances heatmap showing the similarities and
differences between all of the samples. The darker the color, the more
similar the samples are. The diagonal line of very dark blue represents
the comparissons between the same
samples.](DE_analysis_SNI_files/figure-gfm/sample-dist-1.png)

Here we can see that the highest similarities (aside from same samples
comparisons) is between samples from the same organ. Again, S34 is
showing unexpected patterns - I may end up removing all organs related
to S34.

### PCA plots

Next we can look at PCA plots for our samples to see how our different
samples are grouping. The size of the points for these plots will be
scaled by the mating success of the individual to see if that may be at
all correlated with the sample grouping.

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

![Principal components analysis pairsplot showing PCA axis 1 -
4.](DE_analysis_SNI_files/figure-gfm/pca-pairs-plot-1.png)

From the pairs plot we can see that PC1 appears to be explaining
differences between the liver and the other tissues (we can see the one
problematic liver sample from above sticking out from all the rest), PC2
appears to be explaining differences between gonadal and somatic
tissues, PC3 seems to account for differences between male and female
samples, especially in the gonads and in PC4 we can see the one
problematic liver sample again sticking out from all the rest. I am
going to plot all of PC2 - 4 against PC1 below for a clearer picture.

    ## using ntop=100796 top features by variance

![Principal components analysis reflects that most variation corresponds
to differences in expression between organs (green: gill; yellow: gonad;
pink: liver), rather than variation due to sex (circle = female;
triangle = male) or mating status (the size of points reflects an
individual’s mating success). The first axis explains 32% of variation
in the dataset, with gills and liver, two types of somatic tissues,
sitting on opposite ends of that axis. The second axis explains 23% of
the variation in gene expression, with somatic tissues on the opposite
side of the axis from the
gonads.](DE_analysis_SNI_files/figure-gfm/pca-scaled-organ-1.png)

![The same Principal Components Analysis as above but with PC3. This
third axis explains NA% of the variation, largely attributed to
differences between the sexes, particularily in the
gonads.](DE_analysis_SNI_files/figure-gfm/pca1v3-1.png)

![The same Principal Components Analysis as above but with PC4. This
fourth axis explains NA% of the
variation.](DE_analysis_SNI_files/figure-gfm/pca1v4-1.png)

#### Saving the PCA plots to use in figure creation

I will be including the PCAs describing the relationship between PC1 - 3
in the main text, as they account for the majority of the differences
between our samples. The plot with PC1 and PC4 I still want to export,
and will likely include as a supplemental to highlight why samples
related to S34 will be removed from the subsequent analyses.

``` r
# Setting the shapes I want to show up for the different
# organs
organ_shapes <- c(Gill = 16, Liver = 17, Gonad = 15)

# Setting the colors I want for each sex
sex_cols <- c(F = "#7fc97f", M = "#beaed4")
```

``` r
# Create the blank pdf to store the plot in
pdf("figs/Fig_PCA1v2.pdf", height = 6, width = 6)

# Set the plotting parameters
par(mar = c(4, 5, 4, 1), oma = c(2, 2, 2, 2))

# Create the plot for PC1 v PC2
plot(pca_plotting_data$PC1, pca_plotting_data$PC2, col = paste0(sex_cols[pca_plotting_data$sex],
    "75"), pch = organ_shapes[pca_plotting_data$organ], cex = 2,
    cex.lab = 2, cex.axis = 1.75, xlab = paste0("PC1: ", percentVar[1],
        "% variance"), ylab = paste0("PC2: ", percentVar[2],
        "% variance"), bty = "l", xpd = TRUE)

# Add a legend describing the sex
outer_legend("top", c("Female", "Male"), pch = 18, bty = "n",
    col = paste0(sex_cols, "75"), cex = 1.75, ncol = 2, pt.cex = 3)

dev.off()
```

``` r
# Create the blank pdf to store the plot in
pdf("figs/Fig_PCA1v3.pdf", height = 6, width = 6)

# Set the plotting parameters
par(mar = c(4, 5, 4, 1), oma = c(2, 2, 2, 2))

# Create the plot for PC1 v PC3
plot(pca_plotting_data$PC1, pca_plotting_data$PC3, col = paste0(sex_cols[pca_plotting_data$sex],
    "75"), pch = organ_shapes[pca_plotting_data$organ], cex = 2,
    cex.lab = 2, cex.axis = 1.75, xlab = paste0("PC1: ", percentVar[1],
        "% variance"), ylab = paste0("PC3: ", percentVar[3],
        "% variance"), bty = "l", xpd = TRUE)

# Add legend describing the organs
outer_legend("top", c("Gill", "Gonad", "Liver"), pch = c(16,
    15, 17), bty = "n", col = "darkgrey", cex = 1.75, ncol = 3,
    pt.cex = 3)

dev.off()
```

``` r
# Create the blank pdf to store the plot in
pdf("figs/Fig_PCA1v4.pdf", height = 6, width = 8)

# Set the plotting parameters
par(mar = c(4, 5, 1, 5), oma = c(2, 2, 2, 2))

# Create the plot for PC1 v PC4
plot(pca_plotting_data$PC1, pca_plotting_data$PC4, col = paste0(sex_cols[pca_plotting_data$sex],
    "75"), pch = organ_shapes[pca_plotting_data$organ], cex = 2,
    cex.lab = 2, cex.axis = 1.75, xlab = paste0("PC1: ", percentVar[1],
        "% variance"), ylab = paste0("PC4: ", percentVar[4],
        "% variance"), bty = "l", xpd = TRUE)

text(pca_plotting_data$PC1[pca_plotting_data$ID == "S34"] + 50,
    pca_plotting_data$PC4[pca_plotting_data$ID == "S34"] - 25,
    labels = pca_plotting_data$ID[pca_plotting_data$ID == "S34"])

# Add legend describing the organs and sex
outer_legend("right", c("Gill", "Gonad", "Liver"), pch = c(16,
    15, 17), bty = "n", col = "darkgrey", cex = 1.5, ncol = 1,
    pt.cex = 2.25)

outer_legend("topright", c("Female", "Male"), pch = 18, bty = "n",
    col = paste0(sex_cols, "75"), cex = 1.5, ncol = 1, pt.cex = 2.25)
dev.off()
```

#### Creating heatmaps based on the PCA axes

I now want to generate heatmaps based on the genes that have the highest
loadings for each of the PCA axes. I will be using a cut off of 0.02 to
determine which genes will be included in the heatmaps. Once again I am
doing this for PCA1-4, but will likely only inlcude PCA1-3 in a figure
for the maintext and PC4 will go in the supplemental.

``` r
df <- as.data.frame(samples[, c("sex", "organ")])
rownames(df) <- samples$ID
df$organ <- as.character(df$organ)

organ_cols <- c(Gill = "#20B2AA", Gonad = "#EEB422", Liver = "#EE8262")

ann_colors = list(sex = sex_cols, organ = organ_cols)

col_order <- c(rownames(df[df$sex == "F" & df$organ == "Gill",
    ]), rownames(df[df$sex == "M" & df$organ == "Gill", ]), rownames(df[df$sex ==
    "F" & df$organ == "Gonad", ]), rownames(df[df$sex == "M" &
    df$organ == "Gonad", ]), rownames(df[df$sex == "F" & df$organ ==
    "Liver", ]), rownames(df[df$sex == "M" & df$organ == "Liver",
    ]))

pca_rotation <- pca$rotation[, 1:4]
```

``` r
# Function to use to export the heatmaps to pdfs from
# https://stackoverflow.com/questions/43051525/how-to-draw-pheatmap-plot-to-screen-and-also-save-to-file
save_pheatmap_pdf <- function(x, filename, width = 7, height = 7) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    pdf(filename, width = width, height = height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
}
```

``` r
pc1 <- pheatmap(assay(vsd)[which(abs(pca_rotation[, 1]) >= 0.02),
    col_order], cluster_rows = FALSE, show_rownames = FALSE,
    cluster_cols = FALSE, show_colnames = FALSE, annotation_col = df,
    annotation_colors = ann_colors, cellwidth = 9, fontsize = 16,
    annotation_legend = FALSE, main = "Top loading genes on PC1")

save_pheatmap_pdf(pc1, "figs/Fig_pc1_heatmap.pdf", width = 6,
    height = 6)
```

``` r
pc2 <- pheatmap(assay(vsd)[which(abs(pca_rotation[, 2]) >= 0.02),
    col_order], cluster_rows = FALSE, show_rownames = FALSE,
    cluster_cols = FALSE, show_colnames = FALSE, annotation_col = df,
    annotation_colors = ann_colors, cellwidth = 9, fontsize = 16,
    annotation_legend = FALSE, main = "Top loading genes on PC2")

save_pheatmap_pdf(pc2, "figs/Fig_pc2_heatmap.pdf", width = 6,
    height = 6)
```

``` r
pc3 <- pheatmap(assay(vsd)[which(abs(pca_rotation[, 3]) >= 0.02),
    col_order], cluster_rows = FALSE, show_rownames = FALSE,
    cluster_cols = FALSE, show_colnames = FALSE, annotation_col = df,
    annotation_colors = ann_colors, cellwidth = 9, fontsize = 16,
    border_color = NA, main = "Top loading genes on PC3")

save_pheatmap_pdf(pc3, "figs/Fig_pc3_heatmap.pdf", width = 6,
    height = 6)
```

``` r
pc4 <- pheatmap(assay(vsd)[which(abs(pca_rotation[, 4]) >= 0.02),
    col_order], cluster_rows = FALSE, show_rownames = FALSE,
    cluster_cols = FALSE, show_colnames = FALSE, annotation_col = df,
    annotation_colors = ann_colors, cellwidth = 9, fontsize = 16,
    border_color = NA, main = "Top loading genes on PC4")

save_pheatmap_pdf(pc4, "figs/Fig_pc4_heatmap.pdf", width = 6,
    height = 6)
```

#### Making the combined figure

``` r
figPCAa <- image_ggplot(image_read_pdf("figs/Fig_PCA1v2.pdf"),
    interpolate = TRUE)
figPCAb <- image_ggplot(image_read_pdf("figs/Fig_PCA1v3.pdf"),
    interpolate = TRUE)
figPCAc <- image_ggplot(image_read_pdf("figs/Fig_pc1_heatmap.pdf"),
    interpolate = TRUE)
figPCAd <- image_ggplot(image_read_pdf("figs/Fig_pc2_heatmap.pdf"),
    interpolate = TRUE)
figPCAe <- image_ggplot(image_read_pdf("figs/Fig_pc3_heatmap.pdf"),
    interpolate = TRUE)

# make two patchworks
pcas <- figPCAa + figPCAb
hms <- figPCAc + figPCAd + figPCAe

# put the patchworks together
figPCA <- wrap_plots(pcas, hms, ncol = 1) + plot_annotation(tag_levels = "A")

ggsave("figs/FigPCA.pdf", figPCA, height = 4, width = 5)
ggsave("figs/FigPCA.png", figPCA, height = 4, width = 5)  # also save as a png
```

``` r
figPCA4 <- image_ggplot(image_read_pdf("figs/Fig_PCA1v4.pdf"),
    interpolate = TRUE)
figPCA4_hm <- image_ggplot(image_read_pdf("figs/Fig_pc4_heatmap.pdf"),
    interpolate = TRUE)

# put the patchworks together
figPCA_supp <- wrap_plots(figPCA4, figPCA4_hm, ncol = 2) + plot_annotation(tag_levels = "A")

ggsave("docs/figs/FigPCA_supp.pdf", figPCA_supp, height = 4,
    width = 8)
ggsave("docs/figs/FigPCA_supp.png", figPCA_supp, height = 4,
    width = 8)  # also save as a png
```

``` r
figPCA
```

![](DE_analysis_SNI_files/figure-gfm/showAssembledFigPCA-1.png)<!-- -->

In the single-factor analyses we have established that sample S34 is an
outlier. Here, we will remove it and re-run the PCA plots.

``` r
## Remove S34 from the dataset
dds_SN <- dds_SN[, dds_SN$ID != "S34"]

# only keeping rows that have at lead 10 reads total
keep <- rowSums(counts(dds_SN)) >= 10
dds_SN <- dds_SN[keep, ]


# Generate the expression values
dds_SN_exp <- DESeq(dds_SN)
```

    ## estimating size factors

    ## using 'avgTxLength' from assays(dds), correcting for library size

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 5422 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

``` r
# Transform the data
vsd <- vst(dds_SN_exp, blind = FALSE)
```

    ## using ntop=99267 top features by variance

![Principal components analysis reflects that most variation corresponds
to differences in expression between organs (green: gill; yellow: gonad;
pink: liver), rather than variation due to sex (circle = female;
triangle = male) or mating status (the size of points reflects an
individual’s mating success). The first axis explains 31% of variation
in the dataset, with gills and liver, two types of somatic tissues,
sitting on opposite ends of that axis. The second axis explains 20% of
the variation in gene expression, with somatic tissues on the opposite
side of the axis from the
gonads.](DE_analysis_SNI_files/figure-gfm/pca-scaled-organ2-1.png)

    ## Warning: Using size for a discrete variable is not advised.

![The same Principal Components Analysis as above but with colour now
representing sex and shape representing organ type (left) as well as the
comparison of PC1 and PC3 which shows variation in expression patterns
between males and females, particularly in the gonads, explains
approximately 10% of the variance in gene expression data, as reflected
by the third principal components
axis.](DE_analysis_SNI_files/figure-gfm/pca-scaled-sex2-1.png)

    ## Warning: Using size for a discrete variable is not advised.

![The same Principal Components Analysis as above but with colour now
representing sex and shape representing organ type (left) as well as the
comparison of PC1 and PC3 which shows variation in expression patterns
between males and females, particularly in the gonads, explains
approximately 10% of the variance in gene expression data, as reflected
by the third principal components
axis.](DE_analysis_SNI_files/figure-gfm/pca-scaled-sex2-2.png)

# Multifactor design - Comparing M v F across the diff tissue types

If we investigate the column data of our DESeq dataset we can see that
each sample has both a sex and organ type attached to it, we will be
using these two factors in our multi-factor analysis. In this
multi-factor analysis the model was run as counts \~ group, where group
included both the sex and the organ type (i.e. MLiver, FLiver, etc.).
The sample “S34” appeared to a bit of an outlier in the above analysis,
particularly clear in the sample-dist heatmap and the PCA, so it will be
removed from the analysis going forward.

``` r
# Create an additional column the has the two conditions
# combined(sex and organ type)
samples$group <- factor(paste0(samples$sex, samples$organ))

## Create a copy of the DESeq dataset
ddsMF_SN <- DESeqDataSetFromTximport(txi.nigra, colData = samples,
    design = ~group)

## Remove S34 from the dataset
ddsMF_SN <- ddsMF_SN[, ddsMF_SN$ID != "S34"]

## Filter the dataset, only keeping rows that have at least
## 10 reads total
keep <- rowSums(counts(ddsMF_SN)) >= 10  #& rowSums(counts(ddsMF_SN)) < 1e6
ddsMF_SN <- ddsMF_SN[keep, ]

# Run the differential expression analysis
ddsMF_SN_exp <- DESeq(ddsMF_SN)
resultsNames(ddsMF_SN_exp)
```

    ## [1] "Intercept"             "group_FGonad_vs_FGill" "group_FLiver_vs_FGill"
    ## [4] "group_MGill_vs_FGill"  "group_MGonad_vs_FGill" "group_MLiver_vs_FGill"

## Invesitgate the results of the differential expression

Thanks to the multi-factor analysis we can now explore differential
expression between all of the different combinations:

-   Male Liver v. Female Liver
-   Male Gill v. Female Gill
-   Male Gonad v. Female Gonad
-   All of the within sex tissue comparisons (e.g. Male Liver v. Male
    Gill, etc.)

### M-F Liver Comparisson

``` r
## Pulling out the liver M-F results with an alpha of 0.05
liver_con_res <- results(ddsMF_SN_exp, contrast = c("group",
    "MLiver", "FLiver"), alpha = 0.05)
liver_con_res$trin_geneid <- rownames(liver_con_res)
head(liver_con_res)
```

    ## log2 fold change (MLE): group MLiver vs FLiver 
    ## Wald test p-value: group MLiver vs FLiver 
    ## DataFrame with 6 rows and 7 columns
    ##                     baseMean log2FoldChange     lfcSE      stat    pvalue
    ##                    <numeric>      <numeric> <numeric> <numeric> <numeric>
    ## TRINITY_DN0_c0_g1 1935.95781      -0.459818  0.178766 -2.572181 0.0101060
    ## TRINITY_DN0_c0_g2   38.77091       0.387396  0.590144  0.656443 0.5115393
    ## TRINITY_DN0_c0_g3    3.24598       0.000000  1.878160  0.000000 1.0000000
    ## TRINITY_DN0_c1_g1   83.67742       0.207767  1.352079  0.153665 0.8778742
    ## TRINITY_DN0_c2_g1    2.14203      -0.784240  1.793596 -0.437245 0.6619339
    ## TRINITY_DN1_c0_g1 6380.51079      -0.511313  0.246805 -2.071729 0.0382907
    ##                        padj       trin_geneid
    ##                   <numeric>       <character>
    ## TRINITY_DN0_c0_g1 0.0634899 TRINITY_DN0_c0_g1
    ## TRINITY_DN0_c0_g2 0.7905284 TRINITY_DN0_c0_g2
    ## TRINITY_DN0_c0_g3        NA TRINITY_DN0_c0_g3
    ## TRINITY_DN0_c1_g1 1.0000000 TRINITY_DN0_c1_g1
    ## TRINITY_DN0_c2_g1        NA TRINITY_DN0_c2_g1
    ## TRINITY_DN1_c0_g1 0.1587244 TRINITY_DN1_c0_g1

``` r
summary(liver_con_res)
```

    ## 
    ## out of 99246 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 2096, 2.1%
    ## LFC < 0 (down)     : 1645, 1.7%
    ## outliers [1]       : 185, 0.19%
    ## low counts [2]     : 73054, 74%
    ## (mean count < 4)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

There are many criteria that have been employed previously to denote sex
bias. For this study we are classifying sex-biased genes as **genes with
a p-value \< 0.05 AND a logFC $\ge$ \|2\|**. With that criteria in the
liver there are 882 male-biased genes and 665 female-biased genes.

I have then pulled out all of the male-biased, female-biased, and
non-biased genes based on the criteria outlined above. I determined
non-biased genes as a p-vaule \> 0.05.

``` r
# Removing the rows where padj. is NA in results
liver_con_res_noNA <- liver_con_res[!is.na(liver_con_res$padj),
    ]
summary(liver_con_res_noNA)  #We can now see that there are no outliers or low counts since the NAs have been removed
```

    ## 
    ## out of 26028 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 2096, 8.1%
    ## LFC < 0 (down)     : 1645, 6.3%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 0, 0%
    ## (mean count < 4)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
# write.csv(as.data.frame(liver_con_res_noNA),
# 'data/liver_res.csv', row.names = TRUE)

# Creating a vector that contains all of the male-biased
# and female-biased genes in the liver
liver_mal_biased <- liver_con_res_noNA[which(liver_con_res_noNA$log2FoldChange >=
    2 & liver_con_res_noNA$padj <= 0.05), ]
liver_fem_biased <- liver_con_res_noNA[which(liver_con_res_noNA$log2FoldChange <=
    -2 & liver_con_res_noNA$padj <= 0.05), ]

# Creating an object that contains all of the non-biased
# genes in the liver
liver_non_biased <- liver_con_res_noNA[which(liver_con_res_noNA$padj >
    0.05), ]
```

I will be generating a table that outlines the genes that are the most
male or female biased in each organ type. To do this I need to pull the
top 50 differentially expressed genes in males or females, get the
Trinity gene IDs and then BLAST the corresponding sequences. Here I am
pulling out the Trinity gene IDs for those genes.

``` r
# Creating a subset of the results where the p-value is
# less than 0.05
liver_con_res_p05 <- liver_con_res_noNA[liver_con_res_noNA$padj <=
    0.05, ]


# Pulling the top 50 diff expressed genes in Males
top50_MLiver_trin_genes <- cbind(rownames(head(liver_con_res_p05[order(liver_con_res_p05$log2FoldChange,
    decreasing = TRUE), ], n = 50)))
# Run once to create the file and then commented out
# write.table(top50_MLiver_trin_genes,
# 'SN_maleL_top50TRgenes.txt', sep = '', quote=FALSE,
# row.names = FALSE, col.names = FALSE)

# Pulling the top 50 diff expressed genes in Females
top50_Fliver_trin_genes <- cbind(rownames(head(liver_con_res_p05[order(liver_con_res_p05$log2FoldChange,
    decreasing = FALSE), ], n = 50)))

# Run once to create the file and then commented out
# write.table(top50_Fliver_trin_genes,
# 'SN_femL_top50TRgenes.txt', sep = '', quote=FALSE,
# row.names = FALSE, col.names = FALSE)
```

``` r
# Creating a subset of the results where the p-value is
# less than 0.05
liver_con_res_p05 <- liver_con_res_noNA[liver_con_res_noNA$padj <=
    0.05, ]


# Pulling the all diff expressed genes in Males
all_MLiver_trin_genes <- cbind(rownames(head(liver_con_res_p05[order(liver_con_res_p05$log2FoldChange,
    decreasing = TRUE), ], n = 882)))
# Run once to create the file and then commented out
# write.table(all_MLiver_trin_genes,
# 'SN_maleL_allTRgenes.txt', sep = '', quote=FALSE,
# row.names = FALSE, col.names = FALSE)

# Pulling all diff expressed genes in Females
all_Fliver_trin_genes <- cbind(rownames(head(liver_con_res_p05[order(liver_con_res_p05$log2FoldChange,
    decreasing = FALSE), ], n = 665)))

# Run once to create the file and then commented out
# write.table(all_Fliver_trin_genes,
# 'SN_femL_allgenes.txt', sep = '', quote=FALSE, row.names
# = FALSE, col.names = FALSE)
```

### M-F Gill Comparisson

``` r
## Pulling out the gill M-F results
gill_con_res <- results(ddsMF_SN_exp, contrast = c("group", "MGill",
    "FGill"), alpha = 0.05)
gill_con_res$trin_geneid <- rownames(gill_con_res)
head(gill_con_res)
```

    ## log2 fold change (MLE): group MGill vs FGill 
    ## Wald test p-value: group MGill vs FGill 
    ## DataFrame with 6 rows and 7 columns
    ##                     baseMean log2FoldChange     lfcSE      stat    pvalue
    ##                    <numeric>      <numeric> <numeric> <numeric> <numeric>
    ## TRINITY_DN0_c0_g1 1935.95781     -0.0707121  0.164949 -0.428691 0.6681480
    ## TRINITY_DN0_c0_g2   38.77091      1.1627156  0.504883  2.302942 0.0212821
    ## TRINITY_DN0_c0_g3    3.24598      0.0000000  1.770737  0.000000 1.0000000
    ## TRINITY_DN0_c1_g1   83.67742     -2.0239289  1.288330 -1.570971 0.1161894
    ## TRINITY_DN0_c2_g1    2.14203     -2.2071925  1.580612 -1.396417 0.1625890
    ## TRINITY_DN1_c0_g1 6380.51079      0.1937796  0.232346  0.834012 0.4042744
    ##                        padj       trin_geneid
    ##                   <numeric>       <character>
    ## TRINITY_DN0_c0_g1  0.981874 TRINITY_DN0_c0_g1
    ## TRINITY_DN0_c0_g2  0.369603 TRINITY_DN0_c0_g2
    ## TRINITY_DN0_c0_g3        NA TRINITY_DN0_c0_g3
    ## TRINITY_DN0_c1_g1  0.682457 TRINITY_DN0_c1_g1
    ## TRINITY_DN0_c2_g1        NA TRINITY_DN0_c2_g1
    ## TRINITY_DN1_c0_g1  0.908731 TRINITY_DN1_c0_g1

``` r
summary(gill_con_res)
```

    ## 
    ## out of 99246 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 38, 0.038%
    ## LFC < 0 (down)     : 269, 0.27%
    ## outliers [1]       : 185, 0.19%
    ## low counts [2]     : 73054, 74%
    ## (mean count < 4)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

In the gills there were 22 genes that we can consider male-biased and
197 female-biased genes based on our criteria for sex-bias.

I have then pulled out all of the male-biased, female-biased, and
non-biased genes to save them into their own objects.

``` r
# Removing the rows where padj. is NA in results
gill_con_res_noNA <- gill_con_res[!is.na(gill_con_res$padj),
    ]
summary(gill_con_res_noNA)  #We can now see that there are no outliers or 
```

    ## 
    ## out of 26028 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 38, 0.15%
    ## LFC < 0 (down)     : 269, 1%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 0, 0%
    ## (mean count < 4)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
# low counts since the NAs have been removed
# write.csv(as.data.frame(gill_con_res_noNA),
# 'data/gill_res.csv', row.names = TRUE)

# Creating a vector that contains all of the male-biased
# and female-biased genes in the gills
gill_mal_biased <- gill_con_res_noNA[which(gill_con_res_noNA$log2FoldChange >=
    2 & gill_con_res_noNA$padj <= 0.05), ]
gill_fem_biased <- gill_con_res_noNA[which(gill_con_res_noNA$log2FoldChange <=
    -2 & gill_con_res_noNA$padj <= 0.05), ]

# Creating an object that contains all of the non-biased
# genes in the gills, p>0.05
gill_non_biased <- gill_con_res_noNA[which(gill_con_res_noNA$padj >
    0.05), ]
```

Here I am getting the trinity gene IDs that correspond to the top 50
differentially expressed genes for males and females. Because there are
only 27 female-biased and male-biased genes I will only use those.

``` r
# Creating a subset of results where p-value is less than
# 0.05
gill_con_res_p05 <- gill_con_res_noNA[gill_con_res_noNA$padj <=
    0.05, ]

# Pulling the top 27 diff expressed genes in Males
top27_MGill_trin_genes <- cbind(rownames(head(gill_con_res_p05[order(gill_con_res_p05$log2FoldChange,
    decreasing = TRUE), ], n = 27)))

# Run once to create the file and then commented out
# write.table(top27_MGill_trin_genes,
# 'SN_maleG_top27TRgenes.txt', sep = '', quote=FALSE,
# row.names = FALSE, col.names = FALSE)

# Pulling the top 27 diff expressed genes in Females
top27_FGill_trin_genes <- cbind(rownames(head(gill_con_res_p05[order(gill_con_res_p05$log2FoldChange,
    decreasing = FALSE), ], n = 27)))

# Run once to create the file and then commented out
# write.table(top27_FGill_trin_genes,
# 'SN_femG_top27TRgenes.txt', sep = '', quote=FALSE,
# row.names = FALSE, col.names = FALSE)
```

``` r
# Creating a subset of results where p-value is less than
# 0.05
gill_con_res_p05 <- gill_con_res_noNA[gill_con_res_noNA$padj <=
    0.05, ]

# Pulling all diff expressed genes in Males
all_MGill_trin_genes <- cbind(rownames(head(gill_con_res_p05[order(gill_con_res_p05$log2FoldChange,
    decreasing = TRUE), ], n = 22)))

# Run once to create the file and then commented out
# write.table(all_MGill_trin_genes,
# 'SN_maleG_allTRgenes.txt', sep = '', quote=FALSE,
# row.names = FALSE, col.names = FALSE)

# Pulling all diff expressed genes in Females
all_FGill_trin_genes <- cbind(rownames(head(gill_con_res_p05[order(gill_con_res_p05$log2FoldChange,
    decreasing = FALSE), ], n = 197)))

# Run once to create the file and then commented out
# write.table(all_FGill_trin_genes,
# 'SN_femG_allTRgenes.txt', sep = '', quote=FALSE,
# row.names = FALSE, col.names = FALSE)
```

### M-F Gonad Comparisson

``` r
## Pulling out the gonad M-F results
gonad_con_res <- results(ddsMF_SN_exp, contrast = c("group",
    "MGonad", "FGonad"), alpha = 0.05)
gonad_con_res$trin_geneid <- rownames(gonad_con_res)
head(gonad_con_res)
```

    ## log2 fold change (MLE): group MGonad vs FGonad 
    ## Wald test p-value: group MGonad vs FGonad 
    ## DataFrame with 6 rows and 7 columns
    ##                     baseMean log2FoldChange     lfcSE      stat      pvalue
    ##                    <numeric>      <numeric> <numeric> <numeric>   <numeric>
    ## TRINITY_DN0_c0_g1 1935.95781       0.562145  0.134873  4.167964 3.07332e-05
    ## TRINITY_DN0_c0_g2   38.77091      -0.689036  0.364354 -1.891117 5.86087e-02
    ## TRINITY_DN0_c0_g3    3.24598      -2.341826  1.034205 -2.264373 2.35512e-02
    ## TRINITY_DN0_c1_g1   83.67742       0.618702  1.055623  0.586101 5.57807e-01
    ## TRINITY_DN0_c2_g1    2.14203       3.726896  1.276757  2.919032 3.51120e-03
    ## TRINITY_DN1_c0_g1 6380.51079       0.182533  0.190367  0.958847 3.37636e-01
    ##                          padj       trin_geneid
    ##                     <numeric>       <character>
    ## TRINITY_DN0_c0_g1 0.000220537 TRINITY_DN0_c0_g1
    ## TRINITY_DN0_c0_g2 0.134441590 TRINITY_DN0_c0_g2
    ## TRINITY_DN0_c0_g3 0.066673111 TRINITY_DN0_c0_g3
    ## TRINITY_DN0_c1_g1 0.704044301 TRINITY_DN0_c1_g1
    ## TRINITY_DN0_c2_g1 0.014119614 TRINITY_DN0_c2_g1
    ## TRINITY_DN1_c0_g1 0.489693378 TRINITY_DN1_c0_g1

``` r
summary(gonad_con_res)
```

    ## 
    ## out of 99246 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 10944, 11%
    ## LFC < 0 (down)     : 5755, 5.8%
    ## outliers [1]       : 185, 0.19%
    ## low counts [2]     : 48112, 48%
    ## (mean count < 1)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

Between the ovaries and testes (gonads) there were 6705 genes that we
can consider male-biased and 2532 female-biased genes based on our
criteria for sex-bias.

From the gonad M-F results I have then pulled out all of the
male-biased, female-biased, and non-biased genes and saved them to their
own objects.

``` r
# Removing the rows where padj. is NA in results
gonad_con_res_noNA <- gonad_con_res[!is.na(gonad_con_res$padj),
    ]
summary(gonad_con_res_noNA)  #We can now see that there are no outliers or 
```

    ## 
    ## out of 50970 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 10944, 21%
    ## LFC < 0 (down)     : 5755, 11%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 0, 0%
    ## (mean count < 1)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
# low counts since the NAs have been removed
# write.csv(as.data.frame(gonad_con_res_noNA),
# 'data/gonad_res.csv', row.names = TRUE)

# Creating an object that contains all of the male-biased
# and female-biased genes in the gonads
gonad_mal_biased <- gonad_con_res_noNA[which(gonad_con_res_noNA$log2FoldChange >=
    2 & gonad_con_res_noNA$padj <= 0.05), ]
gonad_fem_biased <- gonad_con_res_noNA[which(gonad_con_res_noNA$log2FoldChange <=
    -2 & gonad_con_res_noNA$padj <= 0.05), ]

# Creating an object that contains all of the non-biased
# genes in the gonads, p>0.05
gonad_non_biased <- gonad_con_res_noNA[which(gonad_con_res_noNA$padj >
    0.05), ]
```

I then pulled out the trinity gene IDs that corresponded to the top 50
differentially expressed genes in males and females so that I could
BLAST the sequences.

``` r
# Creating a subset of the results where the p-value is
# less than 0.05
gonad_con_res_p05 <- gonad_con_res_noNA[gonad_con_res_noNA$padj <=
    0.05, ]

# Pulling the top 50 diff expressed genes in Males
top50_MGonad_trin_genes <- cbind(rownames(head(gonad_con_res_p05[order(gonad_con_res_p05$log2FoldChange,
    decreasing = TRUE), ], n = 50)))

# Run once to create the file and then commented out
# write.table(top50_MGonad_trin_genes,
# 'SN_maleGon_top50TRgenes.txt', sep = '', quote=FALSE,
# row.names = FALSE, col.names = FALSE)

# Pulling the top 50 diff expressed genes in Females
top50_FGonad_trin_genes <- cbind(rownames(head(gonad_con_res_p05[order(gonad_con_res_p05$log2FoldChange,
    decreasing = FALSE), ], n = 50)))

# Run once to create the file and then commented out
# write.table(top50_FGonad_trin_genes,
# 'SN_femGon_top50TRgenes.txt', sep = '', quote=FALSE,
# row.names = FALSE, col.names = FALSE)
```

``` r
# Creating a subset of the results where the p-value is
# less than 0.05
gonad_con_res_p05 <- gonad_con_res_noNA[gonad_con_res_noNA$padj <=
    0.05, ]

# Pulling all diff expressed genes in Males
all_MGonad_trin_genes <- cbind(rownames(head(gonad_con_res_p05[order(gonad_con_res_p05$log2FoldChange,
    decreasing = TRUE), ], n = 6705)))

# Run once to create the file and then commented out
# write.table(all_MGonad_trin_genes,
# 'SN_maleGon_allTRgenes.txt', sep = '', quote=FALSE,
# row.names = FALSE, col.names = FALSE)

# Pulling all diff expressed genes in Females
all_FGonad_trin_genes <- cbind(rownames(head(gonad_con_res_p05[order(gonad_con_res_p05$log2FoldChange,
    decreasing = FALSE), ], n = 2532)))

# Run once to create the file and then commented out
# write.table(all_FGonad_trin_genes,
# 'SN_femGon_allTRgenes.txt', sep = '', quote=FALSE,
# row.names = FALSE, col.names = FALSE)
```

### MA Plots

![MA-plots generated for each organ type that compares logFC to the mean
expression. Female-biased and Male-biased genes are represented by color
(green and purple correspondingly) and are determing with a fc cutoff of
2 and a p-value cutoff of
0.05.](DE_analysis_SNI_files/figure-gfm/MA-plot-MF-1.png)

From the MA plots we can see the highest number of biased genes are
found in the gonads (Fig. @ref(fig:MA-plot-MF)) and then the liver and
lastly the gills. The liver appears to be the only organ where there are
more female-biased genes than male-biased genes, which is initially
surprising due to the sex-role-reversed nature of these fish.

## Creating an Upset Plot

![Upset plots to show the number of shared sex-biased genes across the
organs in males.](DE_analysis_SNI_files/figure-gfm/upset-plot-mal-1.png)

![Upset plots to show the number of shared sex-biased genes across the
organs in
females.](DE_analysis_SNI_files/figure-gfm/upset-plot-fem-1.png)

![Upset plots to show the number of shared sex-biased genes across the
organs across both
sexes.](DE_analysis_SNI_files/figure-gfm/upset-plot-all-1.png)

In both males and females the highest number of shared sex-biased genes
is found between the gonads and the liver (Fig.
@ref(fig:upset-plot-mal), @ref(fig:upset-plot-fem)). There are no genes
that are shared across the same organ between males and females (e.g.,
biased in males and females in the liver) which is good (Fig.
@ref(fig:upset-plot-all)). Two sex-biased genes are female-biased in all
organs and three genes are male-biased in all three organs (Fig.
@ref(fig:upset-plot-mal), @ref(fig:upset-plot-fem)).

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
logFC_long<- data.frame(
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
         rep("UB", nrow(gill_non_biased)),
         rep("FB", nrow(gonad_fem_biased)),
         rep("MB", nrow(gonad_mal_biased)),
         rep("UB", nrow(gonad_non_biased)),
         rep("FB", nrow(liver_fem_biased)),
         rep("MB", nrow(liver_mal_biased)),
         rep("UB", nrow(liver_non_biased))
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

With the dataset in the proper format now, we can generate the plot

![Absolute value of the logFC for genes that are female-biased,
male-biased, and unbiased across each tissue type. Raw fold-change
values are added on top of the boxplot as
jitters.](DE_analysis_SNI_files/figure-gfm/FC-var-plot-1.png)

There appears to be differences in the logFC across both the organs and
the sex-bias groups (male-biased, fem-biased, and non-biased) (Fig.
@ref(fig:FC-var-plot)). Let’s explore some analyses to see if we can put
some statistical power behind those differences.

``` r
# Looking at some summary statistics for logFC between the
# different groups
tapply(abs(logFC_long$logFC), list(logFC_long$bias, logFC_long$tissue),
    mean)
tapply(abs(logFC_long$logFC), list(logFC_long$bias, logFC_long$tissue),
    sd)

par(mfrow = c(1, 1))
interaction.plot(logFC_long$tissue, logFC_long$bias, logFC_long$logFC)

logFC_long$tissue <- as.factor(logFC_long$tissue)
logFC_long$bias <- as.factor(logFC_long$bias)

fc_var_aov <- aov(abs(logFC_long$logFC) ~ logFC_long$bias * logFC_long$tissue)
anova(fc_var_aov)

par(mfrow = c(2, 2))
plot(fc_var_aov)

library(car)
leveneTest(fc_var_aov)
par(mfrow = c(1, 1))
hist(resid(fc_var_aov))
# shapiro.test(resid(fc_var_aov))

plot(TukeyHSD(fc_var_aov))
library(multcompView)
multcompLetters4(fc_var_aov, TukeyHSD(fc_var_aov))

## Fitting a negative binomial model
library(MASS)
model <- glm.nb(abs(logFC_long$logFC) ~ logFC_long$tissue * logFC_long$bias)

linear_model <- lm(abs(logFC_long$logFC) ~ logFC_long$tissue *
    logFC_long$bias)
```

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
# Pulling out the geneIDs for genes that were categorized
# as 'outliers' by DESeq2 Calculating the Cooks threshold
# that would have been used
np <- length(resultsNames(ddsMF_SN_exp))
nsamp <- ncol(ddsMF_SN_exp)
cooks_thresh <- qf(0.99, df1 = np, df2 = nsamp - np)

out_ids <- names(mcols(ddsMF_SN_exp)$maxCooks[mcols(ddsMF_SN_exp)$maxCooks >
    cooks_thresh])

# Checking to make sure we have correct number of outliers
length(out_ids)
```

    ## [1] 185

``` r
# Filtering out the dds dataset to remove the outliers
# determined by DESeq
ddsMF_SN_exp_filtered <- ddsMF_SN_exp[!(rownames(ddsMF_SN_exp) %in%
    out_ids), ]

# Create a vector with the different organ types
organs <- levels(colData(ddsMF_SN_exp_filtered)$organ)

# Create an empty list to store my datasets in
SN_sex_specific_genes <- list()

# Generate the for loop to identify MSpecific and FSpecific
# genes in each organ based on Medians
for (organ in organs) {

    # Male-Specific Genes Pull out all of the rows where
    # fem count <=10 in every female sample
    fem0_organ_names <- which(rowSums(t(apply(counts(ddsMF_SN_exp_filtered,
        normalized = TRUE)[, ddsMF_SN_exp_filtered$sex == "F" &
        ddsMF_SN_exp_filtered$organ == organ], 1, function(x) x <=
        10))) == ncol(counts(ddsMF_SN_exp_filtered, normalized = TRUE)[,
        ddsMF_SN_exp_filtered$sex == "F" & ddsMF_SN_exp_filtered$organ ==
            organ]))

    fem0_organ <- counts(ddsMF_SN_exp_filtered, normalized = TRUE)[rownames(counts(ddsMF_SN_exp_filtered,
        normalized = TRUE)) %in% names(fem0_organ_names), ddsMF_SN_exp_filtered$organ ==
        organ]

    ## Pull out rows where median of male count >=20 for
    ## that organ
    mal10_organ <- apply(counts(ddsMF_SN_exp_filtered, normalized = TRUE)[,
        ddsMF_SN_exp_filtered$sex == "M" & ddsMF_SN_exp_filtered$organ ==
            organ], 1, function(x) median(x) >= 20)

    ## Keep only the rows where all fem samples <=10 and
    ## the Male median>=20
    fem0_mal10_organ <- fem0_organ[rownames(fem0_organ) %in%
        names(mal10_organ[mal10_organ == TRUE]), ]

    ## Create a new object with a name based on the organ
    ## type
    organ_malsp <- sub("$", "_male_specific", organ)
    SN_sex_specific_genes[[organ_malsp]] <- fem0_mal10_organ

    # Female-Specific Genes Pull out all of the rows where
    # male count <=10 in every male sample
    mal0_organ_names <- which(rowSums(t(apply(counts(ddsMF_SN_exp_filtered,
        normalized = TRUE)[, ddsMF_SN_exp_filtered$sex == "M" &
        ddsMF_SN_exp_filtered$organ == organ], 1, function(x) x <=
        10))) == ncol(counts(ddsMF_SN_exp_filtered, normalized = TRUE)[,
        ddsMF_SN_exp_filtered$sex == "M" & ddsMF_SN_exp_filtered$organ ==
            organ]))

    mal0_organ <- counts(ddsMF_SN_exp_filtered, normalized = TRUE)[rownames(counts(ddsMF_SN_exp_filtered,
        normalized = TRUE)) %in% names(mal0_organ_names), ddsMF_SN_exp_filtered$organ ==
        organ]

    ## Pull out rows where median of female count >=10 for
    ## that organ
    fem10_organ <- apply(counts(ddsMF_SN_exp_filtered, normalized = TRUE)[,
        ddsMF_SN_exp_filtered$sex == "F" & ddsMF_SN_exp_filtered$organ ==
            organ], 1, function(x) median(x) >= 20)

    # Keep only the rows where male=0 and the fem
    # median>=10
    mal0_fem10_organ <- mal0_organ[rownames(mal0_organ) %in%
        names(fem10_organ[fem10_organ == TRUE]), ]

    # Create a new object with a name based on the organ
    # type
    organ_femsp <- sub("$", "_female_specific", organ)
    SN_sex_specific_genes[[organ_femsp]] <- mal0_fem10_organ

}
```

### Investigating the results of the sex-specific subsetting

Let’s now take a look at how many sex-specific genes we have in each
tissue type:
<table class="table table-striped table-hover table-condensed table-responsive" style="color: black; width: auto !important; ">
<caption>
Sex-specific Genes
</caption>
<thead>
<tr>
<th style="text-align:left;">
Tissue Type
</th>
<th style="text-align:left;">
Sex
</th>
<th style="text-align:right;">
Number of Sex-specific Genes
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Gill
</td>
<td style="text-align:left;">
F
</td>
<td style="text-align:right;">
19
</td>
</tr>
<tr>
<td style="text-align:left;">
Gill
</td>
<td style="text-align:left;">
M
</td>
<td style="text-align:right;">
7
</td>
</tr>
<tr>
<td style="text-align:left;">
Gonad
</td>
<td style="text-align:left;">
F
</td>
<td style="text-align:right;">
318
</td>
</tr>
<tr>
<td style="text-align:left;">
Gonad
</td>
<td style="text-align:left;">
M
</td>
<td style="text-align:right;">
949
</td>
</tr>
<tr>
<td style="text-align:left;">
Liver
</td>
<td style="text-align:left;">
F
</td>
<td style="text-align:right;">
64
</td>
</tr>
<tr>
<td style="text-align:left;">
Liver
</td>
<td style="text-align:left;">
M
</td>
<td style="text-align:right;">
187
</td>
</tr>
</tbody>
</table>

To get a better idea what is going on between the sex-specific genes vs
the sex-biased genes I first determined how many overlaps there were
between our sex-biased and sex-specific genes and then used plotCounts
to plot the counts for each gene separately.

In the **Male Gills** there are 0 genes that overlap between sex-biased
and sex-specific, there are 837 overlapping in the **gonads** and 122
overlapping in the **liver**. For females we have 13 overlapping
sex-specific and sex-biased genes in the **gills**, 315 genes in the
**gonads**, and 51 genes in the **liver**.

``` r
## Female sex-specific
fem_sp_geneIDs <- rownames(SN_sex_specific_genes$Gill_female_specific)

## Set the plotting window
par(mfrow = c(1, 1))

for (gene in fem_sp_geneIDs) {

    plotCounts(ddsMF_SN[, ddsMF_SN$organ == "Gill"], gene = gene,
        intgroup = "sex")

}

## Female sex-biased
gill_fem_biased_geneIDS <- rownames(gill_fem_biased)

## Set the plotting window
par(mfrow = c(1, 1))

for (gene in gill_fem_biased_geneIDS) {

    plotCounts(ddsMF_SN[, ddsMF_SN$organ == "Gill"], gene = gene,
        intgroup = "sex", col = "salmon", pch = 19)

}
```

### Pulling out the gene IDs for all of the sex-specific genes

``` r
# Run once to create the file and then commented out
# write.table(cbind(rownames(SN_sex_specific_genes$Gill_male_specific)),
# 'SN_malG_specific_TRgenes.txt', sep = '', quote=FALSE,
# row.names = FALSE, col.names = FALSE)
# write.table(cbind(rownames(SN_sex_specific_genes$Gill_female_specific)),
# 'SN_femG_specific_TRgenes.txt', sep = '', quote=FALSE,
# row.names = FALSE, col.names = FALSE)
# write.table(cbind(rownames(SN_sex_specific_genes$Gonad_male_specific)),
# 'SN_malGon_specific_TRgenes.txt', sep = '', quote=FALSE,
# row.names = FALSE, col.names = FALSE)
# write.table(cbind(rownames(SN_sex_specific_genes$Gonad_female_specific)),
# 'SN_femGon_specific_TRgenes.txt', sep = '', quote=FALSE,
# row.names = FALSE, col.names = FALSE)
# write.table(cbind(rownames(SN_sex_specific_genes$Liver_male_specific)),
# 'SN_malL_specific_TRgenes.txt', sep = '', quote=FALSE,
# row.names = FALSE, col.names = FALSE)
# write.table(cbind(rownames(SN_sex_specific_genes$Liver_female_specific)),
# 'SN_femL_specific_TRgenes.txt', sep = '', quote=FALSE,
# row.names = FALSE, col.names = FALSE)
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
sex and an absence in the other sex.

``` r
# Make a vector that contains all of the groupings
biased_bins <- c("Unbiased", "Low", "Med", "High", "Extreme",
    "Sex-specific")

# Create a new column in the dataset and use ifelse
# statements to set the category limits abs(logFC) was used
# to account for the fem-biased genes possessing negative
# values
logFC_long$bias_cat <- ifelse(logFC_long$bias == "UB", biased_bins[1],
    ifelse(abs(logFC_long$logFC) >= 2 & abs(logFC_long$logFC) <
        3, biased_bins[2], ifelse(abs(logFC_long$logFC) >= 3 &
        abs(logFC_long$logFC) < 5, biased_bins[3], ifelse(abs(logFC_long$logFC) >=
        5 & abs(logFC_long$logFC) < 10, biased_bins[4], biased_bins[5]))))

# Making sure the genes we categorized as sex-specific are
# labeled as sex-specific for their bias cat. in the
# dataset
organs <- c("Gill", "Gill", "Gonad", "Gonad", "Liver", "Liver")
bias <- c("MB", "FB", "MB", "FB", "MB", "FB")

for (i in 1:length(SN_sex_specific_genes)) {

    tmp <- SN_sex_specific_genes[[i]]
    tmp <- as.data.frame(tmp)
    tmp$geneID <- rownames(tmp)

    for (j in 1:nrow(tmp)) {

        one_gene <- tmp[j, ]

        if (one_gene[["geneID"]] %in% logFC_long[logFC_long$tissue ==
            organs[[i]] & logFC_long$bias == bias[[i]], "geneID"]) {

            logFC_long[logFC_long$geneID == one_gene[["geneID"]] &
                logFC_long$tissue == organs[[i]] & logFC_long$bias ==
                bias[[i]], "bias_cat"] <- "Sex-specific"
        } else {

            one_gene_dat <- data.frame(matrix(ncol = ncol(logFC_long),
                nrow = 1))
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

# Create subset of our long dataset that does not include
# the non-biased genes
logFC_long_noNB <- logFC_long[logFC_long$bias_cat != "Unbiased",
    ]

# Make a table to count the number of genes in each
# category for each organ
bias_cat_table <- table(logFC_long_noNB$bias, logFC_long_noNB$bias_cat,
    logFC_long_noNB$tissue)

# Add the counts of sex-specific genes to the table
bias_cat_gill <- bias_cat_table[, , "Gill"]
# write.table(bias_cat_gill, 'data/bias_cat_gills.txt',
# row.names = TRUE)

bias_cat_gonad <- bias_cat_table[, , "Gonad"]
# write.table(bias_cat_gonad, 'data/bias_cat_gonad.txt',
# row.names = TRUE)

bias_cat_liver <- bias_cat_table[, , "Liver"]
# write.table(bias_cat_liver, 'data/bias_cat_liver.txt',
# row.names = TRUE)

# Plot the counts for each tissue type
par(mfrow = c(1, 3))
barplot(bias_cat_gill[, biased_bins[biased_bins != "Unbiased"]],
    beside = TRUE, ylim = c(0, max(bias_cat_table) + 1), col = c("#7fc97f",
        "#beaed4"), main = "Gill", ylab = "Number of Genes",
    xlab = "Bias Category")
barplot(bias_cat_gonad[, biased_bins[biased_bins != "Unbiased"]],
    beside = TRUE, ylim = c(0, max(bias_cat_table) + 1), col = c("#7fc97f",
        "#beaed4"), main = "Gonad")
barplot(bias_cat_liver[, biased_bins[biased_bins != "Unbiased"]],
    beside = TRUE, ylim = c(0, max(bias_cat_table) + 1), col = c("#7fc97f",
        "#beaed4"), main = "Liver")
legend("topright", legend = c("Female-biased", "Male-biased"),
    fill = c("#7fc97f", "#beaed4"))
```

![](DE_analysis_SNI_files/figure-gfm/bin-sex-bias-1.png)<!-- -->

## Gene Ontology Analysis

We wanted to figure out two things:

1.  What are the genes that are sex-biased within each tissue type  
2.  What are the GO terms for those genes.

To accomplish this I first had to pull out the sequences that
corresponded to the Trinity gene IDS pulled out above. Once those
sequences were fetched they were blasted against the *Corythoichthys
intestinalis* genome. This was completed with the following bash script
on the RCC.

### Using BLAST in command line

``` bash
#Create a BLAST database - Used Corythoichthys intestinalis for our BLAST database
makeblastdb -in /home/rccuser/shared/emily_files/S_nigra_2024/c_intestinalis_data/data/GCF_030265065.1/GCF_030265065.1_ASM3026506v1_genomic.fna -out c_intestinalis_genome -dbtype nucl
```

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
        blastn -db $blast_database -query ${sample}.fasta -out ${sample}_corythoicthys_blast.txt \
                -evalue 0.001 \
                -num_threads 12 \
                -outfmt "6 qseqid qstart qend stitle salltitles sseqid sallseqid sacc sblastnames sstart send evalue bitscore length pident gaps"       

        #Move the outputs into desired output directory
        mv ${sample}* $output_dir

done
```

This script was run as
`bash ../scripts/blast_pipeline.sh nigra_high_exp_genes/TRgenes/ ../scripts/subset_fasta_file nigra_supertranscript.fasta c_intestinalis_genome/c_intestinalis_genome nigra_high_exp_genes/high_exp_results`

Once the BLAST results were obtained in the RCC I read them into R for
further filtering (keeping only one hit per gene)

``` r
# Specify the directory where BLAST results are located
blast_path <- "data/cory_BLAST/"

# Create a list of the files I want
SN_blast_files <- list.files(blast_path, pattern = "cory")

# Create an empty list to store my datasets in
SN_blast_list <- list()

# Create a loop to read in all of the blast results
for (file_name in SN_blast_files) {
    # Read the file and store it in a new object
    file_data <- read.delim(file.path("data/cory_BLAST", file_name),
        stringsAsFactors = FALSE, quote = "", sep = "\t", header = FALSE)

    # Add column names to the dataset
    colnames(file_data) <- c("trin_geneid", "trin_gene_start",
        "trin_gene_end", "cory_gene_info", "cory_gene_start
                           cory_gene_end",
        "evalue", "bit_score", "length", "pident", "gaps")

    # Create a new object with a name based on the file
    # name
    blast_name <- sub(".txt$", "", file_name)  #Removes the file extension
    SN_blast_list[[blast_name]] <- file_data
}
```

### Filtering Blast results

Because of the way BLAST was run in the RCC, we can have multiple hits
for one gene, I want to filter that out so that there is only one hit
kept for each gene.

``` r
# Use lapply to apply the function to each dataset stored
# in the list created above
blast_output_filtered <- lapply(SN_blast_list, function(data_frame) {

    # Pull out the Unique Trinity gene IDs
    uniqueID <- unique(data_frame[1])

    # Create an empty dataframe to store intermediate
    # results in
    output <- data.frame(matrix(data = NA, ncol = ncol(data_frame),
        nrow = 0))
    colnames(output) <- c("trin_geneid", "trin_gene_start", "trin_gene_end",
        "cory_gene_info", "cory_gene_start", "cory_gene_end",
        "evalue", "bit_score", "length", "pident", "gaps")

    # Generate a for loop that pulls out the lowest e-value
    # + highest % identity for each gene
    for (gene in uniqueID$trin_geneid) {

        # Subset the dataset for each Trinity gene
        this_trin_gene <- subset(data_frame, trin_geneid == gene)
        # Pull out gene with smallest e-value
        uni_gene_subset_lowe <- this_trin_gene[which(this_trin_gene$evalue ==
            min(this_trin_gene$evalue)), ]
        # In case mult. genes have same e-value, pull out
        # gene with highest % identity
        uni_gene_subset_lowe_highpid <- uni_gene_subset_lowe[which(uni_gene_subset_lowe$pident ==
            max(uni_gene_subset_lowe$pident)), ]

        # keep only one of multiple identical rows
        uni_gene_subset_lowe_highpid <- unique(uni_gene_subset_lowe_highpid)

        # check that there is a single gene ID in cory, if
        # not, only first row is kept
        if (length(gsub("^.*(GeneID:\\d+)\\].*$", "\\1", uni_gene_subset_lowe_highpid$cory_gene_info)) >
            1) {
            # browser() } else{
            uni_gene_subset_lowe_highpid <- uni_gene_subset_lowe_highpid[1,
                ]
        }

        # Add the final gene into the empty dataframe from
        # above
        output <- rbind(output, uni_gene_subset_lowe_highpid)
    }

    return(output)
})
```

### Merging the BLAST filtered Data frames with the sex-biased information

Not all of the trinity genes had a BLAST hit come back, but I still want
to include what those genes are in the dataset. To do this I am merging
all of the top 50 genes with the filtered BLAST dataset and if there are
genes that did not have any hits (i.e. are not in the filtered BLAST
database) then they will receive an “NA”.

All genes that had an “NA” come back from blasting against the *C.
intestinalis* database will be blasted against the non-redundant BLAST
database as well.

``` r
# Getting all diff. expressed genes in Males
all_MGonad <- head(gonad_con_res_p05[order(gonad_con_res_p05$log2FoldChange,
    decreasing = TRUE), ], n = 6705)
all_MGill <- head(gill_con_res_p05[order(gill_con_res_p05$log2FoldChange,
    decreasing = TRUE), ], n = 22)
all_MLiver <- head(liver_con_res_p05[order(liver_con_res_p05$log2FoldChange,
    decreasing = TRUE), ], n = 882)

# Pulling all diff expressed genes in Females
all_FGonad <- head(gonad_con_res_p05[order(gonad_con_res_p05$log2FoldChange,
    decreasing = FALSE), ], n = 2532)
all_FGill <- head(gill_con_res_p05[order(gill_con_res_p05$log2FoldChange,
    decreasing = FALSE), ], n = 197)
all_FLiver <- head(liver_con_res_p05[order(liver_con_res_p05$log2FoldChange,
    decreasing = FALSE), ], n = 665)
# Create dataframes for all of the above that contains the
# trinity gene ID and the logFoldChange
all_MGonad_df <- as.data.frame(cbind(trin_geneid = all_MGonad$trin_geneid,
    LogFC = all_MGonad$log2FoldChange))
all_MGill_df <- as.data.frame(cbind(trin_geneid = all_MGill$trin_geneid,
    LogFC = all_MGill$log2FoldChange))
all_MLiver_df <- as.data.frame(cbind(trin_geneid = all_MLiver$trin_geneid,
    LogFC = all_MLiver$log2FoldChange))
all_FGonad_df <- as.data.frame(cbind(trin_geneid = all_FGonad$trin_geneid,
    LogFC = all_FGonad$log2FoldChange))
all_FGill_df <- as.data.frame(cbind(trin_geneid = all_FGill$trin_geneid,
    LogFC = all_FGill$log2FoldChange))
all_FLiver_df <- as.data.frame(cbind(trin_geneid = all_FLiver$trin_geneid,
    LogFC = all_FLiver$log2FoldChange))

# Merge the above dataframes with the corresponding
# blast_output_filtered dataframe
blast_output_filtered$SN_femG_allTRgenes_blast <- merge(all_FGill_df,
    blast_output_filtered$SN_femG_allTRgenes_corythoicthys_blast,
    by = "trin_geneid", all.x = TRUE, all.y = TRUE)
blast_output_filtered$SN_femGon_allTRgenes_blast <- merge(all_FGonad_df,
    blast_output_filtered$SN_femL_allTRgenes_corythoicthys_blast,
    by = "trin_geneid", all.x = TRUE, all.y = TRUE)
blast_output_filtered$SN_femL_allTRgenes_blast <- merge(all_FLiver_df,
    blast_output_filtered$SN_femL_allTRgenes_corythoicthys_blast,
    by = "trin_geneid", all.x = TRUE, all.y = TRUE)
blast_output_filtered$SN_maleG_allTRgenes_blast <- merge(all_MGill_df,
    blast_output_filtered$SN_maleG_allTRgenes_corythoicthys_blast,
    by = "trin_geneid", all.x = TRUE, all.y = TRUE)
blast_output_filtered$SN_maleGon_allTRgenes_blast <- merge(all_MGonad_df,
    blast_output_filtered$SN_maleGon_allTRgenes_corythoicthys_blast,
    by = "trin_geneid", all.x = TRUE, all.y = TRUE)
blast_output_filtered$SN_maleL_allTRgenes_blast <- merge(all_MLiver_df,
    blast_output_filtered$SN_maleL_allTRgenes_corythoicthys_blast,
    by = "trin_geneid", all.x = TRUE, all.y = TRUE)
```

### Blasting the sex-specific genes against the *C. intestinalis* genome

The trinity genes that were pulled out for the sex-specific genes above
were blasted with `blastn` against the *C. intestinalis* genome in the
RCC. I then need to load those BLAST results into R, filter them, and
pull out what the *C. intestinalis* gene names are that correspond to
them.

First, I read in all of the BLAST results for the sex-specific genes in
each organ.

``` r
# Specify the directory where BLAST results are located
blast_path <- "data/cory_BLAST/"

# Create a list of the files I want
SN_blast_files <- list.files(blast_path, pattern = "cory")

# Create an empty list to store my datasets in
SN_blast_list_sex_specific <- list()

# Create a loop to read in all of the blast results
for (file_name in SN_blast_files) {
    # Read the file and store it in a new object
    file_data <- read.delim(file.path(blast_path, file_name),
        stringsAsFactors = FALSE, quote = "", sep = "\t", header = FALSE)

    file_data <- file_data[, c(1:4, 10:16)]

    # Add column names to the dataset
    colnames(file_data) <- c("trin_geneid", "trin_gene_start",
        "trin_gene_end", "cory_gene_info", "cory_gene_start",
        "cory_gene_end", "evalue", "bit_score", "length", "pident",
        "gaps")

    # Create a new object with a name based on the file
    # name
    blast_name <- sub(".txt$", "", file_name)  #Removes the file extension
    SN_blast_list_sex_specific[[blast_name]] <- file_data
}
```

I then need to filter to make sure each trinity gene only has one hit in
the same way I did above.

``` r
# Use lapply to apply the function to each dataset stored
# in the list created above
blast_output_filtered_SS <- lapply(SN_blast_list_sex_specific,
    function(data_frame) {

        # Pull out the Unique Trinity gene IDs
        uniqueID <- unique(data_frame[1])

        # Create an empty dataframe to store intermediate
        # results in
        output <- data.frame(matrix(data = NA, ncol = ncol(data_frame),
            nrow = 0))
        colnames(output) <- c("trin_geneid", "trin_gene_start",
            "trin_gene_end", "cory_gene_info", "cory_gene_start",
            "cory_gene_end", "evalue", "bit_score", "length",
            "pident", "gaps")

        # Generate a for loop that pulls out the lowest
        # e-value + highest % identity for each gene
        for (gene in uniqueID$trin_geneid) {

            # Subset the dataset for each Trinity gene
            this_trin_gene <- subset(data_frame, trin_geneid ==
                gene)
            # Pull out gene with smallest e-value
            uni_gene_subset_lowe <- this_trin_gene[which(this_trin_gene$evalue ==
                min(this_trin_gene$evalue)), ]
            # In case mult. genes have same e-value, pull
            # out gene with highest % identity
            uni_gene_subset_lowe_highpid <- uni_gene_subset_lowe[which(uni_gene_subset_lowe$pident ==
                max(uni_gene_subset_lowe$pident)), ]

            # keep only one of multiple identical rows
            uni_gene_subset_lowe_highpid <- unique(uni_gene_subset_lowe_highpid)

            # check that there is a single gene ID in
            # scovelli, if not, only first row is kept
            if (length(gsub("^.*(GeneID:\\d+)\\].*$", "\\1",
                uni_gene_subset_lowe_highpid$scov_gene_info)) >
                1) {

                uni_gene_subset_lowe_highpid <- uni_gene_subset_lowe_highpid[1,
                  ]
            }

            # Add the final gene into the empty dataframe
            # from above
            output <- rbind(output, uni_gene_subset_lowe_highpid)
        }

        return(output)
    })
```

Now I am going to pull out the gene name and protein ID from the
“cory_gene_info” column.

``` r
for (i in 1:length(blast_output_filtered_SS)) {

    # Pulling out the gene name and adding it to a new
    # column
    blast_output_filtered_SS[[i]]$gene <- gsub("^(.*)(gene=)(\\w+\\d*)(.*$)",
        "\\3", blast_output_filtered_SS[[i]]$cory_gene_info)

    # Pulling out the protein name and adding it to a new
    # column
    blast_output_filtered_SS[[i]]$protein <- gsub("^(.*)(protein=)(.*)([[:graph:]]\\s[[:graph:]])(protein_id=)(.*$)",
        "\\3", blast_output_filtered_SS[[i]]$cory_gene_info)


}
```

``` r
# Create a long format dataset that has the tissue
# information, BLAST gene and protein, geneID, and sex
blast_long_SS <- data.frame(tissue = c(rep("Gill", sum(nrow(blast_output_filtered_SS$SN_femG_specific_TRgenes_corythoicthys_blast),
    nrow(blast_output_filtered_SS$SN_malG_specific_TRgenes_corythoicthys_blast))),
    rep("Gonad", sum(nrow(blast_output_filtered_SS$SN_femGon_specific_TRgenes_corythoicthys_blast),
        nrow(blast_output_filtered_SS$SN_malGon_specific_TRgenes_corythoicthys_blast))),
    rep("Liver", sum(nrow(blast_output_filtered_SS$SN_femL_specific_TRgenes_corythoicthys_blast),
        nrow(blast_output_filtered_SS$SN_malL_specific_TRgenes_corythoicthys_blast)))),
    sex = c(rep("Female", nrow(blast_output_filtered_SS$SN_femG_specific_TRgenes_corythoicthys_blast)),
        rep("Male", nrow(blast_output_filtered_SS$SN_malG_specific_TRgenes_corythoicthys_blast)),
        rep("Female", nrow(blast_output_filtered_SS$SN_femGon_specific_TRgenes_corythoicthys_blast)),
        rep("Male", nrow(blast_output_filtered_SS$SN_malGon_specific_TRgenes_corythoicthys_blast)),
        rep("Female", nrow(blast_output_filtered_SS$SN_femL_specific_TRgenes_corythoicthys_blast)),
        rep("Male", nrow(blast_output_filtered_SS$SN_malL_specific_TRgenes_corythoicthys_blast))),
    gene = c(blast_output_filtered_SS$SN_femG_specific_TRgenes_corythoicthys_blast$gene,
        blast_output_filtered_SS$SN_malG_specific_TRgenes_corythoicthys_blast$gene,
        blast_output_filtered_SS$SN_femGon_specific_TRgenes_corythoicthys_blast$gene,
        blast_output_filtered_SS$SN_malGon_specific_TRgenes_corythoicthys_blast$gene,
        blast_output_filtered_SS$SN_femL_specific_TRgenes_corythoicthys_blast$gene,
        blast_output_filtered_SS$SN_malL_specific_TRgenes_corythoicthys_blast$gene),
    protein = c(blast_output_filtered_SS$SN_femG_specific_TRgenes_corythoicthys_blast$protein,
        blast_output_filtered_SS$SN_malG_specific_TRgenes_corythoicthys_blast$protein,
        blast_output_filtered_SS$SN_femGon_specific_TRgenes_corythoicthys_blast$protein,
        blast_output_filtered_SS$SN_malGon_specific_TRgenes_corythoicthys_blast$protein,
        blast_output_filtered_SS$SN_femL_specific_TRgenes_corythoicthys_blast$protein,
        blast_output_filtered_SS$SN_malL_specific_TRgenes_corythoicthys_blast$protein),
    trinityID = c(blast_output_filtered_SS$SN_femG_specific_TRgenes_corythoicthys_blast$trin_geneid,
        blast_output_filtered_SS$SN_malG_specific_TRgenes_corythoicthys_blast$trin_geneid,
        blast_output_filtered_SS$SN_femGon_specific_TRgenes_corythoicthys_blast$trin_geneid,
        blast_output_filtered_SS$SN_malGon_specific_TRgenes_corythoicthys_blast$trin_geneid,
        blast_output_filtered_SS$SN_femL_specific_TRgenes_corythoicthys_blast$trin_geneid,
        blast_output_filtered_SS$SN_malL_specific_TRgenes_corythoicthys_blast$trin_geneid))
# write csv once, then comment out write.csv(blast_long_SS,
# 'SN_BLAST_results_sex_specific.csv', row.names = FALSE)
```

### Blasting against the Zebrafish genome

In order for us to be able to do the GO analysis, we need to BLAST the
sequences against one of the PANTHER organisms. The zebrafish, *Danio
reiro* was the only fish in that list.

A BLAST database was generated with the *D. reiro* proteome as:
`makeblastdb -in ../ncbi_dataset/data/GCF_000002035.6/protein.faa -out d_rerio_prot -dbtyp prot`

The trinity gene IDs that correspond to ALL of the female- or
male-biased genes were pulled out and saved to `.txt` files.

``` r
write.table(cbind(rownames(gill_fem_biased)), "SN_femG_biased_TRgenes.txt",
    sep = "", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(cbind(rownames(gill_mal_biased)), "SN_malG_biased_TRgenes.txt",
    sep = "", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(cbind(rownames(gonad_fem_biased)), "SN_femGon_biased_TRgenes.txt",
    sep = "", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(cbind(rownames(gonad_mal_biased)), "SN_malGon_biased_TRgenes.txt",
    sep = "", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(cbind(rownames(liver_fem_biased)), "SN_femL_biased_TRgenes.txt",
    sep = "", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(cbind(rownames(liver_mal_biased)), "SN_malL_biased_TRgenes.txt",
    sep = "", quote = FALSE, row.names = FALSE, col.names = FALSE)
```

The `blast_pipeline.sh` script was then used like above to pull out all
of the corresponding sequences for each Trinity gene and then to get the
BLAST results. The big difference here was that `blastx` was used
instead of `blastn` because we are searching nucleotides against a
protein database. Once we have the BLAST results we can read and filter
them into R as we did before.

In order to run PANTHER we need IDs that follow one of their supported
formats. From NCBI, the gene name works (e.g. rpl27). That information
is not automatically supplied with the BLAST output, but we can use the
.gff file from NCBI to pull out the gene names that correspond to the
proteinID (e.g. NP_956018.1).

``` r
# Read in the GFF file for zebrafish that has info about
# the geneID
reiro_gff <- read.delim("data/danio_BLAST/d_reiro_genomic.gff",
    header = FALSE, comment.char = "#")

# Keep only the columns I want and rename them
reiro_gff <- reiro_gff[, c(1:5, 9)]
colnames(reiro_gff) <- c("seqname", "source", "feature", "start",
    "end", "gene_info")

# Subset the dataset to only include the rows we're
# interested in
genes <- reiro_gff[reiro_gff$feature == "CDS", ]

# Pull out the gene name that corresponds to each protein
# ID
genes$gene_name <- gsub("^(.*;)(gene=)(\\w+\\d*);(.*$)", "\\3",
    genes$gene_info)
genes$prot_id <- gsub("^(.*;)(protein_id=)(.*)$", "\\3", genes$gene_info)

# Merge the gene names pulled out above with the rest of
# the BLAST datasets based on the proteinID
blast_output_merged <- lapply(blast_output_filtered, function(dataframe) {

    output <- merge(dataframe, unique(genes[, 7:8]), by = "prot_id")

    return(output)
})
```

I am now going to make sure the sex bias information is included with
the BLAST results so we can perform statistical overrepresentation and
enrichment tests with PANTHER.

``` r
blast_output_merged$SN_femG_allTRgenes_danio_blast <- merge(as.data.frame(gill_fem_biased),
    blast_output_merged$SN_femG_allTRgenes_danio_blast, by = "trin_geneid",
    all.y = TRUE)

blast_output_merged$SN_femGon_allTRgenes_danio_blast <- merge(as.data.frame(gonad_fem_biased),
    blast_output_merged$SN_femGon_allTRgenes_danio_blast, by = "trin_geneid",
    all.y = TRUE)

blast_output_merged$SN_femL_allTRgenes_danio_blast <- merge(as.data.frame(liver_fem_biased),
    blast_output_merged$SN_femL_allTRgenes_danio_blast, by = "trin_geneid",
    all.y = TRUE)

blast_output_merged$SN_maleG_allTRgenes_danio_blast <- merge(as.data.frame(gill_mal_biased),
    blast_output_merged$SN_maleG_allTRgenes_danio_blast, by = "trin_geneid",
    all.y = TRUE)

blast_output_merged$SN_maleGon_allTRgenes_danio_blast <- merge(as.data.frame(gonad_mal_biased),
    blast_output_merged$SN_maleGon_allTRgenes_danio_blast, by = "trin_geneid",
    all.y = TRUE)

blast_output_merged$SN_maleL_allTRgenes_danio_blast <- merge(as.data.frame(liver_mal_biased),
    blast_output_merged$SN_maleL_allTRgenes_danio_blast, by = "trin_geneid",
    all.y = TRUE)
```

``` r
# Write once then comment out
write.table(c(blast_output_merged$SN_femG_allTRgenes_danio_blast$gene_name,
    blast_output_merged$SN_femGon_allTRgenes_danio_blast$gene_name,
    blast_output_merged$SN_femL_allTRgenes_danio_blast$gene_name,
    blast_output_merged$SN_maleG_allTRgenes_danio_blast$gene_name,
    blast_output_merged$SN_maleGon_allTRgenes_danio_blast$gene_name,
    blast_output_merged$SN_maleL_allTRgenes_danio_blast$gene_name),
    "SN_GOnames.txt", sep = " ", quote = FALSE, row.names = FALSE,
    col.names = FALSE)

# write once then comment out
write.table(rbind(blast_output_merged$SN_femG_allTRgenes_danio_blast[,
    c("gene_name", "log2FoldChange")], blast_output_merged$SN_femGon_allTRgenes_danio_blast[,
    c("gene_name", "log2FoldChange")], blast_output_merged$SN_femL_allTRgenes_danio_blast[,
    c("gene_name", "log2FoldChange")], blast_output_merged$SN_maleG_allTRgenes_danio_blast[,
    c("gene_name", "log2FoldChange")], blast_output_merged$SN_maleGon_allTRgenes_danio_blast[,
    c("gene_name", "log2FoldChange")], blast_output_merged$SN_maleL_allTRgenes_danio_blast[,
    c("gene_name", "log2FoldChange")]), "SN_GOnames_FC.txt",
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
```

``` r
panther <- read.delim("data/pantherGeneList.txt", header = FALSE)
colnames(panther) <- c("GeneID", "MappedID", "GeneName", "pantherFAM",
    "panther_prot_class", "species")

blast_output_merged$SN_femG_allTRgenes_danio_blast <- merge(blast_output_merged$SN_femG_allTRgenes_danio_blast,
    panther[, c(2, 5)], by.x = "gene_name", by.y = "MappedID",
    all.x = TRUE)
blast_output_merged$SN_femGon_allTRgenes_danio_blast <- merge(blast_output_merged$SN_femGon_allTRgenes_danio_blast,
    panther[, c(2, 5)], by.x = "gene_name", by.y = "MappedID",
    all.x = TRUE)
blast_output_merged$SN_femL_allTRgenes_danio_blast <- merge(blast_output_merged$SN_femL_allTRgenes_danio_blast,
    panther[, c(2, 5)], by.x = "gene_name", by.y = "MappedID",
    all.x = TRUE)
blast_output_merged$SN_maleG_allTRgenes_danio_blast <- merge(blast_output_merged$SN_maleG_allTRgenes_danio_blast,
    panther[, c(2, 5)], by.x = "gene_name", by.y = "MappedID",
    all.x = TRUE)
blast_output_merged$SN_maleGon_allTRgenes_danio_blast <- merge(blast_output_merged$SN_maleGon_allTRgenes_danio_blast,
    panther[, c(2, 5)], by.x = "gene_name", by.y = "MappedID",
    all.x = TRUE)
blast_output_merged$SN_maleL_allTRgenes_danio_blast <- merge(blast_output_merged$SN_maleL_allTRgenes_danio_blast,
    panther[, c(2, 5)], by.x = "gene_name", by.y = "MappedID",
    all.x = TRUE)
```

``` r
go_tabs <- lapply(blast_output_merged, function(dat) {

    tab <- data.frame(table(dat$panther_prot_class))
    tab$Var1 <- as.character(tab$Var1)
    tab$Var1[tab$Var1 == ""] <- "Unclassified"
    tab$prop <- tab$Freq/sum(tab$Freq)

    return(tab)
})


all_go_dat <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(all_go_dat) <- c(colnames(go_tabs$SN_femG_allTRgenes_danio_blast),
    "bias_cat", "tissue")
tissues <- c("Gill", "Gonad", "Liver", "Gill", "Gonad", "Liver")

for (i in 1:length(go_tabs)) {
    tmp <- go_tabs[[i]]
    tmp$cat <- names(go_tabs)[i]
    tmp$tissue <- tissues[[i]]
    all_go_dat <- rbind(all_go_dat, tmp)
}

all_go_sums <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(all_go_sums) <- c("Freq", "prot_class", "tissue", "prop")

for (organ in unique(all_go_dat$tissue)) {

    tmp <- all_go_dat[all_go_dat$tissue == organ, ]
    go_sums <- as.data.frame(tapply(tmp$Freq, tmp$Var1, sum))
    go_sums$prot_class <- rownames(go_sums)
    colnames(go_sums)[1] <- "Freq"
    go_sums$tissue <- organ
    rownames(go_sums) <- NULL

    go_sums$prop <- go_sums$Freq/sum(go_sums$Freq)

    all_go_sums <- rbind(all_go_sums, go_sums)

}

all_go_sums$prot_class_final <- ifelse(all_go_sums$prop < 0.02,
    paste0("Other Protien Class"), paste0(all_go_sums$prot_class))

all_go_sums <- all_go_sums[all_go_sums$prot_class != "Unclassified",
    ]


# make a barplot
my_colors <- c(Gonad = "#EEB422", Liver = "#EE8262", Gill = "#20B2AA")
ggplot(all_go_sums, aes(prot_class, prop, fill = tissue)) + geom_bar(position = "dodge",
    stat = "identity") + scale_fill_manual(values = my_colors) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
```

![](DE_analysis_SNI_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

## Looking at GO groups for the male- and female-biased genes

``` r
# Calculate the frequency at which each protein class
# occurs within the different tissue types Create an empty
# dataset to store the values in
all_go_sums_liver <- data.frame(matrix(ncol = 6, nrow = 0))

## Add appropriate column names
colnames(all_go_sums_liver) <- c(colnames(all_go_dat), "prop")

## Loop through each bias type to calculate frequencies
for (bias in c("SN_femL_biased_allTRgenes_danio_blast", "SN_malL_biased_allTRgenes_danio_blast")) {

    tmp <- all_go_dat[all_go_dat$bias_cat == bias, ]

    # Calculate the respective proportion for each protein
    # class
    tmp$prop <- tmp$Freq/sum(tmp$Freq)

    # Export the data to the previously created dataframe
    all_go_sums_liver <- rbind(tmp, all_go_sums_liver)

}


# Change all of the protein classes that have a frequency
# of less than 0.01 to 'Other'
for (class in unique(all_go_sums_liver$Var1)) {

    # Pull out the rows that contain that protein class
    match_rows <- which(all_go_sums_liver$Var1 == class)

    # If any of the proportions associated with those rows
    # is less than 0.01, classify as other
    if (any(all_go_sums_liver$prop[match_rows] > 0.01) == FALSE) {

        all_go_sums_liver$Var1[match_rows] <- " other"

    }
}

# Change format of 'Unclassified so it will come first'
all_go_sums_liver$Var1[all_go_sums_liver$Var1 == "Unclassified"] <- " unclassified"

# Generate the GO plot
bias_col <- c(SN_femL_biased_TRgenes_danio_blast = "#7fc97f",
    SN_malL_biased_TRgenes_danio_blast = "#beaed4")


pdf("figs/FigGO_SBGliver.pdf", height = 14, width = 12)

ggplot(all_go_sums_liver, aes(fct_rev(Var1), prop, fill = "bias_cat")) +
    geom_bar(position = "dodge", stat = "identity") + scale_fill_manual(values = bias_col) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
            size = 20), axis.text.y = element_text(size = 15,
            hjust = 1, vjust = 0.3), text = element_text(size = 20),
        legend.position = "bottom", axis.ticks.length.y = unit(0.25,
            "cm")) + coord_flip() + labs(y = "Proportion of sex-biased genes",
    x = "")

dev.off()
```

    ## png 
    ##   2

``` r
# Calculate the frequency at which each protein class
# occurs within the different tissue types Create an empty
# dataset to store the values in all_go_sums_gonad <-
# data.frame(matrix(ncol = 6, nrow = 0))

## Add appropriate column names colnames(all_go_sums_gonad)
## <- c(colnames(all_go_dat), 'prop')

## Loop through each bias type to calculate frequencies for
## (bias in c('SN_femGon_biased_allTRgenes_danio_blast',
## 'SN_malGon_biased_allTRgenes_danio_blast')) {

# tmp <- all_go_dat[all_go_dat$bias_cat == bias, ]

# Calculate the respective proportion for each protein
# class tmp$prop <- tmp$Freq/sum(tmp$Freq)

# Export the data to the previously created dataframe
# all_go_sums_gonad <- rbind(tmp, all_go_sums_gonad)

# }


# Change all of the protein classes that have a frequency
# of less than 0.01 to 'Other' for(class in
# unique(all_go_sums_gonad$Var1)){

# Pull out the rows that contain that protein class
# match_rows <- which(all_go_sums_gonad$Var1==class)

# If any of the proportions associated with those rows is
# less than 0.01, classify as other
# if(any(all_go_sums_gonad$prop[match_rows]>0.01)==FALSE){

# all_go_sums_gonad$Var1[match_rows] <- ' other'

# } }

# Change format of 'Unclassified so it will come first'
# all_go_sums_gonad$Var1[all_go_sums_gonad$Var1=='Unclassified']
# <- ' unclassified'

# Generate the GO plot bias_col <-
# c('SN_femGon_biased_allTRgenes_danio_blast' = '#7fc97f',
# 'SN_malGon_biased_allTRgenes_danio_blast' = '#beaed4')


# pdf('figs/FigGO_SBGgonad.pdf',height=14, width=12)

# ggplot(all_go_sums_gonad, aes(fct_rev(Var1), prop, fill =
# bias_cat)) + geom_bar(position = 'dodge',
# stat='identity') + scale_fill_manual(values = bias_col) +
# theme(panel.grid.major = element_blank(),
# panel.grid.minor = element_blank(), panel.background =
# element_blank(), axis.line = element_line(colour =
# 'black'), axis.text.x = element_text(angle = 90, vjust =
# 0.5, hjust=1,size=20), axis.text.y =
# element_text(size=15, hjust = 1, vjust = 0.3),
# text=element_text(size=20), legend.position = 'bottom',
# axis.ticks.length.y = unit(.25, 'cm')) + coord_flip() +
# labs(y='Proportion of sex-biased genes', x='')

# dev.off()
```

``` r
# Combine the GO figures
figGOa <- image_ggplot(image_read_pdf("figs/FigGO_SBGliver.pdf"),
    interpolate = TRUE)
```

    ## PDF error: No display font for 'Symbol'

    ## PDF error: No display font for 'ArialUnicode'

``` r
figGOb <- image_ggplot(image_read_pdf("figs/FigGO_SBGgonad.pdf"),
    interpolate = TRUE)

# put the patchworks together
figGO <- wrap_plots(figGOa, figGOb, ncol = 2) + plot_annotation(tag_levels = "A")

ggsave("figs/FigGO_indv.pdf", figGO, height = 5, width = 8)
ggsave("figs/FigGO_indv.png", figGO, height = 5, width = 8)
```

``` r
write.table(blast_output_merged$SN_femGon_allTRgenes_danio_blast$gene_name,
    "SN_GOnames_SBG_ovary.txt", sep = " ", quote = FALSE, row.names = FALSE,
    col.names = FALSE)

write.table(blast_output_merged$SN_maleGon_allTRgenes_danio_blast$gene_name,
    "SN_GOnames_SBG_testes.txt", sep = " ", quote = FALSE, row.names = FALSE,
    col.names = FALSE)

write.table(blast_output_merged$SN_femL_allTRgenes_danio_blast$gene_name,
    "SN_GOnames_SBG_fLiver.txt", sep = " ", quote = FALSE, row.names = FALSE,
    col.names = FALSE)

write.table(blast_output_merged$SN_maleL_allTRgenes_danio_blast$gene_name,
    "SN_GOnames_SBG_mLiver.txt", sep = " ", quote = FALSE, row.names = FALSE,
    col.names = FALSE)
```