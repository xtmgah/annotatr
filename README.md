[![Travis-CI Build Status](https://travis-ci.org/rcavalcante/annotatr.svg?branch=master)](https://travis-ci.org/rcavalcante/annotatr) [![Coverage Status](https://coveralls.io/repos/rcavalcante/annotatr/badge.svg?branch=master&service=github)](https://coveralls.io/github/rcavalcante/annotatr?branch=master)

# Contents
* [Introduction](https://github.com/rcavalcante/annotatr#introduction)
* [Installation](https://github.com/rcavalcante/annotatr#installation)
* [Annotations](https://github.com/rcavalcante/annotatr#annotations)
    * [CpG Annotations](https://github.com/rcavalcante/annotatr#cpg-annotations)
    * [UCSC knownGenes](https://github.com/rcavalcante/annotatr#ucsc-knowngenes)
* [Usage](https://github.com/rcavalcante/annotatr#usage)
    * [Reading A File](https://github.com/rcavalcante/annotatr#reading-a-file)
    * [Annotating Regions](https://github.com/rcavalcante/annotatr#annotating-regions)
    * [Summarizing Over Annotations](https://github.com/rcavalcante/annotatr#summarizing-over-annotations)
    * [Visualizing](https://github.com/rcavalcante/annotatr#visualizing)
        * [Visualizing Over Regions](https://github.com/rcavalcante/annotatr#visualizing-over-regions)
        * [Visualizing Over Summaries](https://github.com/rcavalcante/annotatr#visualizing-over-summaries)

# Introduction

Genomic regions resulting from next-generation sequencing experiments and bioinformatics pipelines are made more valuable when annotated to genomic features. A SNP occurring in an exon, or an enhancer, is likely of greater interest than one occurring in an inter-genic region. It may be of interest to find that a particular transcription factor overwhelmingly binds in promoters, while another binds mostly in 3’UTRs. Hyper-methylation at promoters containing a CpG island may indicate different regulatory regimes in one condition compared to another.

`annotatr` provides *pre-computed* genomic annotations and a set of functions to read, intersect, summarize, and visualize genomic regions in the context of genomic annotations.

# Installation

`annotatr` package source is available at [http://www.github.com/rcavalcante/annotatr](http://www.github.com/rcavalcante/annotatr). The package can be installed directly from GitHub with the [`devtools`](https://cran.r-project.org/web/packages/devtools/index.html) package.

```{r, eval=FALSE}
devtools::install_github('rcavalcante/annotatr')
```

# Annotations

We downloaded CpG island and UCSC knownGene tracks (for hg19, hg38, mm9, and mm10) from the UCSC Genome Browser. The package source stores these files in `data-raw/`. The package source also contains scripts (in `data-raw/`) used to transform the raw data into `GenomicRanges` objects. The final annotations are bundled with the package and are located in `data/`. Supported annotations are listed with `supported_annotations()` after loading the `annotatr` package.

## CpG Annotations

The base CpG island (CGI) track serves as our CpG island annotations. CpG shores are defined as 2Kb upstream/downstream from the ends of the CpG islands, less the CpG islands. CpG shelves are defined as another 2Kb upstream/downstream of the farthest upstream/downstream limits of the CpG shores, less the CpG islands and CpG shores. The remaining genomic regions make up the inter-CGI annotation.

## UCSC knownGenes

The UCSC knownGenes annotations include 1-5Kb upstream of the TSS, the promoter (<1Kb upstream of the TSS), 5'UTR, exons, introns, CDS, 3'UTR, and 5'UTR exons, 5'UTR introns, 3'UTR exons, and 3'UTR introns. The schematic below gives an idea of how the location coordinates in the knownGene files can be used to determine the annotations.

![Schematic of knownGene annotations.](annotatr_knownGenes.jpeg)

# Usage

## Reading A File

`annotatr` reads in genomic regions of interest encoded as [BED6](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) using the `read_bed()` function. The `name` column may be a non-unique character vector indicating some classification for each region (e.g. hyper-methylated and hypo-methylated). The `score` column may take continuous or discrete numerical values that represent data associated with each region (e.g. the percent methylation).

`read_bed()` requires the `file` name to be read, a `genome` (one of hg19, hg38, mm9, or mm10), a logical value indicating whether the data is `stranded`, and a logical value indicating whether the `score` column should be used.

```{r, echo=FALSE}
suppressWarnings(devtools::load_all())
```

```{r}
# This file in inst/extdata represents regions tested for differential
# methylation between two conditions. Additionally, there are columns
# reporting the p-value on the test for differential meth., the
# meth. difference between the two groups, and the group meth. rates.
dm_file = system.file(
  'extdata',
  'IDH2mut_v_NBM_multi_data_chr9.txt.gz',
  package='annotatr')
dm_regions = read_bed(
  file = dm_file,
  col.names= c(
    'chr','start','end','DM_status','pval',
    'strand','diff_meth','mu1','mu0'),
  genome = 'hg19',
  stranded = FALSE,
  use.score = TRUE)
print(dm_regions)
```

## Annotating Regions

Users may select annotations a la carte via the accessors listed with `supported_annotations()`, or via shortcuts. The `hg19_cpgs` shortcut annotates regions to CpG islands, CpG shores, CpG shelves, and inter-CGI. The `hg19_basicgenes` shortcut annotates regions to 1-5Kb, promoters, 5'UTRs, exons, introns, and 3'UTRs. The `hg19_detailedgenes` shortcut annotates regions to 1-5Kb, promoters, 5'UTR exons, 5'UTR introns, CDS exons, CDS introns, 3'UTR exons, and 3'UTR introns. The `hg19_basicgenes` and `hg19_detailedgenes` shortcuts may not be used at the same time. Shortcuts for other `supported_genomes()` are accessed in the same way, replacing `hg19` with one of `hg38`, `mm9`, or `mm10`.

`annotate_regions()` requires a `GRanges` object (either the result of `read_bed()` or an object the user already has), a character vector indicating the `annotations` to annotate the regions with, the `genome` as above, a logical value indicating whether to `ignore.strand` when calling `GenomicRanges::findOverlaps()`, and a logical value indicating whether to `use.score`.

```{r}
# Select annotations for intersection with regions
annots = c('hg19_cpgs','hg19_basicgenes')
# Intersect the regions we read in with the annotations
dm_annotated = annotate_regions(
  regions = dm_regions,
  annotations = annots,
  ignore.strand = TRUE,
  use.score = TRUE)
# A dplyr::tbl_df object is returned
# NOTE: tbl_df has a special print method that displays
# a tidy version of the tbl_df, variables not shown are listed below.
print(dm_annotated)
```

## Summarizing Over Annotations

The dataset used in the running example has a class (differential methylation status: hyper, hypo, none) and a score (methylation difference between two groups) associated with each region. We will employ all three summarization functions in preparation for visualization.

The three summarization functions -- `summarize_annotations()`, `summarize_numerical()`, and `summarize_categorical()` -- all take the `tbl_df` objects output from `annotate_regions()` as their input, and output either a `dplyr::tbl_df` (for `summarize_annotations()`) or `dplyr::grouped_df` (for the others).

```{r, eval=FALSE}
# Usage of summarize functions with defaults
summarize_annotations(annotated_regions)

summarize_numerical(
  annotated_regions,
  by = c("annot_type", "annot_id"),
  over = "score")

summarize_categorical(
  annotated_regions,
  by = c("annot_type", "annot_id"))
```

When there is no class (`name` column) or score (`score` column) information associated with the regions, `summarize_annotations()` is the only possible function to use. It gives the counts of regions in each annotation type (see example below). If there is class and/or score information, then `summarize_numerical()` and/or `summarize_categorical()` may be used.

```{r}
# Find the number of regions per annotation type
dm_annsum = summarize_annotations(annotated_regions = dm_annotated)
print(dm_annsum)

# Take the mean of the score column across all regions
# occurring in an annotation.
dm_numsum = summarize_numerical(
  annotated_regions = dm_annotated,
  by = c('annot_type', 'annot_id'),
  over = c('diff_meth'))
print(dm_numsum)

# Count the occurrences of classifications in the name
# column across the annotation types.
dm_catsum = summarize_categorical(
  annotated_regions = dm_annotated,
  by = c('annot_type', 'name'))
print(dm_catsum)
```

## Visualizing

`annotatr` has three visualization functions (`visualize_annotation()`, `visualize_numerical()` and `visualize_categorical()`) that are matched with the summarization functions. The visualization functions return an object of type `ggplot` that can be viewed (`print`) and saved (`ggsave`).

```{r, eval=FALSE}
# Usage of visualization functions with defaults

visualize_annotation(summarized_annotations, annotation_order = NULL,
  plot_title = NULL, x_label = NULL, y_label = NULL)

visualize_numerical(tbl, x, y = NULL, facet = "annot_type",
  facet_order = NULL, bin_width = 10, plot_title = NULL, x_label = NULL,
  y_label = NULL)

visualize_categorical(summarized_cats, x, fill = NULL, x_order = NULL,
  fill_order = NULL, position = "stack", plot_title = NULL,
  legend_title = NULL, x_label = NULL, y_label = NULL)
```

All visualizations use an internal function, `tidy_annotations()`, on `annotation_order`, `x_order`, and `fill_order` parameters, which cleans annotation accessors for visualization. Each visualization function is designed to work with a small number of parameters. Order parameters are not required; when `NULL`, all values are displayed. Parameters for plot titles, legend titles, and axis labels are likewise optional, and will default to object name attributes if `NULL`.

### Visualizing Over Regions

With numerical data, the `visualize_numerical()` function allows for the plotting of a single variable (histogram) or of two variables (scatterplot) at the region level, faceting over the annotations, rather than summarizing over them.

```{r, fig.align='center', fig.cap='Methylation Rates in Regions Over Gene Features in Control Group.', fig.height=6, fig.width=6, fig.show='hold'}
dm_vs_regions_annot = visualize_numerical(
  tbl = dm_annotated,
  x = 'mu0',
  facet = 'annot_type',
  facet_order = c('hg19_knownGenes_1to5kb','hg19_knownGenes_promoters','hg19_knownGenes_5UTRs','hg19_knownGenes_3UTRs'),
  plot_title = 'Group 0 Region Methylation In Genes',
  x_label = 'Group 0')
print(dm_vs_regions_annot)
```

```{r, fig.align='center', fig.cap='Methylation Rates in Regions Over DM Status in Group 0 vs Group 1.', fig.height=5, fig.width=7, fig.show='hold'}
dm_vs_regions_name = visualize_numerical(
  tbl = dm_annotated,
  x = 'mu0',
  y = 'mu1',
  facet = 'name',
  facet_order = c('hyper','hypo','none'),
  plot_title = 'Region Methylation: Group 0 vs Group 1',
  x_label = 'Group 0',
  y_label = 'Group 1')
print(dm_vs_regions_name)
```

### Visualizing Over Summaries

In some instances, it is more helpful to visualized the summarized data.

```{r, fig.align='center', fig.cap='Number of DM regions per annotation.', fig.height=5, fig.width=5, fig.show = 'hold'}
# View the number of regions per annotation. This function
# is useful when there is no classification or data
# associated with the regions.
annots_order = c(
  'hg19_knownGenes_1to5kb',
  'hg19_knownGenes_promoters',
  'hg19_knownGenes_5UTRs',
  'hg19_knownGenes_exons',
  'hg19_knownGenes_introns',
  'hg19_knownGenes_3UTRs')
dm_vs_kg_annotations = visualize_annotation(
  summarized_annotations = dm_annsum,
  annotation_order = annots_order,
  plot_title = 'Number of Sites Tested for DM on chr9',
  x_label = 'knownGene Annotations',
  y_label = 'Count')
print(dm_vs_kg_annotations)
```

```{r, fig.align='center', fig.cap='Methylation difference distributions across the CpG annotations.', fig.height=6, fig.width=6, fig.show = 'hold'}
# View the score summarized over the CpG annotations

# A subset of cpgs_order could be chosen to display
# the data distribution only in those annotations.
cpgs_order = c(
  'hg19_cpg_islands',
  'hg19_cpg_shores',
  'hg19_cpg_shelves',
  'hg19_cpg_inter')
dm_vs_cpg_num = visualize_numerical(
  tbl = dm_numsum,
  x = 'mean',
  facet = 'annot_type',
  facet_order = cpgs_order,
  bin_width = 5,
  plot_title = 'Group Meth. Diffs. Across CpG Annotations',
  x_label = 'Mean Diff. in Meth. Over Annotations')
print(dm_vs_cpg_num)
```

```{r, fig.align='center', fig.cap='Methylation difference distributions across the basic genes annotation.', fig.height=6.5, fig.width=6.5, fig.show='hold'}
# View the score summarized over the knownGene annotations.

# A subset of genes_order could be chosen to display
# the data distribution only in those annotations.
genes_order = c(
  'hg19_knownGenes_1to5kb',
  'hg19_knownGenes_promoters',
  'hg19_knownGenes_5UTRs',
  'hg19_knownGenes_exons',
  'hg19_knownGenes_introns',
  'hg19_knownGenes_3UTRs')
dm_vs_kg_num = visualize_numerical(
  tbl = dm_numsum,
  x = 'mean',
  facet = 'annot_type',
  facet_order = genes_order,
  bin_width = 5,
  plot_title = 'Group Meth. Diffs. Across knownGene Annotations',
  x_label = 'Mean Diff. in Meth. Over Annotations')
print(dm_vs_kg_num)
```

```{r, fig.align='center', fig.cap='Differential methylation classification with counts of CpG annotations.', fig.height=6, fig.width=6, fig.show='hold'}
# View the counts of CpG annotations in data classes

# The orders for the x-axis labels. This is also a subset
# of the labels (hyper, hypo, none).
x_order = c(
  'hyper',
  'hypo')
# The orders for the fill labels. Can also use this
# parameter to subset annotation types to fill.
fill_order = c(
  'hg19_cpg_islands',
  'hg19_cpg_shores',
  'hg19_cpg_shelves',
  'hg19_cpg_inter')
# Make a barplot of the data class where each bar
# is composed of the counts of CpG annotations.
dm_vs_cpg_cat1 = visualize_categorical(
  summarized_cats=dm_catsum, x='name', fill='annot_type',
  x_order = x_order, fill_order = fill_order, position='stack',
  plot_title = 'DM Status by CpG Annotation Counts',
  legend_title = 'Annotations',
  x_label = 'DM status',
  y_label = 'Count')
print(dm_vs_cpg_cat1)
```

```{r, fig.align='center', fig.cap='Differential methylation classification with proportion of CpG annotations.', fig.height=6, fig.width=6, fig.show='hold'}
# Use the same order vectors as the previous code block,
# but use proportional fill instead of counts.

# Make a barplot of the data class where each bar
# is composed of the *proportion* of CpG annotations.
dm_vs_cpg_cat2 = visualize_categorical(
  summarized_cats=dm_catsum, x='name', fill='annot_type',
  x_order = x_order, fill_order = fill_order, position='fill',
  plot_title = 'DM Status by CpG Annotation Proportions',
  legend_title = 'Annotations',
  x_label = 'DM status',
  y_label = 'Proportion')
print(dm_vs_cpg_cat2)
```

```{r, fig.align='center', fig.cap='Basic gene annotations with proportions of DM classification.', fig.height=6, fig.width=6, fig.show='hold'}
# View the proportions of data classes in knownGene annotations

# The orders for the x-axis labels.
x_order = c(
  'hg19_knownGenes_1to5kb',
  'hg19_knownGenes_promoters',
  'hg19_knownGenes_5UTRs',
  'hg19_knownGenes_exons',
  'hg19_knownGenes_introns',
  'hg19_knownGenes_3UTRs')
# The orders for the fill labels.
fill_order = c(
  'hyper',
  'hypo',
  'none')
dm_vs_kg_cat = visualize_categorical(
  summarized_cats=dm_catsum, x='annot_type', fill='name',
  x_order = x_order, fill_order = fill_order, position='fill',
  legend_title = 'DM Status',
  x_label = 'knownGene Annotations',
  y_label = 'Proportion')
print(dm_vs_kg_cat)
```
