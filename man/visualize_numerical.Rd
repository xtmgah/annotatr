% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualize.R
\name{visualize_numerical}
\alias{visualize_numerical}
\title{Visualize numerical data over regions or regions summarized over annotations}
\usage{
visualize_numerical(tbl, x, y = NULL, facet = "annot_type",
  facet_order = NULL, bin_width = 10, plot_title = NULL, x_label = NULL,
  y_label = NULL)
}
\arguments{
\item{tbl}{A \code{dplyr::tbl} returned from \code{annotate_regions()} or \code{summarize_numerical()}. If the data is not summarized, the data is at the region level. If it is summarized, it represents the average or standard deviation of the regions by the character vector used for \code{by} in \code{summarize_numerical()}.}

\item{x}{A string indicating the column of the \code{tbl} to use for the x-axis.}

\item{y}{A string indicating the column of the \code{tbl} to use for the y-axis. Default is \code{NULL}, meaning a histogram over \code{x} will be plotted. If it is not \code{NULL}, a scatterplot is plotted.}

\item{facet}{A string indicating which categorical variable in the \code{tbl} to make \code{ggplot2} facets over. Default is \code{annot_type}.}

\item{facet_order}{A character vector which give the order of the facets, and can be used to subset the column in the \code{tbl} used for the \code{facet}. For example, if \code{facet = 'annot_type'}, then the annotations maybe subsetted to just CpG annotations. Default is \code{NULL}, meaning all annotations in their default order are used.}

\item{bin_width}{An integer indicating the bin width of the histogram used for score. Default 10. Select something appropriate for the data. NOTE: This is only used if \code{y} is \code{NULL}.}

\item{plot_title}{A string used for the title of the plot. Default \code{NULL}, no title displayed.}

\item{x_label}{A string used for the x-axis label. Default \code{NULL}, corresponding variable name used.}

\item{y_label}{A string used for the y-axis label. Default \code{NULL}, corresponding variable name used.}
}
\value{
A \code{ggplot} object which can be viewed by calling it, or saved with \code{ggplot2::ggsave}.
}
\description{
The \code{facet_order} vector also determines the categories used to create the background distribution that appears overlayed in red in each facet.
}
\examples{
# An example with multi-columned data
dm = system.file('extdata', 'IDH2mut_v_NBM_multi_data_chr9.txt.gz', package = 'annotatr')
annotations = c('hg19_basicgenes','hg19_cpgs')

dm_d = read_bed(
  file = dm,
  col.names=c('chr','start','end','DM_status','pval','strand','diff_meth','mu1','mu0'),
  genome = 'hg19',
  stranded = FALSE,
  use.score = TRUE)

# Annotate the regions
dm_r = annotate_regions(
  regions = dm_d,
  annotations = annotations,
  ignore.strand = TRUE,
  use.score = TRUE)

# Summarize a numerical column over the annotations
dm_sn = summarize_numerical(
  annotated_regions = dm_r,
  by = c('annot_type', 'annot_id'),
  over = 'diff_meth')

# Plot histograms of the mean methylation differences
# of regions in annotations and facet on annotations.
# NOTE: Background distribution (everything in \\code{facet_order})
# is plotted in each facet for comparison.
dm_vs_sumnum = visualize_numerical(
  tbl = dm_sn,
  x = 'mean',
  facet = 'annot_type',
  facet_order = c('hg19_cpg_islands','hg19_cpg_shores','hg19_cpg_shelves','hg19_cpg_inter'),
  bin_width = 5,
  plot_title = 'Mean Meth. Diff. over CpG Annots.',
  x_label = 'Methylation Difference')

# Plot histogram of group 1 methylation rates across the CpG annotations.
# NOTE: Background distribution (everything in \\code{facet_order})
# is plotted in each facet for comparison.
dm_vs_regions_mu1 = visualize_numerical(
 tbl = dm_r,
 x = 'mu1',
 facet = 'annot_type',
 facet_order = c('hg19_cpg_islands','hg19_cpg_shores','hg19_cpg_shelves','hg19_cpg_inter'),
 bin_width = 5,
 plot_title = 'Group 1 Methylation over CpG Annotations',
 x_label = 'Group 1 Methylation')

# Can also use the result of annotate_regions() to plot two numerical
# data columns against each other for each region, and facet by annotations.
dm_vs_regions_annot = visualize_numerical(
  tbl = dm_r,
  x = 'mu0',
  y = 'mu1',
  facet = 'annot_type',
  facet_order = c('hg19_knownGenes_1to5kb','hg19_knownGenes_promoters','hg19_knownGenes_5UTRs','hg19_knownGenes_3UTRs'),
  plot_title = 'Region Methylation: Group 0 vs Group 1',
  x_label = 'Group 0',
  y_label = 'Group 1')

# Another example, but using differential methylation status as the facets.
dm_vs_regions_name = visualize_numerical(
  tbl = dm_r,
  x = 'mu0',
  y = 'mu1',
  facet = 'DM_status',
  facet_order = c('hyper','hypo','none'),
  plot_title = 'Region Methylation: Group 0 vs Group 1',
  x_label = 'Group 0',
  y_label = 'Group 1')

}

