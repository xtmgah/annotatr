% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualize.R
\name{visualize_numerical_coannotations}
\alias{visualize_numerical_coannotations}
\title{Visualize numerical data occurring in pairs of annotations}
\usage{
visualize_numerical_coannotations(tbl, x, y = NULL, annot1, annot2,
  bin_width = 10, plot_title = NULL, x_label = NULL, y_label = NULL)
}
\arguments{
\item{tbl}{A \code{dplyr::tbl} returned from \code{annotate_regions()}. If the data is not summarized, the data is at the region level. If it is summarized, it represents the average or standard deviation of the regions by the character vector used for \code{by} in \code{summarize_numerical()}.}

\item{x}{A string indicating the column of the \code{tbl} to use for the x-axis.}

\item{y}{A string indicating the column of the \code{tbl} to use for the y-axis. Default is \code{NULL}, meaning a histogram over \code{x} will be plotted. If it is not \code{NULL}, a scatterplot is plotted.}

\item{annot1}{A string indicating the first annotation type.}

\item{annot2}{A string indicating the second annotation type.}

\item{bin_width}{An integer indicating the bin width of the histogram used for score. Default 10. Select something appropriate for the data. NOTE: This is only used if \code{y} is \code{NULL}.}

\item{plot_title}{A string used for the title of the plot. Default \code{NULL}, no title displayed.}

\item{x_label}{A string used for the x-axis label. Default \code{NULL}, corresponding variable name used.}

\item{y_label}{A string used for the y-axis label. Default \code{NULL}, corresponding variable name used.}
}
\value{
A \code{ggplot} object which can be viewed by calling it, or saved with \code{ggplot2::ggsave}.
}
\description{
Visualize numerical data associated with regions occurring in \code{annot1}, \code{annot2} and in both.
}
\examples{

dm = system.file('extdata', 'IDH2mut_v_NBM_multi_data_chr9.txt.gz', package = 'annotatr')
annotations = c('hg19_basicgenes','hg19_cpgs','hg19_enhancers_fantom')

dm_d = read_bed(
  file = dm,
  col.names=c('chr','start','end','DM_status','pval','strand','diff_meth','mu1','mu0'),
  genome = 'hg19',
  stranded = FALSE,
  use.score = TRUE)

dm_r = annotate_regions(
  regions = dm_d,
  annotations = annotations,
  ignore.strand = TRUE,
  use.score = TRUE)

dm_vs_num_co = visualize_numerical_coannotations(
  tbl = dm_r,
  x = 'mu0',
  annot1 = 'hg19_cpg_islands',
  annot2 = 'hg19_knownGenes_promoters',
  bin_width = 5,
  plot_title = 'Group 0 Perc. Meth. in CpG Islands and Promoters',
  x_label = 'Percent Methylation')

}

