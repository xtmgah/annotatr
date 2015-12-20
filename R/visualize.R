#' Visualize the number of regions per annotation
#'
#' Given a \code{dplyr::tbl_df} of counts of regions per annotation (from \code{summarize_annotations()}), visualize the counts per annotation as a bar graph.
#'
#' @param summarized_annotations The \code{tbl_df} result of \code{summarize_annotations()}.
#' @param annotation_order A character vector which doubles as the subset of annotations desired for visualization as well as the ordering. If \code{NULL}, all annotations are displayed.
#' @param plot_title A string used for the title of the plot. Default \code{NULL}, no title is displayed.
#' @param x_label A string used for the x-axis label. Default \code{NULL}, corresponding variable name used.
#' @param y_label A string used for the y-axis label. Default \code{NULL}, corresponding variable name used.
#'
#' @return A \code{ggplot} object which can be viewed by calling it, or saved with \code{ggplot2::ggsave}.
#'
#' @examples
#' ########################################################################
#' # An example of ChIP-seq peaks with signalValue used for score
#' bed = system.file('extdata', 'Gm12878_Ezh2_sorted_scores.narrowPeak.gz', package = 'annotatr')
#' annotations = c('hg19_basicgenes','hg19_cpgs')
#'
#' d = read_bed(file = bed, genome = 'hg19', stranded = FALSE, use.score = TRUE)
#'
#' i = annotate_regions(
#'   regions = d,
#'   annotations = annotations,
#'   ignore.strand = TRUE,
#'   use.score = TRUE)
#'
#' s = summarize_annotations(i)
#'
#' annots_order = c(
#'   'hg19_cpg_islands',
#'   'hg19_cpg_shores',
#'   'hg19_cpg_shelves',
#'   'hg19_cpg_inter',
#'   'hg19_knownGenes_1to5kb',
#'   'hg19_knownGenes_promoters',
#'   'hg19_knownGenes_5UTRs',
#'   'hg19_knownGenes_exons',
#'   'hg19_knownGenes_introns',
#'   'hg19_knownGenes_3UTRs')
#' v_annots = visualize_annotation(s, annots_order)
#'
#' @export
visualize_annotation = function(summarized_annotations, annotation_order=NULL,
  plot_title=NULL, x_label=NULL, y_label=NULL) {

  ########################################################################
  # Argument parsing and error handling

    if(class(summarized_annotations)[1] != "tbl_df") {
      stop('Error: summarized_annotations must have class tbl_df. The best way to ensure this is to pass the result of summarize_annotations() into this function.')
    }

  ########################################################################
  # Order and subset the annotations if annotation_order is not NULL
  summarized_annotations = order_subset_summary(summary = summarized_annotations, col='annot_type', col_order=annotation_order)

  ########################################################################
  # Construct the plot

    # Make the base ggplot
    # NOTE: binwidth may need to be a parameter
    plot =
      ggplot(summarized_annotations, aes(x=annot_type, y=n)) +
      geom_bar(stat='identity') +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))

    # Add any user defined labels to the plot if their values are not NULL
    # if they are NULL, ggplot() will use defaults
    if(!is.null(plot_title)) {
      plot = plot + ggtitle(plot_title)
    }
    if(!is.null(x_label)) {
      plot = plot + xlab(x_label)
    }
    if(!is.null(y_label)) {
      plot = plot + ylab(y_label)
    }

  return(plot)
}

#' Visualize score distributions over annotations
#'
#' Given a \code{dplyr::grouped_df} of score aggregated annotations (from \code{summarize_numerical()}), visualize the distribution of scores over the annotation types.
#'
#' @param summarized_scores The \code{grouped_df} result of \code{summarize_numerical()}.
#' @param annotation_order A character vector which doubles as the subset of annotations desired for visualization as well as the ordering.
#' @param bin_width An integer indicating the bin width of the histogram used for score. Default 10. Select something appropriate for the data.
#' @param plot_title A string used for the title of the plot. Default \code{NULL}, no title displayed.
#' @param x_label A string used for the x-axis label. Default \code{NULL}, corresponding variable name used.
#' @param y_label A string used for the y-axis label. Default \code{NULL}, corresponding variable name used.
#'
#' @return A \code{ggplot} object which can be viewed by calling it, or saved with \code{ggplot2::ggsave}.
#'
#' @examples
#' ########################################################################
#' # An example of ChIP-seq peaks with signalValue used for score
#' bed = system.file('extdata', 'Gm12878_Ezh2_sorted_scores.narrowPeak.gz', package = 'annotatr')
#' annotations = c('hg19_basicgenes','hg19_cpgs')
#'
#' d = read_bed(file = bed, genome = 'hg19', stranded = FALSE, use.score = TRUE)
#'
#' i = annotate_regions(
#'   regions = d,
#'   annotations = annotations,
#'   ignore.strand = TRUE,
#'   use.score = TRUE)
#'
#' s = summarize_numerical(i)
#'
#' hg19_cpgs_order = c(
#'   'hg19_cpg_islands',
#'   'hg19_cpg_shores',
#'   'hg19_cpg_shelves',
#'   'hg19_cpg_inter')
#' v_hg19_cpgs = visualize_numerical(s, hg19_cpgs_order)
#'
#' genes_order = c(
#'   'hg19_knownGenes_1to5kb',
#'   'hg19_knownGenes_promoters',
#'   'hg19_knownGenes_5UTRs',
#'   'hg19_knownGenes_exons',
#'   'hg19_knownGenes_introns',
#'   'hg19_knownGenes_3UTRs')
#' v_genes = visualize_numerical(s, genes_order, bin_width=30)
#'
#' ########################################################################
#' # An example of percent methylation at CpG sites
#' bed = system.file('extdata', '2607_mc_hmc_perc_meth_chr21.txt.gz', package = 'annotatr')
#' annotations = c('hg19_basicgenes','hg19_cpgs')
#'
#' d = read_bed(file = bed, genome = 'hg19', stranded = FALSE, use.score = TRUE)
#'
#' i = annotate_regions(
#'   regions = d,
#'   annotations = annotations,
#'   ignore.strand = TRUE,
#'   use.score = TRUE)
#'
#' s = summarize_numerical(i)
#'
#' hg19_cpgs_order = c(
#'   'hg19_cpg_islands',
#'   'hg19_cpg_shores',
#'   'hg19_cpg_shelves',
#'   'hg19_cpg_inter')
#' v_hg19_cpgs = visualize_numerical(s, hg19_cpgs_order)
#'
#' genes_order = c(
#'   'hg19_knownGenes_1to5kb',
#'   'hg19_knownGenes_promoters',
#'   'hg19_knownGenes_5UTRs',
#'   'hg19_knownGenes_exons',
#'   'hg19_knownGenes_introns',
#'   'hg19_knownGenes_3UTRs')
#' v_genes = visualize_numerical(s, genes_order, bin_width=5)
#'
#' @export
visualize_numerical = function(summarized_scores, annotation_order=NULL, bin_width=10,
  plot_title=NULL, x_label=NULL, y_label=NULL) {

  ########################################################################
  # Argument parsing and error handling
    if(class(summarized_scores)[1] != "grouped_df") {
      stop('Error: summarized_scores must have class grouped_df. The best way to ensure this is to pass the result of summarize_numerical() into this function.')
    }

  ########################################################################
  # Order and subset the annotations if annotation_order is not NULL
  summarized_scores = order_subset_summary(summary = summarized_scores, col = 'annot_type', col_order = annotation_order)

  ########################################################################
  # Construct the plot

    # Make the base ggplot
    # NOTE: binwidth may need to be a parameter
    plot =
      ggplot(summarized_scores, aes(mean)) +
      geom_histogram(binwidth=bin_width, aes(y=..density..)) +
      facet_wrap(~annot_type) +
      theme_bw()

    # Add any user defined labels to the plot if their values are not NULL
    # if they are NULL, ggplot() will use defaults
    if(!is.null(plot_title)) {
      plot = plot + ggtitle(plot_title)
    }
    if(!is.null(x_label)) {
      plot = plot + xlab(x_label)
    }
    if(!is.null(y_label)) {
      plot = plot + ylab(y_label)
    }

  return(plot)
}

#' Visualize names over annotations
#'
#' Given a \code{dplyr::grouped_df} of name aggregated annotations (from \code{summarize_categorical()}), visualize the the distribution of \code{annot_type} in \code{data_name}.
#'
#' @param summarized_names The \code{grouped_df} result of \code{summarize_categorical()}.
#' @param x One of 'annot_type' or 'data_name', indicating whether annotation classes or data classes will appear on the x-axis.
#' @param fill One of 'annot_type', 'data_name', or \code{NULL}, indicating whether annotation classes or data classes will fill the bars. If \code{NULL} then the bars will be the total counts of the x classes.
#' @param x_order A character vector that subsets and orders the x classes. Default \code{NULL}, uses existing values.
#' @param fill_order A character vector that subsets and orders the fill classes. Default \code{NULL}, uses existing values.
#' @param position A string which has the same possible values as in \code{ggplot2::geom_bar(..., position)}, i.e., 'stack', 'fill', 'dodge', etc.
#' @param plot_title A string used for the title of the plot. Default \code{NULL}, no title displayed.
#' @param legend_title A string used for the legend title to describe fills (if fill is not \code{NULL}). Default \code{NULL}, displays corresponding variable name.
#' @param x_label A string used for the x-axis label. Default \code{NULL}, corresponding variable name used.
#' @param y_label A string used for the y-axis label. Default \code{NULL}, corresponding variable name used.
#'
#' @return A \code{ggplot} object which can be viewed by calling it, or saved with \code{ggplot2::ggsave}.
#'
#' @examples
#' # An example of differentially methylated regions classified as DM up, DM down, or no DM
#' bed = system.file('extdata', 'IDH2mut_v_NBM_names_scores_chr9.txt.gz', package = 'annotatr')
#' annotations = c('hg19_basicgenes','hg19_cpgs')
#'
#' d = read_bed(file = bed, genome = 'hg19', stranded = FALSE, use.score = TRUE)
#'
#' i = annotate_regions(
#'   regions = d,
#'   annotations = annotations,
#'   ignore.strand = TRUE,
#'   use.score = TRUE)
#'
#' s = summarize_categorical(i)
#'
#' fill_order = c(
#'   'hg19_cpg_islands',
#'   'hg19_cpg_shores',
#'   'hg19_cpg_shelves',
#'   'hg19_cpg_inter')
#' x_order = c(
#'   'hyper',
#'   'hypo')
#' v_hg19_cpgs_counts_data_annot = visualize_categorical(summarized_names=s, x='name', fill='annot_type',
#'   x_order = x_order, fill_order = fill_order, position='stack')
#' v_hg19_cpgs_proportions_data_annot = visualize_categorical(summarized_names=s, x='name', fill='annot_type',
#'   x_order = x_order, fill_order = fill_order, position='fill')
#' v_hg19_cpgs_nofill_data = visualize_categorical(summarized_names=s, x='name', fill=NULL,
#'   x_order = x_order, fill_order = fill_order, position='stack')
#'
#' x_order = c(
#'   'hg19_cpg_islands',
#'   'hg19_cpg_shores',
#'   'hg19_cpg_shelves',
#'   'hg19_cpg_inter')
#' fill_order = c(
#'   'hyper',
#'   'hypo',
#'   'none')
#' v_hg19_cpgs_counts_annot_data = visualize_categorical(summarized_names=s, x='annot_type', fill='name',
#'   x_order = x_order, fill_order = fill_order, position='stack')
#' v_hg19_cpgs_proportions_annot_data = visualize_categorical(summarized_names=s, x='annot_type', fill='name',
#'   x_order = x_order, fill_order = fill_order, position='fill')
#' v_hg19_cpgs_nofill_annot = visualize_categorical(summarized_names=s, x='annot_type', fill=NULL,
#'   x_order = x_order, fill_order = fill_order, position='stack')
#'
#' @export
visualize_categorical = function(summarized_names, x, fill=NULL, x_order=NULL, fill_order=NULL,
  position = 'stack', plot_title=NULL, legend_title=NULL, x_label=NULL, y_label=NULL) {

  ########################################################################
  # Argument parsing and error handling

    # Check correct class of input
    if(class(summarized_names)[1] != "grouped_df") {
      stop('Error: summarized_names must have class grouped_df. The best way to ensure this is to pass the result of summarize_categorical() into this function.')
    }

    # Ensure the value of x is a column name in summarized_names
    if( !(x %in% colnames(summarized_names)) ) {
      stop('The column name used for x does not exist in summarized_names. Try using "annot_type" or "data_name".')
    }

    # Ensure the value of fill is a column name in summarized_names if it isn't NULL
    # Also ensure fill != x
    if( !is.null(fill) ) {
      if( !(fill %in% colnames(summarized_names)) ) {
        stop('The column name used for fill does not exist in summarized_names. Try using "annot_type" or "data_name".')
      }
      if( x == fill ) {
        stop('Error: x cannot equal fill')
      }
    }

    # Check valid position argument
    if(position != 'stack' && position != 'fill' && position != 'dodge') {
      stop('Error: position must be one of "stack", "fill", or "dodge"')
    }

  ########################################################################
  # Order and subset x if x_order isn't NULL
  summarized_names = order_subset_summary(summary = summarized_names, col = x, col_order = x_order)

  ########################################################################
  # Order and subset fill if fill and fill_order are not NULL
  summarized_names = order_subset_summary(summary = summarized_names, col = fill, col_order = fill_order)

  ########################################################################
  # Construct the plot

    # Make base ggplot
    plot =
      ggplot(summarized_names, aes_string(x=x, y='n', fill=fill)) +
      geom_bar(stat='identity', position=position, aes_string(order = fill), width=0.5) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))

    # Change the fill scale and name if legend_title isn't null
    if(!is.null(legend_title)) {
      plot = plot + scale_fill_brewer(name=legend_title)
    } else {
      plot = plot + scale_fill_brewer()
    }

    # Add any user defined labels to the plot if their values are not NULL
    # if they are NULL, ggplot() will use defaults
    if(!is.null(plot_title)) {
      plot = plot + ggtitle(plot_title)
    }
    if(!is.null(x_label)) {
      plot = plot + xlab(x_label)
    }
    if(!is.null(y_label)) {
      plot = plot + ylab(y_label)
    }

  return(plot)
}
