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

#' Visualize numerical data over regions or regions summarized over annotations
#'
#' @param tbl A \code{dplyr::tbl} returned from \code{annotate_regions()} or \code{summarize_numerical()}. If the data is not summarized, the data is at the region level. If it is summarized, it represents the average or standard deviation of the regions by the character vector used for \code{by} in \code{summarize_numerical()}.
#' @param x A string indicating the column of the \code{tbl} to use for the x-axis.
#' @param y A string indicating the column of the \code{tbl} to use for the y-axis. Default is \code{NULL}, meaning a histogram over \code{x} will be plotted. If it is not \code{NULL}, a scatterplot is plotted.
#' @param facet A string indicating which categorical variable in the \code{tbl} to make \code{ggplot2} facets over. Default is \code{annot_type}.
#' @param facet_order A character vector which give the order of the facets, and can be used to subset the column in the \code{tbl} used for the \code{facet}. For example, if \code{facet = 'annot_type'}, then the annotations maybe subsetted to just CpG annotations. Default is \code{NULL}, meaning all annotations in their default order are used.
#' @param bin_width An integer indicating the bin width of the histogram used for score. Default 10. Select something appropriate for the data. NOTE: This is only used if \code{y} is \code{NULL}.
#' @param plot_title A string used for the title of the plot. Default \code{NULL}, no title displayed.
#' @param x_label A string used for the x-axis label. Default \code{NULL}, corresponding variable name used.
#' @param y_label A string used for the y-axis label. Default \code{NULL}, corresponding variable name used.
#'
#' @return A \code{ggplot} object which can be viewed by calling it, or saved with \code{ggplot2::ggsave}.
#'
#' @examples
#' # An example with multi-columned data
#' dm = system.file('extdata', 'IDH2mut_v_NBM_multi_data_chr9.txt.gz', package = 'annotatr')
#' annotations = c('hg19_basicgenes','hg19_cpgs')
#'
#' dm_d = read_bed(
#'   file = dm,
#'   col.names=c('chr','start','end','DM_status','pval','strand','diff_meth','mu1','mu0'),
#'   genome = 'hg19',
#'   stranded = FALSE,
#'   use.score = TRUE)
#'
#' # Annotate the regions
#' dm_r = annotate_regions(
#'   regions = dm_d,
#'   annotations = annotations,
#'   ignore.strand = TRUE,
#'   use.score = TRUE)
#'
#' # Summarize a numerical column over the annotations
#' dm_sn = summarize_numerical(
#'   annotated_regions = dm_r,
#'   by = c('annot_type', 'annot_id'),
#'   over = 'diff_meth')
#'
#' # Plot histograms of the mean methylation differences
#' # of regions in annotations and facet on annotations.
#' dm_vs_sumnum = visualize_numerical(
#'   tbl = dm_sn,
#'   x = 'mean',
#'   facet = 'annot_type',
#'   facet_order = c('hg19_cpg_islands','hg19_cpg_shores','hg19_cpg_shelves','hg19_cpg_inter'),
#'   bin_width = 5,
#'   plot_title = 'Mean Meth. Diff. over CpG Annots.',
#'   x_label = 'Methylation Difference')
#'
#' # Can also use the result of annotate_regions() to plot two numerical
#' # data columns against each other for each region, and facet by annotations.
#' dm_vs_regions_annot = visualize_numerical(
#'   tbl = dm_r,
#'   x = 'mu0',
#'   y = 'mu1',
#'   facet = 'annot_type',
#'   facet_order = c('hg19_knownGenes_1to5kb','hg19_knownGenes_promoters','hg19_knownGenes_5UTRs','hg19_knownGenes_3UTRs'),
#'   plot_title = 'Region Methylation: Group 0 vs Group 1',
#'   x_label = 'Group 0',
#'   y_label = 'Group 1')
#'
#' # Another example, but using differential methylation status as the facets.
#' dm_vs_regions_name = visualize_numerical(
#'   tbl = dm_r,
#'   x = 'mu0',
#'   y = 'mu1',
#'   facet = 'DM_status',
#'   facet_order = c('hyper','hypo','none'),
#'   plot_title = 'Region Methylation: Group 0 vs Group 1',
#'   x_label = 'Group 0',
#'   y_label = 'Group 1')
#'
#' @export
visualize_numerical = function(tbl, x, y=NULL, facet = 'annot_type', facet_order = NULL, bin_width=10,
  plot_title=NULL, x_label=NULL, y_label=NULL) {

  ########################################################################
  # Argument parsing and error handling
    if(!(class(tbl)[1] == "tbl_df" || class(tbl)[1] == "grouped_df")) {
      stop('Error: tbl must have class tbl_df or grouped_df. Valid inputs can come from annotate_regions() or summarize_numerical().')
    }

  ########################################################################
  # Order and subset the annotations if facet_order is not NULL
  if(facet == 'annot_type' && is.null(facet_order)) {
    facet_order = unique(tbl[[facet]])
  }
  tbl = order_subset_summary(summary = tbl, col = facet, col_order = facet_order)

  ########################################################################
  # Construct the plot

    if(is.null(y)) {
      # Make the base histogram ggplot
      # NOTE: binwidth may need to be a parameter
      plot =
        ggplot(tbl, aes_string(x=x)) +
        geom_histogram(binwidth=bin_width, aes(y=..density..)) +
        facet_wrap( as.formula(paste("~", facet)) ) +
        theme_bw()
    } else {
      # Make the base scatter ggplot
      plot = ggplot(tbl, aes_string(x=x, y=y)) +
        geom_point(alpha = 1/8, size = 1) +
        facet_wrap( as.formula(paste("~", facet)) ) +
        theme_bw()
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

#' Visualize names over annotations
#'
#' Given a \code{dplyr::grouped_df} of name aggregated annotations (from \code{summarize_categorical()}), visualize the the distribution of \code{annot_type} in \code{data_name}.
#'
#' @param summarized_cats The \code{grouped_df} result of \code{summarize_categorical()}.
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
#' dm = system.file('extdata', 'IDH2mut_v_NBM_multi_data_chr9.txt.gz', package = 'annotatr')
#' annotations = c('hg19_basicgenes','hg19_cpgs')
#'
#' dm_d = read_bed(
#'   file = dm,
#'   col.names=c('chr','start','end','DM_status','pval','strand','diff_meth','mu1','mu0'),
#'   genome = 'hg19',
#'   stranded = FALSE,
#'   use.score = TRUE)
#'
#' dm_r = annotate_regions(
#'   regions = dm_d,
#'   annotations = annotations,
#'   ignore.strand = TRUE,
#'   use.score = TRUE)
#'
#' dm_sc = summarize_categorical(
#'   annotated_regions = dm_r,
#'   by = c('annot_type', 'DM_status'))
#'
#' dm_order = c(
#'   'hyper',
#'   'hypo',
#'   'none')
#' genes_order = c(
#'   'hg19_knownGenes_1to5kb',
#'   'hg19_knownGenes_promoters',
#'   'hg19_knownGenes_5UTRs',
#'   'hg19_knownGenes_exons',
#'   'hg19_knownGenes_introns',
#'   'hg19_knownGenes_3UTRs')
#'
#' dm_vn = visualize_categorical(
#'   summarized_cats = dm_sc,
#'   x = 'DM_status',
#'   fill = 'annot_type',
#'   x_order = dm_order,
#'   fill_order = genes_order,
#'   position = 'fill',
#'   legend_title = 'knownGene Annotations',
#'   x_label = 'DM status',
#'   y_label = 'Proportion')
#'
#' @export
visualize_categorical = function(summarized_cats, x, fill=NULL, x_order=NULL, fill_order=NULL,
  position = 'stack', plot_title=NULL, legend_title=NULL, x_label=NULL, y_label=NULL) {

  ########################################################################
  # Argument parsing and error handling

    # Check correct class of input
    if(class(summarized_cats)[1] != "grouped_df") {
      stop('Error: summarized_cats must have class grouped_df. The best way to ensure this is to pass the result of summarize_categorical() into this function.')
    }

    # Ensure the value of x is a column name in summarized_cats
    if( !(x %in% colnames(summarized_cats)) ) {
      stop('The column name used for x does not exist in summarized_cats. Try using "annot_type" or "data_name".')
    }

    # Ensure the value of fill is a column name in summarized_cats if it isn't NULL
    # Also ensure fill != x
    if( !is.null(fill) ) {
      if( !(fill %in% colnames(summarized_cats)) ) {
        stop('The column name used for fill does not exist in summarized_cats. Try using "annot_type" or "data_name".')
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
  if(x == 'annot_type' && is.null(x_order)) {
    x_order = unique(summarized_cats[[x]])
  }
  summarized_cats = order_subset_summary(summary = summarized_cats, col = x, col_order = x_order)

  ########################################################################
  # Order and subset fill if fill and fill_order are not NULL
  if(!is.null(fill) && fill == 'annot_type' && is.null(fill_order)) {
    fill_order = unique(summarized_cats[[fill]])
  }
  summarized_cats = order_subset_summary(summary = summarized_cats, col = fill, col_order = fill_order)

  ########################################################################
  # Construct the plot

    # Make base ggplot
    plot =
      ggplot(summarized_cats, aes_string(x=x, y='n', fill=fill)) +
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
