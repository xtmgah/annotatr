#' Visualize the number of regions per annotation
#'
#' Given a \code{dplyr::tbl_df} of counts of regions per annotation (from \code{summarize_annotations()}), visualize the counts per annotation as a bar graph.
#'
#' @param annotated_regions The \code{tbl_df} result of \code{annotate_regions()}.
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
#' v_annots = visualize_annotation(i, annots_order)
#'
#' @export
visualize_annotation = function(annotated_regions, annotation_order=NULL,
  plot_title=NULL, x_label=NULL, y_label=NULL) {

  ########################################################################
  # Argument parsing and error handling

    if(class(annotated_regions)[1] != "tbl_df") {
      stop('Error: annotated_regions must have class tbl_df. The best way to ensure this is to pass the result of annotate_regions() into this function.')
    }

  ########################################################################
  # Order and subset the annotations
  annotated_regions = subset_order_tbl(tbl = annotated_regions, col='annot_type', col_order=annotation_order)

  ########################################################################
  # Construct the plot

    # Make the base ggplot
    # NOTE: binwidth may need to be a parameter
    plot =
      ggplot(annotated_regions, aes(x=annot_type)) +
      geom_bar() +
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

#' Visualize pair-wise annotations across regions
#'
#' All co-occurring annotations associated with a region are computed and displayed as a heatmap.
#'
#' @param annotated_regions The \code{tbl_df} result of \code{annotate_regions()}.
#' @param annotation_order A character vector which doubles as the subset of annotations desired for visualization as well as the ordering. If \code{NULL}, all annotations are displayed.
#' @param plot_title A string used for the title of the plot. Default \code{NULL}, no title is displayed.
#' @param axes_label A string used for the axis labels. Default \code{NULL}, corresponding variable name used.
#'
#' @return A \code{ggplot} object which can be viewed by calling it, or saved with \code{ggplot2::ggsave}.
#'
#' dm = system.file('extdata', 'IDH2mut_v_NBM_multi_data_chr9.txt.gz', package = 'annotatr')
#' annotations = c('hg19_basicgenes','hg19_cpgs','hg19_enhancers_fantom')
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
#' all_order = c(
#'   'hg19_cpg_islands',
#'   'hg19_cpg_shores',
#'   'hg19_cpg_shelves',
#'   'hg19_cpg_inter',
#'   'hg19_enhancers_fantom',
#'   'hg19_knownGenes_1to5kb',
#'   'hg19_knownGenes_promoters',
#'   'hg19_knownGenes_5UTRs',
#'   'hg19_knownGenes_exons',
#'   'hg19_knownGenes_introns',
#'   'hg19_knownGenes_3UTRs')
#'
#' chip_vcas = visualize_coannotations(
#'   annotated_regions = dm_r,
#'   annotation_order = all_order,
#'   axes_label = 'Annotations',
#'   plot_title = 'Co-occurrence of Annotations')
#'
#' @export
visualize_coannotations = function(annotated_regions, annotation_order=NULL,
  plot_title=NULL, axes_label=NULL) {

  ########################################################################
  # Argument parsing and error handling

    if(class(annotated_regions)[1] != "tbl_df") {
      stop('Error: annotated_regions must have class tbl_df. The best way to ensure this is to pass the result of annotate_regions() into this function.')
    }

  ########################################################################
  # Order and subset the annotations
  annotated_regions = subset_order_tbl(tbl = annotated_regions, col='annot_type', col_order=annotation_order)

  ########################################################################
  # Find the co-annotations

  annotation_pairs_by_region = dplyr::do(
    dplyr::group_by(annotated_regions, data_chrom, data_start, data_end),
    expand.grid(.$annot_type, .$annot_type, stringsAsFactors=F))

  pairwise_annotation_counts = table(annotation_pairs_by_region[['Var1']], annotation_pairs_by_region[['Var2']])

  pac_m = reshape2::melt(pairwise_annotation_counts, value.name = 'Counts')

  ########################################################################
  # Construct the plot

    # Make the base ggplot
    # NOTE: binwidth may need to be a parameter
    plot = ggplot(pac_m, aes(Var1, Var2)) +
      geom_raster(aes(fill = Counts)) +
      geom_text(aes(fill = Counts, label = Counts)) +
      scale_fill_gradient(low = "white", high = "steelblue") +
      theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.text.y = element_text(angle = 30, hjust = 1))

    # Add any user defined labels to the plot if their values are not NULL
    # if they are NULL, ggplot() will use defaults
    if(!is.null(plot_title)) {
      plot = plot + ggtitle(plot_title)
    }
    if(!is.null(axes_label)) {
      plot = plot + xlab(axes_label)
      plot = plot + ylab(axes_label)
    }

  return(plot)
}

#' Visualize numerical data over regions or regions summarized over annotations
#'
#' The \code{facet_order} vector also determines the categories used to create the background distribution that appears overlayed in red in each facet.
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
#' # NOTE: Background distribution (everything in \code{facet_order})
#' # is plotted in each facet for comparison.
#' dm_vs_sumnum = visualize_numerical(
#'   tbl = dm_sn,
#'   x = 'mean',
#'   facet = 'annot_type',
#'   facet_order = c('hg19_cpg_islands','hg19_cpg_shores','hg19_cpg_shelves','hg19_cpg_inter'),
#'   bin_width = 5,
#'   plot_title = 'Mean Meth. Diff. over CpG Annots.',
#'   x_label = 'Methylation Difference')
#'
#' # Plot histogram of group 1 methylation rates across the CpG annotations.
#' # NOTE: Background distribution (everything in \code{facet_order})
#' # is plotted in each facet for comparison.
#' dm_vs_regions_mu1 = visualize_numerical(
#'  tbl = dm_r,
#'  x = 'mu1',
#'  facet = 'annot_type',
#'  facet_order = c('hg19_cpg_islands','hg19_cpg_shores','hg19_cpg_shelves','hg19_cpg_inter'),
#'  bin_width = 5,
#'  plot_title = 'Group 1 Methylation over CpG Annotations',
#'  x_label = 'Group 1 Methylation')
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
  # Order and subset the annotations
  sub_tbl = subset_order_tbl(tbl = tbl, col = facet, col_order = facet_order)

  ########################################################################
  # Construct the plot
  # Note, data must be dplyr::ungroup()-ed before hand for the proper
  # display of the background distribution.

    if(is.null(y)) {
      # Make the base histogram ggplot
      plot =
        ggplot(data = ungroup(sub_tbl), aes_string(x=x, y='..density..')) +
        geom_histogram(binwidth=bin_width, fill = 'black', alpha = 0.75) +
        facet_wrap( as.formula(paste("~", facet)) ) + # Over the facets
        geom_histogram(data = select(ungroup(tbl), -matches(facet)),
          binwidth=bin_width, fill = 'red', alpha = 0.5) + # All the data
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    } else {
      # Make the base scatter ggplot
      plot = ggplot(sub_tbl, aes_string(x=x, y=y)) +
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

#' Visualize numerical data occurring in pairs of annotations
#'
#' Visualize numerical data associated with regions occurring in \code{annot1}, \code{annot2} and in both.
#'
#' @param tbl A \code{dplyr::tbl} returned from \code{annotate_regions()}. If the data is not summarized, the data is at the region level. If it is summarized, it represents the average or standard deviation of the regions by the character vector used for \code{by} in \code{summarize_numerical()}.
#' @param x A string indicating the column of the \code{tbl} to use for the x-axis.
#' @param y A string indicating the column of the \code{tbl} to use for the y-axis. Default is \code{NULL}, meaning a histogram over \code{x} will be plotted. If it is not \code{NULL}, a scatterplot is plotted.
#' @param annot1 A string indicating the first annotation type.
#' @param annot2 A string indicating the second annotation type.
#' @param bin_width An integer indicating the bin width of the histogram used for score. Default 10. Select something appropriate for the data. NOTE: This is only used if \code{y} is \code{NULL}.
#' @param plot_title A string used for the title of the plot. Default \code{NULL}, no title displayed.
#' @param x_label A string used for the x-axis label. Default \code{NULL}, corresponding variable name used.
#' @param y_label A string used for the y-axis label. Default \code{NULL}, corresponding variable name used.
#'
#' @return A \code{ggplot} object which can be viewed by calling it, or saved with \code{ggplot2::ggsave}.
#'
#' @examples
#'
#' dm = system.file('extdata', 'IDH2mut_v_NBM_multi_data_chr9.txt.gz', package = 'annotatr')
#' annotations = c('hg19_basicgenes','hg19_cpgs','hg19_enhancers_fantom')
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
#' dm_vs_num_co = visualize_numerical_coannotations(
#'   tbl = dm_r,
#'   x = 'mu0',
#'   annot1 = 'hg19_cpg_islands',
#'   annot2 = 'hg19_knownGenes_promoters',
#'   bin_width = 5,
#'   plot_title = 'Group 0 Perc. Meth. in CpG Islands and Promoters',
#'   x_label = 'Percent Methylation')
#'
#' @export
visualize_numerical_coannotations = function(tbl, x, y=NULL, annot1, annot2, bin_width=10,
  plot_title=NULL, x_label=NULL, y_label=NULL) {

  ########################################################################
  # Argument parsing and error handling

    if(class(tbl)[1] != "tbl_df") {
      stop('Error: tbl must have class tbl_df. The best way to ensure this is to pass the result of annotate_regions() into this function.')
    }

  ########################################################################
  # Order and subset the annotations
  annotation_order = c(annot1,annot2)
  sub_tbl = subset_order_tbl(tbl = tbl, col='annot_type', col_order=annotation_order)

  ########################################################################
  # Find the co-annotations

  pairs_by_region = dplyr::do(
    dplyr::group_by(sub_tbl, data_chrom, data_start, data_end),
    expand.grid(annot1 = .$annot_type, annot2 = .$annot_type, stringsAsFactors=F))

  # Join on the data chromosome locations
  pairs_by_region = dplyr::inner_join(x = pairs_by_region, y = sub_tbl, by = c('data_chrom','data_start','data_end'))

  # pairs_by_region = dplyr::mutate(
  #   pairs_by_region,
  #   pair_combo = paste(unique(annot1, annot2), collapse=' and '))

  ########################################################################
  # Construct the plot
  # Note, data must be dplyr::ungroup()-ed before hand for the proper
  # display of the background distribution.

    if(is.null(y)) {
      # Make the base histogram ggplot
      plot =
        ggplot(data = ungroup(pairs_by_region), aes_string(x=x, y='..density..')) +
        geom_histogram(binwidth=bin_width, fill = 'black', alpha = 0.75) +
        facet_wrap( annot1 ~ annot2 ) + # Over the facets
        geom_histogram(data = ungroup(tbl),
          binwidth=bin_width, fill = 'red', alpha = 0.5) + # All the data
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    } else {
      # Make the base scatter ggplot
      plot = ggplot(pairs_by_region, aes_string(x=x, y=y)) +
        geom_point(alpha = 1/8, size = 1) +
        facet_wrap( annot1 ~ annot2 ) +
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
#' Given a \code{dplyr::tbl_df} of annotated regions (from \code{annotate_regions()}), visualize the the distribution of \code{fill} in \code{x}. A bar representing the distribution of all \code{fill} in \code{x} will be added according to the contents of \code{fill}. This serves as a background to compare against.
#'
#' @param annotated_regions The \code{grouped_df} result of \code{annotated_regions()}.
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
#'   annotated_regions = dm_r,
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
visualize_categorical = function(annotated_regions, x, fill=NULL, x_order=NULL, fill_order=NULL,
  position = 'stack', plot_title=NULL, legend_title=NULL, x_label=NULL, y_label=NULL) {

  ########################################################################
  # Argument parsing and error handling

    # Check correct class of input
    if(class(annotated_regions)[1] != "tbl_df") {
      stop('Error: annotated_regions must have class tbl_df. The best way to ensure this is to pass the result of summarize_categorical() into this function.')
    }

    # Ensure the value of x is a column name in summarized_cats
    if( !(x %in% colnames(annotated_regions)) ) {
      stop('The column name used for x does not exist in annotated_regions.')
    }

    # Ensure the value of fill is a column name in summarized_cats if it isn't NULL
    # Also ensure fill != x
    if( !is.null(fill) ) {
      if( !(fill %in% colnames(annotated_regions)) ) {
        stop('The column name used for fill does not exist in annotated_regions.')
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
  # Order and subset based on fill_order
  annotated_regions = subset_order_tbl(tbl = annotated_regions, col = fill, col_order = fill_order)

  ########################################################################
  # Order and subset based on x_order
  if(is.null(x_order)) {
    x_order = unique(annotated_regions[[x]])
  }
  sub_annot_regions = subset_order_tbl(tbl = annotated_regions, col = x, col_order = x_order)

  ########################################################################
  # Construct the plot

    # Make base ggplot
    plot =
      ggplot(annotated_regions, aes(x='All')) +
      geom_bar(aes_string(fill=fill), position=position, width=0.5) + # The All bar
      geom_bar(data = sub_annot_regions, aes_string(x=x, fill=fill), position=position, width=0.5) + # The subsets bars
      theme(axis.text.x = element_text(angle = 30, hjust = 1))

    # Change the fill scale and name if legend_title isn't null
    if(!is.null(legend_title)) {
      plot = plot + scale_fill_brewer(name=legend_title)
    } else {
      plot = plot + scale_fill_brewer()
    }

    # Deal with the x-axis labels to make sure the order is correct
    if(x == 'annot_type') {
      plot = plot + scale_x_discrete(limits = c('All', names(tidy_annotations(x_order))))
    } else {
      plot = plot + scale_x_discrete(limits = c('All', x_order))
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
