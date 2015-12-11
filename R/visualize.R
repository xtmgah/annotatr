#' Visualize the number of regions per annotation
#'
#' Given a \code{dplyr::tbl_df} of counts of regions per annotation (from \code{summarize_annotation()}), visualize the counts per annotation as a bar graph.
#'
#' @param summarized_annotations The \code{tbl_df} result of \code{summarize_annotation()}.
#' @param annotation_order A character vector which doubles as the subset of annotations desired for visualization as well as the ordering.
#'
#' @return A \code{ggplot} object which can be viewed by calling it, or saved with \code{ggplot2::ggsave}.
#'
#' @examples
#' ########################################################################
#' # An example of ChIP-seq peaks with signalValue used for score
#' bed = system.file('extdata', 'Gm12878_Ezh2_sorted_scores.narrowPeak.gz', package = 'annotatr')
#' annotations = c('basic_genes','cpgs')
#'
#' d = read_bed(filename = bed, genome = 'hg19', stranded = FALSE, use.score = TRUE)
#'
#' i = annotate_regions(
#'   regions = d,
#'   annotations = annotations,
#'   genome = 'hg19',
#'   ignore.strand = TRUE,
#'   use.score = TRUE)
#'
#' s = summarize_annotation(i)
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
visualize_annotation = function(summarized_annotations, annotation_order) {
  if(class(summarized_annotations)[1] != "tbl_df") {
    stop('Error: summarized_annotations must have class tbl_df. The best way to ensure this is to pass the result of summarize_annotations() into this function.')
  }

  # Collect all annotation types
  all_annot_types = unique(summarized_annotations[['annot_type']])

  # Check set equality of annot_types in the summarized_scores and the annotation_order
  if( !dplyr::setequal(all_annot_types, annotation_order) ) {
    if( all(annotation_order %in% all_annot_types) ) {
      summarized_annotations = subset(summarized_annotations, summarized_annotations[['annot_type']] %in% annotation_order)
    } else {
      stop('There are annotations in annotation_order that are not present in annot_type column of summarized_annotations.')
    }
  }

  # Convert annot_type to a vector with the levels in the desired order
  summarized_annotations[['annot_type']] = factor(
    summarized_annotations[['annot_type']],
    levels = annotation_order)

  # Make the ggplot
  # NOTE: binwidth may need to be a parameter
  plot =
    ggplot(summarized_annotations, aes(x=annot_type, y=n)) +
    geom_bar(stat='identity') +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  return(plot)
}

#' Visualize score distributions over annotations
#'
#' Given a \code{dplyr::grouped_df} of score aggregated annotations (from \code{summarize_score()}), visualize the distribution of scores over the annotation types.
#'
#' @param summarized_scores The \code{grouped_df} result of \code{summarize_score()}.
#' @param annotation_order A character vector which doubles as the subset of annotations desired for visualization as well as the ordering.
#'
#' @return A \code{ggplot} object which can be viewed by calling it, or saved with \code{ggplot2::ggsave}.
#'
#' @examples
#' ########################################################################
#' # An example of ChIP-seq peaks with signalValue used for score
#' bed = system.file('extdata', 'Gm12878_Ezh2_sorted_scores.narrowPeak.gz', package = 'annotatr')
#' annotations = c('basic_genes','cpgs')
#'
#' d = read_bed(filename = bed, genome = 'hg19', stranded = FALSE, use.score = TRUE)
#'
#' i = annotate_regions(
#'   regions = d,
#'   annotations = annotations,
#'   genome = 'hg19',
#'   ignore.strand = TRUE,
#'   use.score = TRUE)
#'
#' s = summarize_score(i)
#'
#' cpgs_order = c(
#'   'hg19_cpg_islands',
#'   'hg19_cpg_shores',
#'   'hg19_cpg_shelves',
#'   'hg19_cpg_inter')
#' v_cpgs = visualize_score(s, cpgs_order)
#'
#' genes_order = c(
#'   'hg19_knownGenes_1to5kb',
#'   'hg19_knownGenes_promoters',
#'   'hg19_knownGenes_5UTRs',
#'   'hg19_knownGenes_exons',
#'   'hg19_knownGenes_introns',
#'   'hg19_knownGenes_3UTRs')
#' v_genes = visualize_score(s, genes_order)
#'
#' ########################################################################
#' # An example of percent methylation at CpG sites
#' bed = system.file('extdata', '2607_mc_hmc_perc_meth_chr21.txt.gz', package = 'annotatr')
#' annotations = c('basic_genes','cpgs')
#'
#' d = read_bed(filename = bed, genome = 'hg19', stranded = FALSE, use.score = TRUE)
#'
#' i = annotate_regions(
#'   regions = d,
#'   annotations = annotations,
#'   genome = 'hg19',
#'   ignore.strand = TRUE,
#'   use.score = TRUE)
#'
#' s = summarize_score(i)
#'
#' cpgs_order = c(
#'   'hg19_cpg_islands',
#'   'hg19_cpg_shores',
#'   'hg19_cpg_shelves',
#'   'hg19_cpg_inter')
#' v_cpgs = visualize_score(s, cpgs_order)
#'
#' genes_order = c(
#'   'hg19_knownGenes_1to5kb',
#'   'hg19_knownGenes_promoters',
#'   'hg19_knownGenes_5UTRs',
#'   'hg19_knownGenes_exons',
#'   'hg19_knownGenes_introns',
#'   'hg19_knownGenes_3UTRs')
#' v_genes = visualize_score(s, genes_order)
#'
#' @export
visualize_score = function(summarized_scores, annotation_order) {
  if(class(summarized_scores)[1] != "grouped_df") {
    stop('Error: summarized_scores must have class grouped_df. The best way to ensure this is to pass the result of summarize_score() into this function.')
  }

  # Collect all annotation types
  all_annot_types = unique(summarized_scores[['annot_type']])

  # Check set equality of annot_types in the summarized_scores and the annotation_order
  if( !dplyr::setequal(all_annot_types, annotation_order) ) {
    if( all(annotation_order %in% all_annot_types) ) {
      summarized_scores = subset(summarized_scores, summarized_scores[['annot_type']] %in% annotation_order)
    } else {
      stop('There are annotations in annotation_order that are not present in annot_type column of summarized_scores.')
    }
  }

  # Convert annot_type to a vector with the levels in the desired order
  summarized_scores[['annot_type']] = factor(
    summarized_scores[['annot_type']],
    levels = annotation_order)

  # Make the ggplot
  # NOTE: binwidth may need to be a parameter
  plot =
    ggplot(summarized_scores, aes(mean)) +
    geom_histogram(binwidth=5, aes(y=..density..)) +
    facet_wrap(~annot_type) +
    theme_bw()

  return(plot)
}

#' Visualize scores over annotations
#'
#' Given a \code{dplyr::grouped_df} of name aggregated annotations (from \code{summarize_name()}), visualize the the distribution of \code{annot_type} in \code{data_name}.
#'
#' @param summarized_names The \code{grouped_df} result of \code{summarize_name()}.
#' @param annotation_order A character vector which doubles as the subset of annotations desired for visualization as well as the ordering.
#' @param data_order A character vector which doubles as the subset of data_name's desired for visualization as well as the ordering.
#' @param fill A boolean indicating whether \code{ggplot::geom_bar(..., position='fill')} or not. When TRUE, the barplot is normalized across the categories on the x-axis, and when FALSE, the barplot reflects counts of the categories on the x-axis. Default is FALSE.
#'
#' @return A \code{ggplot} object which can be viewed by calling it, or saved with \code{ggplot2::ggsave}.
#'
#' @examples
#' # An example of differentially methylated regions classified as DM up, DM down, or no DM
#' bed = system.file('extdata', 'IDH2mut_v_NBM_names_scores_chr9.txt.gz', package = 'annotatr')
#' annotations = c('basic_genes','cpgs')
#'
#' d = read_bed(filename = bed, genome = 'hg19', stranded = FALSE, use.score = TRUE)
#'
#' i = annotate_regions(
#'   regions = d,
#'   annotations = annotations,
#'   genome = 'hg19',
#'   ignore.strand = TRUE,
#'   use.score = TRUE)
#'
#' s = summarize_name(i)
#'
#' cpgs_order = c(
#'   'hg19_cpg_islands',
#'   'hg19_cpg_shores',
#'   'hg19_cpg_shelves',
#'   'hg19_cpg_inter')
#' data_order = c(
#'   'DMup',
#'   'DMdown')
#' v_cpgs_counts = visualize_name(s, cpgs_order, data_order, fill=FALSE)
#' v_cpgs_proportions = visualize_name(s, cpgs_order, data_order, fill=TRUE)
#'
#' genes_order = c(
#'   'hg19_knownGenes_1to5kb',
#'   'hg19_knownGenes_promoters',
#'   'hg19_knownGenes_5UTRs',
#'   'hg19_knownGenes_exons',
#'   'hg19_knownGenes_introns',
#'   'hg19_knownGenes_3UTRs')
#' data_order = c(
#'   'DMup',
#'   'DMdown')
#' v_genes_counts = visualize_name(s, genes_order, data_order, fill=FALSE)
#' v_genes_proportions = visualize_name(s, genes_order, data_order, fill=TRUE)
#'
#' @export
visualize_name = function(summarized_names, annotation_order, data_order, fill = FALSE) {
  if(class(summarized_names)[1] != "grouped_df") {
    stop('Error: summarized_names must have class grouped_df. The best way to ensure this is to pass the result of summarize_name() into this function.')
  }

  # Collect all annotation types
  all_annot_types = unique(summarized_names[['annot_type']])
  all_data_names = unique(summarized_names[['data_name']])

  # Check set equality of annot_types in the summarized_scores and the annotation_order
  if( !dplyr::setequal(all_annot_types, annotation_order) ) {
    if( all(annotation_order %in% all_annot_types) ) {
      summarized_names = subset(summarized_names, summarized_names[['annot_type']] %in% annotation_order)
    } else {
      stop('There are annotations in annotation_order that are not present in the annot_type column of summarized_scores.')
    }
  }

  # Check set equality of data_names in the summarized_scores and the data_order
  if( !dplyr::setequal(all_data_names, data_order) ) {
    if( all(data_order %in% all_data_names) ) {
      summarized_names = subset(summarized_names, summarized_names[['data_name']] %in% data_order)
    } else {
      stop('There are data_names in data_order that are not present in the data_names column of summarized_names.')
    }
  }

  # Convert annot_type to a vector with the levels in the desired order
  summarized_names[['annot_type']] = factor(
    summarized_names[['annot_type']],
    levels = annotation_order)

  summarized_names[['data_name']] = factor(
    summarized_names[['data_name']],
    levels = data_order)

  # Make the ggplot
  if(fill) {
    plot =
      ggplot(summarized_names, aes(x=data_name, y=n, fill=annot_type)) +
      geom_bar(stat='identity', position='fill', aes(order = annot_type)) +
      scale_fill_brewer('Blues') +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
  } else {
    plot =
      ggplot(summarized_names, aes(x=data_name, y=n, fill=annot_type)) +
      geom_bar(stat='identity', aes(order = annot_type)) +
      scale_fill_brewer('Blues') +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
  }

  return(plot)
}
