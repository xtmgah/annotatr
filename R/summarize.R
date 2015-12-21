#' Summarize annotation counts
#'
#' Given a \code{dplyr::tbl_df} of annotated regions, count the number of regions in each annotation.
#'
#' @param annotated_regions The \code{tbl_df} result of \code{annotate_regions()}.
#'
#' @return A \code{tbl_df} of the number of regions per annotation.
#'
#' @examples
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
#' @export
summarize_annotations = function(annotated_regions) {
  if(class(annotated_regions)[1] != "tbl_df") {
    stop('Error: annotated_regions must have class tbl_df. The best way to ensure this is to pass the result of annotate_regions() into this function.')
  }

  message('Summarizing scores over annotation types')
  agg = dplyr::tally(
    dplyr::group_by(annotated_regions, annot_type)
  )

  return(agg)
}

#' Summarize numerical data over groupings of annotated regions
#'
#' Given a \code{dplyr::tbl_df} of annotated regions, summarize numerical data columns based on a grouping.
#'
#' @param annotated_regions The \code{tbl_df} result of \code{annotate_regions()}.
#' @param by A character vector of the columns to group over. Default is \code{c(annot_type, annot_id)}.
#' @param over A character vector of which columns to \code{count}, take the \code{mean}, and take the \code{sd} over after grouping according to the \code{by} column. Default is \code{score}. NOTE: If more than one value is used, the naming scheme for the resuling \code{dplyr::tbl} summary columns are \code{COLNAME_n}, \code{COLNAME_mean}, \code{COLNAME_sd}. If \code{over} has length one, then the column names are \code{n}, \code{mean}, \code{sd}.
#'
#' @return A \code{dplyr::tbl} with groups, and the \code{count}, \code{mean}, and \code{sd} of the \code{cols} \code{by} the groupings.
#'
#' @examples
#' # Test on a very simple bed file to demonstrate different options
#' bed = system.file('extdata', 'test_read_multiple_data_head.bed', package = 'annotatr')
#' annotations = c('hg19_cpg_islands', 'hg19_knownGenes_promoters')
#'
#' d = read_bed(
#'   file = bed,
#'   col.names=TRUE,
#'   genome = 'hg19',
#'   stranded = FALSE,
#'   use.score = TRUE)
#'
#' i = annotate_regions(
#'   regions = d,
#'   annotations = annotations,
#'   ignore.strand = TRUE,
#'   use.score = TRUE)
#'
#' # Testing over normal by
#' sn1 = summarize_numerical(
#'   annotated_regions = i,
#'   by = c('annot_type', 'annot_id'),
#'   over = c('score', 'mu1', 'mu0'))
#'
#' # Testing over a different by
#' sn2 = summarize_numerical(
#'   annotated_regions = i,
#'   by = c('name'),
#'   over = c('score', 'mu1', 'mu0'))
#'
#' @export
summarize_numerical = function(annotated_regions, by = c('annot_type', 'annot_id'), over = 'score') {
  if(class(annotated_regions)[1] != "tbl_df") {
    stop('Error: annotated_regions must have class tbl_df. The best way to ensure this is to pass the result of annotate_regions() into this function.')
  }

  message(sprintf('Grouping regions by %s, and summarizing numerical data over %s',
    paste(by, collapse=' & '), paste(over, collapse=' & ')))
  agg = dplyr::summarize_each_(
    dplyr::group_by_(annotated_regions, .dots = by),
    dplyr::funs(n(), 'mean', 'sd'),
    over)

  return(agg)
}

#' Summarize categoricals over groupings of annotated regions
#'
#' Given a \code{dplyr::tbl_df} of annotated regions, summarize names over (annot_type, name) pairs. The assumption is that the name column in the input BED file consists of a set of non-unique names effectively serving as a label for the region.
#'
#' @param annotated_regions The \code{tbl_df} result of \code{annotate_regions()}.
#' @param by A character vector to group the data by and tally over. Default is \code{c(annot_type, annot_id)}.
#'
#' @return A \code{tbl_df} of the counts of groupings according to the \code{by} vector.
#'
#' @examples
#' bed = system.file('extdata', 'test_read_multiple_data_head.bed', package = 'annotatr')
#' annotations = c('hg19_cpg_islands', 'hg19_knownGenes_promoters')
#'
#' d = read_bed(
#'   file = bed,
#'   col.names=TRUE,
#'   genome = 'hg19',
#'   stranded = FALSE,
#'   use.score = TRUE)
#'
#' i = annotate_regions(
#'   regions = d,
#'   annotations = annotations,
#'   ignore.strand = TRUE,
#'   use.score = TRUE)
#'
#' sc1 = summarize_categorical(
#'   annotated_regions = i,
#'   by = c('annot_type', 'name'))
#' sc2 = summarize_categorical(
#'   annotated_regions = i,
#'   by = c('annot_type', 'diff_exp'))
#' sc3 = summarize_categorical(
#'   annotated_regions = i,
#'   by = c('name','diff_exp'))
#'
#' @export
summarize_categorical = function(annotated_regions, by = c('annot_type', 'annot_id')) {
  if(class(annotated_regions)[1] != "tbl_df") {
    stop('Error: annotated_regions must have class tbl_df. The best way to ensure this is to pass the result of annotate_regions() into this function.')
  }

  message(sprintf('Grouping regions by %s, and tallying',
    paste(by, collapse=' & ')))
  agg = dplyr::tally(
    dplyr::group_by_(annotated_regions, .dots = by))

  return(agg)
}
