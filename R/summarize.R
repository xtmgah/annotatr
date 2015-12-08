#' Summarize scores over annotations
#'
#' Given a dplyr::tbl_df of annotated regions, summarize scores over (annot_type, annot_id) pairs.
#'
#' @param annotated_regions The tbl_df result of \code{annotate_intersections()}.
#'
#' @return A tbl_df of the mean of scores over (annot_type, annot_id) pairs.
#'
#' @examples
#' # An example of ChIP-seq peaks with signalValue used for score
#' bed = system.file('extdata', 'Gm12878_Ezh2_sorted_scores.narrowPeak.gz', package = 'annotatr')
#' annotations = c('basic_genes','cpgs')
#'
#' d = read_bed(filename = bed, genome = 'hg19', stranded = F, use.score = TRUE)
#'
#' i = intersect_annotations(
#'   regions = d,
#'   annotations = annotations,
#'   genome = 'hg19',
#'   ignore.strand = T)
#'
#' t = annotate_intersections(
#'   regions = d,
#'   intersections = i,
#'   use.score = TRUE)
#'
#' s = summarize_score(t)
#'
#' @export
summarize_score = function(annotated_regions) {
  if(class(annotated_regions)[1] != "tbl_df") {
    stop('Error: annotated_regions must have class tbl_df. The best way to ensure this is to pass the result of annotate_intersections() into this function.')
  }

  message('Summarizing scores over (annotation type, annotation ID) pairs')
  agg = dplyr::summarize(
    dplyr::group_by(annotated_regions, annot_type, annot_id),
    count = length(data_score),
    mean = mean(data_score),
    sd = sd(data_score)
  )

  return(agg)
}

#' Summarize scores over annotations
#'
#' Given a dplyr::tbl_df of annotated regions, summarize names over (annot_type, data_name) pairs. The assumption is that the name column in the input BED file consists of a set of non-unique names effectively serving as a label for the region.
#'
#' @param annotated_regions The tbl_df result of \code{annotate_intersections()}.
#'
#' @return A tbl_df of the counts of (annot_type, data_name) pairs.
#'
#' @examples
#' # An example of differentially methylated regions classified as DM up, DM down, or no DM
#' bed = system.file('extdata', 'IDH2mut_v_NBM_names_scores_chr9.txt.gz', package = 'annotatr')
#' annotations = c('basic_genes','cpgs')
#'
#' d = read_bed(filename = bed, genome = 'hg19', stranded = F, use.score = TRUE)
#'
#' i = intersect_annotations(
#'   regions = d,
#'   annotations = annotations,
#'   genome = 'hg19',
#'   ignore.strand = T)
#'
#' t = annotate_intersections(
#'   regions = d,
#'   intersections = i,
#'   use.score = TRUE)
#'
#' s = summarize_name(t)
#'
#' @export
summarize_name = function(annotated_regions) {
  if(class(annotated_regions)[1] != "tbl_df") {
    stop('Error: annotated_regions must have class tbl_df. The best way to ensure this is to pass the result of annotate_intersections() into this function.')
  }

  message('Summarizing names over annotation type')
  agg = dplyr::tally(
    dplyr::group_by(annotated_regions, annot_type, data_name)
  )

  return(agg)
}
