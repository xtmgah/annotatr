# Functions relating to using the intersections to annotate regions to annotations

#' A function to create tabular output of regions intersected with annotations
#'
#' Given output from \code{intersect_annotations()}, create a table for output which shows peak information and corresponding annotation information.
#'
#' @param intersections A list of Hits objects from \code{intersect_annotations()}.
#'
#' @return A table with columns from the GenomicRanges object for the regions and corresponding annotations
#'
#' @examples
#' bed = system.file('extdata', 'Gm12878_Pol2.narrowPeak.gz', package = 'annotatr')
#' annotations = c('basic_genes')
#'
#' d = read_bed(filename = bed, genome = 'hg19', stranded = FALSE)
#'
#' i = intersect_annotations(
#'   regions = d,
#'   annotations = annotations,
#'   genome = 'hg19',
#'   ignore.strand = TRUE)
#'
#' t = tabulate_intersections(
#'   regions = d,
#'   intersections = i)
#'
#' @export
tabulate_intersections = function(regions, intersections) {
  # Checks before moving forward
  if(class(regions)[1] != "GRanges") {
    stop('Error in tabulate_intersections(...): regions object is not GRanges.')
  }

  if(class(i) != "list") {
    stop('Error in tabulate_intersections(...): intersections must be a list.')
  }

  if(unique(sapply(i, class)) != "Hits") {
    stop('Error in tabulate_intersections(...): intersections must be a list of Hits objects.')
  }

  # In the future, we could think about changing this to mclapply to speed it up
  tab_list = lapply(names(intersections), function(n){
    message(sprintf('Tabulating %s', n))

    # Get subsets of
    r_sub = regions[intersections[[n]]@queryHits]
    a_sub = get(n)[intersections[[n]]@subjectHits]

    # Create the data.frame
    df = data.frame(
      data_chrom = as.character(seqnames(r_sub)),
      data_start = start(r_sub),
      data_end = end(r_sub),
      data_strand = as.character(strand(r_sub)),
      annot_chrom = as.character(seqnames(a_sub)),
      annot_start = start(a_sub),
      annot_end = end(a_sub),
      annot_strand = as.character(strand(a_sub)),
      annot_type = n,
      stringsAsFactors=F
    )
  })

  tab = Reduce(rbind, tab_list)

  return(tab)
}
