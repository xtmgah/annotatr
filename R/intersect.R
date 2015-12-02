#' A function to intersect user data with annotation data
#'
#' Given a GRanges object constructed from user supplied data with read.bed(),
#' and desired annotation overlaps, use GenomicRanges::findOverlaps() to return a
#' Hits object indicating which elements of the user-data (query) intersect which
#' elements of the annotation data (subject).
#'
#' @param regions The GRanges object of the user-data returned by read.bed().
#' @param annotations A character vector of annotations to overlap with user-data. Valid annotation codes can be found with supported_annotations(). The "basic_genes" shortcut annotates regions to the 1-5Kb, promoter, 5UTR, CDS, and 3UTR knownGene regions. The "detailed_genes" shortcut annotates regions to the 1-5Kb, promoter, 5UTR exon/intron, CDS exon/intron, and 3UTR exon/intron knownGene regions. The "cpgs" shortcut annotates regions to the CpG islands, shores, shelves, and interCGI regions.
#' @param genome One of the genomes in supported_genomes(). Should match argument used in read.bed().
#' @param ignore.strand A boolean indicating whether strandedness should be respected in findOverlaps().
#'
#' @return A list of Hits objects giving the indices of objects overlapping in query and subject for each annotation.
#'
#' @examples
#' bed = system.file('extdata', 'test_intersect.bed', package = 'annotatr')
#' annotations = c('hg19_cpg_islands','hg19_cpg_shores','hg19_knownGenes_promoters')
#'
#' d = read_bed(filename = bed, genome = 'hg19', stranded = FALSE)
#'
#' i = intersect_annotations(
#'   regions = d,
#'   annotations = annotations,
#'   genome = 'hg19',
#'   ignore.strand = TRUE)
#'
#' @export
intersect_annotations = function(regions, annotations, genome, ignore.strand) {
  # Checks before moving forward
  if(class(regions)[1] != "GRanges") {
    stop('Error in intersect_annotations(...): regions object is not GRanges.')
  }

  if(!(genome %in% supported_genomes())) {
    stop('Error in intersect_annotations(...): unsupported genome given. See supported_genomes().')
  }

  # Check that the annotations are supported, tell the user which are unsupported
  if(!all(annotations %in% c(supported_annotations(),'cpgs','basic_genes','detailed_genes'))) {
    unsupported = setdiff(annotations, c(supported_annotations(),'cpgs','basic_genes','detailed_genes'))

    stop(sprintf('Error in intersect_annotations(...): "%s" is(are) not supported. See supported_annotations().',
      paste(unsupported, collapse=', ')))
  }

  # Collect the annotation objects into a GRangesList
  data(list = annotations, package = 'annotatr')
  annot_list = lapply(annotations, get)
  names(annot_list) = annotations

  # Perform the intersections in an lapply (consider using mclapply)
  overlaps = lapply(annot_list, function(annot){
    GenomicRanges::findOverlaps(regions, annot, ignore.strand = ignore.strand)
  })

  return(overlaps)
}
