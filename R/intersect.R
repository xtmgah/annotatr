#' A function to intersect user data with annotation data
#'
#' Given a GRanges object constructed from user supplied data with read_bed()(),
#' and desired annotation overlaps, use GenomicRanges::findOverlaps() to return a
#' Hits object indicating which elements of the user-data (query) intersect which
#' elements of the annotation data (subject).
#'
#' @param regions The GRanges object of the user-data returned by read_bed().
#' @param annotations A character vector of annotations to overlap with user-data.
#' Valid annotation codes can be found with supported_annotations(). The "basicgenes" shortcut annotates regions to the 1-5Kb, promoter, 5UTR, exon, intron, and 3UTR knownGene regions. The "detailedgenes" shortcut annotates regions to the 1-5Kb, promoter, 5UTR exon/intron, CDS exon/intron, and 3UTR exon/intron knownGene regions. The "cpgs" shortcut annotates regions to the CpG islands, shores, shelves, and interCGI regions. NOTE: basicgenes and detailedgenes annotations cannot be done at the same time, and shortcuts need to be appended by the genome.
#' Custom annotations whose names are of the form [genome]_custom_[annotation_name] can also be used. Custom annotations should be read in with \code{read_annotations()}. They can be for a \code{supported_genome()}, or for an unsupported genome.
#' @param ignore.strand Logical variable indicating whether strandedness should be respected in findOverlaps(). Default FALSE.
#' @param use.score Logical variable. Include the "score" for each genomic region in the tabulated results. Score can mean a variety of things, e.g. percent methylation of a CpG/region or fold-change of a ChIP-seq peak.
#'
#' @return A \code{dplyr::tbl_df} with columns from the GenomicRanges object for the regions and corresponding annotations.
#'
#' @examples
#' # A very simple example with only 3 genomic regions
#' bed = system.file('extdata', 'test_intersect.bed', package = 'annotatr')
#' annotations = c('hg19_cpg_islands','hg19_cpg_shores','hg19_knownGenes_promoters')
#'
#' d = read_bed(file = bed, genome = 'hg19', stranded = FALSE)
#'
#' i = annotate_regions(
#'   regions = d,
#'   annotations = annotations,
#'   ignore.strand = TRUE,
#'   use.score = FALSE)
#'
#' # A more complicated example using Gm12878 Pol2 ChIP-seq from ENCODE and an annotation shortcut
#' bed = system.file('extdata', 'Gm12878_Pol2.narrowPeak.gz', package = 'annotatr')
#' annotations = c('hg19_basicgenes')
#'
#' d = read_bed(file = bed, genome = 'hg19', stranded = FALSE)
#'
#' i = annotate_regions(
#'   regions = d,
#'   annotations = annotations,
#'   ignore.strand = TRUE,
#'   use.score = FALSE)
#'
#' # An example with a combination of custom and built-in annotaitons
#' a_file = system.file('extdata', 'test_annotations_3.bed', package='annotatr')
#' r_file = system.file('extdata', 'test_read_multiple_data_head.bed', package='annotatr')
#'
#' r = read_bed(file=r_file, col.names=TRUE, genome='hg19', stranded=FALSE, use.score=FALSE)
#' hg19_custom_TFBS = read_annotations(file=a_file, genome='hg19', annotation_name='TFBS')
#'
#' annotations = c('hg19_custom_TFBS', 'hg19_cpgs')
#'
#' a = annotate_regions(
#'   regions = r,
#'   annotations = annotations,
#'   ignore.strand = TRUE,
#'   use.score = TRUE
#'   )
#'
#' @export
annotate_regions = function(regions, annotations, ignore.strand = TRUE, use.score) {
  # Checks before moving forward
  if(class(regions)[1] != "GRanges") {
    stop('Error in annotate_regions(...): regions object is not GRanges.')
  }

  # Check annotations and expand any shortcuts
  check_annotations(annotations)
  annotations = expand_annotations(annotations)

  # Split the custom and non custom annotations
  custom_annotations = grep('custom', annotations, value=TRUE)
  non_custom_annotations = setdiff(annotations, custom_annotations)

  # Collect the non_custom_annotation objects into a GRangesList
  if(length(non_custom_annotations) > 0) {
    data(list = non_custom_annotations, package = 'annotatr')
  }

  # Perform the intersections in an lapply (consider using mclapply)
  intersections = lapply(annotations, function(annot){
    message(sprintf('Intersecting %s', annot))

    GenomicRanges::findOverlaps(regions, get(annot), ignore.strand = ignore.strand)
  })
  names(intersections) = annotations

  # Pull in data from the regions and the annotations
  # Think about using mclapply instead of lapply to speed up. May cause some
  # memory issues or overhead.
  if(use.score) {
    tab_list = lapply(names(intersections), function(n){
      message(sprintf('Annotating %s with score', n))

      if(length(intersections[[n]]) > 0) {

        # Get subsets of
        r_sub = regions[intersections[[n]]@queryHits]
        a_sub = get(n)[intersections[[n]]@subjectHits]

        # Create the dplyr::tbl_df with score
        df_ra = dplyr::tbl_df(data.frame(
          annot_start = start(a_sub),
          annot_end = end(a_sub),
          annot_strand = as.character(strand(a_sub)),
          annot_type = n,
          annot_id = a_sub$ID,
          data_chrom = as.character(seqnames(r_sub)),
          data_start = start(r_sub),
          data_end = end(r_sub),
          data_strand = as.character(strand(r_sub)),
          stringsAsFactors = F))

        df_d = dplyr::tbl_df(as.data.frame(GenomicRanges::mcols(r_sub)))

        df = dplyr::bind_cols(df_ra, df_d)

        return(df)
      } else {
        return(NULL)
      }
    })
  } else {
    tab_list = lapply(names(intersections), function(n){
      message(sprintf('Annotating %s without score', n))

      if(length(intersections[[n]]) > 0) {

        # Get subsets of
        r_sub = regions[intersections[[n]]@queryHits]
        a_sub = get(n)[intersections[[n]]@subjectHits]

        # Create the dplyr::tbl_df without score
        df = dplyr::tbl_df(data.frame(
          annot_start = start(a_sub),
          annot_end = end(a_sub),
          annot_strand = as.character(strand(a_sub)),
          annot_type = n,
          annot_id = a_sub$ID,
          data_chrom = as.character(seqnames(r_sub)),
          data_start = start(r_sub),
          data_end = end(r_sub),
          data_strand = as.character(strand(r_sub)),
          name = r_sub$name,
          stringsAsFactors=F
        ))

        return(df)
      } else {
        return(NULL)
      }
    })
  }

  # Combine and sort the list of data.frames into a single data.frame
  message('Combining annotations')
  tab = dplyr::bind_rows(tab_list)
  message('Sorting annotations')
  tab = dplyr::arrange(tab, data_chrom, data_start, data_end, annot_start)

  return(tab)
}
